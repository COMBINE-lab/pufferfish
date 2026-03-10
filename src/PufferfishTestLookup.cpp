#include "FastxParser.hpp"
#include <cmath>
#include <chrono>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>

#include "CLI/Timer.hpp"
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferFS.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "jellyfish/mer_dna.hpp"
#include "spdlog/spdlog.h"

#include "ProgOpts.hpp"
#include "PufferfishIndex.hpp"
#include "Util.hpp"

namespace kmers = combinelib::kmers;

std::vector<CanonicalKmer> get_kmers(const std::string& fasta_file, uint32_t k) {
  (void)k;

  std::vector <CanonicalKmer> kmers;
  std::vector<std::string> read_file = {fasta_file};
  auto cfg = fastx_parser::ParserConfigBuilder{}
                 .with_consumers(1)
                 .with_parsers(1)
                 .with_chunk_size(1000)
                 .within_set_parallelism(false)
                 .build();
  fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(cfg, read_file);
  parser.start();

  CLI::AutoTimer timer{"parsing kmers", CLI::Timer::Big};
   
  pufferfish::CanonicalKmerIterator kit_end;
  auto rg = parser.getReadGroup();
  while (parser.refill(rg)) {
    // Here, rg will contain a chunk of read pairs
    // we can process.

    for (auto& rp : rg) {
      // kmer_pos = 0;
      auto& r1 = rp.first().seq;
      pufferfish::CanonicalKmerIterator kit1(r1);
      for (; kit1 != kit_end; ++kit1) {
        kmers.push_back(kit1->first);
      }
    }
  }

  return kmers;
}


template <typename IndexT>
int doPufferfishTestLookup(IndexT& pi, pufferfish::ValidateOptions& validateOpts) {
  CanonicalKmer::k(pi.k());
  uint32_t k = pi.k();
  size_t found = 0;
  size_t notFound = 0;
  size_t totalHits = 0;
  size_t totalKmers = 0;
  size_t posErrors = 0;
  size_t contigPosErrors = 0;
  size_t contigLenErrors = 0;
  size_t refPosErrors = 0;
  size_t orientationErrors = 0;

  std::vector<std::string> read_file = {validateOpts.refFile};
  auto cfg = fastx_parser::ParserConfigBuilder{}
                 .with_consumers(1)
                 .with_parsers(1)
                 .with_chunk_size(1000)
                 .within_set_parallelism(false)
                 .build();
  fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(cfg, read_file);
  parser.start();
  pufferfish::util::QueryCache qc;
  pufferfish::CanonicalKmerIterator kit_end;

  // Conditionally create SSHash streaming query
  std::unique_ptr<PufferfishIndex::sshash_streaming_query_t> sq_ptr;
  if constexpr (std::is_same_v<IndexT, PufferfishIndex>) {
    if (pi.isSSHashBacked()) {
      sq_ptr = std::make_unique<PufferfishIndex::sshash_streaming_query_t>(
          pi.createStreamingQuery());
    }
  }

  auto start = std::chrono::high_resolution_clock::now();
  auto rg = parser.getReadGroup();
  while (parser.refill(rg)) {
    for (auto& rp : rg) {
      auto& r1 = rp.first().seq;
      pufferfish::CanonicalKmerIterator kit1(r1);
      int32_t prevPos = -2;

      if constexpr (std::is_same_v<IndexT, PufferfishIndex>) {
        if (sq_ptr) { sq_ptr->reset(); }
      }

      for (; kit1 != kit_end; ++kit1) {
        pufferfish::util::ProjectedHits phits;
        ++totalKmers;

        if constexpr (std::is_same_v<IndexT, PufferfishIndex>) {
          if (sq_ptr) {
            int32_t curPos = kit1->second;
            if (curPos != prevPos + 1) { sq_ptr->reset(); }
            prevPos = curPos;
            phits = pi.getRefPos(kit1->first, qc, *sq_ptr, r1.c_str() + curPos);
          } else {
            phits = pi.getRefPos(kit1->first, qc);
          }
        } else {
          phits = pi.getRefPos(kit1->first, qc);
        }

        if (phits.empty()) {
          ++notFound;
        } else if (phits.refRange.size() <= 200) {
          ++found;
          totalHits += phits.refRange.size();

          // --- Positional correctness checks ---
          auto& mer = kit1->first;

          // 1. Verify k-mer at globalPos matches the query k-mer
          uint64_t fk = pi.getSeqBits(2 * phits.globalPos_, 2 * k);
          auto keq = mer.isEquivalent(fk);
          if (keq == KmerMatchType::NO_MATCH) {
            if (posErrors < 10) {
              std::cerr << "POS ERROR: kmer=" << mer.to_str()
                        << " globalPos=" << phits.globalPos_
                        << " contigIdx=" << phits.contigIdx_
                        << " fk=0x" << std::hex << fk << std::dec << "\n";
            }
            ++posErrors;
          } else {
            // 2. Verify orientation consistency
            bool expectedFW = (keq == KmerMatchType::IDENTITY_MATCH);
            if (expectedFW != phits.contigOrientation_) {
              ++orientationErrors;
            }
          }

          // 3. Verify contigPos < contigLen
          if (phits.contigPos_ >= phits.contigLen_) {
            if (contigPosErrors < 10) {
              std::cerr << "CONTIG_POS ERROR: contigPos=" << phits.contigPos_
                        << " >= contigLen=" << phits.contigLen_
                        << " contigIdx=" << phits.contigIdx_ << "\n";
            }
            ++contigPosErrors;
          }

          // 4. Verify contigLen >= k
          if (phits.contigLen_ < k) {
            if (contigLenErrors < 10) {
              std::cerr << "CONTIG_LEN ERROR: contigLen=" << phits.contigLen_
                        << " < k=" << k << " contigIdx=" << phits.contigIdx_ << "\n";
            }
            ++contigLenErrors;
          }

          // 5. For each reference hit, verify decoded position is within reference bounds
          for (auto& rpos : phits.refRange) {
            auto decoded = phits.decodeHit(rpos);
            uint32_t tid = rpos.transcript_id();
            uint32_t refLen = pi.refLength(tid);
            if (decoded.pos + k > refLen) {
              if (refPosErrors < 10) {
                std::cerr << "REF_POS ERROR: tid=" << tid
                          << " decoded.pos=" << decoded.pos
                          << " + k=" << k << " > refLen=" << refLen << "\n";
              }
              ++refPosErrors;
            }
          }
        }
      }
    }
  }
  parser.stop();

  auto finish = std::chrono::high_resolution_clock::now();
  auto time_in_nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count();
  double ns_per_kmer = time_in_nanoseconds / (static_cast<double>(totalKmers));

  std::cerr << "found = " << found << ", not found = " << notFound << "\n";
  std::cerr << "total hits = " << totalHits << "\n";
  std::cerr << "average query time is " << ns_per_kmer << " ns / query\n";

  // Report validation results
  std::cerr << "\n=== Positional Validation ===\n";
  std::cerr << "globalPos k-mer mismatches: " << posErrors << " / " << found << "\n";
  std::cerr << "orientation mismatches:     " << orientationErrors << " / " << found << "\n";
  std::cerr << "contigPos >= contigLen:     " << contigPosErrors << " / " << found << "\n";
  std::cerr << "contigLen < k:              " << contigLenErrors << " / " << found << "\n";
  std::cerr << "ref pos out of bounds:      " << refPosErrors << " / " << totalHits << "\n";

  if (posErrors == 0 && orientationErrors == 0 && contigPosErrors == 0 &&
      contigLenErrors == 0 && refPosErrors == 0) {
    std::cerr << "ALL POSITIONAL CHECKS PASSED\n";
  } else {
    std::cerr << "VALIDATION ERRORS DETECTED\n";
  }

  if constexpr (std::is_same_v<IndexT, PufferfishIndex>) {
    if (sq_ptr) {
      std::cerr << "\nstreaming query stats: "
                << sq_ptr->num_searches() << " full searches, "
                << sq_ptr->num_extensions() << " extensions ("
                << (100.0 * sq_ptr->num_extensions() /
                    (sq_ptr->num_searches() + sq_ptr->num_extensions() + sq_ptr->num_negative_lookups() + sq_ptr->num_invalid_lookups()))
                << "% extended)\n";
    }
  }

  return 0;
}




template <typename IndexT>
int doPufferfishTestLookupStreamingParse(IndexT& pi, pufferfish::ValidateOptions& validateOpts) {
  CanonicalKmer::k(pi.k());
  int k = pi.k();
  (void)k;
  size_t found = 0;
  size_t notFound = 0;
  size_t totalHits = 0;
  {
    CLI::AutoTimer timer{"searching kmers", CLI::Timer::Big};

    std::vector<std::string> read_file = {validateOpts.refFile};
    auto cfg = fastx_parser::ParserConfigBuilder{}
                   .with_consumers(1)
                   .with_parsers(1)
                   .with_chunk_size(1000)
                   .within_set_parallelism(false)
                   .build();
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(cfg, read_file);
    parser.start();
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    size_t rn{0};
    pufferfish::util::QueryCache qc;


    pufferfish::CanonicalKmerIterator kit_end;
    auto rg = parser.getReadGroup();
    while (parser.refill(rg)) {
      // Here, rg will contain a chunk of read pairs
      // we can process.
      for (auto& rp : rg) {
        // kmer_pos = 0;
        if (rn % 500000 == 0) {
          std::cerr << "rn : " << rn << "\n";
          std::cerr << "found = " << found << ", notFound = " << notFound
                    << ", total hits = " << totalHits << "\n";
        }
        ++rn;
        auto& r1 = rp.first().seq;
        /*
        CanonicalKmer mer;
        bool valid = true;
        int lastinvalid = -1;
          for (size_t i = 0; i < r1.length(); ++i) {
          int c = kmers::codeForChar(r1[i]);
          if (c != -1) {
            mer.shiftFw(r1[i]);
            valid = (i - lastinvalid >= k);
          } else {
            lastinvalid = i;
            valid = false;
          }
          if (i >= k - 1 and valid) {
            auto phits = pi.getRefPos(mer);
            if (phits.empty()) {
            ++notFound;
            } else {
            ++found;
            bool cor = false;
            uint32_t clen = 0;
            std::vector<uint32_t> wrongPos;
            for (auto& rpos : phits.refRange) {
            ++totalHits;
            }
            }
           }

           }
          */

        pufferfish::CanonicalKmerIterator kit1(r1);
        for (; kit1 != kit_end; ++kit1) {
          auto phits = pi.getRefPos(kit1->first, qc);
          if (phits.empty()) {
            ++notFound;
          } else {
            ++found;
            totalHits += phits.refRange.size();
          }
        }
      }
    }
  }

  std::cerr << "found = " << found << ", not found = " << notFound << "\n";
  std::cerr << "total hits = " << totalHits << "\n";
  return 0;
}

int pufferfishTestLookup(pufferfish::ValidateOptions& validateOpts) {
  auto indexDir = validateOpts.indexDir;
  std::string indexType;
  {
    std::ifstream infoStream(indexDir + "/info.json");
    cereal::JSONInputArchive infoArchive(infoStream);
    infoArchive(cereal::make_nvp("sampling_type", indexType));
    std::cerr << "Index type = " << indexType << '\n';
    infoStream.close();
  }

  PufferfishIndex pi(validateOpts.indexDir);
  return doPufferfishTestLookup(pi, validateOpts);
}
