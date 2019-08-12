#include "FastxParser.hpp"
#include <cmath>
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
#include "PufferfishSparseIndex.hpp"
#include "Util.hpp"

namespace kmers = combinelib::kmers;

template <typename IndexT>
int doPufferfishTestLookup(IndexT& pi, pufferfish::ValidateOptions& validateOpts) {
  CanonicalKmer::k(pi.k());
  int k = pi.k();
  (void)k;
  size_t found = 0;
  size_t notFound = 0;
  size_t totalHits = 0;
  {
    CLI::AutoTimer timer{"searching kmers", CLI::Timer::Big};
    std::vector<std::string> read_file = {validateOpts.refFile};
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(read_file, 1, 1);
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
        auto& r1 = rp.seq;
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

  if (indexType == "sparse") { 
    PufferfishSparseIndex pi(validateOpts.indexDir);
    return doPufferfishTestLookup(pi, validateOpts);
  } else if (indexType == "dense") {
    PufferfishIndex pi(validateOpts.indexDir);
    return doPufferfishTestLookup(pi, validateOpts);
  }
  return 0;
}

