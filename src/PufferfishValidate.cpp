#include "FastxParser.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <iterator>
#include <type_traits>

#include "spdlog/spdlog.h"
#include "jellyfish/mer_dna.hpp"
#include "PufferFS.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "CanonicalKmer.hpp"

#include "Util.hpp"
#include "PufferfishIndex.hpp"

int pufferfishValidate(util::ValidateOptions& validateOpts) {
  PufferfishIndex pi(validateOpts.indexDir);
  CanonicalKmer::k(pi.k());
  size_t k = pi.k();
  size_t found = 0;
  size_t notFound = 0;
  size_t correctPosCntr = 0;
  size_t incorrectPosCntr = 0;
  size_t numTrueTxp = 0;
  size_t numLostTxp = 0;
  {
    ScopedTimer st;
    std::vector<std::string> read_file = {validateOpts.refFile};
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(read_file, 1, 1);
    parser.start();
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    size_t rn{0};
    //size_t kmer_pos{0};
    auto rg = parser.getReadGroup();
    while (parser.refill(rg)) {
      // Here, rg will contain a chunk of read pairs
      // we can process.
      for (auto& rp : rg) {
        //kmer_pos = 0;
        if (rn % 500000 == 0) {
          std::cerr << "rn : " << rn << "\n";
          std::cerr << "found = " << found << ", notFound = " << notFound
                    << "\n";
        }
        ++rn;
        auto& r1 = rp.seq;
        CanonicalKmer mer;
        bool merOK = mer.fromStr(r1); // mer.from_chars(r1.begin());
        if (!merOK) {
          std::cerr << "contig too short!";
          std::exit(1);
        }

        auto phits = pi.getRefPos(mer);
        if (phits.empty()) {
          ++notFound;
        } else {
          ++found;
          bool correctPos = false;
          bool foundTxp = false;
          for (auto& rpos : phits.refRange) {
            if (pi.refName(rpos.transcript_id()) == rp.name) {
              foundTxp = true;
              auto refInfo = phits.decodeHit(rpos);
              if (refInfo.pos == 0) {
                correctPos = true;
              }
            } else {
            }
          }
          if (foundTxp) {
            ++numTrueTxp;
            if (correctPos) { ++correctPosCntr; } else { ++incorrectPosCntr; }
          } else {
            ++numLostTxp;
          }
        }
        
        for (size_t i = k; i < r1.length(); ++i) {
          mer.shiftFw(r1[i]);
          auto phits = pi.getRefPos(mer);
          if (phits.empty()) {
            ++notFound;
          } else {
            ++found;
            bool correctPos = false;
            bool foundTxp = false;
            bool cor = false;
            uint32_t clen = 0;
            std::vector<uint32_t> wrongPos;
            for (auto& rpos : phits.refRange) {
              if (pi.refName(rpos.transcript_id()) == rp.name) {
                foundTxp = true;
                auto refInfo = phits.decodeHit(rpos);
                if (refInfo.pos == i - k + 1) {
                correctPos = true;
                } else {
                  cor = phits.contigOrientation_;
                  clen = phits.contigLen_;
                  wrongPos.push_back(refInfo.pos);
                }
              } else {
                }
                }
                if (foundTxp) {
                ++numTrueTxp;
                if (correctPos) {
                  ++correctPosCntr;
                } else {
                  std::cerr << "correct pos = " << i - k + 1 << ", found " << util::str(wrongPos)
                            << ", contig orientation = " << cor << ", contig len = " << clen
                            << ", contig rank = " << pi.contigID(mer) << '\n';
                  ++incorrectPosCntr;
                }
            } else {
              ++numLostTxp;
            }
        }
      }
    }
  }
  }
  std::cerr << "found = " << found << ", not found = " << notFound << "\n";
  std::cerr << "correctPos = " << correctPosCntr
            << ", incorrectPos = " << incorrectPosCntr << "\n";
  std::cerr << "corrextTxp = " << numTrueTxp << ", lostTxp = " << numLostTxp << "\n";




  return 0;
}
