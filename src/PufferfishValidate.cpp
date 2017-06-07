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
  size_t pos{0};
  size_t k = pi.k();
  size_t found = 0;
  size_t notFound = 0;
  size_t correctPosCntr = 0;
  size_t incorrectPosCntr = 0;
  size_t numTrueTxp = 0;
  {
    ScopedTimer st;
    std::vector<std::string> read_file = {validateOpts.refFile};
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(read_file, 1, 1);
    parser.start();
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    size_t rn{0};
    size_t kmer_pos{0};
    auto rg = parser.getReadGroup();
    while (parser.refill(rg)) {
      // Here, rg will contain a chunk of read pairs
      // we can process.
      for (auto& rp : rg) {
        kmer_pos = 0;
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
        if (phits.refRange.empty()) {
          ++notFound;
        } else {
          ++found;
          bool correctPos = false;
          bool foundTxp = false;
          for (auto& rpos : phits.refRange) {
            if (pi.refName(rpos.transcript_id()) == rp.name) {
              foundTxp = true;
            } else {
            }
          }
          if (foundTxp) {
            ++numTrueTxp;
          }
        }
        /*
        if (pos <= S - k) {
          
          uint64_t fk = seqVec.get_int(2 * pos, 2 * k);
          my_mer fkm;
          fkm.word__(0) = fk;
          if (mer.fwWord() == fkm.word(0) or mer.rcWord() == fkm.word(0)) {
          
            found += 1;
            bool correctPos = false;
            bool foundTxp = false;
            auto rank = realRank(pos);
            for (auto& tr : pi.refList(rank)) { //}[rankOrderedContigIDs[rank]]) {
              if (pi.refName(tr.transcript_id()) == rp.name) {
                //if (trRefIDs[tr.transcript_id()] == rp.name) {
                foundTxp = true;
                uint64_t sp = (uint64_t)realSelect(rank);
                auto relPos = pos - sp;
                auto clen = (uint64_t)realSelect(rank + 1) - sp;
                if (relPos + tr.pos() == kmer_pos or
                    clen - (relPos - tr.pos() - k - 1) == kmer_pos) {
                  correctPos = true;
                  break;
                } else {
                  // std::cerr << "ours : " << relPos + tr.pos << " , true : "
                  // << kmer_pos << "\n";
                }
              }
            }
            if (correctPos) {
              correctPosCntr++;
            } else {
              incorrectPosCntr++;
            }
            if (foundTxp) {
              numTrueTxp++;
            }
            // if (!foundTxp) { std::cerr << "failed to find true txp [" <<
            // rp.name << "]\n"; }
          } else {
            // std::cerr << "rn = " << rn - 1 << ", fk = " <<
            // fkm.get_canonical().to_str() << ", km = " <<
            // mer.get_canonical().to_str() << ", pkmer = " << pkmer <<" \n";
            // std::cerr << "found = " << found << ", not found = " << notFound
            // << "\n";
            notFound += 1;
          }
        } else {
          // std::cerr << "pos = " << pos << ", shouldn't happen\n";
          notFound += 1;
        }
        */

        for (size_t i = k; i < r1.length(); ++i) {
          mer.shiftFw(r1[i]);
          pos = pi.getRawPos(mer);
          //km = mer.getCanonicalWord();
          //res = bphf->lookup(km);
          //pos = (res < N) ? posVec[res] : std::numeric_limits<uint64_t>::max();
          if (pi.isValidPos(pos)){//}pos <= S - k) {
            ++found;
          } else {
            ++notFound;
          }
          /*
            uint64_t fk = seqVec.get_int(2 * pos, 62);
            my_mer fkm;
            fkm.word__(0) = fk; // fkm.canonicalize();
            if (mer.fwWord() == fkm.word(0) or mer.rcWord() == fkm.word(0)) {
              found += 1;

            } else {
              notFound += 1;
            }
          } else {
            // std::cerr << "pos = " << pos << ", shouldn't happen\n";
            notFound += 1;
          }
        }
        */
      }
    }
  }
  }
  std::cerr << "found = " << found << ", not found = " << notFound << "\n";
  std::cerr << "correctPos = " << correctPosCntr
            << ", incorrectPos = " << incorrectPosCntr << "\n";
  std::cerr << "corrextTxp = " << numTrueTxp << "\n";




  return 0;
}
