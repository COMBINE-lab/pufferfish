#ifndef HIT_COLLECTOR_HPP
#define HIT_COLLECTOR_HPP

#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "Util.hpp"
#include "edlib.h"
#include "jellyfish/mer_dna.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <sparsepp/spp.h>
// using spp:sparse_hash_map;

#define JUMPSIZE 10
template <typename PufferfishIndexT> class MemCollector {

public:
  MemCollector(PufferfishIndexT* pfi) : pfi_(pfi) { k = pfi_->k(); }

  bool clusterMems(std::vector<std::pair<int, util::ProjectedHits>>& hits,
                   spp::sparse_hash_map<size_t, util::TrClusters>& memClusters,
                   uint32_t maxSpliceGap, bool verbose = false) {
    (void)verbose;

    if (hits.empty())
      return false;

    // spp::sparse_hash_map<std::pair<size_t, bool>, std::vector<util::MemInfo>>
    // trMemMap;

    // Map from (reference id, orientation) pair to a cluster of
    // MEMs.
    std::map<std::pair<size_t, bool>, std::vector<util::MemInfo>,
             util::cmpByPair>
        trMemMap;
    for (auto hIter = hits.begin(); hIter != hits.end(); ++hIter) {
      auto& readPos = hIter->first;
      auto& projHits = hIter->second;
      for (auto& posIt : projHits.refRange) {
        auto refPosOri = projHits.decodeHit(posIt);
        trMemMap[std::make_pair(posIt.transcript_id(), refPosOri.isFW)]
            .emplace_back(refPosOri.pos, readPos, projHits.k_);
      }
    }

    for (auto trMemIt = trMemMap.begin(); trMemIt != trMemMap.end();
         ++trMemIt) {
      auto& trOri = trMemIt->first;
      auto& tid = trOri.first;
      auto& isFw = trOri.second;
      auto& memList = trMemIt->second;
      // sort memList according to mem read positions

      // pull the orientation to a top-level condition to avoid
      // checking repeatedly during the sort.
      if (isFw) {
        std::sort(memList.begin(), memList.end(),
                  [](util::MemInfo& q1, util::MemInfo& q2) -> bool {
                      return q1.rpos < q2.rpos;
                  });

      } else {
        std::sort(memList.begin(), memList.end(),
                  [](util::MemInfo& q1, util::MemInfo& q2) -> bool {
                    return q1.rpos > q2.rpos;
                  });
      }

      std::vector<util::MemCluster> currMemClusters;
      // cluster mems so that all the mems in one cluster are concordant.
      for (auto hitIt = memList.begin(); hitIt != memList.end(); hitIt++) {
        bool foundAtLeastOneCluster = false;
        for (auto prevClus = currMemClusters.begin();
             prevClus != currMemClusters.end(); prevClus++) {
          auto& mems = prevClus->mems;
          if (hitIt->tpos - mems[mems.size() - 1].tpos < maxSpliceGap) {
            mems.emplace_back(hitIt->tpos, hitIt->rpos, hitIt->memlen);
            foundAtLeastOneCluster = true;
          }
        }
        if (!foundAtLeastOneCluster) {
          currMemClusters.emplace_back(hitIt->tpos, hitIt->rpos, hitIt->memlen);
        }
      }
      memClusters[tid].addBatchCluster(isFw, currMemClusters);
    }
    return true;
  }

  size_t expandHitEfficient(util::ProjectedHits& hit,
                            pufferfish::CanonicalKmerIterator& kit) {
    auto& allContigs = pfi_->getSeq();
    // startPos points to the next kmer in contig (which can be the left or
    // right based on the orientation of match)
    size_t cStartPos =
        hit.globalPos_ - hit.contigPos_ - 1; // next kmer in the read
    size_t cEndPos = cStartPos + hit.contigLen_;
    size_t cCurrPos = hit.globalPos_; // start from next character if fw match
    if (hit.contigOrientation_) { // if match is fw, go to the next k-mer in the
                                  // contig
      cCurrPos += k;
    }
    kit++; // go to the next kmer in the read (In read we always move fw.)
    size_t readStart = kit->second;
    pufferfish::CanonicalKmerIterator kit_end;
    bool stillMatch = true;
    CanonicalKmer cmer;
    while (stillMatch && cCurrPos < cEndPos && cCurrPos > cStartPos &&
           kit != kit_end) {        // over kmers in contig
      if (hit.contigOrientation_) { // if fw match, compare read last base with
                                    // contig first base and move fw in the
                                    // contig
        auto baseCnt = k < cEndPos - cCurrPos ? k : cEndPos - cCurrPos;
        uint64_t fk = allContigs.get_int(2 * (cCurrPos), 2 * baseCnt);
        cmer.fromNum(fk);
        cCurrPos += baseCnt;
        for (size_t i = 0; i < baseCnt && kit != kit_end; i++) {
          auto readCode = (kit->first.fwWord()) >> (2 * (k - 1)) & 0x3;
          auto contigCode = (fk >> (2 * i)) & 0x3;
          if (readCode != contigCode) {
            stillMatch = false;
            break;
          }
          hit.k_++;
          kit++;
        }
      } else { // if rc match, compare read last base with contig last base and
               // move backward in the contig
        auto baseCnt = k < cCurrPos - cStartPos ? k : cCurrPos - cStartPos;
        uint64_t fk = allContigs.get_int(2 * (cCurrPos - baseCnt), 2 * baseCnt);
        cmer.fromNum(fk);
        cCurrPos -= baseCnt;
        for (size_t i = baseCnt - 1; i >= 0 && kit != kit_end; i--) {
          auto readCode = kit->first.rcWord() & 0x3;
          auto contigCode = (fk >> (2 * i)) & 0x3;
          if (readCode != contigCode) {
            stillMatch = false;
            break;
          }
          hit.k_++;
          kit++;
        }
      }
    }
    return readStart;
  }

  bool operator()(std::string& read,
                  spp::sparse_hash_map<size_t, util::TrClusters>& memClusters,
                  uint32_t maxSpliceGap
                  /*,
                  util::MateStatus mateStatus,
                  bool consistentHits,
                  std::vector<std::string>& refBlocks*/) {
    // currently unused:
    // uint32_t readLen = static_cast<uint32_t>(read.length()) ;

    /*if(refBlocks.size() != readLen)
        refBlocks.resize(readLen) ;
    */
    util::ProjectedHits phits;
    std::vector<std::pair<int, util::ProjectedHits>> rawHits;

    CanonicalKmer::k(k);
    pufferfish::CanonicalKmerIterator kit_end;
    pufferfish::CanonicalKmerIterator kit1(read);
    util::QueryCache qc;

    // string block iterator
    // decltype(refBlocks.begin()) bl ;
    // initialize
    // bl = refBlocks.begin() ;
    while (kit1 != kit_end) {
      auto phits = pfi_->getRefPos(kit1->first, qc);
      if (!phits.empty()) {
        // kit1 gets updated inside expandHitEfficient function
        size_t readPos = expandHitEfficient(phits, kit1);
        // std::cerr<<"Final match length: " << phits.k_ << "\n";
        rawHits.push_back(std::make_pair(readPos, phits));
      } else
        ++kit1;
    }

    if (rawHits.size() > 0) {
      clusterMems(rawHits, memClusters, maxSpliceGap);
      return true;
    }
    return false;
  }

private:
  PufferfishIndexT* pfi_;
  size_t k;
  AlignerEngine ae_;
};
#endif
