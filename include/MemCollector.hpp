#ifndef HIT_COLLECTOR_HPP
#define HIT_COLLECTOR_HPP

#include "CommonTypes.hpp"
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

template <typename PufferfishIndexT> class MemCollector {

public:
  MemCollector(PufferfishIndexT* pfi) : pfi_(pfi) { k = pfi_->k(); }

  bool clusterMems(std::vector<std::pair<int, util::ProjectedHits>>& hits,
                   spp::sparse_hash_map<pufferfish::common_types::ReferenceID, std::vector<util::MemCluster>>& memClusters,
                   uint32_t maxSpliceGap, std::vector<util::UniMemInfo>& memCollection, bool verbose = false) {
    using namespace pufferfish::common_types;
    (void)verbose;

    if (hits.empty()) {
      return false;
    }

    // Map from (reference id, orientation) pair to a cluster of MEMs.
    std::map<std::pair<ReferenceID, bool>, std::vector<util::MemInfo>>
        trMemMap;
    memCollection.reserve(hits.size());
    for (auto& hit : core::range<decltype(hits.begin())>(hits.begin(), hits.end())) {
      auto& readPos = hit.first;
      auto& projHits = hit.second;
      // NOTE: here we rely on internal members of the ProjectedHit (i.e., member variables ending in "_").
      // Maybe we want to change the interface (make these members public or provide accessors)?
      memCollection.emplace_back(projHits.contigIdx_, projHits.contigOrientation_, readPos, projHits.k_, projHits.contigPos_);
      auto memItr = std::prev(memCollection.end());
      for (auto& posIt : projHits.refRange) {
        auto refPosOri = projHits.decodeHit(posIt);
        trMemMap[std::make_pair(posIt.transcript_id(), refPosOri.isFW)]
          .emplace_back(memItr, refPosOri.pos);
      }
    }

    for (auto& trMem : core::range<decltype(trMemMap.begin())>(trMemMap.begin(), trMemMap.end())) {
      auto& trOri = trMem.first;
      auto& tid = trOri.first;
      auto& isFw = trOri.second;
      auto& memList = trMem.second;
      // sort memList according to mem reference positions
      std::sort(memList.begin(), memList.end(),
                [](util::MemInfo& q1, util::MemInfo& q2) -> bool {
                  return q1.tpos < q2.tpos; // sort based on tpos
                  });
      std::vector<util::MemCluster> currMemClusters;
      // cluster MEMs so that all the MEMs in one cluster are concordant.
      for (auto& hit : core::range<decltype(memList.begin())>(memList.begin(), memList.end())) {
        bool foundAtLeastOneCluster = false;
        bool gapIsSmall = false;
        for (auto prevClus = currMemClusters.rbegin();
             prevClus != currMemClusters.rend(); prevClus++) {
          if (hit.tpos - prevClus->getTrLastHitPos() < maxSpliceGap) { // if the distance between last mem and the new one is NOT longer than maxSpliceGap
            gapIsSmall = true;
            if ( (isFw && hit.memInfo->rpos >= prevClus->getReadLastHitPos()) ||
                 (!isFw && hit.memInfo->rpos <= prevClus->getReadLastHitPos())) {
              foundAtLeastOneCluster = true;
              prevClus->mems.emplace_back(hit.memInfo, hit.tpos);
            }
          } else {
            break;
          }
        }

        if (!foundAtLeastOneCluster) {
          auto prevLastIndex = static_cast<int32_t>(currMemClusters.size()) - 1;
          // Create the new clusters on the end of the currMemClusters vector
          currMemClusters.emplace_back(isFw);
          auto& newClus = currMemClusters.back();
          if ((prevLastIndex > 0) and gapIsSmall) {
            auto& lastClus = currMemClusters[prevLastIndex];
            // add all previous compatable mems before this last one that was crossed
            for (auto mem = lastClus.mems.begin(); mem != lastClus.mems.end() && mem->memInfo->rpos < hit.memInfo->rpos; mem++) {
              newClus.mems.emplace_back(mem->memInfo, mem->tpos);
            }
          }
          newClus.mems.emplace_back(hit.memInfo, hit.tpos);
        }
      }
      /*
      if (verbose) {
        std::cout << "t" << tid << " , isFw:" << isFw << " cluster size:" << currMemClusters.size() << "\n";
        for (auto& clus : currMemClusters) {
          for (auto& mem : clus.mems) {
            std::cout << "\t t" << mem.tpos << " r" << mem.memInfo->rpos << " len" << mem.memInfo->memlen << "\n";
          }
          std::cout << "\n";
        }
      }
      */
      // This is kind of inefficient (copying the currMemClusters while probably we can build it on the fly
      memClusters[tid].insert(memClusters[tid].end(), std::make_move_iterator(currMemClusters.begin()),
                              std::make_move_iterator(currMemClusters.end()));
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

    int currReadStart = kit->second + 1;
    auto readSeqView = kit.seq();
    auto readSeqLen = readSeqView.size();
    auto readSeqStart = currReadStart;
    auto readSeqOffset = currReadStart + k - 1;
    int fastNextReadCode{0};
    bool stillMatch = true;

    while (stillMatch and
           (cCurrPos < cEndPos) and
           (cCurrPos > cStartPos) and
           readSeqOffset < readSeqLen) { // over kmers in contig

      if (hit.contigOrientation_) { // if fw match, compare read last base with
                                    // contig first base and move fw in the
                                    // contig
        auto baseCnt = k < cEndPos - cCurrPos ? k : cEndPos - cCurrPos;
        uint64_t fk = allContigs.get_int(2 * (cCurrPos), 2 * baseCnt);
        cCurrPos += baseCnt;
        for (size_t i = 0; i < baseCnt && readSeqOffset < readSeqLen; i++) {
          // be dirty and peek into the underlying read
          fastNextReadCode = my_mer::code(readSeqView[readSeqOffset]);
          int contigCode = (fk >> (2 * i)) & 0x3;
          if (fastNextReadCode != contigCode) {
            stillMatch = false;
            break;
          }
          hit.k_++;
          readSeqOffset++;
          readSeqStart++;
        }
      } else { // if rc match, compare read last base with contig last base and
               // move backward in the contig
        auto baseCnt = k < cCurrPos - cStartPos ? k : cCurrPos - cStartPos;
        uint64_t fk = allContigs.get_int(2 * (cCurrPos - baseCnt), 2 * baseCnt);
        cCurrPos -= baseCnt;
        for (size_t i = baseCnt - 1; i >= 0 && readSeqOffset < readSeqLen; i--) {
          // be dirty and peek into the underlying read
          fastNextReadCode = my_mer::code(my_mer::complement(readSeqView[readSeqOffset]));
          int contigCode = (fk >> (2 * i)) & 0x3;
          if (fastNextReadCode != contigCode) {
            stillMatch = false;
            break;
          }
          hit.k_++;
          readSeqOffset++;
          readSeqStart++;
        }
      }
    }
    if (!hit.contigOrientation_) {
      hit.contigPos_ -= (hit.k_ - k);
      hit.globalPos_ -= (hit.k_ - k);
    }
    kit.jumpTo(readSeqStart);
    return currReadStart;
  }

  bool operator()(std::string& read,
                  spp::sparse_hash_map<size_t, std::vector<util::MemCluster>>& memClusters,
                  uint32_t maxSpliceGap,
                  util::MateStatus mateStatus) {
    // currently unused:
    // uint32_t readLen = static_cast<uint32_t>(read.length()) ;

    util::ProjectedHits phits;
    std::vector<std::pair<int, util::ProjectedHits>> rawHits;

    CanonicalKmer::k(k);
    pufferfish::CanonicalKmerIterator kit_end;
    pufferfish::CanonicalKmerIterator kit1(read);
    util::QueryCache qc;

    while (kit1 != kit_end) {
      auto phits = pfi_->getRefPos(kit1->first, qc);
      if (!phits.empty()) {
        // kit1 gets updated inside expandHitEfficient function
        size_t readPos = expandHitEfficient(phits, kit1);
        rawHits.push_back(std::make_pair(readPos, phits));
      } else
        ++kit1;
    }

    // if this is the right end of a paired-end read, use memCollectionRight,
    // otherwise (left end or single end) use memCollectionLeft.
    auto* memCollection = (mateStatus == util::MateStatus::PAIRED_END_RIGHT) ?
      &memCollectionRight : &memCollectionLeft;
    if (rawHits.size() > 0) {
      clusterMems(rawHits, memClusters, maxSpliceGap, *memCollection);
      return true;
    }
    return false;
  }

  void clear() {
    memCollectionLeft.clear();
    memCollectionRight.clear();
  }

private:
  PufferfishIndexT* pfi_;
  size_t k;
  AlignerEngine ae_;
  std::vector<util::UniMemInfo> memCollectionLeft;
  std::vector<util::UniMemInfo> memCollectionRight;
};
#endif
