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

#define JUMPSIZE 10
template <typename PufferfishIndexT> class MemCollector {

public:
  MemCollector(PufferfishIndexT* pfi) : pfi_(pfi) { k = pfi_->k(); }

  bool clusterMems(std::vector<std::pair<int, util::ProjectedHits>>& hits,
                   spp::sparse_hash_map<pufferfish::common_types::ReferenceID, std::vector<util::MemCluster>>& memClusters,
                   uint32_t maxSpliceGap, std::vector<util::UniMemInfo>& memCollection, bool verbose = false) {
    using namespace pufferfish::common_types;
    //(void)verbose;

    if (hits.empty()) {
      return false;
    }

    // spp::sparse_hash_map<std::pair<size_t, bool>, std::vector<util::MemInfo>>
    // trMemMap;

    // Map from (reference id, orientation) pair to a cluster of
    // MEMs.
    std::map<std::pair<ReferenceID, bool>, std::vector<util::MemInfo>>
        trMemMap;
    memCollection.reserve(hits.size());
    for (auto& hit : core::range<decltype(hits.begin())>(hits.begin(), hits.end())) {
      auto& readPos = hit.first;
      auto& projHits = hit.second;
      // NOTE: here we rely on internal members of the ProjectedHit (i.e., member variables ending in "_").
      // Maybe we want to change the interface (make these members public or provide accessors)?
      memCollection.emplace_back(projHits.contigIdx_, projHits.contigOrientation_, readPos, projHits.k_);
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
          }
          else break;
        }
        if (!foundAtLeastOneCluster) {
          auto& lastClus = currMemClusters.back();
          util::MemCluster newClus(isFw);
          if (!currMemClusters.empty() && gapIsSmall) {
          // add all previous compatable mems before this last one that was crossed
            for (auto mem = lastClus.mems.begin(); mem != lastClus.mems.end() && mem->memInfo->rpos < hit.memInfo->rpos; mem++) {
              newClus.mems.emplace_back(mem->memInfo, mem->tpos);
            }
          }
          newClus.mems.emplace_back(hit.memInfo, hit.tpos);
          currMemClusters.push_back(newClus);
        }
      }
      if (verbose) {
        std::cout << "t" << tid << " , isFw:" << isFw << " cluster size:" << currMemClusters.size() << "\n";
        for (auto& clus : currMemClusters) {
          for (auto& mem : clus.mems) {
            std::cout << "\t t" << mem.tpos << " r" << mem.memInfo->rpos << " len" << mem.memInfo->memlen << "\n";
          }
          std::cout << "\n";
        }
      }
      // This is kind of inefficient (copying the currMemClusters while probably we can build it on the fly
      memClusters[tid].insert(memClusters[tid].end(), currMemClusters.begin(), currMemClusters.end());
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
    if (!hit.contigOrientation_) {
      hit.contigPos_ -= (hit.k_ - k);
      hit.globalPos_ -= (hit.k_ - k);
    }
    return readStart;
  }

  bool operator()(std::string& read,
                  spp::sparse_hash_map<size_t, std::vector<util::MemCluster>>& memClusters,
                  uint32_t maxSpliceGap,
                  util::MateStatus mateStatus
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
