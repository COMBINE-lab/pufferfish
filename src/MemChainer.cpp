#include "MemChainer.hpp"

static inline float fastlog2(float x) {
  union {
      float f;
      uint32_t i;
  } vx = {x};
  union {
      uint32_t i;
      float f;
  } mx = {(vx.i & 0x007FFFFF) | 0x3f000000};
  float y = vx.i;
  y *= 1.1920928955078125e-7f;

  return y - 124.22551499f
         - 1.498030302f * mx.f
         - 1.72587999f / (0.3520887068f + mx.f);
}

static inline float fasterlog2(float x) {
  union {
      float f;
      uint32_t i;
  } vx = {x};
  float y = vx.i;
  y *= 1.1920928955078125e-7f;
  return y - 126.94269504f;
}

void MemClusterer::setMaxAllowedRefsPerHit(uint32_t max){
  maxAllowedRefsPerHit = max;
}

uint32_t MemClusterer::getMaxAllowedRefsPerHit() {
  return maxAllowedRefsPerHit;
}

bool MemClusterer::fillMemCollection(std::vector<std::pair<int, util::ProjectedHits>> &hits,
                                     std::map<std::pair<pufferfish::common_types::ReferenceID, bool>, std::vector<util::MemInfo>> &trMemMap,
                                     std::vector<util::UniMemInfo> &memCollection, util::ReadEnd re,
                                     spp::sparse_hash_map<pufferfish::common_types::ReferenceID, bool> & other_end_refs, bool verbose) {
  if (verbose)
    std::cerr << "\n[FIND_OPT_CHAIN]\n";

  using namespace pufferfish::common_types;
  //(void)verbose;

  if (hits.empty()) {
    return false;
  }

  // here we guarantee that even if later we fill up
  // every gap between the hits and before the first and after the last hit
  // still the memCollection size doesn't change and hence the pointers are valid
  memCollection.reserve(maxAllowedRefsPerHit * 2 * hits.size() + 1);
  if (verbose)
    std::cerr << "\nreserved memCollection size: " << maxAllowedRefsPerHit * 2 * hits.size() + 1 << "\n";
  for (auto &hit : core::range<decltype(hits.begin())>(hits.begin(), hits.end())) {
    auto &readPos = hit.first;
    auto &projHits = hit.second;
    // NOTE: here we rely on internal members of the ProjectedHit (i.e., member variables ending in "_").
    // Maybe we want to change the interface (make these members public or provide accessors)?
    auto &refs = projHits.refRange;
    if (verbose)
      std::cerr << "total number of references found: " << refs.size() << "\n";
    uint32_t mappings{0};
    if (static_cast<uint64_t>(refs.size()) < maxAllowedRefsPerHit) {
      memCollection.emplace_back(projHits.contigIdx_, projHits.contigOrientation_,
                                 readPos, projHits.k_, projHits.contigPos_,
                                 projHits.globalPos_ - projHits.contigPos_, projHits.contigLen_, re);
      auto memItr = std::prev(memCollection.end());
      for (auto &posIt : refs) {
      //If we want to let the the hits to the references also found by the other end to be accepted
      //if (static_cast<uint64_t>(refs.size()) < maxAllowedRefsPerHit or other_end_refs.find(posIt.transcript_id()) != other_end_refs.end() ) {
        auto refPosOri = projHits.decodeHit(posIt);
        trMemMap[std::make_pair(posIt.transcript_id(), refPosOri.isFW)]
                .emplace_back(memItr, refPosOri.pos, refPosOri.isFW);
        mappings++;
      //}
      }

      if (verbose)
        std::cerr << "total number of mappings found: " << mappings << "\n";
    }
  }
  return true;
}

bool MemClusterer::findOptChain(std::vector<std::pair<int, util::ProjectedHits>> &hits,
                                spp::sparse_hash_map<pufferfish::common_types::ReferenceID, std::vector<util::MemCluster>> &memClusters,
                                uint32_t maxSpliceGap, std::vector<util::UniMemInfo> &memCollection, uint32_t readLen,
                                spp::sparse_hash_map<pufferfish::common_types::ReferenceID, bool>& other_end_refs,
                                bool hChain, bool verbose) {
  using namespace pufferfish::common_types;
  //(void)verbose;

  // Map from (reference id, orientation) pair to a cluster of MEMs.
  std::map<std::pair<ReferenceID, bool>, std::vector<util::MemInfo>>
          trMemMap;
  if (!fillMemCollection(hits, trMemMap, memCollection, util::ReadEnd::LEFT, other_end_refs, verbose))
    return false;

  std::vector<double> f;
  std::vector<int32_t> p;
  for (auto &trMem : core::range<decltype(trMemMap.begin())>(trMemMap.begin(), trMemMap.end())) {
    auto &trOri = trMem.first;
    auto &tid = trOri.first;
    auto &isFw = trOri.second;
    auto &memList = trMem.second;
    // sort memList according to mem reference positions
    std::sort(memList.begin(), memList.end(),
              [isFw](util::MemInfo &q1, util::MemInfo &q2) -> bool {
                  auto q1ref = q1.tpos + q1.memInfo->memlen;
                  auto q2ref = q2.tpos + q2.memInfo->memlen;
                  auto q1read = q1.memInfo->rpos + q1.memInfo->memlen;
                  auto q2read = q2.memInfo->rpos + q2.memInfo->memlen;
                  return q1ref != q2ref ? q1ref < q2ref :
                         (isFw ? q1read < q2read : q1read > q2read);// sort based on tpos
              });
    /*if (verbose) {
      std::cerr << "\ntid" << tid << " , isFw:" << isFw << "\n";
      for (auto &m : memList) {
        std::cerr << "\ttpos:" << m.tpos << " rpos:" << m.memInfo->rpos << " len:" << m.memInfo->memlen
                  << "\n";
      }
    }*/

    //auto minPosIt = memList.begin();
    // find the valid chains
    // Use variant of minimap2 scoring (Li 2018)
    // https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty191/4994778
    auto alpha = [](int32_t qdiff, int32_t rdiff, int32_t ilen) -> double {
        double score = ilen;
        double mindiff = (qdiff < rdiff) ? qdiff : rdiff;
        return (score < mindiff) ? score : mindiff;
    };

    auto beta = [maxSpliceGap](int32_t qdiff, int32_t rdiff, double avgseed) -> double {
        if (qdiff < 0 or ((uint32_t) std::max(qdiff, rdiff) > maxSpliceGap)) {
          return std::numeric_limits<double>::infinity();
        }
        double l = qdiff - rdiff;
        int32_t al = std::abs(l);
        // To penalize cases with organized gaps for reads such as
        // CTCCTCATCCTCCTCATCCTCCTCCTCCTCCTCCTCCTCCGCTGCCGCCGCCGACCGACTGAACCGCACCCGCCGCGCCGCACCGCCTCCAAGTCCCGGC
        // polyester simulated on human transcriptome. 0.01 -> 0.05
        return (l == 0) ? 0.0 : (0.05 * avgseed * al + 0.5 * fastlog2(static_cast<float>(al)));
    };
    double bottomScore = std::numeric_limits<double>::lowest();
    double bestScore = bottomScore;
    int32_t bestChainEnd = -1;
    std::vector<int32_t> bestChainEndList;
    double avgseed = 31.0;
    f.clear();
    p.clear();
    //auto lastHitId = static_cast<int32_t>(memList.size() - 1);
    std::vector<util::MemInfo> newMemList;
    uint32_t prev_qposi = 0;
    uint32_t prev_rposi = 0;
    for (int32_t i = 0; i < static_cast<int32_t>(memList.size()); ++i) {
      auto &hi = memList[i];
      auto qposi = hi.isFw ? hi.memInfo->rpos + hi.extendedlen : readLen - hi.rpos;
      auto rposi = hi.tpos + hi.extendedlen;
      if (i > 0 and std::labs(qposi - prev_qposi) == std::labs(rposi - prev_rposi) and hi.tpos < prev_rposi) {
        auto &lastMem = newMemList.back();
        uint32_t extension = rposi - prev_rposi;
        //lastMem.memInfo->memlen += extension;
        lastMem.extendedlen += extension;
        if (!isFw) {
          lastMem.rpos = hi.memInfo->rpos;
        }
        //lastMem.rpos -= extension;
      } else {
        newMemList.emplace_back(hi.memInfo, hi.tpos, hi.isFw);
      }
      prev_qposi = qposi;
      prev_rposi = rposi;
    }
    auto before = memList.size();
    if (memList.size() != newMemList.size()) {
      trMem.second = newMemList;
      memList = trMem.second;
    }
    auto after = newMemList.size();
    if (verbose and after != before) {
      std::cerr << before << "\t" << after << "\t" << memList.size() << "\n";
    }
    if (verbose) {
      std::cerr << "\ntid" << tid << " , isFw:" << isFw << "\n";
      for (auto &m : memList) {
        std::cerr << "\ttpos:" << m.tpos << " rpos:" << m.rpos << " len:" << m.extendedlen
                  << "\n";
      }
    }

    for (int32_t i = 0; i < static_cast<int32_t>(memList.size()); ++i) {
      auto &hi = memList[i];
      //if (hi.extendedlen != hi.memInfo->memlen)
      //  std::cerr<< hi.extendedlen << "  " <<  hi.memInfo->memlen << "\n";

      auto qposi = hi.rpos + hi.extendedlen;
      auto rposi = hi.tpos + hi.extendedlen;

      double baseScore = static_cast<double>(hi.extendedlen);
      p.push_back(i);
      f.push_back(baseScore);

      // possible predecessors in the chain
      int32_t numRounds{2};
      (void) numRounds;
      for (int32_t j = i - 1; j >= 0; --j) {
        auto &hj = memList[j];

        auto qposj = hj.rpos + hj.extendedlen;
        auto rposj = hj.tpos + hj.extendedlen;

        int32_t qdiff = isFw ? qposi - qposj :
                        (qposj - hj.extendedlen) - (qposi - hi.extendedlen);
        int32_t rdiff = rposi - rposj;

        auto extensionScore = f[j] + alpha(qdiff, rdiff, hi.extendedlen) - beta(qdiff, rdiff, avgseed);
        //To fix cases where there are repetting sequences in the read or reference
        int32_t rdiff_mem = hi.tpos - (hj.tpos + hj.extendedlen);
        int32_t qdiff_mem = isFw ? hi.rpos - (hj.rpos + hj.extendedlen) : hj.rpos - (hi.rpos + hi.extendedlen);
        if (rdiff == 0 or qdiff == 0 or rdiff * qdiff < 0 or rdiff_mem * qdiff_mem < 0 or hi.rpos == hj.rpos or
            hi.tpos == hj.tpos)
          extensionScore = -std::numeric_limits<double>::infinity();
        if (verbose) {
          std::cerr << i << " " << j <<
                    " extendedleni:" << hi.extendedlen << " extendedlenj:" << hj.extendedlen <<
                    " f[i]:" << f[i] << " f[j]:" << f[j] <<
                    " readDiff:" << qdiff << " refDiff:" << rdiff <<
                    " alpha:" << alpha(qdiff, rdiff, hi.extendedlen) <<
                    " beta:" << beta(qdiff, rdiff, avgseed) <<
                    " extensionScore: " << extensionScore << "\n";
        }
        bool extendWithJ = (extensionScore > f[i]);
        p[i] = extendWithJ ? j : p[i];
        f[i] = extendWithJ ? extensionScore : f[i];
        // HEURISTIC : if we connected this match to an earlier one
        // i.e. if we extended the chain.
        // This implements Heng Li's heuristic ---
        // "
        // We note that if anchor i is chained to j, chaining i to a predecessor of j
        // is likely to yield a lower score.
        // "
        // here we take this to the extreme, and stop at the first j to which we chain.
        // we can add a parameter "h" as in the minimap paper.  But here we expect the
        // chains of matches in short reads to be short enough that this may not be worth it.
        if (hChain and p[i] < i) {
          numRounds--;
          if (numRounds <= 0) { break; }
        }
        // If the last two hits are too far from each other, we are sure that 
        // every other hit will be even further since the mems are sorted
        if (rdiff > readLen * 2) { break; }
        // Mohsen: This heuristic hurts the accuracy of the chain in the case of this read:
        // TGAACGCTCTATGATGTCAGCCTACGAGCGCTCTATGATGTTAGCCTACGAGCGCTCTATGATGTCCCCTATGGCTGAGCGCTCTATGATGTCAGCTTAT
        // from Polyester simalted sample aligning to the human transcriptome
      }
      if (f[i] > bestScore) {
        bestScore = f[i];
        bestChainEnd = i;
        bestChainEndList.clear();
        bestChainEndList.push_back(bestChainEnd);
      } else if (f[i] == bestScore) {
        bestChainEndList.push_back(i);
      }
    }
    // Do backtracking
    std::vector<bool> seen(f.size());
    for (uint64_t i = 0; i < seen.size(); i++) seen[i] = false;
    for (auto bestChainEnd : bestChainEndList) {
      bool shouldBeAdded = true;
      if (bestChainEnd >= 0) {
        std::vector<uint64_t> memIndicesInReverse;
        auto lastPtr = p[bestChainEnd];
        while (lastPtr < bestChainEnd) {
          if (seen[bestChainEnd]) {
            shouldBeAdded = false;
            //break;
          }
          memIndicesInReverse.push_back(bestChainEnd);
          seen[bestChainEnd] = true;
          bestChainEnd = lastPtr;
          lastPtr = p[bestChainEnd];
          //lastPtr = bestChainEnd;
          //bestChainEnd = p[bestChainEnd];
        }
        if (seen[bestChainEnd]) {
          shouldBeAdded = false;
        }
        memIndicesInReverse.push_back(bestChainEnd);
        if (shouldBeAdded) {
          memClusters[tid].insert(memClusters[tid].begin(), util::MemCluster(isFw, readLen));
          for (auto it = memIndicesInReverse.rbegin(); it != memIndicesInReverse.rend(); it++) {
            memClusters[tid][0].addMem(memList[*it].memInfo, memList[*it].tpos,
                                       memList[*it].extendedlen, memList[*it].rpos, isFw);
          }
          memClusters[tid][0].coverage = bestScore;
          if (memClusters[tid][0].coverage == readLen)
            memClusters[tid][0].perfectChain = true;
          if (verbose)
            std::cerr<<"Added position: " << memClusters[tid][0].coverage << " " << memClusters[tid][0].mems[0].tpos << "\n";
        }
        //minPosIt += lastPtr;
      } else {
        // should not happen
        std::cerr << "[FATAL] : Cannot find any valid chain for quasi-mapping\n";
        std::cerr << "num hits = " << memList.size() << "\n";
        std::cerr << "bestChainEnd = " << bestChainEnd << "\n";
        std::cerr << "bestChainScore = " << bestScore << "\n";
        std::exit(1);
      }
    }

  }
  if (verbose)
    std::cerr << "\n[END OF FIND_OPT_CHAIN]\n";

  return true;
}
