//
// Created by Fatemeh Almodaresi on 8/28/18.
//

#ifndef PUFFERFISH_CHAINFINDER_H
#define PUFFERFISH_CHAINFINDER_H

#include "PufferfishIndex.hpp"
#include "CommonTypes.hpp"
#include "Util.hpp"

#include <sparsepp/spp.h>
#include "ankerl/unordered_dense.h"
#include "metro/metrohash64.h"

class MemClusterer {
private:
  uint32_t maxAllowedRefsPerHit_ = 1000;
  double consensusFraction_ = 0.65;
  double pre_merge_chain_sub_thresh_{0.9};
  double inv_pre_merge_chain_sub_thresh_{1.0/0.9};

  // Encode (tid, isFW) as a single uint64_t: (tid << 1) | isFW
  // This avoids pair construction/hashing overhead in the hot loop.
  using RefMemMap = pufferfish::util::CachedVectorMap<uint64_t, std::vector<pufferfish::util::MemInfo>, std::hash<uint64_t>>;
  static inline uint64_t encodeRefKey(uint64_t tid, bool isFW) { return (tid << 1) | static_cast<uint64_t>(isFW); }
  static inline uint64_t decodeRefTid(uint64_t key) { return key >> 1; }
  static inline bool decodeRefFW(uint64_t key) { return key & 1; }

public:
  void set_chain_sub_opt_thresh_(double pre_merge_chain_sub_thresh, double inv_pre_merge_chain_sub_thresh);
  double chainSubOptThresh() const;
  void setConsensusFraction(double cf);
  double getConsensusFraction() const;
  void setMaxAllowedRefsPerHit(uint32_t max);
  uint32_t getMaxAllowedRefsPerHit();
  void setHitFilterPolicy(pufferfish::util::HitFilterPolicy hfp);
  pufferfish::util::HitFilterPolicy getHitFilterPolicy() const;

  size_t fillMemCollection(std::vector<std::pair<int, pufferfish::util::ProjectedHits>> &hits,
                         //pufferfish::common_types::RefMemMapT& trMemMap,
                         RefMemMap& trMemMap,
                           std::vector<pufferfish::util::UniMemInfo> &memCollection, uint64_t firstDecoyIndex,
                           ankerl::unordered_dense::map<pufferfish::common_types::ReferenceID, bool>& other_end_refs);

  bool findOptChain(std::vector<std::pair<int, pufferfish::util::ProjectedHits>> &hits,
                    pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& memClusters,
                    //ankerl::unordered_dense::map<pufferfish::common_types::ReferenceID, std::vector<pufferfish::util::MemCluster>> &memClusters,
                      uint32_t maxSpliceGap, std::vector<pufferfish::util::UniMemInfo> &memCollection, uint32_t readLen,
                    ankerl::unordered_dense::map<pufferfish::common_types::ReferenceID, bool>& other_end_refs,
                    bool hChain,
                    RefMemMap& trMemMap,
                    uint64_t firstDecoyIndex,
                    //pufferfish::common_types::RefMemMapT& trMemMap,
                    bool verbose = false);

private:
  itlib::small_vector<double> f;
  itlib::small_vector<int32_t> p;
  itlib::small_vector<uint8_t> keepMem;
  itlib::small_vector<uint64_t> memIndicesInReverse;
  itlib::small_vector<int32_t> bestChainEndList;
  itlib::small_vector<int32_t> chainQuerySig;
  pufferfish::util::HitFilterPolicy hitFilterPolicy_{pufferfish::util::HitFilterPolicy::FILTER_AFTER_CHAINING};
};

#endif //PUFFERFISH_CHAINFINDER_H
