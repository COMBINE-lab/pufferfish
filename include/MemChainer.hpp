//
// Created by Fatemeh Almodaresi on 8/28/18.
//

#ifndef PUFFERFISH_CHAINFINDER_H
#define PUFFERFISH_CHAINFINDER_H

#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "PufferfishLossyIndex.hpp"
#include "CommonTypes.hpp"

#include <sparsepp/spp.h>
#include "parallel_hashmap/phmap.h"

class MemClusterer {
private:
    uint32_t maxAllowedRefsPerHit = 1000;

public:

  void setMaxAllowedRefsPerHit(uint32_t max);
  uint32_t getMaxAllowedRefsPerHit();


  bool fillMemCollection(std::vector<std::pair<int, pufferfish::util::ProjectedHits>> &hits,
                         pufferfish::common_types::RefMemMapT& trMemMap,
                           std::vector<pufferfish::util::UniMemInfo> &memCollection, pufferfish::util::ReadEnd re,
                           phmap::flat_hash_map<pufferfish::common_types::ReferenceID, bool>& other_end_refs,
                           bool verbose = false);

  bool findOptChain(std::vector<std::pair<int, pufferfish::util::ProjectedHits>> &hits,
                      phmap::flat_hash_map<pufferfish::common_types::ReferenceID, std::vector<pufferfish::util::MemCluster>> &memClusters,
                      uint32_t maxSpliceGap, std::vector<pufferfish::util::UniMemInfo> &memCollection, uint32_t readLen,
                    phmap::flat_hash_map<pufferfish::common_types::ReferenceID, bool>& other_end_refs,
                    bool hChain,
                    pufferfish::common_types::RefMemMapT& trMemMap,
                    bool verbose = false);
};

#endif //PUFFERFISH_CHAINFINDER_H
