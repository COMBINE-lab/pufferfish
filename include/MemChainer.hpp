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


class MemClusterer {
private:
    uint64_t maxAllowedRefsPerHit = 1000;

public:

  bool fillMemCollection(std::vector<std::pair<int, util::ProjectedHits>> &hits,
                           std::map<std::pair<pufferfish::common_types::ReferenceID, bool>, std::vector<util::MemInfo>> &trMemMap,
                           std::vector<util::UniMemInfo> &memCollection, util::ReadEnd re, bool verbose = false);

  bool findOptChain(std::vector<std::pair<int, util::ProjectedHits>> &hits,
                      spp::sparse_hash_map<pufferfish::common_types::ReferenceID, std::vector<util::MemCluster>> &memClusters,
                      uint32_t maxSpliceGap, std::vector<util::UniMemInfo> &memCollection, uint32_t readLen, bool hChain, bool verbose = false);
};

#endif //PUFFERFISH_CHAINFINDER_H
