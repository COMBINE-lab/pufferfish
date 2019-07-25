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
    uint32_t maxAllowedRefsPerHit = 1000;

public:

  void setMaxAllowedRefsPerHit(uint32_t max);
  uint32_t getMaxAllowedRefsPerHit();

  bool fillMemCollection(std::vector<std::pair<int, pufferfish::util::ProjectedHits>> &hits,
                           std::map<std::pair<pufferfish::common_types::ReferenceID, bool>, std::vector<pufferfish::util::MemInfo>> &trMemMap,
                           std::vector<pufferfish::util::UniMemInfo> &memCollection, pufferfish::util::ReadEnd re,
                           spp::sparse_hash_map<pufferfish::common_types::ReferenceID, bool>& other_end_refs,
                           bool verbose = false);

  bool findOptChain(std::vector<std::pair<int, pufferfish::util::ProjectedHits>> &hits,
                      spp::sparse_hash_map<pufferfish::common_types::ReferenceID, std::vector<pufferfish::util::MemCluster>> &memClusters,
                      uint32_t maxSpliceGap, std::vector<pufferfish::util::UniMemInfo> &memCollection, uint32_t readLen,
                      spp::sparse_hash_map<pufferfish::common_types::ReferenceID, bool>& other_end_refs,
                      bool hChain, bool verbose = false);
};

#endif //PUFFERFISH_CHAINFINDER_H
