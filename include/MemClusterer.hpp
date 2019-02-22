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
  bool clusterMems(std::vector<std::pair<int, util::ProjectedHits>> &hits,
                     spp::sparse_hash_map<pufferfish::common_types::ReferenceID, std::vector<util::MemCluster>> &memClusters,
                     uint32_t maxSpliceGap, std::vector<util::UniMemInfo> &memCollection, bool verbose = false);

  bool fillMemCollection(std::vector<std::pair<int, util::ProjectedHits>> &hits,
                           std::map<std::pair<pufferfish::common_types::ReferenceID, bool>, std::vector<util::MemInfo>> &trMemMap,
                           std::vector<util::UniMemInfo> &memCollection, util::ReadEnd re, bool verbose = false);


  // We have both left and right end mems
  bool findOptChain(std::map<std::pair<pufferfish::common_types::ReferenceID, bool>, std::vector<util::MemInfo>> &trMemMap,
                 std::vector<util::JointMems> &jointMemsList,
                 std::vector<util::MemCluster> &all,
                 uint32_t maxSpliceGap,
                 uint32_t maxFragmentLength,
                 bool verbose = false);

	template<typename PufferfishIndexT>
  bool findOptChainAllowingOneJumpBetweenTheReadEnds(
            std::map<std::pair<pufferfish::common_types::ReferenceID, bool>, std::vector<util::MemInfo>> &trMemMap,
            std::vector<util::JointMems> &jointMemsList,
            std::vector<util::MemCluster> &all,
            uint32_t maxSpliceGap,
            uint32_t maxFragmentLength,
						PufferfishIndexT* pfi_,
						uint32_t readLenLeft,
            uint32_t readLenRight,
            bool verbose = false);

  bool findOptChain(std::vector<std::pair<int, util::ProjectedHits>> &hits,
                      spp::sparse_hash_map<pufferfish::common_types::ReferenceID, std::vector<util::MemCluster>> &memClusters,
                      uint32_t maxSpliceGap, std::vector<util::UniMemInfo> &memCollection, bool verbose = false);
};

#endif //PUFFERFISH_CHAINFINDER_H
