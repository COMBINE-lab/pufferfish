#ifndef HIT_COLLECTOR_HPP
#define HIT_COLLECTOR_HPP

#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "PufferfishLossyIndex.hpp"
#include "Util.hpp"
#include "edlib.h"
#include "jellyfish/mer_dna.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <sparsepp/spp.h>

#include "MemClusterer.hpp"


template <typename PufferfishIndexT> class MemCollector {
enum class ExpansionTerminationType : uint8_t { MISMATCH = 0, CONTIG_END, READ_END };  

public:
  MemCollector(PufferfishIndexT* pfi) : pfi_(pfi) { k = pfi_->k(); }

  void recoverGaps(spp::sparse_hash_map<pufferfish::common_types::ReferenceID, 
                                        std::vector<util::MemCluster>>& memClustersMap,
                                        std::vector<util::UniMemInfo>& memCollection,
                                        size_t rlen, 
                                        bool verbose=false);

  size_t expandHitEfficient(util::ProjectedHits& hit,
                          pufferfish::CanonicalKmerIterator& kit,
                          ExpansionTerminationType& et);
    
  bool operator()(std::string &read,
                  //spp::sparse_hash_map<size_t, std::vector<util::MemCluster>>& memClusters,
                  uint32_t maxSpliceGap,
                  util::MateStatus mateStatus,
                  util::QueryCache& qc,
                  bool verbose=false);

  void findBestChain(std::vector<util::JointMems> &jointHits,
                       std::vector<util::MemCluster> &all,
                       uint32_t maxSpliceGap, uint32_t maxFragmentLength, bool verbose); 

  void findOptChainAllowingOneJumpBetweenTheReadEnds(std::vector<util::JointMems> &jointHits,
                                                       std::vector<util::MemCluster> &all,
                                                       uint32_t maxSpliceGap, uint32_t maxFragmentLength,
                                                       uint32_t readLenLeft,
                                                       uint32_t readLenRight,
                                                       bool verbose);

  void clear();

private:
  PufferfishIndexT* pfi_;
  size_t k;
  //AlignerEngine ae_;
  std::vector<util::UniMemInfo> memCollectionLeft;
  std::vector<util::UniMemInfo> memCollectionRight;
  bool isSingleEnd = false;
  MemClusterer mc;
  std::map<std::pair<pufferfish::common_types::ReferenceID, bool>, std::vector<util::MemInfo>> trMemMap;
};
#endif
