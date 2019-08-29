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
#include "parallel_hashmap/phmap.h"
#include "chobo/small_vector.hpp"
#include "MemChainer.hpp"


template <typename PufferfishIndexT> class MemCollector {
enum class ExpansionTerminationType : uint8_t { MISMATCH = 0, CONTIG_END, READ_END };  

public:
  explicit MemCollector(PufferfishIndexT* pfi) : pfi_(pfi) { k = pfi_->k(); }

  size_t expandHitEfficient(pufferfish::util::ProjectedHits& hit,
                          pufferfish::CanonicalKmerIterator& kit,
                          ExpansionTerminationType& et);

  bool operator()(std::string &read,
                  pufferfish::util::QueryCache& qc,
                  bool isLeft=false,
                  bool verbose=false);

  bool findChains(std::string &read,
                  pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& memClusters,
                  //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>>& memClusters,
                  uint32_t maxSpliceGap,
                  pufferfish::util::MateStatus mateStatus,
                  bool hChain=false,
                  bool isLeft=false,
                  bool verbose=false);
  void clear();

  void configureMemClusterer(uint32_t max);

  void setConsensusFraction(double cf);

  double getConsensusFraction() const;

  void setHitFilterPolicy(pufferfish::util::HitFilterPolicy hfp);

  pufferfish::util::HitFilterPolicy getHitFilterPolicy() const;

private:
  PufferfishIndexT* pfi_;
  size_t k;
  //AlignerEngine ae_;
  std::vector<pufferfish::util::UniMemInfo> memCollectionLeft;
  std::vector<pufferfish::util::UniMemInfo> memCollectionRight;
  bool isSingleEnd = false;
  MemClusterer mc;
  //pufferfish::common_types::RefMemMapT trMemMap;
  pufferfish::util::CachedVectorMap<std::pair<pufferfish::common_types::ReferenceID, bool>, std::vector<pufferfish::util::MemInfo>, pufferfish::util::pair_hash> trMemMap;


  phmap::flat_hash_map<pufferfish::common_types::ReferenceID, bool> left_refs;
  phmap::flat_hash_map<pufferfish::common_types::ReferenceID, bool> right_refs;

  std::vector<std::pair<int, pufferfish::util::ProjectedHits>> left_rawHits;
  std::vector<std::pair<int, pufferfish::util::ProjectedHits>> right_rawHits;
};
#endif
