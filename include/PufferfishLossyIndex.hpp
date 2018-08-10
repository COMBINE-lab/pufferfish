#ifndef _PUFFERFISHLOSSY_INDEX_HPP_
#define _PUFFERFISHLOSSY_INDEX_HPP_

#include <vector>

#include "core/range.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"
#include "cereal/archives/json.hpp"

#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "BooPHF.h"
#include "Util.hpp"
#include "PufferfishBaseIndex.hpp"

class PufferfishLossyIndex : public PufferfishBaseIndex<PufferfishLossyIndex> {
  friend PufferfishBaseIndex;
  using hasher_t = pufferfish::types::hasher_t;
  using boophf_t = pufferfish::types::boophf_t;
  using pos_vector_t = PufferfishBaseIndex<PufferfishLossyIndex>::pos_vector_t;
  using seq_vector_t = PufferfishBaseIndex<PufferfishLossyIndex>::seq_vector_t;
  using edge_vector_t = PufferfishBaseIndex<PufferfishLossyIndex>::edge_vector_t;
  using bit_vector_t = PufferfishBaseIndex<PufferfishLossyIndex>::bit_vector_t;

private:
  uint32_t k_{0};
  uint32_t twok_{0};
  uint64_t numKmers_{0};
  uint64_t lastSeqPos_{0};
  uint64_t numSampledKmers_{0};
  bool haveEdges_{false};
  std::vector<uint32_t> eqClassIDs_;
  std::vector<std::vector<uint32_t>> eqLabels_;
  std::vector<std::string> refNames_;
  std::vector<uint32_t> refLengths_;
  std::vector<util::Position> contigTable_;
  std::vector<uint64_t> contigOffsets_;
  uint64_t numContigs_{0};
  bit_vector_t contigBoundary_;
  sdsl::bit_vector::rank_1_type contigRank_;
  sdsl::bit_vector::select_1_type contigSelect_;
  seq_vector_t seq_;
  edge_vector_t edge_;

  //for lossy representation
  bit_vector_t presenceVec_;
  sdsl::bit_vector::rank_1_type presenceRank_;
  sdsl::bit_vector::select_1_type presenceSelect_;
  pos_vector_t sampledPos_;

  std::unique_ptr<boophf_t> hash_{nullptr};
  boophf_t* hash_raw_{nullptr};

public:
  PufferfishLossyIndex();
  PufferfishLossyIndex(const std::string& indexPath);
  // Returns a ProjectedHits object that contains all of the
  // projected reference hits for the given kmer.
  auto getRefPos(CanonicalKmer& mer) -> util::ProjectedHits;
  auto getRefPos(CanonicalKmer& mer, util::QueryCache& qc) -> util::ProjectedHits;
};

#endif // _PUFFERFISH_INDEX_HPP_
