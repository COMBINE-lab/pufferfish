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

private:
  uint32_t k_{0};
  uint32_t twok_{0};
  uint64_t numKmers_{0};
  uint64_t lastSeqPos_{0};
  uint64_t numSampledKmers_{0};
  std::vector<uint32_t> eqClassIDs_;
  std::vector<std::vector<uint32_t>> eqLabels_;
  std::vector<std::string> refNames_;
  std::vector<uint32_t> refLengths_;
  std::vector<std::vector<util::Position>> contigTable_;
  uint64_t numContigs_{0};
  sdsl::bit_vector contigBoundary_;
  sdsl::bit_vector::rank_1_type contigRank_;
  sdsl::bit_vector::select_1_type contigSelect_;
  sdsl::int_vector<2> seq_;
  sdsl::int_vector<8> edge_;

  //for lossy representation
  sdsl::bit_vector presenceVec_;
  sdsl::bit_vector::rank_1_type presenceRank_;
  sdsl::bit_vector::select_1_type presenceSelect_;
  sdsl::int_vector<> sampledPos_;

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
