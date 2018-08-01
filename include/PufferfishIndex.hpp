#ifndef _PUFFERFISH_INDEX_HPP_
#define _PUFFERFISH_INDEX_HPP_

#include <vector>

#include "cereal/archives/json.hpp"
#include "core/range.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"

#include "BooPHF.h"
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "Util.hpp"
#include "PufferfishBaseIndex.hpp"

class PufferfishIndex : public PufferfishBaseIndex<PufferfishIndex> {
  friend PufferfishBaseIndex;
  using hasher_t = pufferfish::types::hasher_t;
  using boophf_t = pufferfish::types::boophf_t;

private:
  uint32_t k_{0};
  uint32_t twok_{0};
  uint64_t numKmers_{0};
  bool haveEdges_{false};
  std::vector<uint32_t> eqClassIDs_;
  std::vector<std::vector<uint32_t>> eqLabels_;
  std::vector<std::string> refNames_;
  std::vector<uint32_t> refLengths_;
  std::vector<util::Position> contigTable_;
  std::vector<uint64_t> contigOffsets_;

  uint64_t numContigs_{0};
  sdsl::bit_vector contigBoundary_;
  sdsl::bit_vector::rank_1_type contigRank_;
  sdsl::bit_vector::select_1_type contigSelect_;
  sdsl::int_vector<2> seq_;
  sdsl::int_vector<8> edge_;
  sdsl::int_vector<> pos_;

  std::unique_ptr<boophf_t> hash_{nullptr};
  boophf_t* hash_raw_{nullptr};
  size_t lastSeqPos_{std::numeric_limits<size_t>::max()};

public:
  PufferfishIndex();
  PufferfishIndex(const std::string& indexPath);

  // Returns a ProjectedHits object that contains all of the
  // projected reference hits for the given kmer.
  auto getRefPos(CanonicalKmer& mer) -> util::ProjectedHits;
  // Returns a ProjectedHits object that contains all of the
  // projected reference hits for the given kmer.  Uses the results
  // of the previous contig info (start, end) from qc if the same
  // contig contains the match.  For correlated searches (e.g., from a read)
  // this can considerably speed up querying.
  auto getRefPos(CanonicalKmer& mer, util::QueryCache& qc) -> util::ProjectedHits;
};

#endif // _PUFFERFISH_INDEX_HPP_
