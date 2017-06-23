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
#include "Util.hpp"

class PufferfishIndex {
  using hasher_t = boomphf::SingleHashFunctor<uint64_t>;
  using boophf_t = boomphf::mphf<uint64_t, hasher_t>;
  using EqClassID = uint32_t;
  using EqClassLabel = std::vector<uint32_t>;

private:
  uint32_t k_{0};
  uint32_t twok_{0};
  uint64_t numKmers_{0};
  std::vector<uint32_t> eqClassIDs_;
  std::vector<std::vector<uint32_t>> eqLabels_;
  std::vector<std::string> refNames_;
  // std::vector<util::ContigPosInfo> cPosInfo_;
  std::vector<std::vector<util::Position>> contigTable_;
  uint64_t numContigs_{0};
  sdsl::bit_vector contigBoundary_;
  sdsl::bit_vector::rank_1_type contigRank_;
  sdsl::bit_vector::select_1_type contigSelect_;
  sdsl::int_vector<2> seq_;
  sdsl::int_vector<> pos_;
  std::unique_ptr<boophf_t> hash_{nullptr};
  boophf_t* hash_raw_{nullptr};
  size_t lastSeqPos_{std::numeric_limits<size_t>::max()};

public:
  PufferfishIndex();
  PufferfishIndex(const std::string& indexPath);

  // Get the equivalence class ID (i.e., rank of the equivalence class)
  // for a given contig.
  EqClassID getEqClassID(uint32_t contigID);

  // Get the equivalence class label for a contig (i.e., the set of reference
  // sequences containing
  // the contig).
  const EqClassLabel& getEqClassLabel(uint32_t contigID);

  // Get the k value with which this index was built.
  uint32_t k();
  // Get the list of reference sequences & positiosn corresponding to a contig
  const std::vector<util::Position>& refList(uint64_t contigRank);
  // Get the name of a given reference sequence
  const std::string& refName(uint64_t refRank);

  // Returns true if the given k-mer appears in the dBG, false otherwise
  bool contains(CanonicalKmer& mer);

  uint32_t contigID(CanonicalKmer& mer);

  // Returns the position in the compacted dBG sequence vector where the
  // given k-mer occurs, or std::numeric_limits<uint32_t>::max() otherwise.
  uint64_t getRawPos(CanonicalKmer& mer);
  // Returns true if pos is a valid position in the compacted sequence array
  // and false otherwise.
  bool isValidPos(uint64_t pos);
  // Returns a ProjectedHits object that contains all of the
  // projected reference hits for the given kmer.
  auto getRefPos(CanonicalKmer& mer) -> util::ProjectedHits;
};

#endif // _PUFFERFISH_INDEX_HPP_
