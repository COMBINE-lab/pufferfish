#ifndef _PUFFERFISHSPARSE_INDEX_HPP_
#define _PUFFERFISHSPARSE_INDEX_HPP_

#include <vector>

#include "rank9b.hpp"
#include "core/range.hpp"
#include "cereal/archives/json.hpp"

#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "BooPHF.h"
#include "Util.hpp"
#include "PufferfishBaseIndex.hpp"
#include "rank9sel.hpp"

class PufferfishSparseIndex : public PufferfishBaseIndex<PufferfishSparseIndex> {
  friend PufferfishBaseIndex;
  using hasher_t = pufferfish::types::hasher_t;
  using boophf_t = pufferfish::types::boophf_t;
  using pos_vector_t = PufferfishBaseIndex<PufferfishSparseIndex>::pos_vector_t;
  using seq_vector_t = PufferfishBaseIndex<PufferfishSparseIndex>::seq_vector_t;
  using edge_vector_t = PufferfishBaseIndex<PufferfishSparseIndex>::edge_vector_t;
  using bit_vector_t = PufferfishBaseIndex<PufferfishSparseIndex>::bit_vector_t;

private:
  uint32_t k_{0};
  uint32_t twok_{0};
  int32_t extensionSize_{0};
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
  compact::vector<uint64_t, 1> contigBoundary_;
  //bit_vector_t contigBoundary_;
  rank9sel rankSelDict;

  seq_vector_t seq_;
  edge_vector_t edge_;
  //for sparse representation
  compact::vector<uint64_t,1> presenceVec_;
  compact::vector<uint64_t,1> canonicalNess_;
  compact::vector<uint64_t,1> directionVec_ ;
  rank9b presenceRank_;
  compact::vector<uint64_t> extSize_{16};
  compact::vector<uint64_t> auxInfo_{16};
  pos_vector_t sampledPos_{16};

  std::unique_ptr<boophf_t> hash_{nullptr};

  static const constexpr uint64_t shiftTable_[] = {
    0x0, 0x7, 0x38, 0x1c0, 0xe00, 0x7000, 0x38000, 0x1c0000,
    0xe00000, 0x7000000, 0x38000000, 0x1c0000000, 0xe00000000,
    0x7000000000, 0x38000000000, 0x1c0000000000, 0xe00000000000,
    0x7000000000000, 0x38000000000000, 0x1c0000000000000,
    0xe00000000000000, 0x7000000000000000};

public:
  compact::vector<uint64_t, 2> refseq_;
  std::vector<uint64_t> refAccumLengths_;


  PufferfishSparseIndex();
  PufferfishSparseIndex(const std::string& indexPath);
  //Returns the position in the compacted bBG sequence from the sparse
  //index the above routine can be replaced by this code in
  //future versions
  // uint64_t getSparseRawPos(CanonicalKmer& mer);

  // Returns a ProjectedHits object that contains all of the
  // projected reference hits for the given kmer.
  auto getRefPos(CanonicalKmer mer) -> util::ProjectedHits;
  auto getRefPos(CanonicalKmer mer, util::QueryCache& qc) -> util::ProjectedHits;

private:
  auto getRefPosHelper_(CanonicalKmer& mer, uint64_t pos, bool didWalk = false) -> util::ProjectedHits;
  auto getRefPosHelper_(CanonicalKmer& mer, uint64_t pos, util::QueryCache& qc, bool didWalk = false) -> util::ProjectedHits;

};

#endif // _PUFFERFISH_INDEX_HPP_
