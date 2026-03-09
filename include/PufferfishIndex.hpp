#ifndef _PUFFERFISH_INDEX_HPP_
#define _PUFFERFISH_INDEX_HPP_

#include <vector>

#include "cereal/archives/json.hpp"
#include "core/range.hpp"

#include "BooPHF.hpp"
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "Util.hpp"
#include "PufferfishBaseIndex.hpp"
#include "compact_vector/compact_vector.hpp"
#include "rank9sel.hpp"
#include "dictionary_types.hpp"
#include "streaming_query.hpp"

class PufferfishIndex : public PufferfishBaseIndex<PufferfishIndex> {
  friend PufferfishBaseIndex;
  using hasher_t = pufferfish::types::hasher_t;
  using boophf_t = pufferfish::types::boophf_t;
  using pos_vector_t = PufferfishBaseIndex<PufferfishIndex>::pos_vector_t;
  using seq_vector_t = PufferfishBaseIndex<PufferfishIndex>::seq_vector_t;
  using edge_vector_t = PufferfishBaseIndex<PufferfishIndex>::edge_vector_t;
  using bit_vector_t = PufferfishBaseIndex<PufferfishIndex>::bit_vector_t;


private:
  // SSHash dictionary (replaces hash_, pos_, seq_ for k-mer lookups)
  sshash::dictionary_type dict_;
  bool useSSHash_{false};

  uint32_t k_{0};
  uint32_t twok_{0};
  uint64_t numKmers_{0};
  bool haveEdges_{false};
  bool haveRefSeq_{false};
  bool haveEqClasses_{true};

  std::vector<uint32_t> eqClassIDs_;
  std::vector<std::vector<uint32_t>> eqLabels_;
  std::vector<std::string> refNames_;
  std::vector<uint32_t> refLengths_;
  std::vector<uint32_t> completeRefLengths_;
  std::vector<uint32_t> refExt_;
  pufferfish::ContigTable contigTable_;

  uint64_t numContigs_{0};
  // Legacy members kept for CRTP base class compilation compatibility.
  // Not loaded for SSHash-backed indexes.
  bit_vector_t  contigBoundary_;
  rank9sel rankSelDict;
  seq_vector_t seq_;
  edge_vector_t edge_;
  pos_vector_t pos_{16};

  std::unique_ptr<boophf_t> hash_{nullptr};
  boophf_t* hash_raw_{nullptr};
  size_t lastSeqPos_{std::numeric_limits<size_t>::max()};
  uint64_t numDecoys_{0};
  uint64_t firstDecoyIndex_{0};
  uint64_t firstDecoyEncodedIndex_{0};

public:
  PufferfishIndex();
  PufferfishIndex(const std::string& indexPath, pufferfish::util::IndexLoadingOpts opts = pufferfish::util::IndexLoadingOpts());
  ~PufferfishIndex();

  compact::vector<uint64_t, 2> refseq_;
  std::vector<uint64_t> refAccumLengths_;

  // Returns a ProjectedHits object that contains all of the
  // projected reference hits for the given kmer.
  auto getRefPos(CanonicalKmer& mer) -> pufferfish::util::ProjectedHits;
  // Returns a ProjectedHits object that contains all of the
  // projected reference hits for the given kmer.  Uses the results
  // of the previous contig info (start, end) from qc if the same
  // contig contains the match.  For correlated searches (e.g., from a read)
  // this can considerably speed up querying.
  auto getRefPos(CanonicalKmer& mer, pufferfish::util::QueryCache& qc) -> pufferfish::util::ProjectedHits;

  // Read 2-bit encoded sequence bits from the underlying sequence store.
  // For SSHash-backed indexes, reads from SSHash's internal bit_vector.
  // For legacy indexes, reads from compact seq_ vector.
  // bitPos and bitLen are in BITS (multiply base positions by 2).
  inline uint64_t getSeqBits(uint64_t bitPos, uint64_t bitLen) {
    if (useSSHash_) {
      return dict_.strings().get_bits(bitPos, bitLen);
    }
    return seq_.get_int(bitPos, bitLen);
  }

  // SSHash-specific accessors
  bool isSSHashBacked() const { return useSSHash_; }
  const sshash::dictionary_type& getDict() const { return dict_; }

  // SSHash streaming query support
  using sshash_streaming_query_t = sshash::streaming_query<sshash::dictionary_type, true>;

  // Create a streaming query object for sequential k-mer lookups on a read.
  // The caller should reset() between non-consecutive k-mers.
  sshash_streaming_query_t createStreamingQuery() const {
    return sshash_streaming_query_t(&dict_);
  }

  // Streaming query variant of getRefPos.  Uses the streaming query to avoid
  // full dictionary lookups for consecutive k-mers on the same SSHash string.
  // kmer_str must point to a k-length substring of the read.
  auto getRefPos(CanonicalKmer& mer, pufferfish::util::QueryCache& qc,
                 sshash_streaming_query_t& sq, char const* kmer_str)
      -> pufferfish::util::ProjectedHits;

  // Shadow base class methods to use SSHash when available
  CanonicalKmer getStartKmer(uint64_t cid);
  CanonicalKmer getEndKmer(uint64_t cid);
  uint32_t getContigLen(uint64_t rank);
  uint64_t getGlobalPos(uint64_t rank);
  std::string getSeqStr(size_t globalPos, int64_t length, bool isFw = true);
};

#endif // _PUFFERFISH_INDEX_HPP_
