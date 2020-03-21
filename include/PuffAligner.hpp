#ifndef _PUFFALIGNER_
#define _PUFFALIGNER_

#include "tsl/hopscotch_map.h"
#include "metro/metrohash64.h"

#include "ProgOpts.hpp"
#include "Util.hpp"
#include "compact_vector/compact_vector.hpp"
#include "ksw2pp/KSW2Aligner.hpp"
#include "edlib.h"

#include "parallel_hashmap/phmap.h"

struct PassthroughHash {
	std::size_t operator()(uint64_t const& u) const { return u; }
};

/*enum class PuffAlignmentMode : uint8_t { SCORE_ONLY, CIGAR };

struct PuffAlignmentOptions {
  PuffAlignmentMode mode;
};*/

using HitCounters = pufferfish::util::HitCounters;
using AlignmentResult = pufferfish::util::AlignmentResult;
using AlnCacheMap = phmap::flat_hash_map<uint64_t, AlignmentResult, PassthroughHash>;

class PuffAligner {
public:
  PuffAligner(compact::vector<uint64_t, 2>& ar, std::vector<uint64_t>& ral, uint32_t k_, 
              pufferfish::util::AlignmentConfig& m, ksw2pp::KSW2Aligner& a) : 
    allRefSeq(ar), refAccumLengths(ral), k(k_), 
    mopts(m), aligner(a) {
    ksw_reset_extz(&ez);
		alnCacheLeft.reserve(32);
		alnCacheRight.reserve(32);
  }

/*
  PuffAligner(compact::vector<uint64_t, 2>& ar, std::vector<uint64_t>& ral, uint32_t k_,
              std::string r1, std::string r2, AlignmentOpts* m, ksw2pp::KSW2Aligner& a, bool mult) :
              allRefSeq(ar), refAccumLengths(ral), k(k_), read_left_(r1), read_right_(r2),
              mopts(m), aligner(a), multiMapping(mult) {
		memset(&ez, 0, sizeof(ksw_extz_t));

		alnCacheLeft.reserve(32);
		alnCacheRight.reserve(32);
  }
*/
  PuffAligner(const PuffAligner& other) = delete;
  PuffAligner& operator=(const PuffAligner& other) = delete;
  PuffAligner(PuffAligner&& other) = delete;
  PuffAligner& operator=(PuffAligner&& other) = delete;

  int32_t calculateAlignments(std::string& rl, std::string& rr, pufferfish::util::JointMems& jointHit, HitCounters& hctr, bool isMultimapping, bool verbose);
  int32_t calculateAlignments(std::string& read, pufferfish::util::JointMems& jointHit, HitCounters& hctr, bool isMultimapping, bool verbose);

  bool alignRead(std::string& read, std::string& read_rc, const std::vector<pufferfish::util::MemInfo>& mems, uint64_t queryChainHash, bool perfectChain, bool isFw, size_t tid, AlnCacheMap& alnCache, HitCounters& hctr, AlignmentResult& arOut, bool verbose);

  bool recoverSingleOrphan(std::string& rl, std::string& rr, pufferfish::util::MemCluster& clust, std::vector<pufferfish::util::MemCluster> &recoveredMemClusters, uint32_t tid, bool anchorIsLeft, bool verbose);

  void clearAlnCaches() {alnCacheLeft.clear(); alnCacheRight.clear();}
  void clear() {clearAlnCaches(); orphanRecoveryMemCollection.clear();  read_left_rc_.clear(); read_right_rc_.clear(); ksw_reset_extz(&ez); }

  std::vector<pufferfish::util::UniMemInfo> orphanRecoveryMemCollection;
private:
  compact::vector<uint64_t, 2>& allRefSeq;
  std::vector<uint64_t>& refAccumLengths;
  uint32_t k;
  pufferfish::util::AlignmentConfig mopts;
  ksw2pp::KSW2Aligner& aligner;
  ksw_extz_t ez;

  pufferfish::util::CIGARGenerator cigarGen_;
  std::string rc1_;
  std::string rc2_;
  std::string read_left_rc_;
  std::string read_right_rc_;
  std::string refSeqBuffer_;
  AlignmentResult ar_left;
  AlignmentResult ar_right;

  bool isMultimapping_;
  AlnCacheMap alnCacheLeft;
  AlnCacheMap alnCacheRight;
};


#endif // _PUFFALIGNER_
