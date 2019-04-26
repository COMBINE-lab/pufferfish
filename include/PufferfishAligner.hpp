#ifndef _PUFFERFISH_ALIGNER_
#define _PUFFERFISH_ALIGNER_


#include <vector>
#include <string>

#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"

#include "tsl/hopscotch_map.h"

#include "Util.hpp"
#include "KSW2Aligner.hpp"
#include "ProgOpts.hpp"
#include "compact_vector/compact_vector.hpp"

using HitCounters = util::HitCounters;
 
struct PassthroughHash {
	std::size_t operator()(uint64_t const& u) const { return u; }
};
using AlnCacheMap = tsl::hopscotch_map<uint64_t, int32_t, PassthroughHash>;

class PufferfishAligner {
public:
	PufferfishAligner(compact::vector<uint64_t, 2>& ar, std::vector<uint64_t>& ral, uint32_t k_, AlignmentOpts* m, ksw2pp::KSW2Aligner& a, bool mult) : allRefSeq(ar), refAccumLengths(ral), k(k_), mopts(m), aligner(a), multiMapping(mult) {
		ksw2pp::KSW2Config config;
		config.dropoff = -1;
		config.gapo = mopts->gapOpenPenalty;
		config.gape = mopts->gapExtendPenalty;
		config.bandwidth = 10;
		config.flag = 0;
		aligner.config() = config;
		config.flag |= KSW_EZ_SCORE_ONLY;
		memset(&ez, 0, sizeof(ksw_extz_t));

		alnCacheLeft.reserve(32);
		alnCacheRight.reserve(32);
	};
	int32_t calculateAlignments(std::string& read_left, std::string& read_right, util::JointMems& jointHit, HitCounters& hctr, bool verbose);
	int32_t alignRead(std::string read, std::vector<util::MemInfo>& mems, bool perfectChain, bool isFw, size_t tid, AlnCacheMap& alnCache, HitCounters& hctr, bool verbose);
	void clearAlnCaches() {alnCacheLeft.clear(); alnCacheRight.clear();}
private:
	compact::vector<uint64_t, 2>& allRefSeq;
	std::vector<uint64_t>& refAccumLengths;
	uint32_t k;
	AlignmentOpts* mopts;
	ksw2pp::KSW2Aligner& aligner;
	ksw_extz_t ez;

	bool multiMapping;
	AlnCacheMap alnCacheLeft;
	AlnCacheMap alnCacheRight;
};


#endif // _PUFFERFISH_ALIGNER_
