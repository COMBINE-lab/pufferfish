#ifndef _PUFFALIGNER_
#define _PUFFALIGNER_

#include "tsl/hopscotch_map.h"
#include "metro/metrohash64.h"

#include "Util.hpp"
#include "KSW2Aligner.hpp"
#include "ProgOpts.hpp"
#include "compact_vector/compact_vector.hpp"

struct PassthroughHash {
	std::size_t operator()(uint64_t const& u) const { return u; }
};

using HitCounters = util::HitCounters;
using AlignmentResult = util::AlignmentResult;
using AlnCacheMap = tsl::hopscotch_map<uint64_t, AlignmentResult, PassthroughHash>;

class PuffAligner {
public:
	PuffAligner(compact::vector<uint64_t, 2>& ar, std::vector<uint64_t>& ral, uint32_t k_, AlignmentOpts* m, ksw2pp::KSW2Aligner& a, bool mult) : allRefSeq(ar), refAccumLengths(ral), k(k_), mopts(m), aligner(a), multiMapping(mult) {
		ksw2pp::KSW2Config config;
		config.dropoff = -1;
		config.gapo = mopts->gapOpenPenalty;
		config.gape = mopts->gapExtendPenalty;
		config.bandwidth = 10;
		config.flag = 0;
		aligner.config() = config;
		//config.flag |= KSW_EZ_SCORE_ONLY;
  	config.flag |= KSW_EZ_RIGHT;
    aligner.setConfig(config);
		memset(&ez, 0, sizeof(ksw_extz_t));

		alnCacheLeft.reserve(32);
		alnCacheRight.reserve(32);
	};
	int32_t calculateAlignments(std::string& read_left, std::string& read_right, util::JointMems& jointHit, HitCounters& hctr, bool verbose);
	AlignmentResult alignRead(std::string read, std::vector<util::MemInfo>& mems, bool perfectChain, bool isFw, size_t tid, AlnCacheMap& alnCache, HitCounters& hctr, bool verbose);
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


#endif // _PUFFALIGNER_
