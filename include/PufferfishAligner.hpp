#ifndef _PUFFERFISH_ALIGNER_
#define _PUFFERFISH_ALIGNER_


#include <vector>
#include <string>

#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"


#include "Util.hpp"
#include "KSW2Aligner.hpp"
#include "ProgOpts.hpp"

class PufferfishAligner {
public:
	PufferfishAligner(sdsl::int_vector<2>& ar, std::vector<uint64_t>& ral, uint32_t k_, AlignmentOpts* m, ksw2pp::KSW2Aligner& a) : allRefSeq(ar), refAccumLengths(ral), k(k_), mopts(m), aligner(a) {
		ksw2pp::KSW2Config config;
		config.dropoff = -1;
		config.gapo = mopts->gapOpenPenalty;
		config.gape = mopts->gapExtendPenalty;
		config.bandwidth = 10;
		config.flag = 0;
		aligner.config() = config;
		config.flag |= KSW_EZ_SCORE_ONLY;
		memset(&ez, 0, sizeof(ksw_extz_t));
	};
	int32_t calculateAlignments(std::string& read_left, std::string& read_right, util::JointMems& jointHit, bool verbose);
	int32_t alignmentScore(std::string& read, std::vector<util::MemInfo> mems, bool isFw, size_t tid, bool verbose);
private:
	sdsl::int_vector<2>& allRefSeq;
	std::vector<uint64_t>& refAccumLengths;
	uint32_t k;
	AlignmentOpts* mopts;
	ksw2pp::KSW2Aligner& aligner;
	ksw_extz_t ez;
};


#endif // _PUFFERFISH_ALIGNER_
