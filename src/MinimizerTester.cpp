
#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"

#include "PufferfishIndex.hpp"

void findMinimizer(uint64_t kmer, uint8_t k, uint8_t m, uint8_t* minimizer, uint64_t* minimizerPos) {
  uint8_t mask = (0x01 << 2*m)-1;
  for (uint8_t cntr = 0; cntr < k; cntr+=m) {
    if ((*minimizer) > ((kmer >> 2*cntr) & mask)) {
      (*minimizer) = (kmer >> 2*cntr) & mask;
      (*minimizerPos) = cntr;
    }
  }
}

int main(int argc, char* argv[]) {

  PufferfishIndex pfi(argv[1]);
  std::cout << "Successfully read.\n";
  uint8_t k = pfi.k();
  sdsl::int_vector<2>& seq = pfi.getSeq();
  sdsl::bit_vector::select_1_type& unitigS = pfi.getBoundarySelect();

  for (uint8_t m = 2; m < 3; m++) {
    uint64_t contigCntr = 1;
    uint64_t start = 0;
    uint64_t cur = 0;
    uint64_t next = unitigS(contigCntr);
    uint8_t minimizer = 0x01 << 2*m;
    uint64_t minimizerPos = 0;

    // bucket values:
    // first: sum of unitig lengths
    // second: unitig cnt
    std::vector<std::pair<uint64_t, uint64_t>> buckets(0x01 << 2*m);
    // initial all pairs to zero
    for (auto i = 0; i < buckets.size(); i++) buckets[i] = std::make_pair(0, 0);
    uint64_t unitigSplitCnt = 0; // total number of splits happening across all unitigs
    uint64_t splittedUnitigCnt = 0; // total number of unitigs that we had a split in
    uint64_t tmpUnitigSplitCnt = 0; // number of splits in current unitig
    uint64_t prevPos = 0;
    while (cur < seq.size()) {
      // if we've already passed the last kmer in current unitig, go to next unitig
      if (cur >= next - (k-1)) {
        start = next;
        cur = next;
        next = unitigS(contigCntr++);
        if (tmpUnitigSplitCnt > 0) {
          splittedUnitigCnt++;
          unitigSplitCnt += tmpUnitigSplitCnt;
          // reset current unitig split cnt to zero
          tmpUnitigSplitCnt = 0;
        }
        if (cur >= seq.size()-(k-1)) break;
      }

      // if current kmer is the first kmer in unitig or we've surpassed current minimizer
      if (cur == start or cur > minimizerPos) {
        uint64_t kmer = seq.get_int(2*cur, 2*k);
        uint8_t tmpMinimizer = 0x01 << 2*m;
        uint64_t relMinimizerPos = 0;
        findMinimizer(kmer, k, m, &tmpMinimizer, &relMinimizerPos);
        minimizerPos = cur + relMinimizerPos;
        // If we are not at a new unitig and the minimizer is different from previous one, this means a split in unitig
        if (cur != start and minimizer != tmpMinimizer) {
          tmpUnitigSplitCnt++;
          prevPos = cur;
        }
      }
      else {
        uint64_t nextPotentialMinPos = cur+(k-(m-1));
        uint8_t nextPotentialMinimizer = seq.get_int(2*nextPotentialMinPos, 2*m);
        if (nextPotentialMinimizer < minimizer) {
          tmpUnitigSplitCnt++;
          minimizer = nextPotentialMinimizer;
          minimizerPos = nextPotentialMinPos;
          buckets[minimizer].first += (nextPotentialMinimizer-prevPos+1);
          buckets[minimizer].second++;
          prevPos = cur+1;
        }
      }
      cur++;
    }
  }

}
