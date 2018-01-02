//#include <limits>

#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"

#include "PufferfishIndex.hpp"
#include "xxhash.h"

uint64_t hashIt(uint64_t val, std::map<uint64_t, uint64_t>& val2hash) {
  // return val;
  // return val2hash[val];
  return XXH64(&val, sizeof(uint64_t), 0);
}

uint64_t decodeHash(uint64_t hash, std::map<uint64_t, uint64_t>& hash2bucketIdx) {
  //return hash;
  return hash2bucketIdx[hash];
}

struct Bucket {
  uint64_t numOfKmers = 0;
  uint64_t numOfUnitigs = 0;
  uint64_t seqLength = 0;
};
void findMinimizer(uint64_t kmer, uint8_t k, uint8_t m, uint64_t* minimizer, uint64_t* minimizerPos, std::map<uint64_t, uint64_t>& val2hash) {
  uint8_t mask = (0x01 << 2*m)-1;
  //std::cout << " -- sw " << *minimizer;
  for (uint64_t cntr = 0; cntr <= (uint64_t)k-m; cntr++) {
    //std::cout << "c" << cntr << " " << ((kmer >> 2*cntr) & mask) << " ";
    uint64_t val = (kmer >> 2*cntr) & mask;
    uint64_t minHash = hashIt(val, val2hash);
    //std::cout << val << " " << minHash << " ";
    if ((*minimizer) > minHash) {
      (*minimizer) = minHash;
      (*minimizerPos) = cntr;
    }
  }
  //std::cout << " --- newMin : " << *minimizer << " ";
}

int main(int argc, char* argv[]) {

  PufferfishIndex pfi(argv[1]);
  std::cout << "\nloaded!\n\n";
  uint8_t k = pfi.k();
  sdsl::int_vector<2>& seq = pfi.getSeq();
  sdsl::bit_vector::select_1_type& unitigS = pfi.getBoundarySelect();

  for (uint8_t m = 2; m < 13; m++) {
    std::cout << "\n\nMINIMIZER LENGTH : "<< static_cast<size_t>(m) << "\n\n";
    std::map<uint64_t, uint64_t> hash2bucketIdx;
    std::map<uint64_t, uint64_t> val2hash;
    size_t bucketCnt = (0x01 << 2*m);
    for (uint64_t i = 0; i < bucketCnt; i++) {
      uint64_t val = i;
      size_t a = XXH64(&val, sizeof(uint64_t), 0);
      hash2bucketIdx[a] = i;
      val2hash[i] = a;
      //std::cout << i << " " << a << "\n";
    }
    uint64_t contigCntr = 1;
    uint64_t start = 0;
    uint64_t cur = 0;
    uint64_t next = unitigS(contigCntr) + 1;//need +1 to get to the start of next contig
    /*sdsl::bit_vector & boundaryBV = pfi.getBoundary();
    std::cout << " A small test\n next : " << next << "\n";
    uint64_t wrdCnt = std::min(next, (uint64_t)64);
    for (size_t i = 0; i < next; i+=wrdCnt) {
      wrdCnt = std::min(next-i, (uint64_t)64);
      uint64_t wrd = boundaryBV.get_int(i, wrdCnt);
      for (size_t wrdCntr = 0; wrdCntr < wrdCnt; wrdCntr++)
        if (wrd >> wrdCntr & 0x01) std::cout << "ours: " << i+wrdCntr;
    }
    std::cout << "\nEnd of test.\n";
    std::exit(1);
    */
    uint64_t minimizer = -1;//std::numeric_limits<uint64_t>::max; //0x01 << 2*m;
    uint64_t minimizerPos = 0;
    uint64_t totalNumOfKmers = (next - cur) -k+1;
    uint64_t totalSeqSize = next-cur;
    //std::cout << totalNumOfKmers << " " ;

    // findMinimizer test:
    /*x
    uint8_t tmpMinimizer = 0x01 << 2*m;
    uint64_t relMinimizerPos = 0;
    findMinimizer(0xFFFFFFFEFFFFFFFF, k, m, &tmpMinimizer, &relMinimizerPos);
    std::cout << "IT SHOULD BE 11 AND IS " << (size_t)tmpMinimizer << "\n";
    */

    std::vector<Bucket> buckets(0x01 << 2*m);
    // initial all pairs to zero
    uint64_t unitigSplitCnt = 0; // total number of splits happening across all unitigs
    uint64_t splittedUnitigCnt = 0; // total number of unitigs that we had a split in
    uint64_t tmpUnitigSplitCnt = 0; // number of splits in current unitig
    uint64_t prevPos = 0;
    while (cur <= seq.size()-k) {
      // if current kmer is the first kmer in unitig or we've surpassed current minimizer
      if (cur == start or cur > minimizerPos) {
        //std::cout << "2\n ";
        uint64_t kmer = seq.get_int(2*cur, 2*k);
        uint64_t tmpMinimizer = -1;//std::numeric_limits<uint64_t>::max; //0x01 << 2*m; // reset minimizer to max
        uint64_t relMinimizerPos = 0;
        findMinimizer(kmer, k, m, &tmpMinimizer, &relMinimizerPos, val2hash);
        minimizerPos = cur + relMinimizerPos;
        // If we are not at a new unitig and the minimizer is different from previous one, this means a split in unitig
        if (cur != start and minimizer != tmpMinimizer) {
          // store the info about the new unitig up until this minimizer
          //std::cout << "2 " << minimizer << " " << hash2bucketIdx[minimizer] << " -- ";
          buckets[decodeHash(minimizer, hash2bucketIdx)].seqLength += (cur+(k-1)-prevPos);
          buckets[decodeHash(minimizer, hash2bucketIdx)].numOfKmers += cur-prevPos;
          buckets[decodeHash(minimizer, hash2bucketIdx)].numOfUnitigs++;
          tmpUnitigSplitCnt++; // increase unitig split count
          prevPos = cur;
        }
        minimizer = tmpMinimizer;
        //std::cout << (size_t)minimizer << "\n";
      }
      else {
        //std::cout << "3 \n";
        uint64_t val = seq.get_int(2*(cur+(k-m)), 2*m);
        uint64_t nextPotentialMinimizer = hashIt(val, val2hash); //val2hash[val];//XXH64(&val, 64, 0);
        if (nextPotentialMinimizer < minimizer) {
          tmpUnitigSplitCnt++;
          minimizer = nextPotentialMinimizer;
          minimizerPos = cur+(k-m);
          //std::cout << "3 " << minimizer << " " << hash2bucketIdx[minimizer] << " -- ";
          buckets[decodeHash(minimizer, hash2bucketIdx)].seqLength += (cur+(k-1)-prevPos);
          buckets[decodeHash(minimizer, hash2bucketIdx)].numOfKmers += (cur-prevPos);
          buckets[decodeHash(minimizer, hash2bucketIdx)].numOfUnitigs++;
          prevPos = cur;
        }
      }
      cur++;
      // if we've already passed the last kmer in current unitig,
      // first add the last part of the contig not added (at least one kmer)
      // then go to the next unitig
      if (cur >= next - (k-1)) {
        ///std::cout << "4 \n";
        // got to the end of the unitig. Add the new contig to the right minimizer bucket
        //std::cout << "4 " << minimizer << " " << hash2bucketIdx[minimizer] << " -- ";
        buckets[decodeHash(minimizer, hash2bucketIdx)].seqLength += (cur+(k-1)-prevPos);
        buckets[decodeHash(minimizer, hash2bucketIdx)].numOfKmers += (cur-prevPos);
        buckets[decodeHash(minimizer, hash2bucketIdx)].numOfUnitigs++;
        if (tmpUnitigSplitCnt > 0) {
          splittedUnitigCnt++;
          unitigSplitCnt += tmpUnitigSplitCnt;
          // reset current unitig split cnt to zero
          tmpUnitigSplitCnt = 0;
        }
        minimizer = -1;//std::numeric_limits<uint64_t>::max; //0x01 << 2*m; // reset minimizer to max for next contig
        start = next;
        cur = next;
        prevPos = cur;
        next = unitigS(++contigCntr)+1;
        if (next == 1) {
          //next = seq.size();
          //std::cout << next << " " << cur << "\n";
          //std::cout << " -- break -- ";
          break;
        }
        totalNumOfKmers += next-cur-k+1;
        totalSeqSize += next-cur;
      }


    }
    std::cout << "\n";
    uint64_t totalBits = totalNumOfKmers*ceil(log2(totalNumOfKmers+1)) + (seq.size()*2);
    std::cout << "# of unitigs before splitting: " << contigCntr << "\n"
              << "# of unitigs after splitting: " << unitigSplitCnt << "\n"
              << "# of unitigs that were splitted: " << splittedUnitigCnt << "\n";
    std::cout << "\nBefore applying minimizers:\n"
              << "# of kmers: " << totalNumOfKmers << " , "
              << "bits per kmer pos "<< ceil(log2(totalNumOfKmers)) << "\n"
              << "total sequence length: " << totalSeqSize << "\n"
              <<"total bits: " << totalBits << " or " << totalBits/(1024*1024*8) << "MB\n";


    size_t seqSize = 0;
    totalBits = 0;
    uint64_t sumNumKmers = 0;
    std::cout << "\nAfter applying minimizers:\n";
    uint8_t bCntr = 0;
    for (auto bIt = buckets.begin(); bIt != buckets.end(); bIt++) {
      if (bIt->numOfKmers != 0) {
        sumNumKmers += bIt->numOfKmers;
        totalBits += bIt->seqLength*2 + bIt->numOfKmers*ceil(log2(bIt->numOfKmers+1));
        seqSize += bIt->seqLength;
        /*std::cout << "b" << static_cast<size_t>(bCntr++) << ":"
                << "u" << bIt->numOfUnitigs << ","
                << "k" << bIt->numOfKmers << ","
                << ceil(log2(bIt->numOfKmers+1)) << "\t";
        */
      }
    }
    std::cout <<"\ntotal number of kmers: " << sumNumKmers << "\n"
              <<"Sum of all sequence lengths: " << seqSize << "\n"
              <<"total bits: " << totalBits << " or " << totalBits/(1024*1024*8) << "MB\n";
  }

}
