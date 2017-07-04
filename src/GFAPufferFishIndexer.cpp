#include "FastxParser.hpp"
#include <cmath>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>
#include <sstream>
#include <bitset>

#include "CanonicalKmer.hpp"
#include "OurGFAReader.hpp"
#include "PufferFS.hpp"
#include "PufferfishIndex.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "cereal/archives/json.hpp"
#include "jellyfish/mer_dna.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"
#include "spdlog/spdlog.h"
//#include "gfakluge.hpp"

uint64_t swap_uint64(uint64_t val) {
  val = ((val << 8) & 0xFF00FF00FF00FF00ULL) |
        ((val >> 8) & 0x00FF00FF00FF00FFULL);
  val = ((val << 16) & 0xFFFF0000FFFF0000ULL) |
        ((val >> 16) & 0x0000FFFF0000FFFFULL);
  return (val << 32) | (val >> 32);
}

// adapted from :
// http://stackoverflow.com/questions/34875315/implementation-my-own-list-and-iterator-stl-c
class ContigKmerIterator {
public:
  typedef ContigKmerIterator self_type;
  typedef uint64_t value_type;
  typedef value_type& reference;
  typedef value_type* pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef int64_t difference_type;

  ContigKmerIterator(sdsl::int_vector<2>* storage, sdsl::bit_vector* rank,
                     uint8_t k, uint64_t startAt)
      : storage_(storage), rank_(rank), k_(k), curr_(startAt) {
    if (curr_ + k_ <= rank_->size()) {
      mer_.fromNum(storage_->get_int(2 * curr_, 2 * k_));
      // mer_.word__(0) = storage_->get_int(2 * curr_, 2 * k_);
    }
    // rcMer_ = mer_.get_reverse_complement();
  }

  ContigKmerIterator&
  operator=(ContigKmerIterator& other) { //}= default;//(sdsl::int_vector<>*
                                         // storage, sdsl::bit_vector* rank,
    // uint8_t k, uint64_t startAt) :
    storage_ = other.storage_;
    rank_ = other.rank_;
    k_ = other.k_;
    curr_ = other.curr_;
    mer_ = other.mer_;
    // rcMer_ = other.rcMer_;
    word_ = other.word_;
    return *this;
  }

  ContigKmerIterator operator++() {
    ContigKmerIterator i = *this;
    advance_();
    return i;
  }

  ContigKmerIterator operator++(int) {
    advance_();
    return *this;
  }

  reference operator*() {
    // word_ = (mer_.word(0) < rcMer_.word(0)) ? mer_.word(0) : rcMer_.word(0);
    word_ = mer_.getCanonicalWord();
    return word_;
  }

  difference_type pos() { return curr_; }

  bool isCanonical(){
	  return mer_.fwWord() == mer_.getCanonicalWord() ;
  }

  bool isEndKmer() {
	  size_t endPos = curr_ + k_ - 1;
	  return ((*rank_)[endPos] == 1);
  }

  pointer operator->() {
    word_ = mer_.getCanonicalWord(); //(mer_.word(0) < rcMer_.word(0)) ?
                                     // mer_.word(0) : rcMer_.word(0);
    return &word_;
  }
  bool operator==(const self_type& rhs) { return curr_ == rhs.curr_; }

  bool operator!=(const self_type& rhs) { return curr_ != rhs.curr_; }

  bool operator<(const self_type& rhs) { return curr_ < rhs.curr_; }

  bool operator<=(const self_type& rhs) { return curr_ <= rhs.curr_; }

private:
  void advance_() {
    size_t endPos = curr_ + k_ - 1;
    if (endPos + 1 < rank_->size() and (*rank_)[endPos] == 1) {
      curr_ += k_;
      mer_.fromNum(storage_->get_int(2 * curr_, 2 * k_));
    } else {
      if (curr_ + k_ < rank_->size()) {
        int c = (*storage_)[curr_ + k_];
        mer_.shiftFw(c);
      } else {
        mer_.fromNum(storage_->get_int(2 * (rank_->size() - k_), 2 * k_));
      }
      ++curr_;
    }
  }
  sdsl::int_vector<2>* storage_{nullptr};
  sdsl::bit_vector* rank_{nullptr};
  uint8_t k_{0};
  uint64_t curr_{0};
  CanonicalKmer mer_;
  uint64_t word_{0};
};

int pufferfishTest(util::TestOptions& testOpts) {
  (void)testOpts;
  std::cerr << "this command is not yet implemented\n";
  return 1;
}

void computeSampledPositions(size_t tlen, uint32_t k, int sampleSize, std::vector<size_t>& sampledInds){
  sampledInds.clear() ;
  auto numOfKmers = tlen - k;
  size_t lastCovered = 0 ;
  for(size_t j = 0 ; j <= numOfKmers; j++){
    if(j > lastCovered){
      auto next_samp = std::min(j + sampleSize/2 - 1 ,numOfKmers) ;
      sampledInds.push_back(next_samp) ;
      lastCovered = next_samp + sampleSize/2 + 1;
    }
  }
  if(lastCovered == numOfKmers)
    sampledInds.push_back(lastCovered) ;
}

std::string packedToString(sdsl::int_vector<2>& seqVec, uint64_t offset, uint32_t len) {
  std::stringstream s;
  for (size_t i = offset; i < offset + len; ++i) {
    auto c = seqVec[i];
    s << my_mer::rev_code(c);
  }
  auto st = s.str();
  return st;
}

enum class NextSampleDirection : uint8_t { FORWARD = 0, REVERSE=1 };

uint32_t getEncodedExtension(sdsl::int_vector<2>& seqVec, uint64_t firstSampPos, uint64_t distToSamplePos,
                             uint32_t maxExt, NextSampleDirection dir, bool print=false) {
  uint32_t encodedNucs{0};
  uint32_t bitsPerCode{2};
  std::vector<uint32_t> charsToExtend;
  size_t i = 0;
  for (; i < distToSamplePos; ++i) {
    if ( firstSampPos + i >= seqVec.size() ) {
      std::cerr << "seqVec.size() " << seqVec.size() << ", looking for index " << firstSampPos + i << "\n";
      std::cerr << "dist to sample is " << distToSamplePos << ", i = " << i << "\n";
    }
    auto c = seqVec[firstSampPos + i];
    charsToExtend.push_back(c);
  }
  if (dir == NextSampleDirection::REVERSE) {
    std::reverse(charsToExtend.begin(), charsToExtend.end());
  }

  for (size_t j = 0; j < charsToExtend.size(); ++j) {
    auto c = charsToExtend[j];
    encodedNucs |= (c << (bitsPerCode  * (maxExt - j - 1)));
  }
  return encodedNucs;
}

int pufferfishIndex(util::IndexOptions& indexOpts) {
  uint32_t k = indexOpts.k;
  std::string gfa_file = indexOpts.gfa_file;
  std::string rfile = indexOpts.rfile;
  std::string outdir = indexOpts.outdir;

  // If the user included the '/' in the output directory path, remove
  // it here
  if (outdir.back() == '/') {
    outdir.pop_back();
  }
  std::vector<std::string> read_file = {rfile};

  auto console = spdlog::stderr_color_mt("console");

  size_t tlen{0};
  size_t numKmers{0};
  size_t nread{0};
  CanonicalKmer::k(k);

  puffer::fs::MakeDir(outdir.c_str());

  PosFinder pf(gfa_file.c_str(), k - 1);
  pf.parseFile();
  // std::exit(1);
  pf.mapContig2Pos();
  pf.serializeContigTable(outdir);
  pf.clearContigTable();

  {
    auto& cnmap = pf.getContigNameMap();
    for (auto& kv : cnmap) {
      auto& r1 = kv.second;
      tlen += r1.length;
      numKmers += r1.length - k + 1;
      ++nread;
    }
    console->info("# segments = {}", nread);
    console->info("total length = {}", tlen);
  }

  // now we know the size we need --- create our bitvectors and pack!
  size_t w = std::log2(tlen) + 1;
  console->info("positional integer width = {}", w);
  // sdsl::int_vector<> seqVec(tlen, 0, 2);
  auto& seqVec = pf.getContigSeqVec();
  sdsl::bit_vector rankVec(tlen);
  auto& cnmap = pf.getContigNameMap();
  for (auto& kv : cnmap) {
    rankVec[kv.second.offset + kv.second.length - 1] = 1;
  }
  size_t nkeys{numKmers};
  size_t numContigs{cnmap.size()};

  std::cerr << "seqSize = " << sdsl::size_in_mega_bytes(seqVec) << "\n";
  std::cerr << "rankSize = " << sdsl::size_in_mega_bytes(rankVec) << "\n";
  // std::cerr << "posSize = " << sdsl::size_in_mega_bytes(posVec) << "\n";
  std::cerr << "num keys = " << nkeys << "\n";
  ContigKmerIterator kb(&seqVec, &rankVec, k, 0);
  ContigKmerIterator ke(&seqVec, &rankVec, k, seqVec.size() - k + 1);

#ifdef PUFFER_DEBUG
  auto ks = kb;
  size_t nkeyIt{0};
  for (; ks < ke; ++ks) {
    nkeyIt++;
  }
  std::cerr << "num keys (iterator)= " << nkeyIt << "\n";
#endif // PUFFER_DEBUG
 
  typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
  typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

  auto keyIt = boomphf::range(kb, ke);
  boophf_t* bphf =
      new boophf_t(nkeys, keyIt, 16, 3.5); // keys.size(), keys, 16);
  std::cerr << "mphf size = " << (bphf->totalBitSize() / 8) / std::pow(2, 20)
            << "\n";

  sdsl::store_to_file(seqVec, outdir + "/seq.bin");
  sdsl::store_to_file(rankVec, outdir + "/rank.bin");

  // size_t slen = seqVec.size();
  //#ifndef PUFFER_DEBUG
  // seqVec.resize(0);
  // rankVec.resize(0);
  //#endif

  if (!indexOpts.isSparse) {
    sdsl::int_vector<> posVec(nkeys, 0, w);
    {
      size_t i = 0;
      ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
      ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
      for (; kb1 != ke1; ++kb1) {
        auto idx = bphf->lookup(*kb1); // fkm.word(0));
        if (idx >= posVec.size()) {
          std::cerr << "i =  " << i << ", size = " << seqVec.size()
                    << ", idx = " << idx << ", size = " << posVec.size() << "\n";
        }
        posVec[idx] = kb1.pos();

        // validate
#ifdef PUFFER_DEBUG
        uint64_t kn = seqVec.get_int(2 * kb1.pos(), 2 * k);
        CanonicalKmer sk;
        sk.fromNum(kn);
        if (sk.isEquivalent(*kb1) == KmerMatchType::NO_MATCH) {
          my_mer r;
          r.word__(0) = *kb1;
          std::cerr << "I thought I saw " << sk.to_str() << ", but I saw "
                    << r.to_str() << "\n";
        }
#endif
      }
    }

    /** Write the index **/
    std::ofstream descStream(outdir + "/info.json");
    {
      cereal::JSONOutputArchive indexDesc(descStream);
      std::string sampStr = "dense";
      indexDesc(cereal::make_nvp("sampling_type", sampStr));
      indexDesc(cereal::make_nvp("k", k));
      indexDesc(cereal::make_nvp("num_kmers", nkeys));
      indexDesc(cereal::make_nvp("num_contigs", numContigs));
      indexDesc(cereal::make_nvp("seq_length", tlen));
    }
    descStream.close();

    sdsl::store_to_file(posVec, outdir + "/pos.bin");
    std::ofstream hstream(outdir + "/mphf.bin");
    bphf->save(hstream);
    hstream.close();



  } else { // sparse index; it's GO time!
    int extensionSize = indexOpts.extensionSize;
    int sampleSize = 2 * extensionSize + 1;
    sdsl::bit_vector presenceVec(nkeys);
    size_t sampledKmers{0};
    std::vector<size_t> sampledInds;
    std::vector<size_t> contigLengths;
    //fill up optimal positions
    {
      auto& cnmap = pf.getContigNameMap() ;
      std::vector<size_t> sampledInds ;
      for(auto& kv : cnmap){
        auto& r1 = kv.second ;
        sampledInds.clear();
        computeSampledPositions(r1.length, k, sampleSize, sampledInds) ;
        sampledKmers += sampledInds.size() ;
        contigLengths.push_back(r1.length) ;
      }
      console->info("# sampled kmers = {}", sampledKmers) ;
      console->info("# skipped kmers = {}", numKmers - sampledKmers) ;
   }

    //fill up the vectors
    uint32_t extSymbolWidth = 2;
    uint32_t extWidth = std::log2(extensionSize);
    std::cerr << "extWidth = " << extWidth << "\n";
    sdsl::int_vector<> auxInfo((numKmers-sampledKmers), 0, extSymbolWidth*extensionSize) ;
    sdsl::int_vector<> extSize((numKmers-sampledKmers), 0, extWidth) ;
    //extSize[idx - rank] = extensionDist;
    sdsl::bit_vector direction(numKmers - sampledKmers) ;
    sdsl::bit_vector canonicalNess(numKmers - sampledKmers);
    sdsl::int_vector<> samplePosVec(sampledKmers, 0, w);


  // new presence Vec
  {
    std::cerr << "\nFilling presence Vector \n" ;

    size_t i = 0 ;
    ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
    ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
    size_t contigId{0};
    int sampleCounter = 0 ;

    //debug flags
    int loopCounter = 0;
    size_t ourKeys = 0 ;
    while(kb1 != ke1){
        sampledInds.clear();
        auto clen = contigLengths[contigId];
        computeSampledPositions(clen, k, sampleSize, sampledInds) ;
        contigId++;
        loopCounter++ ;

        my_mer r;
        auto zeroPos = kb1.pos();
        auto nextSampIter = sampledInds.begin();
        auto prevSamp = *nextSampIter;
        auto skipLen = kb1.pos() - zeroPos;
        bool didSample = false;
        bool done = false;

        for (size_t j = 0; j < clen - k + 1; ++kb1, ++j) {
          skipLen = kb1.pos() - zeroPos;
          if (!done and skipLen == *nextSampIter) {
            auto idx = bphf->lookup(*kb1);
            presenceVec[idx] = 1 ;
            i++ ;
            didSample = true;
            prevSamp = *nextSampIter;
            ++nextSampIter;
            if (nextSampIter == sampledInds.end()) {
              done = true;
            }
          }
          didSample = false;
        }
        if (nextSampIter != sampledInds.end()) {
          std::cerr << "I didn't sample " << std::distance(nextSampIter, sampledInds.end()) << " samples for contig " << contigId - 1 << "\n";
          std::cerr << "last sample is " << sampledInds.back() << "\n";
          std::cerr << "contig length is " << contigLengths[contigId-1] << "\n";
        }
    }

    std::cerr << " i = " << i
              << " sampledKmers = " << sampledKmers
              << " Loops = "<< loopCounter
              << " Contig array = "<<contigLengths.size()
              << "\n" ;

  }

  sdsl::bit_vector::rank_1_type realPresenceRank(&presenceVec) ;
  std::cerr << " num ones in presenceVec = " << realPresenceRank(presenceVec.size()-1) << "\n" ;

  //bidirectional sampling
  {

    ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
    ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);

    size_t contigId{0} ;
    size_t coveredKeys{0} ;
    size_t totalKmersIshouldSee{0} ;

    // For every valid k-mer (i.e. every contig)
    while(kb1 != ke1){
      sampledInds.clear();
      auto clen = contigLengths[contigId];
      auto thisContigLength = clen;
      computeSampledPositions(clen, k, sampleSize, sampledInds) ;
      totalKmersIshouldSee += (thisContigLength - k + 1);

      contigId++ ;

      size_t skip = 0 ;

      my_mer r;

      auto zeroPos = kb1.pos();
      auto nextSampIter = sampledInds.begin();
      auto prevSampIter = sampledInds.end();
      auto skipLen = kb1.pos() - zeroPos;
      NextSampleDirection sampDir = NextSampleDirection::FORWARD;
      bool done = false;
      for (size_t j = 0; j < clen - k + 1; ++kb1, ++j) {
          int64_t nextSampPos = (nextSampIter != sampledInds.end()) ? *nextSampIter : -1;
          int64_t prevSampPos = (prevSampIter != sampledInds.end()) ? *prevSampIter : -1;
          uint64_t distToNext = (nextSampPos >= 0) ? nextSampPos - j : std::numeric_limits<uint64_t>::max();
          uint64_t distToPrev = (prevSampPos >= 0) ? j - prevSampPos : std::numeric_limits<uint64_t>::max();

          if (distToNext == std::numeric_limits<uint64_t>::max() and
              distToPrev == std::numeric_limits<uint64_t>::max()) {
            std::cerr << "Could not find valid sample position, should not happen!\n";
            std::exit(1);
          }

          sampDir = (distToNext < distToPrev) ? NextSampleDirection::FORWARD : NextSampleDirection::REVERSE;
          skipLen = kb1.pos() - zeroPos;
          // If this is a sampled position
          if (!done and skipLen == *nextSampIter) {
            prevSampIter = nextSampIter;
            ++nextSampIter;
            if (nextSampIter == sampledInds.end()) {
              done = true;
            }
            auto idx = bphf->lookup(*kb1);
            auto rank = (idx == 0) ? 0 : realPresenceRank(idx);
            samplePosVec[rank] = kb1.pos();
          } else { // not a sampled position
            uint32_t ext = 0;
            size_t firstSampPos = 0;
            uint32_t extensionDist = 0;
            if (sampDir == NextSampleDirection::FORWARD) {
              firstSampPos = zeroPos + j + k;
              extensionDist = distToNext - 1;
              ext = getEncodedExtension(seqVec, firstSampPos, distToNext, extensionSize, sampDir);
            } else if (sampDir == NextSampleDirection::REVERSE) {
              firstSampPos = zeroPos + prevSampPos;
              extensionDist = distToPrev - 1;
              ext = getEncodedExtension(seqVec, firstSampPos, distToPrev, extensionSize, sampDir);
            } else {
              std::cerr << "Error during extension encoding, should not happen!\n";
              std::exit(1);
            }
            auto idx = bphf->lookup(*kb1);
            auto rank = (idx == 0) ? 0 : realPresenceRank(idx);
            canonicalNess[idx - rank] = kb1.isCanonical();
            extSize[idx - rank] = extensionDist;
            auxInfo[idx - rank] = ext;
            direction[idx - rank] = (sampDir == NextSampleDirection::FORWARD) ? 1 : 0;
          }
        }
    }


  }


  /** Write the index **/
  std::ofstream descStream(outdir + "/info.json");
  {
    cereal::JSONOutputArchive indexDesc(descStream);
    std::string sampStr = "sparse";
    indexDesc(cereal::make_nvp("sampling_type", sampStr));
    indexDesc(cereal::make_nvp("sample_size", sampleSize));
    indexDesc(cereal::make_nvp("extension_size", extensionSize));
    indexDesc(cereal::make_nvp("k", k));
    indexDesc(cereal::make_nvp("num_kmers", nkeys));
    indexDesc(cereal::make_nvp("num_sampled_kmers",sampledKmers));
    indexDesc(cereal::make_nvp("num_contigs", numContigs));
    indexDesc(cereal::make_nvp("seq_length", tlen));
  }
  descStream.close();

  //sdsl::store_to_file(posVec, outdir + "/pos.bin");
  std::ofstream hstream(outdir + "/mphf.bin");
  sdsl::store_to_file(presenceVec, outdir + "/presence.bin");
  sdsl::store_to_file(samplePosVec, outdir + "/sample_pos.bin");
  sdsl::store_to_file(auxInfo, outdir + "/extension.bin");
  sdsl::store_to_file(extSize, outdir + "/extensionSize.bin");
  sdsl::store_to_file(canonicalNess, outdir + "/canonical.bin");
  sdsl::store_to_file(direction, outdir + "/direction.bin");
  bphf->save(hstream);
  hstream.close();

  } 
  return 0;
}
