#include "FastxParser.hpp"
#include <cmath>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>
#include <sstream>
#include <bitset>
#include <cerrno>
#include <cstring>
#include <cereal/archives/binary.hpp>

#include "ProgOpts.hpp"
#include "CanonicalKmer.hpp"
#include "PufferfishBinaryGFAReader.cpp"
#include "PufferFS.hpp"
#include "PufferfishIndex.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "PufferfishConfig.hpp"
#include "cereal/archives/json.hpp"
#include "jellyfish/mer_dna.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"
#include "spdlog/spdlog.h"
#include "Kmer.hpp" // currently requires k <= 32
#include "compact_vector/compact_vector.hpp"

namespace kmers = combinelib::kmers;

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
  using self_type = ContigKmerIterator;
  using value_type = uint64_t;
  using reference = value_type&;
  using pointer = value_type*;
  using iterator_category = std::forward_iterator_tag;
  using difference_type = int64_t;

  ContigKmerIterator() = delete;

  ContigKmerIterator(compact::vector<uint64_t, 2>* storage, compact::vector<uint64_t, 1>* rank,
                     uint8_t k, uint64_t startAt)
      : storage_(storage), rank_(rank), k_(k), curr_(startAt) {
    if (curr_ + k_ <= rank_->size()) {
      //nextValidPosition_();
      mer_.fromNum(storage_->get_int(2*curr_, 2*k_));
      // mer_.word__(0) = storage_->get_int(2 * curr_, 2 * k_);
    }
    // rcMer_ = mer_.get_reverse_complement();
  }

  bool advanceToValid() {
    nextValidPosition_();
    return (curr_ + k_ <= rank_->size());
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

  inline bool isValid() {
    size_t endPos = curr_ + k_ - 1;
    return !(endPos + 1 >= rank_->size() or (*rank_)[endPos] == 1);
  }

  inline bool isInvalid() {
    return !isValid();
  }

private:

  void nextValidPosition_() {
    size_t endPos = curr_ + k_ - 1;
    if (endPos < rank_->size()) {
      // See if we cross a rank boundary
      bool crossesBoundary{false};
      bool isNextValid{false};
      size_t boundaryIndex{0};
      for (size_t i = curr_; i < endPos; ++i) {
        if ((*rank_)[i] == 1) {
          crossesBoundary = true;
          boundaryIndex = i;
          isNextValid = (i + k_) < rank_->size();
          break;
        }
      }

      // If so, that's the start of the next valid k-mer
      // if that position is valid
      if (crossesBoundary) {
        if (isNextValid) {
          curr_ = boundaryIndex + 1;
        } else {
          // if that position is invalid, then go to the end.
          goto endPos;
        }
      }
      // At this point, either curr_ points to the next valid
      // start position, or we have skipped over to the endPos label.
      mer_.fromNum(storage_->get_int(2*curr_, 2*k_));
      return;
    }

  endPos:
    // Fallthrough if we couldn't find a valid position.
    mer_.fromNum(storage_->get_int(2*(rank_->size() - k_), 2*k_));
    curr_ = storage_->size() - k_ + 1;
  }

  void advance_() {
    size_t endPos = curr_ + k_ - 1;
    if (endPos + 1 < rank_->size() and (*rank_)[endPos] == 1) {
      curr_ += k_;
      mer_.fromNum(storage_->get_int(2*curr_, 2*k_));
    } else {
      if (curr_ + k_ < rank_->size()) {
        int c = (*storage_)[curr_ + k_];
        mer_.shiftFw(c);
        ++curr_;
      } else {
        mer_.fromNum(storage_->get_int(2*(rank_->size() - k_), 2*k_));
        curr_ = storage_->size() - k_ + 1;
        //curr_ = rank_->size();
      }
    }
  }
  compact::vector<uint64_t, 2>* storage_{nullptr};
  compact::vector<uint64_t, 1>* rank_{nullptr};
  uint8_t k_{0};
  uint64_t curr_{0};
  CanonicalKmer mer_;
  uint64_t word_{0};
};

int pufferfishTest(TestOptions& testOpts) {
  (void)testOpts;
  std::cerr << "this command is not yet implemented\n";
  return 1;
}

void computeSampledPositionsLossy(size_t tlen, uint32_t k, int32_t sampleSize, std::vector<size_t>& sampledInds){
  // Let's start out simple, we always keep the first & last k-mers in a unipath.
  // If the unipath is longer than sampleSize, we keep intermediate samples as well.
  sampledInds.clear();
  auto numOfKmers = tlen - k;
  size_t lastSampled = 0;
  while (lastSampled < numOfKmers) {
    sampledInds.push_back(lastSampled) ;
    auto next_samp = std::min(lastSampled + sampleSize, numOfKmers) ;
    lastSampled = next_samp;
  }
  if (lastSampled ==  numOfKmers) {
    sampledInds.push_back(numOfKmers);
  }
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

std::string packedToString(compact::vector<uint64_t, 2>& seqVec, uint64_t offset, uint32_t len) {
  std::stringstream s;
  for (size_t i = offset; i < offset + len; ++i) {
    auto c = seqVec[i];
    s << kmers::charForCode(c);
  }
  auto st = s.str();
  return st;
}

enum class NextSampleDirection : uint8_t { FORWARD = 0, REVERSE=1 };

uint32_t getEncodedExtension(compact::vector<uint64_t, 2>& seqVec, uint64_t firstSampPos, uint64_t distToSamplePos,
                             uint32_t maxExt, NextSampleDirection dir) {
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


int buildGraphMain(std::vector<std::string>& args);
int dumpGraphMain(std::vector<std::string>& args);


int pufferfishIndex(IndexOptions& indexOpts) {
  uint32_t k = indexOpts.k;
  std::string gfa_file = indexOpts.gfa_file;
  std::string rfile = indexOpts.rfile;
  std::string outdir = indexOpts.outdir;
  bool buildEdgeVec = indexOpts.buildEdgeVec;

  // If the user included the '/' in the output directory path, remove
  // it here
  if (outdir.back() == '/') {
    outdir.pop_back();
  }

  auto console = spdlog::stderr_color_mt("console");

  size_t tlen{0};
  size_t numKmers{0};
  size_t nread{0};
  CanonicalKmer::k(k);

  if (puffer::fs::MakePath(outdir.c_str()) != 0) {
    console->error(std::strerror(errno));
    std::exit(1);
  }

  {
    std::vector<std::string> args;
    args.push_back("twopaco");
    args.push_back("-k");
    args.push_back(std::to_string(k));
    args.push_back("-t");
    args.push_back(std::to_string(indexOpts.p));
    args.push_back("-f");
    args.push_back("32");
    args.push_back("--outfile");
    args.push_back(outdir+"/tmp_dbg.bin");
    args.push_back("--tmpdir");
    args.push_back(outdir);
    args.push_back(rfile);
    buildGraphMain(args);
    //std::exit(0);
  }

  {
    std::vector<std::string> args;
    args.push_back("graphdump");
    args.push_back("-k");
    args.push_back(std::to_string(k));
    args.push_back("-s");
    args.push_back(rfile);
    args.push_back("-f");
    args.push_back("binPufferized");
    args.push_back(outdir+"/tmp_dbg.bin");
    args.push_back("-p");
    args.push_back(outdir);
    dumpGraphMain(args);
  }

  pufferfish::BinaryGFAReader pf(gfa_file.c_str(), k - 1, buildEdgeVec, console);
  pf.parseFile();
  pf.mapContig2Pos();
  pf.serializeContigTable(outdir);
  {
    auto& cnmap = pf.getContigNameMap();
    for (auto& kv : cnmap) {
      auto& r1 = kv.second;
      tlen += r1.length;
      numKmers += r1.length - k + 1;
      ++nread;
    }
    console->info("# segments = {:n}", nread);
    console->info("total length = {:n}", tlen);
  }

  // parse the reference list and store the strings in a 2bit-encoded vector
  bool keepRef = (rfile.size() > 0);
  if (keepRef) {
    auto &refIds = pf.getRefIDs();
    auto &refLengths = pf.getRefLengths();
    std::vector<uint64_t> refAccumLengths(refLengths.size());
    std::unordered_map<std::string, uint64_t> refIdMap;
    uint64_t prev = 0;
    for (uint64_t i = 0; i < refIds.size(); i++) {
      if (refIdMap.find(refIds[i]) != refIdMap.end()) {
        console->error("Two references with the same name but different sequences: {}. "
                       "We require that all input records have a unique name "
                       "up to the first whitespace character.", refIds[i]);
        std::exit(1);
      }
      refIdMap[refIds[i]] = i;
      refAccumLengths[i] = prev + refLengths[i];
      prev = refAccumLengths[i];
//      std::cerr << i << ":" << refLengths[i] << "\n";
    }
    //compact 2bit vector
//    std::cerr << "\nrefAccumLengths.size():" << refAccumLengths.size() << "\n";
//    std::cerr << refAccumLengths.back() << "\n";
    compact::vector<uint64_t, 2> refseq(refAccumLengths.back());
    // go over all the reference files
    console->info("Reading the reference files ...");
    std::vector<std::string> ref_files = {rfile};
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(ref_files, 1, 1);
    parser.start();
    auto rg = parser.getReadGroup();
    // read the reference sequences and encode them into the refseq int_vector

    while (parser.refill(rg)) {
      for (auto &rp : rg) {
        stx::string_view seqv(rp.seq);
        auto &refIdx = refIdMap[rp.name];
        auto offset = refIdx == 0 ? 0 : refAccumLengths[refIdx - 1];
        pf.encodeSeq(refseq, offset, seqv);
      }
    }

    parser.stop();
    // store the 2bit-encoded references
    // store reference accumulative lengths
    std::string accumLengthsFilename = outdir + "/refAccumLengths.bin";
    std::ofstream ral(accumLengthsFilename);
    cereal::BinaryOutputArchive ralAr(ral);
    ralAr(refAccumLengths);

    // store reference sequences
    std::ofstream seqFile(outdir + "/refseq.bin", std::ios::binary);
    refseq.serialize(seqFile);
    seqFile.close();

  }
  pf.clearContigTable();

  // now we know the size we need --- create our bitvectors and pack!
  size_t w = std::log2(tlen) + 1;
  console->info("positional integer width = {:n}", w);

  auto& seqVec = pf.getContigSeqVec();
  auto& rankVec = pf.getRankVec();
  auto& edgeVec = pf.getEdgeVec() ;

 /* compact::vector<uint64_t, 1> rankVec(tlen);
  for(uint64_t i=0; i<tlen; i++) rankVec[i]=0;//if(rankVec[i]!=0) {std::cerr<<"Not zero\n"; break;}
  auto& cnmap = pf.getContigNameMap();
  for (auto& kv : cnmap) {
    rankVec[kv.second.offset + kv.second.length - 1] = 1;
  }*/
  size_t nkeys{numKmers};
  size_t numContigs{pf.getContigNameMap().size()};

  console->info("seqSize = {:n}", seqVec.size());
  console->info("rankSize = {:n}", rankVec.size());

  console->info("edgeVecSize = {:n}", edgeVec.size());

  console->info("num keys = {:n}", nkeys);
  ContigKmerIterator kb(&seqVec, &rankVec, k, 0);
  ContigKmerIterator ke(&seqVec, &rankVec, k, seqVec.size() - k + 1);

#ifdef PUFFER_DEBUG
  auto ks = kb;
  size_t nkeyIt{0};
  for (; ks < ke; ++ks) {
    nkeyIt++;
  }
  console->info("num keys (iterator)= {:n}", nkeyIt);
#endif // PUFFER_DEBUG
 
  typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
  typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

  auto keyIt = boomphf::range(kb, ke);
  boophf_t* bphf =
      new boophf_t(nkeys, keyIt, indexOpts.p, 3.5); // keys.size(), keys, 16);
  console->info("mphf size = {} MB", (bphf->totalBitSize() / 8) / std::pow(2, 20));

/*  std::ofstream seqFile(outdir + "/seq.bin", std::ios::binary);
  seqVec.serialize(seqFile);
  seqFile.close();

  std::ofstream rankFile(outdir + "/rank.bin", std::ios::binary);
  rankVec.serialize(rankFile);
  rankFile.close();*/

  bool haveEdgeVec = edgeVec.size() > 0;
  if (haveEdgeVec) {
    std::ofstream edgeFile(outdir + "/edge.bin", std::ios::binary);
    edgeVec.serialize(edgeFile);
    edgeFile.close();
  }

  // if using quasi-dictionary idea (https://arxiv.org/pdf/1703.00667.pdf)
  //uint32_t hashBits = 4;
  if (!indexOpts.isSparse and !indexOpts.lossySampling) {  
    // if using quasi-dictionary idea (https://arxiv.org/pdf/1703.00667.pdf)
    compact::ts_vector<uint64_t> posVec(w, nkeys);
    {

      struct ContigVecChunk {
        uint64_t s;
        uint64_t e;
      };

      // Build position table in parallel. We have up to indexOpts.p threads
      // so divide the contig array into even chunks.
      auto nthread = indexOpts.p;
      double chunkSizeFrac = seqVec.size() / static_cast<double>(nthread);
      // reduce the number of threads until chunks are big enough;
      while (chunkSizeFrac < 8192 and nthread > 1) {
        nthread /= 2;
        chunkSizeFrac = seqVec.size() / static_cast<double>(nthread);
      }

      auto chunkSize = static_cast<uint64_t>(std::ceil(chunkSizeFrac));
      console->info("chunk size = {:n}", chunkSize);

      std::vector<ContigVecChunk> chunks;
      chunks.reserve(nthread);
      size_t itOffset{0};
      for (size_t i = 0; i < nthread; ++i) {
        ContigKmerIterator startIt(&seqVec, &rankVec, k, itOffset);
        startIt.advanceToValid();
        uint64_t s = startIt.pos();
        if (i == nthread - 1) {
          itOffset = seqVec.size() - k + 1;
        } else {
          itOffset = s + chunkSize;
        }
        ContigKmerIterator endIt(&seqVec, &rankVec, k, itOffset);
        endIt.advanceToValid();
        uint64_t e = endIt.pos();
        chunks.push_back({s,e});
        console->info("chunk {} = [{:n}, {:n})", i, s, e);
      }

      auto fillPos = [&seqVec, &rankVec, k, &bphf, &console, &posVec](ContigVecChunk chunk) -> void {

        ContigKmerIterator kb1(&seqVec, &rankVec, k, chunk.s);
        ContigKmerIterator ke1(&seqVec, &rankVec, k, chunk.e);
        for (; kb1 < ke1; ++kb1) {
          auto idx = bphf->lookup(*kb1); // fkm.word(0));
          if (idx >= posVec.size()) {
            std::cerr<<*kb1<<"\n";
            console->info("seq size = {:n}, idx = {:n}, pos size = {:n}",
                          seqVec.size(), idx, posVec.size());
            std::cerr<<*kb1<<"\n";
          }
          // ContigKmerIterator::value_type mer = *kb1;
          // if using quasi-dictionary idea (https://arxiv.org/pdf/1703.00667.pdf)
          //posVec[idx] = (kb1.pos() << hashBits) | (mer & 0xF);
          posVec[idx] = kb1.pos();
          // validate
#ifdef PUFFER_DEBUG
          uint64_t kn = seqVec.get_int(2*kb1.pos(), 2*k);
          CanonicalKmer sk;
          sk.fromNum(kn);
          if (sk.isEquivalent(*kb1) == KmerMatchType::NO_MATCH) {
            my_mer r;
            r.word__(0) = *kb1;
            console->error("I thought I saw {}, but I saw {} --- pos {}", sk.to_str(), r.toStr(), kb1.pos());
          }
#endif
        }
      };

      std::vector<std::thread> workers;
      workers.reserve(chunks.size());
      for (auto chunk : chunks) {
        workers.push_back(std::thread(fillPos, chunk));
      }
      for (auto& w : workers) {
        w.join();
      }
      console->info("finished populating pos vector");

    }

    console->info("writing index components");
    /** Write the index **/
    std::ofstream descStream(outdir + "/info.json");
    {
      cereal::JSONOutputArchive indexDesc(descStream);
      std::string sampStr = "dense";
      std::vector<std::string> refGFA{gfa_file};
      indexDesc(cereal::make_nvp("index_version", pufferfish::indexVersion));
      indexDesc(cereal::make_nvp("reference_gfa", refGFA));
      indexDesc(cereal::make_nvp("sampling_type", sampStr));
      indexDesc(cereal::make_nvp("k", k));
      indexDesc(cereal::make_nvp("num_kmers", nkeys));
      indexDesc(cereal::make_nvp("num_contigs", numContigs));
      indexDesc(cereal::make_nvp("seq_length", tlen));
      indexDesc(cereal::make_nvp("have_ref_seq", keepRef));
      indexDesc(cereal::make_nvp("have_edge_vec", haveEdgeVec));
    }
    descStream.close();

    std::ofstream posFile(outdir + "/pos.bin", std::ios::binary);
    posVec.serialize(posFile);
    posFile.close();
    std::ofstream hstream(outdir + "/mphf.bin");
    bphf->save(hstream);
    hstream.close();
    console->info("finished writing dense pufferfish index");

  } else if (indexOpts.isSparse) { // sparse index; it's GO time!
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
      console->info("# sampled kmers = {:n}", sampledKmers) ;
      console->info("# skipped kmers = {:n}", numKmers - sampledKmers) ;
    }

    //fill up the vectors
    uint32_t extSymbolWidth = 2;
    uint32_t extWidth = std::log2(extensionSize);
    console->info("extWidth = {}", extWidth);
    sdsl::int_vector<> auxInfo((numKmers-sampledKmers), 0, extSymbolWidth*extensionSize) ;
    sdsl::int_vector<> extSize((numKmers-sampledKmers), 0, extWidth) ;
    //extSize[idx - rank] = extensionDist;
    sdsl::bit_vector direction(numKmers - sampledKmers) ;
    sdsl::bit_vector canonicalNess(numKmers - sampledKmers);
    sdsl::int_vector<> samplePosVec(sampledKmers, 0, w);


  // new presence Vec
  {
    console->info("\nFilling presence Vector");

    size_t i = 0 ;
    ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
    ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
    size_t contigId{0};

    //debug flags
    int loopCounter = 0;
    //size_t ourKeys = 0 ;
    while(kb1 != ke1){
        sampledInds.clear();
        auto clen = contigLengths[contigId];
        computeSampledPositions(clen, k, sampleSize, sampledInds) ;
        contigId++;
        loopCounter++ ;

        my_mer r;
        auto zeroPos = kb1.pos();
        auto skipLen = kb1.pos() - zeroPos;
        auto nextSampIter = sampledInds.begin();
        //auto prevSamp = *nextSampIter;
        //bool didSample = false;
        bool done = false;

        for (size_t j = 0; j < clen - k + 1; ++kb1, ++j) {
          skipLen = kb1.pos() - zeroPos;
          if (!done and skipLen == static_cast<decltype(skipLen)>(*nextSampIter)) {
            auto idx = bphf->lookup(*kb1);
            presenceVec[idx] = 1 ;
            i++ ;
            //didSample = true;
            //prevSamp = *nextSampIter;
            ++nextSampIter;
            if (nextSampIter == sampledInds.end()) {
              done = true;
            }
          }
          //didSample = false;
        }
        if (nextSampIter != sampledInds.end()) {
          console->info("I didn't sample {}, samples for contig {}", std::distance(nextSampIter, sampledInds.end()), contigId - 1);
          console->info("last sample is " , sampledInds.back());
          console->info("contig length is " , contigLengths[contigId-1]);
        }
    }

    console->info("i = {:n}, sampled kmers = {:n}, loops = {:n}, contig array = {:n}",
                  i, sampledKmers, loopCounter, contigLengths.size());
  }

  sdsl::bit_vector::rank_1_type realPresenceRank(&presenceVec) ;
  console->info("num ones in presenceVec = {:n}", realPresenceRank(presenceVec.size()-1));

  //bidirectional sampling
  {

    ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
    ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);

    size_t contigId{0} ;
    //size_t coveredKeys{0} ;
    size_t totalKmersIshouldSee{0} ;

    // For every valid k-mer (i.e. every contig)
    while(kb1 != ke1){
      sampledInds.clear();
      auto clen = contigLengths[contigId];
      auto thisContigLength = clen;
      computeSampledPositions(clen, k, sampleSize, sampledInds) ;
      totalKmersIshouldSee += (thisContigLength - k + 1);

      contigId++ ;

      //size_t skip = 0 ;

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
            console->error("Could not find valid sample position, should not happen!");
            std::exit(1);
          }

          sampDir = (distToNext < distToPrev) ? NextSampleDirection::FORWARD : NextSampleDirection::REVERSE;
          skipLen = kb1.pos() - zeroPos;
          // If this is a sampled position
          if (!done and skipLen == static_cast<decltype(skipLen)>(*nextSampIter)) {
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
    std::vector<std::string> refGFA{gfa_file};
    indexDesc(cereal::make_nvp("index_version", pufferfish::indexVersion));
    indexDesc(cereal::make_nvp("reference_gfa", refGFA));
    indexDesc(cereal::make_nvp("sampling_type", sampStr));
    indexDesc(cereal::make_nvp("sample_size", sampleSize));
    indexDesc(cereal::make_nvp("extension_size", extensionSize));
    indexDesc(cereal::make_nvp("k", k));
    indexDesc(cereal::make_nvp("num_kmers", nkeys));
    indexDesc(cereal::make_nvp("num_sampled_kmers",sampledKmers));
    indexDesc(cereal::make_nvp("num_contigs", numContigs));
    indexDesc(cereal::make_nvp("seq_length", tlen));
    indexDesc(cereal::make_nvp("have_ref_seq", keepRef));
    indexDesc(cereal::make_nvp("have_edge_vec", haveEdgeVec));
  }
  descStream.close();

  std::ofstream hstream(outdir + "/mphf.bin");
  sdsl::store_to_file(presenceVec, outdir + "/presence.bin");
  sdsl::store_to_file(samplePosVec, outdir + "/sample_pos.bin");
  sdsl::store_to_file(auxInfo, outdir + "/extension.bin");
  sdsl::store_to_file(extSize, outdir + "/extensionSize.bin");
  sdsl::store_to_file(canonicalNess, outdir + "/canonical.bin");
  sdsl::store_to_file(direction, outdir + "/direction.bin");
  bphf->save(hstream);
  hstream.close();

  } else { // lossy sampling index
    int32_t sampleSize = static_cast<int32_t>(indexOpts.lossy_rate);
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
        computeSampledPositionsLossy(r1.length, k, sampleSize, sampledInds) ;
        sampledKmers += sampledInds.size() ;
        contigLengths.push_back(r1.length) ;
      }
      console->info("# sampled kmers = {:n}", sampledKmers) ;
      console->info("# skipped kmers = {:n}", numKmers - sampledKmers) ;
    }

    sdsl::int_vector<> samplePosVec(sampledKmers, 0, w);

    // new presence Vec
    {
      {
      console->info("\nFilling presence vector");
      size_t i = 0 ;
      ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
      ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
      size_t contigId{0};
      int loopCounter = 0;
      while(kb1 != ke1){
        sampledInds.clear();
        auto clen = contigLengths[contigId];
        computeSampledPositionsLossy(clen, k, sampleSize, sampledInds) ;
        contigId++;
        loopCounter++ ;

        my_mer r;
        auto zeroPos = kb1.pos();
        auto skipLen = kb1.pos() - zeroPos;
        auto nextSampIter = sampledInds.begin();
        bool done = false;
        
        for (size_t j = 0; j < clen - k + 1; ++kb1, ++j) {
          skipLen = kb1.pos() - zeroPos;
          if (!done and skipLen == static_cast<decltype(skipLen)>(*nextSampIter)) {
            auto idx = bphf->lookup(*kb1);
            presenceVec[idx] = 1 ;
            samplePosVec[i] = kb1.pos();
            i++;
            ++nextSampIter;
            if (nextSampIter == sampledInds.end()) {
              done = true;
            }
          }
        }
        if (nextSampIter != sampledInds.end()) {
          console->info("I didn't sample {:n}, samples for contig {:n}", std::distance(nextSampIter, sampledInds.end()), contigId - 1);
          console->info("last sample is {:n}", sampledInds.back());
          console->info("contig length is {:n}", contigLengths[contigId-1]);
        }
      }
      }

      {
        console->info("\nFilling sampled position vector");
      size_t i = 0 ;
      ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
      ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
      size_t contigId{0};
      int loopCounter = 0;
      sdsl::bit_vector::rank_1_type realPresenceRank(&presenceVec) ;
      while(kb1 != ke1){
        sampledInds.clear();
        auto clen = contigLengths[contigId];
        computeSampledPositionsLossy(clen, k, sampleSize, sampledInds) ;
        contigId++;
        loopCounter++ ;

        my_mer r;
        auto zeroPos = kb1.pos();
        auto skipLen = kb1.pos() - zeroPos;
        auto nextSampIter = sampledInds.begin();
        bool done = false;
        
        for (size_t j = 0; j < clen - k + 1; ++kb1, ++j) {
          skipLen = kb1.pos() - zeroPos;
          if (!done and skipLen == static_cast<decltype(skipLen)>(*nextSampIter)) {
            auto idx = bphf->lookup(*kb1);
            auto rank = (idx == 0) ? 0 : realPresenceRank(idx);
            samplePosVec[rank] = kb1.pos();
            ++i;
            ++nextSampIter;
            if (nextSampIter == sampledInds.end()) {
              done = true;
            }
          }
        }
        if (nextSampIter != sampledInds.end()) {
          console->info("I didn't sample {:n}, samples for contig {}", std::distance(nextSampIter, sampledInds.end()), contigId - 1);
          console->info("last sample is {:n}", sampledInds.back());
          console->info("contig length is {:n}", contigLengths[contigId-1]);
        }
      }
      console->info("i = {:n}, sampled kmers = {:n}, loops = {:n}, contig array = {:n}",
                    i, sampledKmers, loopCounter, contigLengths.size());

    }
    }

    /** Write the index **/
    std::ofstream descStream(outdir + "/info.json");
    {
      cereal::JSONOutputArchive indexDesc(descStream);
      std::string sampStr = "lossy";
      std::vector<std::string> refGFA{gfa_file};
      indexDesc(cereal::make_nvp("index_version", pufferfish::indexVersion));
      indexDesc(cereal::make_nvp("reference_gfa", refGFA));
      indexDesc(cereal::make_nvp("sampling_type", sampStr));
      indexDesc(cereal::make_nvp("sample_size", sampleSize));
      indexDesc(cereal::make_nvp("k", k));
      indexDesc(cereal::make_nvp("num_kmers", nkeys));
      indexDesc(cereal::make_nvp("num_sampled_kmers",sampledKmers));
      indexDesc(cereal::make_nvp("num_contigs", numContigs));
      indexDesc(cereal::make_nvp("seq_length", tlen));
      indexDesc(cereal::make_nvp("have_ref_seq", keepRef));
      indexDesc(cereal::make_nvp("have_edge_vec", haveEdgeVec));
    }
    descStream.close();

    std::ofstream hstream(outdir + "/mphf.bin");
    sdsl::store_to_file(presenceVec, outdir + "/presence.bin");
    sdsl::store_to_file(samplePosVec, outdir + "/sample_pos.bin");
    bphf->save(hstream);
    hstream.close();
  }
  return 0;
}
