#include "FastxParser.hpp"
#include <cmath>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>

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
  size_t gpos{0};
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
  /*
  tlen = 0;

  size_t numContigs{0};
  size_t nkeys{0};
  {
    auto& cnmap = pf.getContigNameMap();
    for (auto& kv : cnmap) {
      ++numContigs;
      size_t len{0};
      const auto& r1 = kv.second;
      CanonicalKmer mer;
      for (size_t i = 0; i < r1.length(); ++i) {
        auto offset = i; // r1.length() - i - 1;
        // NOTE: Having to add things in the reverse order here is strange
        // we should make sure that this doesn't mess with the positions we
        // end up storing!
        auto c = my_mer::code(r1[i]);
        seqVec[gpos + offset] = c;
        mer.shiftFw(c);
        ++len;
        if (len >= k) {
#ifdef PUFFER_DEBUG
          my_mer mm;
          uint64_t num = seqVec.get_int(2 * (gpos + len - k), 2 * k);
          mm.word__(0) = num; // mm.canonicalize();
          if (!(mm == mer.fwMer() or
                mm == mer.rcMer())) { //}mm != mer.get_canonical()) {
            std::cerr << "num & 0x3 = " << (num & 0x3) << "\n";
            std::cerr << "i = " << i << "\n";
            std::cerr << mer.to_str() << ", " << mm.to_str() << "\n";
          }
#endif // PUFFER_DEBUG
          ++nkeys;
        }
        //++gpos;
      }
      gpos += r1.length();
      tlen += r1.length();
      rankVec[tlen - 1] = 1;
    }
  }
  */
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

  sdsl::int_vector<> posVec(tlen, 0, w);
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

  return 0;
}
