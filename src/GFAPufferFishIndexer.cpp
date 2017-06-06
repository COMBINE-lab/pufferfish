#include "FastxParser.hpp"
#include <cmath>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>

#include "cereal/archives/json.hpp"
#include "CanonicalKmer.hpp"
#include "OurGFAReader.hpp"
#include "PufferFS.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "jellyfish/mer_dna.hpp"
#include "popl.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"
#include "spdlog/spdlog.h"
#include "PufferfishIndex.hpp"
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

  ContigKmerIterator(sdsl::int_vector<>* storage, sdsl::bit_vector* rank,
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
                                         //storage, sdsl::bit_vector* rank,
                                         //uint8_t k, uint64_t startAt) :
    storage_ = other.storage_;
    rank_ = other.rank_;
    k_ = other.k_;
    curr_ = other.curr_;
    mer_ = other.mer_;
    // rcMer_ = other.rcMer_;
    word_ = other.word_;
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
                                     //mer_.word(0) : rcMer_.word(0);
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
  sdsl::int_vector<>* storage_{nullptr};
  sdsl::bit_vector* rank_{nullptr};
  uint8_t k_{0};
  uint64_t curr_{0};
  CanonicalKmer mer_;
  uint64_t word_{0};
};


int pufferfishTest(util::TestOptions& testOpts) {
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
  if (outdir.back() == '/') { outdir.pop_back(); }
  std::vector<std::string> read_file = {rfile};

  auto console = spdlog::stderr_color_mt("console");

  size_t tlen{0};
  size_t numKmers{0};
  size_t nread{0};
  CanonicalKmer::k(k);

  puffer::fs::MakeDir(outdir.c_str());

  PosFinder pf(gfa_file.c_str(), k - 1);
  pf.parseFile();
  pf.mapContig2Pos();
  pf.serializeContigTable(outdir);

  {
    auto& cnmap = pf.getContigNameMap();
    for (auto& kv : cnmap) {
      auto& r1 = kv.second;
      tlen += r1.length();
      numKmers += r1.length() - k + 1;
      ++nread;
    }
    console->info("# segments = {}", nread);
    console->info("total length = {}", tlen);
  }

  // now we know the size we need --- create our bitvectors and pack!
  size_t gpos{0};
  size_t w = std::log2(tlen) + 1;
  console->info("positional integer width = {}", w);
  sdsl::int_vector<> seqVec(tlen, 0, 2);
  sdsl::int_vector<> posVec(tlen, 0, w);
  sdsl::bit_vector rankVec(tlen);
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
      uint64_t km;
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
          my_mer mm;
          uint64_t num = seqVec.get_int(2 * (gpos + len - k), 2 * k);
          mm.word__(0) = num; // mm.canonicalize();
          if (!(mm == mer.fwMer() or
                mm == mer.rcMer())) { //}mm != mer.get_canonical()) {
            std::cerr << "num & 0x3 = " << (num & 0x3) << "\n";
            std::cerr << "i = " << i << "\n";
            std::cerr << mer.to_str() << ", " << mm.to_str() << "\n";
          }
          // km = (mer.word(0) < rcMer.word(0)) ? mer.word(0) :
          // rcMer.word(0);//mer.get_canonical().word(0);
          km = mer.getCanonicalWord();
          ++nkeys;
        }
        //++gpos;
      }
      gpos += r1.length();
      tlen += r1.length();
      rankVec[tlen - 1] = 1;
    }
  }
  std::cerr << "seqSize = " << sdsl::size_in_mega_bytes(seqVec) << "\n";
  std::cerr << "rankSize = " << sdsl::size_in_mega_bytes(rankVec) << "\n";
  std::cerr << "posSize = " << sdsl::size_in_mega_bytes(posVec) << "\n";
  std::cerr << "num keys = " << nkeys << "\n";
  typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
  typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

  ContigKmerIterator kb(&seqVec, &rankVec, k, 0);
  ContigKmerIterator ke(&seqVec, &rankVec, k, seqVec.size() - k + 1);
  auto ks = kb;
  size_t nkeyIt{0};
  for (; ks < ke; ++ks) {
    nkeyIt++;
  }
  std::cerr << "num keys (iterator)= " << nkeyIt << "\n";
  auto keyIt = boomphf::range(kb, ke);
  boophf_t* bphf = new boophf_t(nkeys, keyIt, 16); // keys.size(), keys, 16);
  std::cerr << "mphf size = " << (bphf->totalBitSize() / 8) / std::pow(2, 20)
            << "\n";
  size_t kstart{0};
  size_t tlen2{0};
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
    }
  }

  /** Write the index **/
  {
    std::ofstream descStream(outdir + "/info.json");
    cereal::JSONOutputArchive indexDesc(descStream);
    std::string sampStr = "dense";
    indexDesc(cereal::make_nvp("sampling_type", sampStr));
    indexDesc(cereal::make_nvp("k", k));
    indexDesc(cereal::make_nvp("num_kmers", nkeyIt));
    indexDesc(cereal::make_nvp("num_contigs", numContigs));
    indexDesc(cereal::make_nvp("seq_length", tlen));
    descStream.close();
  }

  sdsl::store_to_file(seqVec, outdir + "/seq.bin");
  sdsl::store_to_file(rankVec, outdir + "/rank.bin");
  sdsl::store_to_file(posVec, outdir + "/pos.bin");
  std::ofstream hstream(outdir + "/mphf.bin");
  bphf->save(hstream);
  hstream.close();

   return 0;

  size_t N = nkeys; // keys.size();
  size_t S = seqVec.size();
  size_t found = 0;
  size_t notFound = 0;
  size_t correctPosCntr = 0;
  size_t incorrectPosCntr = 0;
  size_t numTrueTxp = 0;
  {
    std::vector<std::string> rankOrderedContigIDs;

    auto& cnmap = pf.getContigNameMap();
    for (auto& kv : cnmap) {
      rankOrderedContigIDs.push_back(kv.first);
    }
    auto& trRefIDs = pf.getRefIDs();
    ScopedTimer st;
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(read_file, 1, 1);
    parser.start();
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    size_t rn{0};
    size_t kmer_pos{0};
    auto rg = parser.getReadGroup();
    sdsl::bit_vector::rank_1_type realRank(&rankVec);
    sdsl::bit_vector::select_1_type realSelect(&rankVec);
    while (parser.refill(rg)) {
      // Here, rg will contain a chunk of read pairs
      // we can process.
      for (auto& rp : rg) {
        kmer_pos = 0;
        if (rn % 500000 == 0) {
          std::cerr << "rn : " << rn << "\n";
          std::cerr << "found = " << found << ", notFound = " << notFound
                    << "\n";
        }
        ++rn;
        auto& r1 = rp.seq;
        CanonicalKmer mer;
        bool merOK = mer.fromStr(r1); // mer.from_chars(r1.begin());
        if (!merOK) {
          std::cerr << "contig too short!";
          std::exit(1);
        }
        auto km = mer.getCanonicalWord();
        size_t res = bphf->lookup(km);
        uint64_t pos =
            (res < N) ? posVec[res] : std::numeric_limits<uint64_t>::max();
        if (pos <= S - k) {
          uint64_t fk = seqVec.get_int(2 * pos, 2 * k);
          my_mer fkm;
          fkm.word__(0) = fk;
          if (mer.fwWord() == fkm.word(0) or mer.rcWord() == fkm.word(0)) {
            found += 1;
            bool correctPos = false;
            bool foundTxp = false;
            auto rank = realRank(pos);
            for (auto& tr : pf.contig2pos[rankOrderedContigIDs[rank]]) {
              if (trRefIDs[tr.transcript_id] == rp.name) {
                foundTxp = true;
                uint64_t sp = (uint64_t)realSelect(rank);
                auto relPos = pos - sp;
                auto clen = (uint64_t)realSelect(rank + 1) - sp;
                if (relPos + tr.pos == kmer_pos or
                    clen - (relPos - tr.pos - k - 1) == kmer_pos) {
                  correctPos = true;
                  break;
                } else {
                  // std::cerr << "ours : " << relPos + tr.pos << " , true : "
                  // << kmer_pos << "\n";
                }
              }
            }
            if (correctPos) {
              correctPosCntr++;
            } else {
              incorrectPosCntr++;
            }
            if (foundTxp) {
              numTrueTxp++;
            }
            // if (!foundTxp) { std::cerr << "failed to find true txp [" <<
            // rp.name << "]\n"; }
          } else {
            // std::cerr << "rn = " << rn - 1 << ", fk = " <<
            // fkm.get_canonical().to_str() << ", km = " <<
            // mer.get_canonical().to_str() << ", pkmer = " << pkmer <<" \n";
            // std::cerr << "found = " << found << ", not found = " << notFound
            // << "\n";
            notFound += 1;
          }
        } else {
          // std::cerr << "pos = " << pos << ", shouldn't happen\n";
          notFound += 1;
        }

        for (size_t i = k; i < r1.length(); ++i) {
          mer.shiftFw(r1[i]);
          km = mer.getCanonicalWord();
          res = bphf->lookup(km);
          pos = (res < N) ? posVec[res] : std::numeric_limits<uint64_t>::max();
          if (pos <= S - k) {
            uint64_t fk = seqVec.get_int(2 * pos, 62);
            my_mer fkm;
            fkm.word__(0) = fk; // fkm.canonicalize();
            if (mer.fwWord() == fkm.word(0) or mer.rcWord() == fkm.word(0)) {
              found += 1;

            } else {
              notFound += 1;
            }
          } else {
            // std::cerr << "pos = " << pos << ", shouldn't happen\n";
            notFound += 1;
          }
        }
      }
    }
  }
  std::cerr << "found = " << found << ", not found = " << notFound << "\n";
  std::cerr << "correctPos = " << correctPosCntr
            << ", incorrectPos = " << incorrectPosCntr << "\n";
  std::cerr << "corrextTxp = " << numTrueTxp << "\n";

  return 0;
}
