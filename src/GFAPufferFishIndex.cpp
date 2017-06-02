#include "FastxParser.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <iterator>
#include <type_traits>
#include "sdsl/int_vector.hpp"
#include "jellyfish/mer_dna.hpp"
#include "PufferFS.hpp"
#include "BooPHF.h"
#include "OurGFAReader.hpp"
#include "popl.hpp"
#include "ScopedTimer.hpp"

//#include "gfakluge.hpp"

uint64_t swap_uint64( uint64_t val )
{
	    val = ((val << 8) & 0xFF00FF00FF00FF00ULL ) | ((val >> 8) & 0x00FF00FF00FF00FFULL );
	        val = ((val << 16) & 0xFFFF0000FFFF0000ULL ) | ((val >> 16) & 0x0000FFFF0000FFFFULL );
		    return (val << 32) | (val >> 32);
}

using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

bool merFromStr(my_mer& mer, const std::string& s) {
  auto k = my_mer::k();
  if (s.length() < k) { return false; }
  for (size_t i = 0; i < k; ++i) {
    mer.shift_right(s[i]);
  }
  return true;
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

		ContigKmerIterator(sdsl::int_vector<>* storage, sdsl::bit_vector* rank, uint8_t k, uint64_t startAt) : 
			storage_(storage), rank_(rank), k_(k), curr_(startAt) {
      /* 
      if (startAt + k <= rank_->size()) {
        mer_.word__(0) = storage_->get_int(2 * curr_, 2 * k_);
      }
      rcMer_ = mer_.get_reverse_complement();
      */
    }

  ContigKmerIterator& operator=(ContigKmerIterator& other) = default;//(sdsl::int_vector<>* storage, sdsl::bit_vector* rank, uint8_t k, uint64_t startAt) : 


		ContigKmerIterator operator++() { ContigKmerIterator i = *this; advance_(); return i; }
		ContigKmerIterator operator++(int) { advance_(); return *this; }
		reference operator*() {
      //word_ = (mer_.word(0) < rcMer_.word(0)) ? mer_.word(0) : rcMer_.word(0); return word_;
      mer_.word__(0) = storage_->get_int(2 * curr_, 2 * k_); mer_.canonicalize(); word_ = mer_.word(0); return word_; 
		}
		difference_type pos() { return curr_; }
  //pointer operator->() { mer_.word__(0) = storage_->get_int(2 * curr_, 2 * k_); mer_.canonicalize(); word_ = mer_.word(0); return &word_; }
		bool operator==(const self_type& rhs) { return curr_ == rhs.curr_; }
		bool operator!=(const self_type& rhs) { return curr_ != rhs.curr_; }
		bool operator<(const self_type& rhs) { return curr_ < rhs.curr_; }
		bool operator<=(const self_type& rhs) { return curr_ <= rhs.curr_; }

	private:
		void advance_() {
			size_t endPos = curr_ + k_;
			if (endPos < rank_->size() and (*rank_)[endPos] == 1) {
				curr_ += k_;
        /*
        mer_.word__(0) = storage_->get_int(2 * curr_, 2 * k_); 
        rcMer_ = mer_.get_reverse_complement();
        //word_ = std::min(mer_.word(0), rcMer_.word(0));
        */
        
			} else {
        ++curr_;
        /* 
        if (curr_ + k_ <= rank_->size()) {
          auto c = my_mer::code((*storage_)[curr_+k_-1]);
          mer_.shift_right(c);
          rcMer_.shift_left(my_mer::complement(c));
          //word_ = std::min(mer_.word(0), rcMer_.word(0));
        }
        */
			}
		}
		sdsl::int_vector<>* storage_;
		sdsl::bit_vector* rank_;
		uint8_t k_;
		uint64_t curr_;
  my_mer mer_;
  my_mer rcMer_;
		uint64_t word_{0};
};

int pufferfishTest(int argc, char* argv[]) { std::cerr << "not yet implemented\n"; return 1; }

int pufferfishIndex(int argc, char* argv[]) {


  uint32_t k;
  std::string gfa_file;
  std::string cfile;
  std::string rfile;
  std::string outdir;
  popl::Switch helpOption("h", "help", "produce help message");
  popl::Value<uint32_t> kOpt("k", "klen", "length of the k-mer with which the compacted dBG was built", 31, &k);
  popl::Value<std::string> gfaOpt("g", "gfa", "path to the GFA file");
  popl::Value<std::string> outOpt("o", "output", "directory where index is written");
  popl::Value<std::string> readOpt("f", "fa", "path to the Fasta file with reads");

  popl::OptionParser op("Allowed options");
	op.add(helpOption)
    .add(kOpt)
    .add(gfaOpt)
    .add(readOpt)
    .add(outOpt);

  op.parse(argc, argv);
  if (helpOption.isSet()) {
    std::cout << op << '\n';
    std::exit(0);
  }

  gfa_file = gfaOpt.getValue();
  rfile = readOpt.getValue();
  outdir = outOpt.getValue();

	//std::vector<std::string> contig_file = {cfile};
	std::vector<std::string> read_file = {rfile};
	size_t tlen{0};	
	size_t numKmers{0};
	size_t nread{0};
	my_mer::k(k);

  puffer::fs::MakeDir(outdir.c_str());

	PosFinder pf(gfa_file.c_str(), k-1);
	pf.parseFile();	
	pf.mapContig2Pos();
  pf.serializeContigTable(outdir + "/ctable.bin");
	
/*
	auto gg = gfak::GFAKluge();
	gg.parse_gfa_file(contig_file[0]); 	
	auto seqs = gg.get_name_to_seq();
	for (auto& kv : seqs) {
		std::cerr << "k = " << kv.first << ", " << kv.second.id << "\n";
	}
	std::exit(1);
*/
  {
    auto& cnmap = pf.getContigNameMap();
    for (auto& kv : cnmap) {
      auto& r1 = kv.second;
      tlen += r1.length();
      numKmers += r1.length() - k + 1;
      ++nread;
    }
    std::cerr << "numread = " << nread << "\n";
    std::cerr << "tlen = " << tlen << "\n";
  }
  // now we know the size we need --- create our bitvectors and pack!
	size_t gpos{0};
	size_t w = std::log2(tlen) + 1;
	std::cerr << "w = " << w << "\n";
	sdsl::int_vector<> seqVec(tlen, 0, 2);
	sdsl::int_vector<> posVec(tlen, 0, w);
	sdsl::bit_vector rankVec(tlen);
	tlen = 0;
	//std::vector<uint64_t> keys;
	//keys.reserve(numKmers);
	size_t nkeys{0};
	{
    auto& cnmap = pf.getContigNameMap();
    for (auto& kv : cnmap) {
      const auto& r1 = kv.second;
      //std::cerr << "read : [" << kv.second << "]\n";
      rankVec[tlen] = 1;
      my_mer mer;
      my_mer rcMer;
      //mer.from_chars(r1.c_str());
      bool merOK = merFromStr(mer, r1);
      if (!merOK) { std::cerr << "contig too short!"; std::exit(1); }
      rcMer = mer.get_reverse_complement();
      uint64_t km = (mer.word(0) < rcMer.word(0)) ? mer.word(0) : rcMer.word(0);//mer.get_canonical().word(0);
      //keys.push_back(km);
      ++nkeys;
      for (size_t i = 0;  i < r1.length(); ++i) {
        auto offset = i;//r1.length() - i - 1;
        // NOTE: Having to add things in the reverse order here is strange 
        // we should make sure that this doesn't mess with the positions we
        // end up storing!
        auto c = my_mer::code(r1[i]);
        seqVec[gpos + offset] = c;
        if (i >= k) { 
          //keys.push_back(km);
          my_mer mm;
          uint64_t num = seqVec.get_int(2*(gpos + offset - k), 2*k);
          mm.word__(0) = num;// mm.canonicalize();
          if (!(mm == mer or mm == rcMer)) {//}mm != mer.get_canonical()) {
            std::cerr << "num & 0x3 = " << (num & 0x3) << "\n";
            std::cerr << "i = " << i << "\n";
            std::cerr << mer.to_str() << ", " << mm.to_str() <<"\n";
          }
          mer.shift_right(c);
          rcMer.shift_left(my_mer::complement(c));
          km = (mer.word(0) < rcMer.word(0)) ? mer.word(0) : rcMer.word(0);//mer.get_canonical().word(0);
          ++nkeys;
        }
        //++gpos;
      }
      gpos += r1.length();
      tlen += r1.length();
    }
  }
	 std::cerr << "seqSize = " << sdsl::size_in_mega_bytes(seqVec) << "\n";
	 std::cerr << "rankSize = " << sdsl::size_in_mega_bytes(rankVec) << "\n";
	 std::cerr << "posSize = " << sdsl::size_in_mega_bytes(posVec) << "\n";
	 std::cerr << "num keys = " << nkeys << "\n";
	 typedef boomphf::SingleHashFunctor<uint64_t>  hasher_t;
	 typedef boomphf::mphf<  uint64_t, hasher_t  > boophf_t;

	 ContigKmerIterator kb(&seqVec, &rankVec, k, 0);
	 ContigKmerIterator ke(&seqVec, &rankVec, k, seqVec.size() - k + 1);
	 auto ks = kb;
	 size_t nkeyIt{0}; for (; ks != ke; ++ks) { nkeyIt++; }
	 std::cerr << "num keys (iterator)= " << nkeyIt << "\n";
	 auto keyIt = boomphf::range(kb, ke);
	 boophf_t * bphf = new boophf_t(nkeys, keyIt, 16);//keys.size(), keys, 16);
	 std::cerr << "mphf size = " << (bphf->totalBitSize() / 8) / std::pow(2,20) << "\n";
	 size_t kstart{0};
	 size_t tlen2{0};
	 {
		 size_t i = 0;
		 my_mer fkm;
		 my_mer fkc;
		 ContigKmerIterator kb1(&seqVec, &rankVec, k, 0);
		 ContigKmerIterator ke1(&seqVec, &rankVec, k, seqVec.size() - k + 1);
		 //while (i <= seqVec.size() - k) {
		 for (; kb1 != ke1; ++kb1) {
			 //i = kb1.pos();
			 //uint64_t fk = seqVec.get_int(2*i, 2*k);
			 //fkm.word__(0) = fk; fkm.canonicalize();
			 auto idx = bphf->lookup(*kb1);//fkm.word(0));
			 if (idx >= posVec.size()) { std::cerr << "i =  " << i << ", size = " << seqVec.size() << ", idx = " << idx << ", size = " << posVec.size() << "\n";}
			 posVec[idx] = kb1.pos();
			/* 
			size_t endPos = i + k;
			if (endPos < rankVec.size() and rankVec[endPos] == 1) {
				i += k;
			} else {
				++i;
			}
			*/
		 }
	 }

   /** Write the index **/
   sdsl::store_to_file(seqVec, outdir + "/seq.bin");
   sdsl::store_to_file(rankVec, outdir + "/rank.bin");
   sdsl::store_to_file(posVec, outdir + "/pos.bin");
   std::ofstream hstream(outdir + "/mphf.bin");
   bphf->save(hstream);
   hstream.close();

   return 0;

	 size_t N = nkeys;//keys.size();
	 size_t S = seqVec.size();
	 size_t found = 0;
	 size_t notFound = 0;
	{
    ScopedTimer st;
	fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(read_file, 1, 1);
	parser.start();
	// Get the read group by which this thread will
	// communicate with the parser (*once per-thread*)
	size_t rn{0};
	 auto rg = parser.getReadGroup();
	 while (parser.refill(rg)) {
		 // Here, rg will contain a chunk of read pairs
		 // we can process.
		 for (auto& rp : rg) {
			 if (rn % 500000 == 0) { 
				 std::cerr << "rn : " << rn << "\n"; 
				 std::cerr << "found = " << found << ", notFound = " << notFound << "\n";
			 }
			 ++rn;
			 auto& r1 = rp.seq;
			 my_mer mer;
       my_mer rcMer;
			 bool merOK = merFromStr(mer, r1);//mer.from_chars(r1.begin());
       if (!merOK) { std::cerr << "contig too short!"; std::exit(1); }
       rcMer = mer.get_reverse_complement();

       //rcMer = mer.get
			 auto km = (mer.word(0) < rcMer.word(0)) ? mer.word(0) : rcMer.word(0);//mer.get_canonical().word(0); 
			 size_t res = bphf->lookup(km);
			 uint64_t pos = (res < N) ? posVec[res] : std::numeric_limits<uint64_t>::max();
			 if (pos <= S - k) { 
				 uint64_t fk = seqVec.get_int(2*pos, 2*k);
				 my_mer fkm;
				 //fkm.word__(0) = fk; fkm.canonicalize();
         fkm.word__(0) = fk;
				 if ( mer.word(0) == fkm.word(0) or rcMer.word(0) == fkm.word(0) ) { found += 1; } else { 
					 //std::cerr << "rn = " << rn - 1 << ", fk = " << fkm.get_canonical().to_str() << ", km = " << mer.get_canonical().to_str() << ", pkmer = " << pkmer <<" \n";
					 //std::cerr << "found = " << found << ", not found = " << notFound << "\n";
					 notFound += 1; }
			 } else {
				 //std::cerr << "pos = " << pos << ", shouldn't happen\n";
				 notFound += 1;
			 }
			 
			 for (size_t i = k;  i < r1.length(); ++i) {
				mer.shift_right(r1[i]);
        rcMer.shift_left(my_mer::code(my_mer::complement(r1[i])));
			 	km = (mer.word(0) < rcMer.word(0)) ? mer.word(0) : rcMer.word(0);//mer.get_canonical().word(0);
			 	res = bphf->lookup(km);
			 	pos = (res < N) ? posVec[res] : std::numeric_limits<uint64_t>::max();
				 if ( pos <= S - k) { 
					 uint64_t fk = seqVec.get_int(2*pos, 62);
					 my_mer fkm;
					 fkm.word__(0) = fk; //fkm.canonicalize();
					 if ( mer.word(0) == fkm.word(0) or rcMer.word(0) == fkm.word(0) ) { found += 1; } else { notFound += 1; }
				 } else {
				 	//std::cerr << "pos = " << pos << ", shouldn't happen\n";
					 notFound += 1;
				 }
			 }
			 
		 }
	 }
	}
	std::cerr << "found = " << found << ", not found = " << notFound << "\n";
	return 0;
}
