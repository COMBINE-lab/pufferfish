#include "FastxParser.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include "sdsl/int_vector.hpp"
#include "jellyfish/mer_dna.hpp"
#include "BooPHF.h"

uint64_t swap_uint64( uint64_t val )
{
	    val = ((val << 8) & 0xFF00FF00FF00FF00ULL ) | ((val >> 8) & 0x00FF00FF00FF00FFULL );
	        val = ((val << 16) & 0xFFFF0000FFFF0000ULL ) | ((val >> 16) & 0x0000FFFF0000FFFFULL );
		    return (val << 32) | (val >> 32);
}

int main(int argc, char* argv[]) {
	std::vector<std::string> contig_file = {argv[1]};
	std::vector<std::string> read_file = {argv[2]};
	size_t tlen{0};	
	size_t numKmers{0};
	size_t k{31};
	size_t nread{0};
	using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;
	my_mer::k(k);
		{
	fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(contig_file, 1, 1);
	parser.start();
	// Get the read group by which this thread will
	// communicate with the parser (*once per-thread*)
	 auto rg = parser.getReadGroup();
	 while (parser.refill(rg)) {
		 // Here, rg will contain a chunk of read pairs
		 // we can process.
		 for (auto& rp : rg) {
			 auto& r1 = rp.seq;
			 tlen += r1.length();
			 numKmers += r1.length() - k + 1;
			 ++nread;
		 }
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
	std::vector<uint64_t> keys;
	keys.reserve(numKmers);
	{
	fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(contig_file, 1, 1);
	parser.start();
	// Get the read group by which this thread will
	// communicate with the parser (*once per-thread*)
	 auto rg = parser.getReadGroup();
	 while (parser.refill(rg)) {
		 // Here, rg will contain a chunk of read pairs
		 // we can process.
		 for (auto& rp : rg) {
			 auto& r1 = rp.seq;
			 rankVec[tlen] = 1;
			 my_mer mer;
			 mer.from_chars(r1.begin());
			 uint64_t km = mer.get_canonical().word(0);
			 keys.push_back(km);
			 for (size_t i = 0;  i < r1.length(); ++i) {
				 auto offset = r1.length() - i - 1;
				 seqVec[gpos + offset] = my_mer::code(r1[i]);
				 if (i >= k) { 
					 mer.shift_left(r1[i]);
					 km = mer.get_canonical().word(0);
					 keys.push_back(km);
					 my_mer mm;
					 uint64_t num = seqVec.get_int(2*(gpos + offset), 2*k);
					 mm.word__(0) = num; mm.canonicalize();
					 if (mm != mer.get_canonical()) {
					 std::cerr << "num & 0x3 = " << (num & 0x3) << "\n";
					 std::cerr << "i = " << i << "\n";
					 std::cerr << mer.to_str() << ", " << mm.to_str() <<"\n";
					 }

				 }
				 //++gpos;
			 }
			 gpos += r1.length();
			 tlen += r1.length();
		 }
	 }
	}
	 std::cerr << "seqSize = " << sdsl::size_in_mega_bytes(seqVec) << "\n";
	 std::cerr << "rankSize = " << sdsl::size_in_mega_bytes(rankVec) << "\n";
	 std::cerr << "posSize = " << sdsl::size_in_mega_bytes(posVec) << "\n";
	 std::cerr << "num keys = " << keys.size() << "\n";
	 typedef boomphf::SingleHashFunctor<uint64_t>  hasher_t;
	 typedef boomphf::mphf<  uint64_t, hasher_t  > boophf_t;
	 boophf_t * bphf = new boophf_t(keys.size(), keys, 16);
	 std::cerr << "mphf size = " << (bphf->totalBitSize() / 8) / std::pow(2,20) << "\n";
	 size_t kstart{0};
	 size_t tlen2{0};
	 {
		 size_t i = 0;
		 my_mer fkm;
		 my_mer fkc;
		 while (i < seqVec.size() - k) {

			 uint64_t fk = seqVec.get_int(2*i, 2*k);
			 fkm.word__(0) = fk; fkm.canonicalize();
			 auto idx = bphf->lookup(fkm.word(0));
			 if (idx > posVec.size()) { std::cerr << "i =  " << i << ", size = " << seqVec.size() << ", idx = " << idx << ", size = " << posVec.size() << "\n";}
			 posVec[idx] = i;
			 if (rankVec[i + k] == 1) {
				 i += k;
			 } else { 
				 ++i;
			 }
		 }
	 }

	 size_t N = keys.size();
	 size_t S = seqVec.size();
	 size_t found = 0;
	 size_t notFound = 0;
	{
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
			 mer.from_chars(r1.begin());
			 auto km = mer.get_canonical().word(0); 
			 size_t res = bphf->lookup(km);
			 uint64_t pos = (res < N) ? posVec[res] : std::numeric_limits<uint64_t>::max();
			 if ( pos < S - k) { 
				 uint64_t fk = seqVec.get_int(2*pos, 62);
				 my_mer fkm;
				 fkm.word__(0) = fk; fkm.canonicalize();
				 if ( km == fkm.word__(0) ) { found += 1; } else { 
					 //std::cerr << "fk = " << fkm.get_canonical().to_str() << ", km = " << mer.get_canonical().to_str() << "\n";
					 //std::cerr << "found = " << found << ", not found = " << notFound << "\n";
					 notFound += 1; }
			 } else {
				 //std::cerr << "pos = " << pos << "\n";
				 notFound += 1;
			 }
			 /*
			 for (size_t i = k;  i < r1.length(); ++i) {
				mer.shift_left(r1[i]);
			 	km = mer.get_canonical().word(0); 
			 	res = bphf->lookup(km);
			 	pos = (res < N) ? posVec[res] : std::numeric_limits<uint64_t>::max();
				 if ( pos < N - k) { 
					 uint64_t fk = seqVec.get_int(2*pos, 62);
					 if ( fk == km ) { found += 1; } else { notFound += 1; }
				 }
			 }
			 */
		 }
	 }
	}
	std::cerr << "found = " << found << ", not found = " << notFound << "\n";
	return 0;
}
