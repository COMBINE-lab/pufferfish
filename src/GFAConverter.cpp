#include <algorithm>
#include "xxhash.h"
#include "GFAConverter.hpp"
#include "cereal/archives/binary.hpp"

GFAConverter::GFAConverter(const char* gfaFileName, size_t input_k) {
	std::cerr << "Reading GFA file " << gfaFileName << "\n";
	file.reset(new zstr::ifstream(gfaFileName));
	k = input_k;
}

void GFAConverter::parseFile() {
	std::string ln;
	std::string tag, id, value;
	size_t contig_cnt{0};
	spp::sparse_hash_map<std::string, std::string> seq2newid;
	while(std::getline(*file, ln)) {
			char firstC = ln[0];
			if (firstC != 'S' and firstC != 'P') continue;
			stx::string_view lnview(ln);
			std::vector<stx::string_view> splited = util::split(lnview, '\t');
			tag = splited[0].to_string();
			id = splited[1].to_string();
			value = splited[2].to_string();
			// A segment line
			if (tag == "S") {						
				if (util::is_number(id)) {
					processContigSeq(id, value, seq2newid);
					//contigid2seq[id] = value;
				}
				contig_cnt++;
			} 
			// A path line
			if (tag == "P") {
				auto pvalue = splited[2];
				std::vector<std::pair<std::string, bool> > contigVec = util::explode(pvalue, ',');
				// parse value and add all conitgs to contigVec
				path[id] = contigVec;
				pathStart[std::make_pair(contigVec[0].first,contigVec[0].second)] = true;
				pathEnd[std::make_pair(contigVec[contigVec.size()-1].first,contigVec[contigVec.size()-1].second)] = true;
			}
	}

	// build the graph
	spp::sparse_hash_map<std::string, std::string> kmer2contigid;
	std::string prefix = "00";
	uint32_t newContigCntr = 0;
	for (auto & kv : path) {
		std::string id = kv.first;
		auto & contigVec = kv.second;
		
		std::pair<std::string, bool> prev;
		bool pathStart = true;
		for (auto & idOri : contigVec){
				std::vector<std::pair<std::string, bool> > newIdList = old2newids[idOri.first];				
				if (!idOri.second) {						
						std::reverse(newIdList.begin(), newIdList.end());	// As newIdList is a copy it won't affect the actual vector in old2newids (tested)
				}
				for (auto & newIdOri : old2newids[idOri.first]) {
						std::string id = newIdOri.first;
						bool ori = (idOri.second?newIdOri.second:!newIdOri.second);
						if (pathStart) {
								prev = std::make_pair(id, ori);
								pathStart = false;
						}
						else {
								semiCG.addEdge(prev.first, prev.second, id, ori) ;						
						}
				}
		}
	}
	std::cerr << " Total # of Contigs : " << contig_cnt << " Total # of numerical Contigs : " << old2newids.size() << "\n\n";
}

void GFAConverter::processContigSeq(std::string & contigId, std::string & contigSeq, spp::sparse_hash_map<std::string, std::string> & seq2newid) {
	/**
	 * Divide every segment into 3 pieces
	 * first k 
	 * middle one without the first and last nucleotide
	 * last k
	 * special case is segment with length = k+1 that doesn't have the middle case
	 * ATTENTION : The order we insert these 3 pieces into the vector matters. So keep it the way it is
	**/
	std::string prefix = "00";
	std::vector<std::string> seqParts;
	seqParts.push_back(contigSeq.substr(0, k));
	if (contigSeq.size() > k + 1) {
		seqParts.push_back(contigSeq.substr(1, contigSeq.size()-1));
	}
	seqParts.push_back(contigSeq.substr(contigSeq.size()-k, contigSeq.size()));

	
	for (std::string seq : seqParts) {
		std::string newId;
		bool plus = true;
		if (seq2newid.find(seq) != seq2newid.end()) {
			newId = seq2newid[seq];
		}
		else if (seq2newid.find(util::revcomp(seq)) != seq2newid.end()) {
			plus = false;
			newId = seq2newid[util::revcomp(seq)];
		}
		else {
			newId = prefix + "0" + contigId;
			seq2newid[seq] = newId;			
		}
		if (new2seqAoldids.find(newId) != new2seqAoldids.end()) {
			new2seqAoldids[newId].second.emplace_back(contigId, plus);
		}
		else {
			std::vector< std::pair<std::string, bool> > ids;
			ids.emplace_back(contigId, plus);
			new2seqAoldids[newId] = std::make_pair(seq, ids);
		}
		old2newids[contigId].emplace_back(newId, plus);
	}
}
