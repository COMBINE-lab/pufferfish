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
	std::map< std::pair<std::string, bool>, bool> oldPathStart;
	std::map< std::pair<std::string, bool>, bool> oldPathEnd;
	{
	std::cerr << "Start reading GFA file... \n";
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
				}
				contig_cnt++;
			} 
			// A path line
			if (tag == "P") {
				auto pvalue = splited[2];
				std::vector<std::pair<std::string, bool> > contigVec = util::explode(pvalue, ',');
				// parse value and add all conitgs to contigVec
				path[id] = contigVec;
				oldPathStart[std::make_pair(contigVec[0].first,contigVec[0].second)] = true;
				oldPathEnd[std::make_pair(contigVec[contigVec.size()-1].first,contigVec[contigVec.size()-1].second)] = true;
			}
	}
	std::cerr << "Done Reading\nStart updating pathStart and pathEnd...\n";
	}
	// set pathStart and pathEnd to newids	
	auto cntr = 0;
	for (auto & kv : oldPathStart) {
		auto s = kv.first; // just existence of key matters
		auto newIdList = old2newids[s.first];	
		cntr++;
		if (s.second) { // the positive orientation of the old id
			// Add the first new id and keep its orientation
			pathStart[newIdList[0]] = true;
		}
		else {
			// Add last element in new is list and flip the orientation
			auto newIdOri = newIdList[newIdList.size()-1];
			pathStart[std::make_pair(newIdOri.first, !newIdOri.second)] = true;
		}		
		oldPathStart.erase(s);
	}
	for (auto & kv : oldPathEnd) {
		auto s = kv.first; // just existence of key matters
		auto newIdList = old2newids[s.first];	
		if (s.second) { // the positive orientation of the old id
			// Add the first new id and keep its orientation
			pathEnd[newIdList[newIdList.size()-1]] = true;
		}
		else {
			// Add last element in new is list and flip the orientation
			auto newIdOri = newIdList[0];
			pathEnd[std::make_pair(newIdOri.first, !newIdOri.second)] = true;
		}
		oldPathEnd.erase(s);
	}
	std::cerr << "Done updating pathStart and pathEnd based on the newIds\n";
	std::cerr << "Total # of Contigs : " << contig_cnt << "\tTotal # of numerical Contigs : " << old2newids.size() << " \tTotal # of 30 overlaped contigs : " << new2seqAoldids.size() << "\n";
}

void GFAConverter::buildGraph() {
		std::cerr << "Start building the graph...\n";
	// build the graph
	for (auto & kv : path) {
		std::string pathId = kv.first;
		auto & contigVec = kv.second;
		
		std::pair<std::string, bool> prev;
		bool pathStart = true;
		for (auto & idOri : contigVec){
				std::vector<std::pair<std::string, bool> > newIdList = old2newids[idOri.first];			
				if (!idOri.second) {						
						std::reverse(newIdList.begin(), newIdList.end());	// As newIdList is a copy it won't affect the actual vector in old2newids (tested)
				}
				if (pathStart) {
						std::string id = newIdList[0].first;						
						bool ori = (idOri.second?newIdList[0].second:!newIdList[0].second);
						prev = std::make_pair(id, ori);
						pathStart = false;
				}
				for (size_t i = 1; i < newIdList.size(); i++) {
						std::string id = newIdList[i].first;
						bool ori = (idOri.second?newIdList[i].second:!newIdList[i].second);
						semiCG.addEdge(prev.first, prev.second, id, ori) ;						
						prev = std::make_pair(id, ori);
				}
		}
	}
	std::cerr << "Done Building the graph\n";

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
	if (util::isRevcomp(contigSeq)) {
		for (size_t i = 0; i <= contigSeq.size()-k; i++) {
			seqParts.push_back(contigSeq.substr(i, k));
		}
	}
	else {
		seqParts.push_back(contigSeq.substr(0, k));
		if (contigSeq.size() > k + 1) {
			seqParts.push_back(contigSeq.substr(1, contigSeq.size()-2)); 
		}
		seqParts.push_back(contigSeq.substr(contigSeq.size()-k));
	}

	// for each part, search in hash_map whether it already existed or not and assign the proper (new/retrieved) contig id to it
	int cntr = 1;
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
			newId = prefix + std::to_string(cntr) + contigId;
			seq2newid[seq] = newId;			
		}

		if (new2seqAoldids.find(newId) != new2seqAoldids.end()) {
			new2seqAoldids[newId].second.emplace_back(contigId, plus);
		}
		else {
			std::vector< std::pair<std::string, bool> > ids;
			ids.emplace_back(contigId, plus);
			new2seqAoldids[newId] = std::make_pair(seq, ids);
			if (seq.size() == k)
				notVisited[newId] = true;
		}
		old2newids[contigId].emplace_back(newId, plus);
		cntr++;
	}
//	if (old2newids[contigId].size() > 3) std::cerr << "RevComp is true: " << contigId << " " << contigSeq << "\n";
}

void GFAConverter::randomWalk() {
		spp::sparse_hash_map<std::string, pufg::Node> & nodes = semiCG.getVertices();
		std::cerr << "# of contigs with length = 31 : " << notVisited.size() << "\n";
		std::cerr << "\nStart merging .. \n";
		for (auto & kv : notVisited) {
				std::string curId = kv.first;
				pufg::Node& curNode = nodes[curId];
				// I strongly rely on TwoPaCo here for not having a case of possible merging for in & out nodes both while none of in/out nodes are of size k!!
				if (curNode.getRealIndeg() == 1 and !isCornerCase(curNode, true)) {
						mergeIn(curNode);
				}
				else if (curNode.getRealOutdeg() == 1 and !isCornerCase(curNode, false)) {
						mergeOut(curNode);
				}
				//otherwise it is complex and you should keep the node and not merge it with any left or right neighbors
		}
		std::cerr << "Done merging \n";
}

void GFAConverter::eraseFromOldList(std::string nodeId) {
		if (new2seqAoldids.contains(nodeId)) {
			auto & seqAoldids = new2seqAoldids[nodeId];
			std::vector< std::pair<std::string, bool> >& oldids = seqAoldids.second;
			for (auto& idOri : oldids) {
				std::string& id = idOri.first;
				auto& newids = old2newids[id];
				for (size_t i = 0; i < newids.size(); i++) {
					if (newids[i].first == nodeId) {
							newids.erase(newids.begin()+i); 
					}
				}
			}
		}
}
void GFAConverter::mergeIn(pufg::Node& n) {
		std::string id = n.getId();
		pufg::edgetuple& edge = n.getOnlyRealIn();
    if (!new2seqAoldids.contains(id)) {
      std::cerr << "[mergeIn] NO; the id " << id << " was not in new2seqAoldids!\n";
    }
    if (!new2seqAoldids.contains(edge.contigId)) {
      std::cerr << "[mergeIn] NO; the edge.contigId " << edge.contigId << " was not in new2seqAoldids!\n";
    }
		std::string& tobeMerged = new2seqAoldids[id].first;		
		std::string& seq = new2seqAoldids[edge.contigId].first;
		if (edge.baseSign != edge.neighborSign) {
				tobeMerged = util::revcomp(tobeMerged);
//				if (tobeMerged.substr(tobeMerged.size()-(k-1)) != seq.substr(0, k-1)) std::cerr << "1 " << id << " " << edge.contigId << " " << seq << "\n" << tobeMerged << "\n";
				seq = tobeMerged.substr(0, tobeMerged.size()-(k-1)) + seq;
		}
		else {
//			if (tobeMerged.substr(0, k-1) != seq.substr(seq.size() - (k-1))) std::cerr << "2 " << seq << "\n" << tobeMerged << "\n";
			seq += tobeMerged.substr(k-1);
		}
        eraseFromOldList(id);
        new2seqAoldids.erase(id);
		semiCG.removeNode(id);
		if (is_start(id)) update_start(edge.contigId, edge.baseSign == edge.neighborSign);
		if (is_end(id)) update_end(edge.contigId, edge.baseSign == edge.neighborSign);
//		return edge.contigId;
}

void GFAConverter::mergeOut(pufg::Node& n) {
		std::string id = n.getId();
        pufg::edgetuple& edge = n.getOnlyRealOut();
        if (!new2seqAoldids.contains(id)) {
          std::cerr << "[mergeOut] NO; the id " << id << " was not in new2seqAoldids!\n";
        }
        if (!new2seqAoldids.contains(edge.contigId)) {
          std::cerr << "[mergeOut] NO; the edge.contigId " << edge.contigId << " was not in new2seqAoldids!\n";
    }
        std::string& tobeMerged = new2seqAoldids[id].first;
        std::string& seq = new2seqAoldids[edge.contigId].first;
        if (edge.baseSign != edge.neighborSign) {
                tobeMerged = util::revcomp(tobeMerged);
//				if (tobeMerged.substr(0, k-1) != seq.substr(seq.size() - (k-1))) std::cerr << id << " " << edge.contigId << " " << "3 " << seq << "\n" << tobeMerged << "\n";
	            seq += tobeMerged.substr(k-1);
		}
        else {
//			if (tobeMerged.substr(tobeMerged.size()-(k-1)) != seq.substr(0, k-1)) std::cerr << "4 " << seq << "\n" << tobeMerged << "\n";
            seq = tobeMerged.substr(0, tobeMerged.size()-(k-1)) + seq;
        }
        eraseFromOldList(id);
        new2seqAoldids.erase(id);
		semiCG.removeNode(id);
		if (is_start(id)) update_start(edge.contigId, edge.baseSign == edge.neighborSign);
		if (is_end(id)) update_end(edge.contigId, edge.baseSign == edge.neighborSign);
//		return edge.contigId;
}

/**
 * Check 3 cases
 * 1. node is not corner node
 * 2. the neighbor we want to merge the node with is not corner node
 * 3. the current node is not the same as the node it will be merged with
 */
bool GFAConverter::isCornerCase(pufg::Node& n, bool mergeIn) {
		if (mergeIn) {				
			if (is_start(n.getId())) return true;			
			pufg::edgetuple& edge = n.getOnlyRealIn();
			if (edge.contigId == n.getId()) return true;
			pufg::Node& neighbor = semiCG.getVertices()[edge.contigId];
			//if ( (edge.baseSign and edge.neighborSign) or (!edge.baseSign and ! edge.neighborSign) )
			if (edge.baseSign == edge.neighborSign) {
					return is_end(edge.contigId) or neighbor.getRealOutdeg() != 1;
			}
			else return is_start(edge.contigId) or neighbor.getRealIndeg() != 1; 
		}
		else { // Merge out case
			if (is_end(n.getId())) return true;
			pufg::edgetuple& edge = n.getOnlyRealOut();
			if (edge.contigId == n.getId()) return true;
			pufg::Node& neighbor = semiCG.getVertices()[edge.contigId];
			//if ( (edge.baseSign and edge.neighborSign) or (!edge.baseSign and ! edge.neighborSign) ) {
			if (edge.baseSign == edge.neighborSign) {
					return is_start(edge.contigId) or neighbor.getRealIndeg() != 1;
			}
			else {
					return is_end(edge.contigId) or neighbor.getRealOutdeg() != 1; 
			}
		}
		return false;
}

bool GFAConverter::is_start(std::string& nodeId) {
	if (pathStart.find(std::make_pair(nodeId, true)) != pathStart.end() or
			pathEnd.find(std::make_pair(nodeId, false)) != pathEnd.end()) return true;
	return false;
}

bool GFAConverter::is_end(std::string& nodeId) {
	if (pathStart.find(std::make_pair(nodeId, false)) != pathStart.end() or
			pathEnd.find(std::make_pair(nodeId, true)) != pathEnd.end()) return true;
	return false;
}

void GFAConverter::update_start(std::string& newId, bool newOri) {
		pathStart[std::make_pair(newId, newOri)] = true;
		pathEnd[std::make_pair(newId, !newOri)] = true;
}

void GFAConverter::update_end(std::string& newId, bool newOri) {
		pathStart[std::make_pair(newId, !newOri)] = true;
		pathEnd[std::make_pair(newId, newOri)] = true;
}

void GFAConverter::writeFile(const char* gfaFileName) {

	std::ofstream gfa_file(gfaFileName) ;
    uint32_t contigCntr = 0;
    for(auto& kv : new2seqAoldids) {		
		gfa_file << "S" << "\t" << kv.first <<"\t" << (kv.second).first << "\n" ;
        contigCntr++;
    }
    for(auto& p : path) {
        auto tid = p.first ;
        gfa_file << "P" << "\t" << tid << "\t"  ;
        auto vec = p.second ;
		bool first = true;
		std::pair<std::string, bool> prev;
		bool firstPiece = true;
        for(size_t i = 0 ; i < vec.size(); i++) {
			auto newidVec = old2newids[vec[i].first];
			if (!vec[i].second)
				std::reverse(newidVec.begin(), newidVec.end());						
			for (auto& idOri : newidVec) {
				std::string id = idOri.first;
				bool ori = vec[i].second?idOri.second:!idOri.second;
				if (!first) {
						if (firstPiece and prev.first == id and prev.second == ori) {
								firstPiece = false;
								if (new2seqAoldids[id].first.size() == k)
									continue;
						}
/*						if (firstPiece and prev.first == id and prev.second != ori) {
							std::cerr << "found the bug " << prev.first << "\n";
						}*/
						gfa_file << ",";				
				}
	            gfa_file << id << (ori?"+":"-");
				prev = std::make_pair(id, ori);
				first = false;
			}
			if (newidVec.size() != 0)
				firstPiece = true;			
			// here's the crazy case.
			// If an old contig's list of new ids is completely empty, this means that we are actually jumping over this contig, 
			// so we shouldn't be worried about the next new contig with previous one, as prev is not actually the real prev.
			// Say we have the case A- B+ A- in old set of contigs
			// A : k-, k'+, k''+
			// B : k-, k''+
			// Now if we merge all 3 pieces of A into one final contig F, all parts of B are removed.
			// Then if we see F+, F+ we shouldn't remove the second F+ b/c F+ is a newly created contig that has overlap of k-1 with itself.
			// The path is actually F+ followed by F+ and not the case that these two contigs are equal kmers with overlap of k.
			// Adding one "else" and a whole story! It's not even counted as a comment anymore!!!!!! :/
			else firstPiece = false;
        }
		gfa_file << "\t*\n";
    }
    std::cerr << "# of contigs written to file : " << contigCntr << "\n";
}
