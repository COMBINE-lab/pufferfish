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
	std::map< std::pair<std::string, bool>, bool> oldPathStart;
	std::map< std::pair<std::string, bool>, bool> oldPathEnd;
	std::cerr << "Start reading GFA file... \n";
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
				oldPathStart[std::make_pair(contigVec[0].first,contigVec[0].second)] = true;
				oldPathEnd[std::make_pair(contigVec[contigVec.size()-1].first,contigVec[contigVec.size()-1].second)] = true;
			}
	}
	std::cerr << "Done Reading\nStart updating pathStart and pathEnd...\n";

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
	std::cerr << "Done updating pathStart and pathEnd based on the newIds\nStart building the graph...\n";

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
				for (auto & newIdOri : newIdList) {
						std::string id = newIdOri.first;
						bool ori = (idOri.second?newIdOri.second:!newIdOri.second);
						if (pathStart) {
								pathStart = false;
						}
						else {
								if (prev.first == id and prev.second == ori) continue;
								semiCG.addEdge(prev.first, prev.second, id, ori) ;						
						}
						prev = std::make_pair(id, ori);
				}
		}
	}
	std::cerr << "Done Building the graph\n\n\n";
	std::cerr << "\nTotal # of Contigs : " << contig_cnt << "\tTotal # of numerical Contigs : " << old2newids.size() << " \tTotal # of 30 overlaped contigs : " << new2seqAoldids.size() << "\n\n";
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
		seqParts.push_back(contigSeq.substr(1, contigSeq.size()-2)); // WTF is going on here????
	}
	seqParts.push_back(contigSeq.substr(contigSeq.size()-k, contigSeq.size()));


	short cntr = 1;
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
//		if (contigId == "486284" ) std::cerr << newId << " " << plus << " " << seq << "\n";
	}
	if (old2newids[contigId].size() > 3) std::cerr << "WTH " << contigId << "\n";
}

void GFAConverter::randomWalk() {
		spp::sparse_hash_map<std::string, pufg::Node> & nodes = semiCG.getVertices();
		std::cerr << "# of contigs with length = 31 : " << notVisited.size() << "\n";
		std::cerr << "Start merging .. \n";
		for (auto & kv : notVisited) {
				std::string curId = kv.first;
				pufg::Node& curNode = nodes[curId];
				// I hardly rely on TwoPaCo for not having a case of possible merging for in & out nodes both while none of in/out nodes are of size k!!
				if (curNode.getRealIndeg() == 1 and !isCornerCase(curNode, true)) {
						mergeIn(curNode);
				}
				else if (curNode.getRealOutdeg() == 1 and !isCornerCase(curNode, false)) {
//						std::cerr << curId << " In : " << (int)curNode.getRealIndeg() << " Out : " << (int)curNode.getRealOutdeg() << "\n";
						mergeOut(curNode);
				}
				//otherwise it is complex and you should keep the node and not merge it with any left or right neighbors
		}
}

void GFAConverter::eraseFromOldList(std::string nodeId) {
		if (new2seqAoldids.contains(nodeId)) {
			auto & seqAoldids = new2seqAoldids[nodeId];
			std::vector< std::pair<std::string, bool> >& oldids = seqAoldids.second;
//			std::cerr <<nodeId << " OldList size: " << oldids.size() << "\n";	
			for (auto& idOri : oldids) {
				std::string& id = idOri.first;
				auto& newids = old2newids[id];
//				std::cerr << "\t" << id << " was " << newids.size() ;
				for (short i = 0; i < newids.size(); i++) {
					if (newids[i].first == nodeId) {
							newids.erase(newids.begin()+i); 
					}
				}
//				std::cerr << " is " << newids.size() << "\n";
			}
		}
}
void GFAConverter::mergeIn(pufg::Node& n) {
		std::string id = n.getId();
		pufg::edgetuple edge;
		if (n.getIndegP() > 0) {
			for (auto& in : n.getIn())
					if (in.baseSign) {
						edge = in;
						break;
					}
		} else {
			for (auto& out : n.getOut())
					if (!out.baseSign) {
						edge = out;
						break;
					}
		}		
		std::string& tobeMerged = new2seqAoldids[id].first;		
		std::string& seq = new2seqAoldids[edge.contigId].first;
		if (edge.baseSign != edge.neighborSign) {
				tobeMerged = util::revcomp(tobeMerged);
				//if (tobeMerged.substr(tobeMerged.size()-(k-1)) != seq.substr(0, k-1)) std::cerr << "1 " << id << " " << edge.contigId << " " << seq << "\n" << tobeMerged << "\n";
				seq = tobeMerged.substr(0, tobeMerged.size()-(k-1)) + seq;
		}
		else {
//			if (tobeMerged.substr(0, k-1) != seq.substr(seq.size() - (k-1))) std::cerr << "2 " << seq << "\n" << tobeMerged << "\n";
			seq += tobeMerged.substr(k-1);
		}
        eraseFromOldList(id);
        new2seqAoldids.erase(id);
		semiCG.removeNode(id);
//		return edge.contigId;
}

void GFAConverter::mergeOut(pufg::Node& n) {
		std::string id = n.getId();
        pufg::edgetuple edge;
        if (n.getIndegM() > 0) { // which is equal to saying if n.getOutdegP() because only one of these two can be one!!
        	for (auto& in : n.getIn())
				if (!in.baseSign) {
					edge = in;
					break;
				}
        } else {
        	for (auto& out : n.getOut())
				if (out.baseSign) {
					edge = out;
					break;
				}
        }
        std::string& tobeMerged = new2seqAoldids[id].first;
        std::string& seq = new2seqAoldids[edge.contigId].first;
        if (edge.baseSign != edge.neighborSign) {
                tobeMerged = util::revcomp(tobeMerged);
//				if (tobeMerged.substr(0, k-1) != seq.substr(seq.size() - (k-1))) std::cerr << "3 " << seq << "\n" << tobeMerged << "\n";
	            seq += tobeMerged.substr(k-1);
		}
        else {
//			if (tobeMerged.substr(tobeMerged.size()-(k-1)) != seq.substr(0, k-1)) std::cerr << "4 " << seq << "\n" << tobeMerged << "\n";
            seq = tobeMerged.substr(0, tobeMerged.size()-(k-1)) + seq;
        }
        eraseFromOldList(id);
        new2seqAoldids.erase(id);
		semiCG.removeNode(id);
//		return edge.contigId;
}

bool GFAConverter::isCornerCase(pufg::Node& n, bool mergeIn) {
		pufg::edgetuple edge;
		if (mergeIn) {
			if (is_start(n.getId())) return true;
			if (n.getIndegP() > 0) {				
				edge = n.getIn()[0];
			} else {
				edge = n.getOut()[0];
			}
			if ( (edge.baseSign and edge.neighborSign) or (!edge.baseSign and ! edge.neighborSign) )
					return is_end(edge.contigId);
			else return is_start(edge.contigId); 
		}
		else { // Merge out case
			if (is_end(n.getId())) return true;
			if (n.getIndegM() > 0) {				
				edge = n.getIn()[0];
			} else {
				edge = n.getOut()[0];
			}
			if ( (edge.baseSign and edge.neighborSign) or (!edge.baseSign and ! edge.neighborSign) )
					return is_start(edge.contigId);
			else return is_end(edge.contigId); 
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
        for(int i = 0 ; i < vec.size(); i++) {
			auto newidVec = old2newids[vec[i].first];
			if (!vec[i].second)
				std::reverse(newidVec.begin(), newidVec.end());			
			for (auto& idOri : newidVec) {
				bool ori = vec[i].second?idOri.second:!idOri.second;
				if (!first) {
						if (prev.first == idOri.first and prev.second == ori) continue;
						gfa_file << ",";				
				}
	            gfa_file << idOri.first << (ori?"+":"-");
				prev = std::make_pair(idOri.first, ori);
				first = false;
			}
        }
		gfa_file << "\t*\n";
    }
    std::cerr << "# of contigs written to file : " << contigCntr << "\n";
}
