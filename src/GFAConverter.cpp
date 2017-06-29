#include "FastxParser.hpp"
#include "GFAConverter.hpp"
#include "cereal/archives/binary.hpp"
#include "xxhash.h"
#include <algorithm>
#include <assert.h>

GFAConverter::GFAConverter(const char* gfaFileName, size_t input_k) {
		gfaFileName_ = gfaFileName;
  std::cerr << "Reading GFA file " << gfaFileName << "\n";
  file.reset(new zstr::ifstream(gfaFileName));
  k = input_k;
}

void GFAConverter::parseFile() {
  std::string ln;
  std::string tag, id, value;
  size_t contig_cnt{0};
    std::cerr << "Start loading segments... \n";
    spp::sparse_hash_map<std::string, uint64_t> seq2id;
    uint64_t idCntr = 0;
    while (std::getline(*file, ln)) {
      char firstC = ln[0];
 		if (firstC != 'S')
			continue;
      	stx::string_view lnview(ln);
	    std::vector<stx::string_view> splited = util::split(lnview, '\t');
	   // tag = splited[0].to_string();
    	id = splited[1].to_string();
	    value = splited[2].to_string();
        // A segment line
        if (util::is_number(id)) {
		  uint64_t contigId = std::stoll(id);
          processContigSeq(contigId, value, seq2id, idCntr);
        }
        contig_cnt++;
	}	
    std::cerr << "Done loading & processing segments ...\n\nStart loading paths\n";

	file.reset(new zstr::ifstream(gfaFileName_));
	while (std::getline(*file, ln)) {
      	char firstC = ln[0];
		if (firstC != 'P')
			continue;
      	stx::string_view lnview(ln);
	    std::vector<stx::string_view> splited = util::split(lnview, '\t');
  	    // tag = splited[0].to_string();
    	id = splited[1].to_string();
        // A path line
        auto pvalue = splited[2];
        // parse value and add all conitgs to contigVec
        std::vector<std::pair<uint64_t, bool>> contigVec =
            util::explode(pvalue, ',');

		contigVec = convertOld2newPath(contigVec);

		// update graph and add the new set of segments to it
		updateGraph(contigVec);

		// set pathStart and pathEnd
		setPathStart(id, contigVec[0]);
		setDonTouchSide(contigVec);
      }

    std::cerr << "Done Reading gfa file\n";

  std::cerr << "Total # of Contigs : " << contig_cnt
            << "\tTotal # of numerical Contigs : " << old2newids.size()
            << " \tTotal # of contigs after spliting :" << id2seq.size()
            << "\n\n";
}

void GFAConverter::processContigSeq(
    uint64_t& contigId, std::string& contigSeq,
    spp::sparse_hash_map<std::string, uint64_t>& seq2id,
    uint64_t& idCntr) {
  /**
   * Divide every segment into 3 pieces
   * first k
   * middle one without the first and last nucleotide
   * last k
   * special case is segment with length = k+1 that doesn't have the middle case
   * ATTENTION : The order we insert these 3 pieces into the vector matters. So
   *keep it the way it is
  **/
  std::vector<std::string> seqParts;
  if (util::isRevcomp(contigSeq)) {
    for (size_t i = 0; i <= contigSeq.size() - k; i++) {
      seqParts.push_back(contigSeq.substr(i, k));
    }
  } else {
    seqParts.push_back(contigSeq.substr(0, k));
    if (contigSeq.size() > k + 1) {
      seqParts.push_back(contigSeq.substr(1, contigSeq.size() - 2));
    }
    seqParts.push_back(contigSeq.substr(contigSeq.size() - k));
  }

  // for each part, search in hash_map whether it already existed or not and
  // assign the proper (new/retrieved) contig id to it
  // int cntr = 1;
  for (std::string seq : seqParts) {
    uint64_t newId;
    bool plus = true;
    if (seq2id.contains(seq)) {
      newId = seq2id[seq];
    } else if (seq2id.contains(util::revcomp(seq))) {
      plus = false;
      newId = seq2id[util::revcomp(seq)];
    } else {
      newId = idCntr;
      idCntr++;
      seq2id[seq] = newId;
	  id2seq[newId] = seq;
//	  if (newId == 70550)
//			  std::cerr <<contigId << " to " << newId << " " << seq << "\n";
      if (seq.size() == k)
        ksizeContig.push_back(newId);
    }
    old2newids[contigId].emplace_back(newId, plus);
  }
}

std::vector<std::pair<uint64_t, bool> >& GFAConverter::convertOld2newPath(std::vector<std::pair<uint64_t, bool> >& oContigVec) {
	std::vector<std::pair<uint64_t, bool> > nContigVec = {};
	bool isPathStart = true;
	for (auto& idOri : oContigVec) {
			auto newIdList = old2newids[idOri.first];	  
			if (!idOri.second) {
				// As newIdList is a copy it won't affect the actual vector in old2newids (tested)
	    	    std::reverse(newIdList.begin(), newIdList.end());
			}			
			if (isPathStart) {
					nContigVec.emplace_back(newIdList[0].first, idOri.second?newIdList[0].second:!newIdList[0].second);
					isPathStart = false;
			}			
			for (size_t i = 1; i < newIdList.size(); i++) {
				nContigVec.emplace_back(newIdList[i].first, idOri.second?newIdList[i].second:!newIdList[i].second);
			}
	}
	oContigVec = nContigVec;
	return oContigVec;
}

void GFAConverter::updateDonTouchSide(uint64_t& id, Side side) {
		if (!donTouch.contains(id))
			donTouch[id] = side; 
		else if (donTouch[id] != side)
			donTouch[id] = Side::LEFT_RIGHT;
//		if (id == 70547 or id == 70550)
//				std::cerr << id << " " << 
//						(donTouch[id] == Side::LEFT?"left":(donTouch[id] == Side::RIGHT?"right":"both")) << "\n";
}

void GFAConverter::setDonTouchSide(std::vector<std::pair<uint64_t, bool> >& contigVec) {
		if (contigVec[0].second) 
				updateDonTouchSide(contigVec[0].first, Side::LEFT);
		else 
				updateDonTouchSide(contigVec[0].first, Side::RIGHT);

		if (contigVec[contigVec.size()-1].second) 
				updateDonTouchSide(contigVec[contigVec.size()-1].first, Side::RIGHT);
		else 
				updateDonTouchSide(contigVec[contigVec.size()-1].first, Side::LEFT);
}

void GFAConverter::mergeDonTouchSides(uint64_t& mergedId, bool mergedOri, uint64_t& mergingId, bool mergingOri) {		
		if (donTouch.contains(mergedId)) {
				Side s;
			if (mergedOri == mergingOri) {
					s = donTouch[mergedId];
			}
			else {	
					switch (donTouch[mergedId]) {
							case Side::LEFT: s = Side::RIGHT;
											break;
							case Side::RIGHT: s = Side::LEFT;
											 break;
							default: s = donTouch[mergedId];
					}
			}
			updateDonTouchSide(mergingId, s);
		}
}

void GFAConverter::setPathStart(std::string& pathId, std::pair<uint64_t, bool>& startNode) {
//		if (startNode.first == 70550) 
//			std::cerr << "path start " << startNode.first << " " << startNode.second << "\n";
		
	path2start[pathId] = startNode;
	if (start2path.find(startNode) == start2path.end())
			start2path[startNode] = {};
	start2path[startNode].push_back(pathId);
}

void GFAConverter::updatePathStart(uint64_t& oStartNodeId, bool oStartNodeOri, uint64_t& nStartNodeId, bool nStartNodeOri) {
	std::pair<uint64_t, bool> oStartNode = std::make_pair(oStartNodeId, oStartNodeOri);
	std::pair<uint64_t, bool> nStartNode = std::make_pair(nStartNodeId, nStartNodeOri);
	if (start2path.find(oStartNode) != start2path.end()) {		
		auto & oldP = start2path[oStartNode];
		if (start2path.find(nStartNode) == start2path.end())
				start2path[nStartNode] = {};
		for (auto& v : oldP) {
//			if (nStartNode.first == 70550) 
//			std::cerr << "1 path start from " <<  oStartNode.first << " " << oStartNode.second << " to " << nStartNode.first << " " << nStartNode.second << "\n";

			path2start[v] = nStartNode;
			start2path[nStartNode].push_back(v);
		}
	}

	oStartNode = std::make_pair(oStartNodeId, !oStartNodeOri);
	nStartNode = std::make_pair(nStartNodeId, !nStartNodeOri);
	if (start2path.find(oStartNode) != start2path.end()) {		
		auto & oldP = start2path[oStartNode];
		if (start2path.find(nStartNode) == start2path.end())
				start2path[nStartNode] = {};
		for (auto& v : oldP) {
//			if (nStartNode.first == 70550) 
//			std::cerr << "2 path start " << nStartNode.first << " " << nStartNode.second << "\n";
		path2start[v] = nStartNode;
			start2path[nStartNode].push_back(v);
		}
	}

}

void GFAConverter::updateGraph(std::vector<std::pair<uint64_t, bool> >& contigVec) {
    for (size_t i = 0; i < contigVec.size()-1; i++) 
        semiCG.addEdge(contigVec[i].first, contigVec[i].second, contigVec[i+1].first, contigVec[i+1].second);
}

void GFAConverter::randomWalk() {
  spp::sparse_hash_map<uint64_t, pufg::Node>& nodes = semiCG.getVertices();
  std::cerr << "# of contigs with length = " << k << " : " << ksizeContig.size() << "\n";
  std::cerr << "\nStart merging .. \n";
  for (auto& v : ksizeContig) {
    uint64_t curId = v;
    pufg::Node& curNode = nodes[curId];
    // I strongly rely on TwoPaCo here for not having a case of possible merging
    // for in & out nodes both while none of in/out nodes are of size k!!
    if (curNode.getRealIndeg() == 1 and !isCornerCase(curNode, true)) {
      mergeLeft(curNode);
    } else if (curNode.getRealOutdeg() == 1 and !isCornerCase(curNode, false)) {
      mergeRight(curNode);
    }
    // otherwise it is complex and you should keep the node and not merge it
    // with any left or right neighbors
  }
  std::cerr << "Done merging \n\n";
}

void GFAConverter::mergeLeft(pufg::Node& n) {
  uint64_t id = n.getId();
  pufg::edgetuple& edge = n.getOnlyRealIn();
  if (!id2seq.contains(id)) {
    std::cerr << "[mergeLeft] NO; the id " << id
              << " was not in id2seq!\n";
  }
  if (!id2seq.contains(edge.contigId)) {
    std::cerr << "[mergeLeft] NO; the edge.contigId " << edge.contigId
              << " was not in id2seq!\n";
  }
//  std::cerr << "mergeLeft " << id << " " << edge.contigId << "\n";
  std::string& tobeMerged = id2seq[id];
  std::string& seq = id2seq[edge.contigId];
  if (edge.baseSign() != edge.neighborSign()) {
    tobeMerged = util::revcomp(tobeMerged);
    //				if (tobeMerged.substr(tobeMerged.size()-(k-1)) != seq.substr(0,
    //k-1)) std::cerr << "1 " << id << " " << edge.contigId << " " << seq <<
    //"\n" << tobeMerged << "\n";
    seq = tobeMerged.substr(0, tobeMerged.size() - (k - 1)) + seq;
  } else {
    //			if (tobeMerged.substr(0, k-1) != seq.substr(seq.size() - (k-1)))
    //std::cerr << "2 " << seq << "\n" << tobeMerged << "\n";
    seq += tobeMerged.substr(k - 1);
  }
//  if (edge.contigId == 70550 or edge.contigId == 70547) {
//  	std::cerr << "MergeLeft " << id << " " << edge.baseSign() << " to " << edge.contigId << " " << edge.neighborSign() << " " << seq << "\n";
//  }
  id2seq.erase(id);
  //semiCG.removeNode(id);  
//  if (edge.contigId == 70550) std::cerr <<"Ajabaa " <<edge.contigId << " " <<edge.neighborSign() << "\n";
  updatePathStart(id, edge.baseSign(), edge.contigId, edge.neighborSign());
  mergeDonTouchSides(id, edge.baseSign(), edge.contigId, edge.neighborSign());
  semiCG.merge(id, edge.baseSign(), edge.contigId, edge.neighborSign(), false);
}

void GFAConverter::mergeRight(pufg::Node& n) {
  uint64_t id = n.getId();
  pufg::edgetuple& edge = n.getOnlyRealOut();
  if (!id2seq.contains(id)) {
    std::cerr << "[mergeRight] NO; the id " << id
              << " was not in id2seq!\n";
  }
  if (!id2seq.contains(edge.contigId)) {
    std::cerr << "[mergeRight] NO; the edge.contigId " << edge.contigId
              << " was not in id2seq!\n";
  }
  //std::cerr << "mergeRight " << id << " " << edge.contigId << " ";
  std::string& tobeMerged = id2seq[id];
  std::string& seq = id2seq[edge.contigId];
  if (edge.baseSign() != edge.neighborSign()) {
//		  std::cerr << " 1 " << tobeMerged.size() <<" ";
    tobeMerged = util::revcomp(tobeMerged);
    //				if (tobeMerged.substr(0, k-1) != seq.substr(seq.size() - (k-1)))
    //std::cerr << id << " " << edge.contigId << " " << "3 " << seq << "\n" <<
    //tobeMerged << "\n";
    seq += tobeMerged.substr(k - 1);
  } else {
//		  std::cerr << " 2 " << tobeMerged.size() <<" ";
    //			if (tobeMerged.substr(tobeMerged.size()-(k-1)) != seq.substr(0,
    //k-1)) std::cerr << "4 " << seq << "\n" << tobeMerged << "\n";
    seq = tobeMerged.substr(0, tobeMerged.size() - (k - 1)) + seq;
  }
//  if (edge.contigId == 70550 or edge.contigId == 70547) {
//  	std::cerr << "MergeRight " << id << " " << edge.baseSign() << " to " << edge.contigId << " " << edge.neighborSign() << " " << seq << "\n";
//  }

  id2seq.erase(id);
//  semiCG.removeNode(id);
  updatePathStart(id, edge.baseSign(), edge.contigId, edge.neighborSign());
  mergeDonTouchSides(id, edge.baseSign(), edge.contigId, edge.neighborSign());
  semiCG.merge(id, edge.baseSign(), edge.contigId, edge.neighborSign(), true);
}

/**
 * Check 3 cases
 * 1. node is not corner node
 * 2. the neighbor we want to merge the node with is not corner node
 * 3. the current node is not the same as the node it will be merged with
 */
bool GFAConverter::isCornerCase(pufg::Node& n, bool mergeLeft) {
  if (mergeLeft) {
    if (dont_touch(n.getId(), Side::LEFT))
      return true;
    pufg::edgetuple& edge = n.getOnlyRealIn();
    if (edge.contigId == n.getId())
      return true;
    pufg::Node& neighbor = semiCG.getVertices()[edge.contigId];
    // if ( (edge.baseSign and edge.neighborSign) or (!edge.baseSign and !
    // edge.neighborSign) )
    if (edge.baseSign() == edge.neighborSign()) {
      return dont_touch(edge.contigId, Side::RIGHT) or neighbor.getRealOutdeg() != 1;
    } else
      return dont_touch(edge.contigId, Side::LEFT) or neighbor.getRealIndeg() != 1;
  } else { // Merge out case
    if (dont_touch(n.getId(), Side::RIGHT))
      return true;
    pufg::edgetuple& edge = n.getOnlyRealOut();
    if (edge.contigId == n.getId())
      return true;
    pufg::Node& neighbor = semiCG.getVertices()[edge.contigId];
    // if ( (edge.baseSign and edge.neighborSign) or (!edge.baseSign and !
    // edge.neighborSign) ) {
    if (edge.baseSign() == edge.neighborSign()) {
      return dont_touch(edge.contigId, Side::LEFT) or neighbor.getRealIndeg() != 1;
    } else {
      return dont_touch(edge.contigId, Side::RIGHT) or neighbor.getRealOutdeg() != 1;
    }
  }
  return false;
}

bool GFAConverter::dont_touch(uint64_t& nodeId, Side side) {
	return donTouch.contains(nodeId) and (donTouch[nodeId] == side or donTouch[nodeId] == Side::LEFT_RIGHT);
}

void GFAConverter::reconstructPathAndWrite(const char* gfaOutFileName, const char* refFileName) {
	//	exit(1);
		std::cerr << "Start writing segments into the file...\n";
	std::ofstream gfa_out(gfaOutFileName);
	std::vector<std::string> ref_in = {refFileName};
	uint64_t contigCntr = 0;
  	for (auto& kv : id2seq) {
    	gfa_out << "S"
             << "\t" << kv.first << "\t" << kv.second << "\n";
    	contigCntr++;
  	}
	std::cerr << "# of segments written to file : " << contigCntr << "\n\n";

	std::cerr << "Start writing paths into the file ...\n";
	spp::sparse_hash_map<uint64_t, pufg::Node>& nodes = semiCG.getVertices();	
	fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(ref_in, 1, 1);
    parser.start();
    auto rg = parser.getReadGroup();
    while (parser.refill(rg)) {
      for (auto& rp : rg) {

    	gfa_out << "P\t" << rp.name << "\t";

		auto& seq = rp.seq;	
		if (!path2start.contains(rp.name)) {
				std::cerr << "THIS IS EVEN WORSE!!\n path: " << rp.name << " not in path2start map\n";
		}
		auto& segment = path2start[rp.name];		
		uint64_t& id = segment.first;
		bool& ori = segment.second;
		if (!id2seq.contains(id)) {
			std::cerr << "THIS IS VERY BAD!! UNACCEPTABLE!!!! couldn't find id " << id << " in id2seq map\n";
			exit(1);
		}		
		gfa_out << id << (ori ? "+" : "-");

	//	std::cerr << "\nStart Node:\n " ;
	//	std::cerr << rp.name << " -->\n" << id << " " << ori << ":" << id2seq[id] << "\n";
		auto i = id2seq[id].size();
		while (i < seq.size()) {
			bool isFound = false;
			gfa_out << ",";
			pufg::Node& curNode = nodes[id];			
			auto& succs = curNode.getSuccessors();
			for (auto& succ : succs) {
					auto& neighborSeq = id2seq[succ.contigId];
	//				std::cerr << ori << " " << succ.baseSign() << " " << succ.contigId << " " << succ.neighborSign() << " at index " << i << " " << seq[i] << " " << neighborSeq << "\n";
				if (succ.baseSign() == ori) {
						if (succ.neighborSign() and neighborSeq[k-1] == seq[i]) {								 
							isFound = true;
						} else if (!succ.neighborSign()) {
							uint64_t idx = neighborSeq.size() - k;
							if (  (neighborSeq[idx] == 'T' and seq[i] == 'A')
								or (neighborSeq[idx] == 'A' and seq[i] == 'T')
								or (neighborSeq[idx] == 'C' and seq[i] == 'G')
								or (neighborSeq[idx] == 'G' and seq[i] == 'C')
							   ) {
								isFound = true;
							}
						}
						
						if (isFound) {
							id = succ.contigId;
							ori = succ.neighborSign();
	//						std::cerr << "Found\n";
							break;						
						}
				}
			}
			if (!isFound) {
				std::cerr << "Couldn't find any successors that match current nucleotide \n";
				std::cerr << "\tPath: " << seq << "\n";
				std::cerr << "\tNucleotide at " << i << " which was " << seq[i] << "\n";
				exit(1);
			}
			gfa_out << id << (ori ? "+" : "-");
			i += (id2seq[id].size()-(k-1));			
		}
		assert(i==seq.size());
	  	gfa_out << "\n";
      }
    }
	std::cerr << "Done writing paths into the file...\n\n";
}


