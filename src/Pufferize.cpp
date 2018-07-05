#include "FastxParser.hpp"
#include "Pufferize.hpp"
#include "cereal/archives/binary.hpp"
#include "xxhash.h"
#include <algorithm>
#include <assert.h>

Pufferize::Pufferize(const char *gfaFileName, size_t input_k) {
    gfaFileName_ = gfaFileName;
    std::cerr << "Reading GFA file " << gfaFileName << "\n";
    file.reset(new zstr::ifstream(gfaFileName));
    k = input_k;
}

void Pufferize::encodeSeq(sdsl::int_vector<2> &contig, std::string &str) {
    for (size_t i = 0; i < str.length(); ++i) {
        auto c = kmers::codeForChar(str[i]);
        contig[i] = c;
    }
}

std::string Pufferize::convertSubview2str(SubseqView ss) {
    std::string res;
    auto &contig = contigs[ss.contigId];
    assert(contig.size() > ss.startp);
    assert(contig.size() >= ss.len - ss.startp);
    for (size_t i = ss.startp; i < ss.len; ++i) {
        auto c = kmers::charForCode(contig[i]);
        res += c;
    }
    return res;
}

void Pufferize::parseFile() {
    std::string ln;
    std::string tag, id, value;
    size_t contig_cnt{0};
    std::cerr << "Start loading segments... \n";
    spp::sparse_hash_map<uint64_t, uint64_t> seq2id;
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
        if (contig_cnt % 1000000 == 0) {
            std::cerr << contig_cnt << " segments processed\n";
        }
    }
    std::cerr << "Done loading & processing segments ..."
              << ksizeContig.size() << " out of " << idCntr << " contigs were of size k\n";
    std::cerr << "\nStart loading paths\n";

    file.reset(new zstr::ifstream(gfaFileName_));
    uint64_t pathCntr{0};
    while (std::getline(*file, ln)) {
        char firstC = ln[0];
        if (firstC != 'P')
            continue;
        pathCntr++;
        stx::string_view lnview(ln);
        std::vector<stx::string_view> splited = util::split(lnview, '\t');
        if (splited.size() != 4) {
            std::cerr << "the path line should have 4 fields (tab separated) --- skipping!\n";
            continue;
        }
        // tag = splited[0].to_string();
        id = splited[1].to_string();
        if (id == "") {
            std::cerr << "interesting; the ID is empty.  Skipping\n";
            continue;
        }
        // A path line
        auto pvalue = splited[2];
        // parse value and add all conitgs to contigVec
        std::vector<std::pair<uint64_t, bool>> contigVec =
                util::explode(pvalue, ',');
        //std::cerr << "p " << pathCntr << " , contig count: " << contigVec.size() << "\n";
        contigVec = convertOld2newPath(contigVec);

        // update graph and add the new set of segments to it
        updateGraph(contigVec);

        //std::cerr << "updateGraph done\n";
        // set pathStart and pathEnd
        setPathStart(id, contigVec[0]);
        //std::cerr << "setpathstart done\n";
        setDonTouchSide(contigVec);
        //std::cerr << "setDon'tTouchSide done\n";
    }

    std::cerr << "Done Reading gfa file\n";

    std::cerr << "Total # of Contigs : " << contig_cnt
              << "\tTotal # of numerical Contigs : " << old2newids.size()
              << " \tTotal # of contigs after spliting :" << id2seq.size()
              << "\tTotal # of paths:" << pathCntr
              << "\t# of path starts:" << start2path.size()
              << "\n\n";
}

void Pufferize::processContigSeq(
        uint64_t &contigId, std::string &contigSeq,
        spp::sparse_hash_map<uint64_t, uint64_t> &seq2id,
        uint64_t &idCntr) {
    // Add the new contig to the list of 2-bit encoded contigs
    sdsl::int_vector<2> newContig(contigSeq.size());
    encodeSeq(newContig, contigSeq);
    contigs.push_back(newContig);
    /**
     * Divide every segment into 3 pieces
     * first k
     * middle one without the first and last nucleotide
     * last k
     * special case is segment with length = k+1 that doesn't have the middle case
     * ATTENTION : The order we insert these 3 pieces into the vector matters. So
     *keep it the way it is
    **/
    //std::cerr << contigSeq.size() << "\n";
    std::vector<SubseqView> seqParts;
    if (util::isRevcomp(contigSeq)) {
        for (size_t i = 0; i <= contigSeq.size() - k; i++) {
            seqParts.push_back(SubseqView(contigId, i, k));
        }
    } else {
        seqParts.push_back(SubseqView(contigId, 0, k));
        if (contigSeq.size() > k + 1) {
            seqParts.push_back(SubseqView(contigId, 1, contigSeq.size() - 2));
        }
        seqParts.push_back(SubseqView(contigId, contigSeq.size() - k, k));
    }

    // for each part, search in hash_map whether it already existed or not and
    // assign the proper (new/retrieved) contig id to it
    // int cntr = 1;
    for (SubseqView seq : seqParts) {
        uint64_t newId;
        bool plus = true;
        CanonicalKmer newKmer;
        newKmer.k(k);
        if (seq.len == k) {
            newKmer.fromStr(contigSeq.substr(seq.startp, seq.len));
            /*std::cerr << idCntr << " size: " << seq2id.size()
                      << " seq:" << contigSeq.substr(seq.startp, seq.len) << " fw:"
                      << newKmer.fwWord() << " rc:" << newKmer.rcWord()
                      << "\n";*/
            if (seq2id.contains(newKmer.fwWord())) {
                newId = seq2id[newKmer.fwWord()];
            } else if (seq2id.contains(newKmer.rcWord())) {
                plus = false;
                newId = seq2id[newKmer.rcWord()];
            } else {
                newId = idCntr;
                seq2id[newKmer.fwWord()] = newId;
                idCntr++;
                ksizeContig.push_back(newId);
            }
        } else {
            newId = idCntr;
            idCntr++;
        }
        id2seq[newId] = seq;
        old2newids[contigId].push_back(ContigOri(newId, plus));
    }
}

std::vector<std::pair<uint64_t, bool> > &
Pufferize::convertOld2newPath(std::vector<std::pair<uint64_t, bool> > &oContigVec) {
    std::vector<std::pair<uint64_t, bool> > nContigVec = {};
    bool isPathStart = true;
    for (auto &idOri : oContigVec) {
        auto newIdList = old2newids[idOri.first];
        //std::cerr << "\tnewlist size:"<<newIdList.size() <<"\n";
        if (!idOri.second) {
            // As newIdList is a copy it won't affect the actual vector in old2newids (tested)
            std::reverse(newIdList.begin(), newIdList.end());
        }
        if (isPathStart) {
            nContigVec.emplace_back(newIdList[0].id(), idOri.second ? newIdList[0].ori() : !newIdList[0].ori());
            isPathStart = false;
        }
        for (size_t i = 1; i < newIdList.size(); i++) {
            nContigVec.emplace_back(newIdList[i].id(), idOri.second ? newIdList[i].ori() : !newIdList[i].ori());
        }
    }
    oContigVec = nContigVec;
    return oContigVec;
}

void Pufferize::updateDonTouchSide(uint64_t &id, Side side) {
    if (!donTouch.contains(id))
        donTouch[id] = side;
    else if (donTouch[id] != side)
        donTouch[id] = Side::LEFT_RIGHT;
//		if (id == 70547 or id == 70550)
//				std::cerr << id << " " << 
//						(donTouch[id] == Side::LEFT?"left":(donTouch[id] == Side::RIGHT?"right":"both")) << "\n";
}

void Pufferize::setDonTouchSide(std::vector<std::pair<uint64_t, bool> > &contigVec) {
    if (contigVec[0].second)
        updateDonTouchSide(contigVec[0].first, Side::LEFT);
    else
        updateDonTouchSide(contigVec[0].first, Side::RIGHT);

    if (contigVec[contigVec.size() - 1].second)
        updateDonTouchSide(contigVec[contigVec.size() - 1].first, Side::RIGHT);
    else
        updateDonTouchSide(contigVec[contigVec.size() - 1].first, Side::LEFT);
}

void Pufferize::mergeDonTouchSides(uint64_t &mergedId, bool mergedOri, uint64_t &mergingId, bool mergingOri) {
    if (donTouch.contains(mergedId)) {
        Side s;
        if (mergedOri == mergingOri) {
            s = donTouch[mergedId];
        } else {
            switch (donTouch[mergedId]) {
                case Side::LEFT:
                    s = Side::RIGHT;
                    break;
                case Side::RIGHT:
                    s = Side::LEFT;
                    break;
                default:
                    s = donTouch[mergedId];
            }
        }
        updateDonTouchSide(mergingId, s);
    }
}

void Pufferize::setPathStart(std::string &pathId, std::pair<uint64_t, bool> &startNode) {
//		if (startNode.first == 70550) 
//			std::cerr << "path start " << startNode.first << " " << startNode.second << "\n";

    ContigOri co(startNode.first, startNode.second);
    path2start[pathId] = co;
    if (start2path.find(co.wrd) == start2path.end())
        start2path[co.wrd] = {};
    start2path[co.wrd].push_back(pathId);
    if (co.id() == 234251) {
        std::cerr << "lets see " << pathId << " " << co.ori() << "\n";
    }
}

void
Pufferize::updatePathStart(uint64_t &oStartNodeId, bool oStartNodeOri, uint64_t &nStartNodeId, bool nStartNodeOri) {
    //std::pair<uint64_t, bool> oStartNode = std::make_pair(oStartNodeId, oStartNodeOri);
    //std::pair<uint64_t, bool> nStartNode = std::make_pair(nStartNodeId, nStartNodeOri);
    {
        ContigOri oStartNode(oStartNodeId, oStartNodeOri);
        ContigOri nStartNode(nStartNodeId, nStartNodeOri);
       /* std::cerr << oStartNode.id() << " " << oStartNode.ori() << " "
                  << nStartNode.id() << " " << nStartNode.ori() << " ";*/
        if (start2path.find(oStartNode.wrd) != start2path.end()) {
            auto &oldP = start2path[oStartNode.wrd];
            //std::cerr << "s" << oldP.size() << " ";
            if (start2path.find(nStartNode.wrd) == start2path.end()) {
                start2path[nStartNode.wrd] = {};

            }
            for (auto &v : oldP) {
//			if (nStartNode.first == 70550) 
//			std::cerr << "1 path start from " <<  oStartNode.first << " " << oStartNode.second << " to " << nStartNode.first << " " << nStartNode.second << "\n";
                //std::cerr << v << " ";
                if (oStartNode.id() == 234251) {
                    std::cerr << "using " << v << "\n";
                }
                path2start[v] = nStartNode;
                start2path[nStartNode.wrd].push_back(v);
                if (nStartNode.id() == 234251) {
                    std::cerr << "lets see " << v << "\n";
                }
            }
            //std::cerr << "e ";
        }
    }

//    oStartNode = std::make_pair(oStartNodeId, !oStartNodeOri);
//    nStartNode = std::make_pair(nStartNodeId, !nStartNodeOri);
    {
        ContigOri oStartNode(oStartNodeId, !oStartNodeOri);
        ContigOri nStartNode(nStartNodeId, !nStartNodeOri);
        /*std::cerr << oStartNode.id() << " " << oStartNode.ori() << " "
                  << nStartNode.id() << " " << nStartNode.ori() << " ";*/

        if (start2path.find(oStartNode.wrd) != start2path.end()) {
            auto &oldP = start2path[oStartNode.wrd];
            //std::cerr << "s" << oldP.size();
            if (start2path.find(nStartNode.wrd) == start2path.end())
                start2path[nStartNode.wrd] = {};
            for (auto &v : oldP) {
                if (oStartNode.id() == 234251) {
                    std::cerr << "2 using " << v << "\n";
                }
//			if (nStartNode.first == 70550) 
//			std::cerr << "2 path start " << nStartNode.first << " " << nStartNode.second << "\n";
                path2start[v] = nStartNode;
                start2path[nStartNode.wrd].push_back(v);
                if (nStartNode.id() == 234251) {
                    std::cerr << "lets see " << v << "\n";
                }
            }
        }
        //std::cerr << "\n";
    }

}

void Pufferize::updateGraph(std::vector<std::pair<uint64_t, bool> > &contigVec) {
    for (size_t i = 0; i < contigVec.size() - 1; i++) {
        /*std::cerr << "\t" << i << " of " << contigVec.size() << " " << contigVec[i].first << ","
                  << contigVec[i + 1].first << "\n";*/
        semiCG.addEdge(contigVec[i].first, contigVec[i].second, contigVec[i + 1].first, contigVec[i + 1].second);
    }
}

void Pufferize::randomWalk() {
    spp::sparse_hash_map<uint64_t, pufg::Node> &nodes = semiCG.getVertices();
    std::cerr << "# of contigs with length = " << k << " : " << ksizeContig.size() << "\n";
    std::cerr << "\nStart merging .. \n";
    for (auto &v : ksizeContig) {
        uint64_t curId = v;
        pufg::Node &curNode = nodes[curId];
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

void Pufferize::mergeLeft(pufg::Node &n) {
    uint64_t id = n.getId();
    pufg::edgetuple &edge = n.getOnlyRealIn();
    if (!id2seq.contains(id)) {
        std::cerr << "[mergeLeft] NO; the id " << id
                  << " was not in id2seq!\n";
    }
    if (!id2seq.contains(edge.contigId)) {
        std::cerr << "[mergeLeft] NO; the edge.contigId " << edge.contigId
                  << " was not in id2seq!\n";
    }
    auto &seqView = id2seq[edge.contigId];
    seqView.len++;
    /*
//  std::cerr << "mergeLeft " << id << " " << edge.contigId << "\n";
    auto &tobeMerged = id2seq[id];
    auto &seq = id2seq[edge.contigId];
    if (edge.baseSign() != edge.neighborSign()) {
        tobeMerged = util::revcomp(tobeMerged);
        //				if (tobeMerged.substr(tobeMerged.size()-(k-1)) != seq.substr(0,
        //k-1)) std::cerr << "1 " << id << " " << edge.contigId << " " << seq <<
        //"\n" << tobeMerged << "\n";
        seq = tobeMerged.substr(0, tobeMerged.size() - (k - 1)) + seq;
    } else {
        //			if (tobeMerged.substr(0, k-1) != seq.substr(seq.size() - (k-1)))
        std::cerr << "2 " << seq << "\n" << tobeMerged << "\n";
        //seq += tobeMerged.substr(k - 1);
    }
//  if (edge.contigId == 70550 or edge.contigId == 70547) {
//  	std::cerr << "MergeLeft " << id << " " << edge.baseSign() << " to " << edge.contigId << " " << edge.neighborSign() << " " << seq << "\n";
//  }
     */
    id2seq.erase(id);
    //semiCG.removeNode(id);
//  if (edge.contigId == 70550) std::cerr <<"Ajabaa " <<edge.contigId << " " <<edge.neighborSign() << "\n";
    updatePathStart(id, edge.baseSign(), edge.contigId, edge.neighborSign());
    mergeDonTouchSides(id, edge.baseSign(), edge.contigId, edge.neighborSign());
    semiCG.merge(id, edge.baseSign(), edge.contigId, edge.neighborSign(), false);
}

void Pufferize::mergeRight(pufg::Node &n) {
    uint64_t id = n.getId();
    pufg::edgetuple &edge = n.getOnlyRealOut();
    if (!id2seq.contains(id)) {
        std::cerr << "[mergeRight] NO; the id " << id
                  << " was not in id2seq!\n";
    }
    if (!id2seq.contains(edge.contigId)) {
        std::cerr << "[mergeRight] NO; the edge.contigId " << edge.contigId
                  << " was not in id2seq!\n";
    }
    auto &seqView = id2seq[edge.contigId];
    seqView.len++;
    seqView.startp--;
    /*
    //std::cerr << "mergeRight " << id << " " << edge.contigId << " ";
    std::string &tobeMerged = id2seq[id];
    std::string &seq = id2seq[edge.contigId];
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
*/
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
bool Pufferize::isCornerCase(pufg::Node &n, bool mergeLeft) {
    if (mergeLeft) {
        if (dont_touch(n.getId(), Side::LEFT))
            return true;
        pufg::edgetuple &edge = n.getOnlyRealIn();
        if (edge.contigId == n.getId())
            return true;
        pufg::Node &neighbor = semiCG.getVertices()[edge.contigId];
        // if ( (edge.baseSign and edge.neighborSign) or (!edge.baseSign and !
        // edge.neighborSign) )
        if (edge.baseSign() == edge.neighborSign()) {
            return dont_touch(edge.contigId, Side::RIGHT) or neighbor.getRealOutdeg() != 1;
        } else
            return dont_touch(edge.contigId, Side::LEFT) or neighbor.getRealIndeg() != 1;
    } else { // Merge out case
        if (dont_touch(n.getId(), Side::RIGHT))
            return true;
        pufg::edgetuple &edge = n.getOnlyRealOut();
        if (edge.contigId == n.getId())
            return true;
        pufg::Node &neighbor = semiCG.getVertices()[edge.contigId];
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

bool Pufferize::dont_touch(uint64_t &nodeId, Side side) {
    return donTouch.contains(nodeId) and (donTouch[nodeId] == side or donTouch[nodeId] == Side::LEFT_RIGHT);
}

void Pufferize::reconstructPathAndWrite(const char *gfaOutFileName, const char *refFileName) {
    //	exit(1);
    std::cerr << "Start writing segments into the file...\n";
    std::ofstream gfa_out(gfaOutFileName);
    std::vector<std::string> ref_in = {refFileName};
    uint64_t contigCntr = 0;
    for (auto &kv : id2seq) {
        SubseqView &ss = kv.second;
        gfa_out << "S"
                << "\t" << kv.first << "\t" << convertSubview2str(ss) << "\n";
        contigCntr++;
    }
    std::cerr << "# of segments written to file : " << contigCntr << "\n\n";

    std::cerr << "Start writing paths into the file ...\n";
    spp::sparse_hash_map<uint64_t, pufg::Node> &nodes = semiCG.getVertices();
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(ref_in, 1, 1);
    parser.start();
    auto rg = parser.getReadGroup();
    while (parser.refill(rg)) {
        for (auto &rp : rg) {

            gfa_out << "P\t" << rp.name << "\t";

            auto &seq = rp.seq;
            if (!path2start.contains(rp.name)) {
                std::cerr << "THIS IS EVEN WORSE!!\n path: " << rp.name << " not in path2start map\n";
            }
            auto &segment = path2start[rp.name];
            uint64_t id = segment.id();
            bool ori = segment.ori();
            if (!id2seq.contains(id)) {
                std::cerr << "THIS IS VERY BAD!! UNACCEPTABLE!!!! couldn't find id " << id << " in id2seq map\n";
                exit(1);
            }
            gfa_out << id << (ori ? "+" : "-");

            //	std::cerr << "\nStart Node:\n " ;
            //	std::cerr << rp.name << " -->\n" << id << " " << ori << ":" << id2seq[id] << "\n";
            auto i = id2seq[id].len;
            while (i < seq.size()) {
                bool isFound = false;
                gfa_out << ",";
                pufg::Node &curNode = nodes[id];
                auto &succs = curNode.getSuccessors();
                for (auto &succ : succs) {
                    auto &neighborSeq = id2seq[succ.contigId];
                    auto &contig = contigs[neighborSeq.contigId];
                    //				std::cerr << ori << " " << succ.baseSign() << " " << succ.contigId << " " << succ.neighborSign() << " at index " << i << " " << seq[i] << " " << neighborSeq << "\n";
                    if (succ.baseSign() == ori) {
                        assert(contig.size() > (neighborSeq.startp + k - 1));
                        if (succ.neighborSign() and
                            contig[neighborSeq.startp + k - 1] == combinelib::kmers::codeForChar(seq[i])) {
                            isFound = true;
                        } else if (!succ.neighborSign()) {
                            uint64_t idx = neighborSeq.startp + neighborSeq.len - k;
                            assert(contig.size() > idx);
                            if (contig[idx] == combinelib::kmers::codeForChar(combinelib::kmers::complement(seq[i])))
                                /*if ((neighborSeq[idx] == 'T' and seq[i] == 'A')
                                    or (neighborSeq[idx] == 'A' and seq[i] == 'T')
                                    or (neighborSeq[idx] == 'C' and seq[i] == 'G')
                                    or (neighborSeq[idx] == 'G' and seq[i] == 'C')
                                        ) */
                            {
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
                i += ((id2seq[id].len - id2seq[id].startp) - (k - 1));
            }
            assert(i == seq.size());
            gfa_out << "\n";
        }
    }
    std::cerr << "Done writing paths into the file...\n\n";
}


