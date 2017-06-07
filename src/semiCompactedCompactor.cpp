#include <map>
#include "sparsepp/spp.h"
#include "FatPufferGraph.hpp"

class SemiCompactedCompactor {

private:
		std::vector<std::pair<std::string, std::string> > newSegments;
		spp::sparse_hash_map<std::string, std::string> contigid2seq;
		spp::sparse_hash_map<std::string, bool > right_clipped;
		spp::sparse_hash_map<std::string, bool > left_clipped;
		spp::sparse_hash_map<std::string, bool > pathStart;
		spp::sparse_hash_map<std::string, bool > pathEnd;

		pufg::Graph semiGraph;

		struct cmpByPair {
				    bool operator() (std::pair<std::string, bool> a, std::pair<std::string, bool> b) const {
									return (a.first != b.first) ? (a.first < b.first) : a.second;
									    }
		};
		//spp::sparse_hash_map<std::pair<std::string, bool>, std::pair<std::string, bool>, cmpByPair > needsNewNext;
		std::map<std::pair<std::string, bool>, std::pair<std::string, bool>, cmpByPair > needsNewNext;
		bool is_boundary(std::vector<pufg::edgetuple> nodes) {
			for (auto & n : nodes) {
				if (n.neighborSign) {
					if (pathStart.find(n.contigId) != pathStart.end()) return true;
				} else if (pathEnd.find(n.contigId) != pathEnd.end()) return true;
			}
			return false;
		} 
public:
		SemiCompactedCompactor(std::vector<std::pair<std::string, std::string> >& newSegmentsIn, 
						spp::sparse_hash_map<std::string, std::string>& contigid2seqIn,
							spp::sparse_hash_map<std::string, bool >& pathStartIn,
								spp::sparse_hash_map<std::string, bool >& pathEndIn,
								pufg::Graph& semiGraphIn) {
			newSegments = newSegmentsIn;
			contigid2seq = contigid2seqIn;
			pathStart = pathStartIn;
			pathEnd = pathEndIn;
			semiGraph = semiGraphIn;
		}

		void compact() {			
			spp::sparse_hash_map<std::string, pufg::Node> vertices = semiGraph.getVertices();
			for (auto & idSeq : newSegments) {
				pufg::Node s = vertices[idSeq.first];
				bool clipOuts = false, clipIns = false;
				// realIndegs : sum of Indegs to s+ and outdegs from s-
				// realOutdegs : sum of outdegs from s+ and indegs to s-
				if ((s.getRealIndeg() > 1 and s.getRealOutdeg() > 1) // it is a complex k-mer seg!! which should be added to the path at the end
								or (s.getRealIndeg() == 1 and is_boundary(s.getOut())) 
								or (s.getRealOutdeg() == 1 and is_boundary(s.getIn()))
				   ) {
					clipOuts = true;
					clipIns = true;			
					// add it to the list of contigIds (segments) as a new segment
					contigid2seq[idSeq.first] = idSeq.second;
					for (auto & n : s.getIn()) // keep all the incoming nodes to this node in a map of the id to the new segment id
						needsNewNext[std::make_pair(n.contigId, n.baseSign)] = std::make_pair(s.getId(), n.neighborSign);				
				} 		
				else if (s.getRealIndeg() == 1) clipOuts= true;
				else if (s.getRealOutdeg() == 1) clipIns = true;
				/*
				 * 0 ---> 1
				 * 0 ---> 2
				 * 0 ---> 3
				 */
				if (clipOuts) { // just clip off from s+ outgoing nodes and s- incoming nodes
					for (auto & n : s.getOut()) { // outgoing nodes
						if (n.baseSign) { // outgoing nodes from s+
							if (n.neighborSign and left_clipped.find(n.contigId) == left_clipped.end()) {// s+ --> n+
									contigid2seq[n.contigId].erase(0,1);
									left_clipped[n.contigId] = true;
							}
							else if (!n.neighborSign and right_clipped.find(n.contigId) == right_clipped.end()) {// s+ --> n-
								contigid2seq[n.contigId].erase(contigid2seq[n.contigId].size()-1, 1); //last nucleotide
								right_clipped[n.contigId] = true;
							}
						}
					}
					for (auto & n : s.getIn()) { // incoming nodes
						if (!n.baseSign) { // incoming nodes to s-
								if (n.neighborSign and right_clipped.find(n.contigId) == right_clipped.end()) {// n+ --> s-
									contigid2seq[n.contigId].erase(contigid2seq[n.contigId].size()-1, 1); //last nucleotide
									right_clipped[n.contigId] = true;
								}
								else if (!n.neighborSign and left_clipped.find(n.contigId) == left_clipped.end()) {// n- --> s-
									contigid2seq[n.contigId].erase(0, 1); //first nucleotide
									left_clipped[n.contigId] = true;
								}
						}
					}
				}
				/*
				 * 1 --> 0
				 * 2 --> 0
				 * 3 --> 0
				 */
				else if (clipIns) { // just clip off from s+ incoming nodes and s- outgoing nodes
					for (auto & n : s.getIn()) { // incoming nodes
						if (n.baseSign) { // incoming nodes to s+
							if (n.neighborSign and right_clipped.find(n.contigId) == right_clipped.end()) {// s+ --> n+
									contigid2seq[n.contigId].erase(contigid2seq[n.contigId].size()-1, 1);//last nucleotide
									right_clipped[n.contigId] = true;
							}
							else if (!n.neighborSign and left_clipped.find(n.contigId) == left_clipped.end()) {// s+ --> n-
								contigid2seq[n.contigId].erase(0, 1); //first nucleotide
								left_clipped[n.contigId] = true;
							}
						}
					}
					for (auto & n : s.getOut()) { // outgoing nodes
						if (!n.baseSign) { // outgoing nodes from s-
								if (n.neighborSign and left_clipped.find(n.contigId) == left_clipped.end()) {// n+ --> s-
									contigid2seq[n.contigId].erase(0, 1); //first nucleotide
									left_clipped[n.contigId] = true;
								}
								else if (!n.neighborSign and right_clipped.find(n.contigId) == right_clipped.end()) {// n- --> s-
									contigid2seq[n.contigId].erase(contigid2seq[n.contigId].size()-1, 1); //last nucleotide
									right_clipped[n.contigId] = true;
								}
						}
					}
				}
			}
		}

		void updatePath(spp::sparse_hash_map<std::string, std::vector<std::pair<std::string, bool> > >  path) {
			
			for (auto & kv : path) {
					std::string id = kv.first;
					std::vector<std::pair<std::string, bool> > & pathOfContigs = kv.second;
					std::vector<std::pair<std::string, bool> > updatedPath;
					for (auto & c : pathOfContigs) {
						updatedPath.push_back(c);
						if (needsNewNext.find(c) != needsNewNext.end())
								updatedPath.push_back(needsNewNext[c]);
					}

					path[id] = updatedPath;
			}
		}
};
