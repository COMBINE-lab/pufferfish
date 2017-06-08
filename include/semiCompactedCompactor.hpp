#include <map>
#include "sparsepp/spp.h"
#include "FatPufferGraph.hpp"

class SemiCompactedCompactor {

private:
		short k_ = 31;
		std::vector<std::pair<std::string, std::string> > newSegments;
		spp::sparse_hash_map<std::string, std::string> contigid2seq;
		spp::sparse_hash_map<std::string, bool > right_clipped;
		spp::sparse_hash_map<std::string, bool > left_clipped;
		spp::sparse_hash_map<std::string, bool > pathStart;
		spp::sparse_hash_map<std::string, bool > pathEnd;
		spp::sparse_hash_map<std::string, std::vector<std::pair<std::string, bool> > > path;
		pufg::Graph semiGraph;

		struct cmpByPair {
				    bool operator() (std::pair<std::string, bool> a, std::pair<std::string, bool> b) const {
									return (a.first != b.first) ? (a.first < b.first) : a.second;
									    }
		};
		//spp::sparse_hash_map<std::pair<std::string, bool>, std::pair<std::string, bool>, cmpByPair > needsNewNext;
		std::map<std::pair<std::string, bool>, std::pair<std::string, bool>, cmpByPair > needsNewNext;
		spp::sparse_hash_map<std::string, bool> shortContigs;

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
								pufg::Graph& semiGraphIn, short k) {
			newSegments = newSegmentsIn;
			contigid2seq = contigid2seqIn;
			pathStart = pathStartIn;
			pathEnd = pathEndIn;
			semiGraph = semiGraphIn;
			k_ = k;
		}

		void compact() {			
			spp::sparse_hash_map<std::string, pufg::Node> vertices = semiGraph.getVertices();
			uint32_t newSegmentCntr = 0;
			uint32_t boundaryCaseCntr = 0;
			for (auto & idSeq : newSegments) {
				pufg::Node s = vertices[idSeq.first];
				bool clipOuts = false, clipIns = false;
				// realIndegs : sum of Indegs to s+ and outdegs from s-
				// realOutdegs : sum of outdegs from s+ and indegs to s-
				if ((s.getRealIndeg() > 1 and s.getRealOutdeg() > 1)) // it is a complex k-mer seg!! which should be added to the path at the end
//								or (s.getRealIndeg() == 1 and is_boundary(s.getOut())) 
//								or (s.getRealOutdeg() == 1 and is_boundary(s.getIn()))
//				   ) 
				   {
				/*		   std::cerr << s.getId() << "\n\tIn:\n";
						   for (auto & n : s.getIn())
								std::cerr << n.contigId << "-" << n.baseSign << "-" << n.neighborSign << ", ";
						   std::cerr << "\n\tOut:\n";
						   for (auto & n : s.getOut())
								std::cerr << n.contigId << "-" << n.baseSign << "-" << n.neighborSign << ", ";
						   std::cerr << "\n";
				*/	
					clipOuts = true;
					clipIns = true;			
					// add it to the list of contigIds (segments) as a new segment
					//std::cerr << idSeq.first << idSeq.second << "\n";
					contigid2seq[idSeq.first] = idSeq.second;
					newSegmentCntr++;
					for (auto & n : s.getIn()) {// keep all the incoming nodes to this node in a map of the id to the new segment id
						needsNewNext[std::make_pair(n.contigId, n.neighborSign)] = std::make_pair(s.getId(), n.baseSign);				
					}

				} 		
				else if (s.getRealIndeg() == 1) clipOuts= true;
				else if (s.getRealOutdeg() == 1) clipIns = true;
				if ( (s.getRealIndeg() == 1 and is_boundary(s.getOut())) or (s.getRealOutdeg() == 1 and is_boundary(s.getIn())) ) {
					boundaryCaseCntr++;
				}
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
							if (contigid2seq[n.contigId].size() < k_) {
									shortContigs[n.contigId] = true;
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
								if (contigid2seq[n.contigId].size() < k_)
										shortContigs[n.contigId] = true;
						}
					}
				}
				/*
				 * 1 --> 0
				 * 2 --> 0
				 * 3 --> 0
				 */
				if (clipIns) { // just clip off from s+ incoming nodes and s- outgoing nodes
					for (auto & n : s.getIn()) { // incoming nodes
						if (n.baseSign) { // incoming nodes to s+
							if (n.neighborSign and right_clipped.find(n.contigId) == right_clipped.end()) {// n+ --> s+
									contigid2seq[n.contigId].erase(contigid2seq[n.contigId].size()-1, 1);//last nucleotide
									right_clipped[n.contigId] = true;
							}
							else if (!n.neighborSign and left_clipped.find(n.contigId) == left_clipped.end()) {// n- --> s+
								contigid2seq[n.contigId].erase(0, 1); //first nucleotide
								left_clipped[n.contigId] = true;
							}
							if (contigid2seq[n.contigId].size() < k_)
									shortContigs[n.contigId] = true;						
						}
					}
					for (auto & n : s.getOut()) { // outgoing nodes
						if (!n.baseSign) { // outgoing nodes from s-
							if (n.neighborSign and left_clipped.find(n.contigId) == left_clipped.end()) {// s- --> n+
								contigid2seq[n.contigId].erase(0, 1); //first nucleotide
								left_clipped[n.contigId] = true;
							}
							else if (!n.neighborSign and right_clipped.find(n.contigId) == right_clipped.end()) {// s- --> n-
								contigid2seq[n.contigId].erase(contigid2seq[n.contigId].size()-1, 1); //last nucleotide
								right_clipped[n.contigId] = true;
							}
							if (contigid2seq[n.contigId].size() < k_)
									shortContigs[n.contigId] = true;
							
						}
					}
				}
			}
			std::cerr << "# of Added Segments : " << newSegmentCntr << "\n# of Short Contigs : " << shortContigs.size() << "\n# of Boundary Cases : " << boundaryCaseCntr << "\n";			
		}

		void updatePath(spp::sparse_hash_map<std::string, std::vector<std::pair<std::string, bool> > >&  pathIn) {
			path = pathIn;
			for (auto & kv : path) {
					std::string id = kv.first;
					std::vector<std::pair<std::string, bool> > & pathOfContigs = kv.second;
					std::vector<std::pair<std::string, bool> > updatedPath;
					for (auto & c : pathOfContigs) {
						if (shortContigs.find(c.first) == shortContigs.end())
							updatedPath.push_back(c);
						if (needsNewNext.find(c) != needsNewNext.end())
								updatedPath.push_back(needsNewNext[c]);
					}

					path[id] = updatedPath;
			}
		}
		
void writeFile(std::string fileName){
	std::ofstream gfa_file(fileName) ;
	for(auto& cseq : contigid2seq){
		gfa_file << "S" << "\t" << cseq.first <<"\t" << cseq.second << "\n" ;
	}
	for(auto& p : path){
		auto tid = p.first ;
		gfa_file << "P" << "\t" << tid << "\t"  ;
		auto vec = p.second ;
		for(int i = 0 ; i < vec.size()-1 ; i++){
			gfa_file << vec[i].first << ((vec[i].second)?"+":"-") << "," ;
		}
		gfa_file << vec[vec.size()-1].first << ((vec[vec.size()-1].second)?"+":"-") << "\t*\n";

	}
}

};
