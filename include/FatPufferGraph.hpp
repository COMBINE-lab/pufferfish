#ifndef PUFFER_GRAPH_H
#define PUFFER_GRAPH_H

#include <map>
#include "Util.hpp"
#include "sparsepp/spp.h"
#include <iostream>

namespace pufg{

  enum class EdgeType : uint8_t { PLUS_PLUS, PLUS_MINUS, MINUS_PLUS, MINUS_MINUS };

	struct edgetuple{
		edgetuple(bool fSign,std::string cId, bool tSign):
      contigId(cId) {
      if (fSign and tSign) {
        t = EdgeType::PLUS_PLUS;
      } else if (fSign and !tSign) {
        t = EdgeType::PLUS_MINUS;
      } else if (!fSign and tSign) {
        t = EdgeType::MINUS_PLUS;
      } else {
        t = EdgeType::MINUS_MINUS;
      }
    }

		edgetuple() {}

    bool baseSign() {
      return (t == EdgeType::PLUS_PLUS or t == EdgeType::PLUS_MINUS);
    }
    bool neighborSign() {
      return (t == EdgeType::PLUS_PLUS or t == EdgeType::MINUS_PLUS);
    }

    EdgeType t;
		std::string contigId;
	};

	class Node {
	public:
		 	Node() {}
			Node(std::string idIn) {
				id = idIn ;
				indegp = 0 ;
				outdegp = 0 ;
				indegm = 0 ;
				outdegm = 0 ;
			}

			uint8_t getIndegP() {return indegp ;}
			uint8_t getOutdegP() {return outdegp ;}
			uint8_t getIndegM() {return indegm ;}
			uint8_t getOutdegM() {return outdegm ;}
			uint8_t getRealOutdeg() {					
					return distinctRealOut.size();}//(outdegp + outdegm) ;}
			uint8_t getRealIndeg() {					
					return distinctRealIn.size();}//(indegp + indegm); }

			edgetuple& getOnlyRealIn() {
				if (indegp > 0) {
					for (auto& e : in)
            if (e.baseSign()) {
									return e;
							}
				} else {
					for (auto& e : out)
            if (!e.baseSign()) {
									return e;
							}
				}
        // should not get here
        return in.front();
			}

			edgetuple& getOnlyRealOut() {
				if (outdegp > 0) {
					for (auto& e : out)
						if (e.baseSign()) {
							return e;
						}			
				}
				else { // The real outgoing edge should be an incoming edge to negative if it's not an outgoing edge from positive
					for (auto& e : in)
						if (!e.baseSign()) {
								return e;
						}
				}
        // should not get here
        return out.front();
			}
			//std::vector<std::string> getPositiveToNodes()
			std::string& getId() {return id;}

			void insertNodeTo(std::string nodeId, bool sign, bool toSign){
/*					if (id == "0032308305" or id == "0011206970")
							std::cerr << id << " insert node to " << nodeId << toSign << sign << " indegp:" << (int)indegp << " outdegp:" << (int)outdegp << " indegm:" << (int)indegm << " outdegm:" << (int)outdegm << " distinctIn: " << distinctRealIn.size() << " distinctOut: " << distinctRealOut.size() << " in:" << in.size() << " out:" << out.size() << "\n";
*/							
				if(sign) {
					outdegp++;
					std::string key = nodeId + (toSign?"+":"-");
					if (distinctRealOut.contains(key))
						distinctRealOut[key] += 1;
					else distinctRealOut[key] = 1;
				}
				else {
					outdegm++;
					std::string key = nodeId + (toSign?"-":"+");
					if (distinctRealIn.contains(key))
							distinctRealIn[key] += 1;
					else distinctRealIn[key] = 1;
				}

				out.emplace_back(sign, nodeId, toSign);
			}
			
			void removeEdgeTo(std::string nodeId) {
/*					if (id == "0032308305" or id == "0011206970")
							std::cerr << id << " remove edge to " << nodeId << " indegp:" << (int)indegp << " outdegp:" << (int)outdegp << " indegm:" << (int)indegm << " outdegm:" << (int)outdegm << " distinctIn: " << distinctRealIn.size() << " distinctOut: " << distinctRealOut.size() << " in:" << in.size() << " out:" << out.size() << "\n";
*/							
			for (std::vector<edgetuple>::iterator it=out.begin(); it!=out.end();) {
					auto& edge = *it;
					if (edge.contigId == nodeId) {
            if (edge.baseSign()) {
								outdegp--;
								std::string key = nodeId + (edge.neighborSign()?"+":"-");
								if (distinctRealOut.contains(key)) {
									distinctRealOut[key] -= 1;
									if (distinctRealOut[key] == 0)
											distinctRealOut.erase(key);
								}
							}
							else {
								outdegm--;
								std::string key = nodeId + (edge.neighborSign()?"-":"+");
								if (distinctRealIn.contains(key)) {
									distinctRealIn[key] -= 1;
									if (distinctRealIn[key] == 0)
											distinctRealIn.erase(key);
								}
						}
						it = out.erase(it);
					}
					else it++;
				}
			}

			void insertNodeFrom(std::string nodeId, bool sign, bool fromSign){
/*					if (id == "0032308305" or id == "0011206970")
							std::cerr << id << " insert node from " << nodeId << fromSign << sign << " indegp:" << (int)indegp << " outdegp:" << (int)outdegp << " indegm:" << (int)indegm << " outdegm:" << (int)outdegm << " distinctIn: " << distinctRealIn.size() << " distinctOut: " << distinctRealOut.size() << " in:" << in.size() << " out:" << out.size() << "\n";
*/							
			if(sign) {
					indegp++ ;
					std::string key = nodeId + (fromSign?"+":"-");
					if (distinctRealIn.contains(key))
						distinctRealIn[key] += 1;
					else distinctRealIn[key] = 1;
				}
				else {
					indegm++ ;
					std::string key = nodeId + (fromSign?"-":"+");
					if (distinctRealOut.contains(key))
							distinctRealOut[key] += 1;
					else distinctRealOut[key] = 1;
				}
				in.emplace_back(sign, nodeId, fromSign) ;

			}
			
			void removeEdgeFrom(std::string nodeId) {
/*					if (id == "0032308305" or id == "0011206970")
							std::cerr << id << " remove edge from " << nodeId << " indegp:" << indegp << " outdegp:" << outdegp << " indegm:" << indegm << " outdegm:" << outdegm << " distinctIn: " << distinctRealIn.size() << " distinctOut: " << distinctRealOut.size() << " in:" << in.size() << " out:" << out.size() << "\n";
*/							
			for (std::vector<edgetuple>::iterator it=in.begin(); it!=in.end();) {
						auto& edge = *it;
					if (edge.contigId == nodeId) {
            if (edge.baseSign()) {
								indegp--;
								std::string key = nodeId + (edge.neighborSign()?"+":"-");
								if (distinctRealIn.contains(key)) {
									distinctRealIn[key] -= 1;
									if (distinctRealIn[key] == 0)
											distinctRealIn.erase(key);
								}
							}
							else {
								indegm--;
								std::string key = nodeId + (edge.neighborSign()?"-":"+");
								if (distinctRealOut.contains(key)) {
									distinctRealOut[key] -= 1;
									if (distinctRealOut[key] == 0)
											distinctRealOut.erase(key);
								}
						}
						it = in.erase(it);
					}
					else it++;
				}					

			}

			bool checkExistence(bool bSign, std::string toId, bool toSign){
				for(auto& i : out){
					if(i.baseSign() == bSign and i.neighborSign() == toSign){
						if(i.contigId == toId)
							return true ;
					}
				}
				return false ;
			}

			std::vector<edgetuple>& getIn() {return in;}
			std::vector<edgetuple>& getOut() {return out;}

	private:
			std::string id ;
			uint8_t indegp ;
			uint8_t outdegp ;
			uint8_t indegm ;
			uint8_t outdegm ;
			std::vector<edgetuple> in;
			std::vector<edgetuple> out ;
			// We use these two just to count total number of incoming edges or outgoing edges and they literaly have no other usecases!!
			spp::sparse_hash_map<std::string, short> distinctRealIn;
			spp::sparse_hash_map<std::string, short> distinctRealOut;

	};

	class Graph{
	private:
		spp::sparse_hash_map<std::string,Node> Vertices ;
		//std::vector<std::pair<Node,Node> > Edges ;

	public:
		//case where I do
		spp::sparse_hash_map<std::string, Node>& getVertices(){
			return Vertices;
		}

		bool getNode(std::string nodeId){
			if(Vertices.find(nodeId) == Vertices.end())
				return true;
			else
				return false;
		}

		bool addEdge(std::string fromId, bool fromSign, std::string toId, bool toSign){
			//case 1 where the from node does not exist
			//None of the nodes exists
			if(Vertices.find(fromId) == Vertices.end()){
				Node newNode(fromId) ;
				Vertices[fromId] = newNode ;
			}
			if(Vertices.find(toId) == Vertices.end()){
				Node newNode(toId) ;
				Vertices[toId] = newNode ;
			}
			auto& fromNode = Vertices[fromId] ;
			auto& toNode = Vertices[toId] ;

			if(!fromNode.checkExistence(fromSign,toNode.getId(),toSign)){
				fromNode.insertNodeTo(toNode.getId(),fromSign,toSign) ;
				toNode.insertNodeFrom(fromNode.getId(),toSign,fromSign) ;
				return true ;
			}

			return false ;


			//Edges.emplace_back(fromNode,toNode) ;
		}

		bool removeNode(std::string id) {
			if (Vertices.find(id) == Vertices.end()) return false;
			Node& n = Vertices[id];
			for (auto & in : n.getIn()) {
				for (auto & out : n.getOut()) {
					addEdge(in.contigId, in.neighborSign(), out.contigId, out.neighborSign());					
				}
			}
			for (auto& in : n.getIn()) {
				Node& from = Vertices[in.contigId];
				from.removeEdgeTo(id);
			}
			for (auto& out : n.getOut()) {
				Node& to = Vertices[out.contigId];
				to.removeEdgeFrom(id);
			}
      return false;
	}

		/*
		void buildGraph(std::string gfa_file){
			std::ifstream file(gfa_file);
			while(std::get_line(file, ln)){
				char firstC = ln[0];
				if(firstC != 'L') continue ;

			}
		}*/

	};
}

#endif //end puffer
