#ifndef PUFFER_GRAPH_H
#define PUFFER_GRAPH_H

#include <map>
#include "Util.hpp"
#include "sparsepp/spp.h"
#include <iostream>

namespace pufg{
	struct edgetuple{
		edgetuple(bool fSign,std::string cId, bool tSign):
			baseSign(fSign),contigId(cId),neighborSign(tSign) {}

		bool baseSign;
		std::string contigId;
		bool neighborSign;
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

			//std::vector<std::string> getPositiveToNodes()
			std::string getId() {return id;}

			void insertNodeTo(std::string nodeId, bool sign, bool toSign){
				if(sign) {
					outdegp++;
					distinctRealOut[nodeId + (toSign?"+":"-")] = true;
				}
				else {
					outdegm++;
					distinctRealIn[nodeId + (toSign?"-":"+")] = true;
				}

				out.emplace_back(sign, nodeId, toSign);
			}


			void insertNodeFrom(std::string nodeId, bool sign, bool fromSign){
				if(sign) {
					indegp++ ;
					distinctRealIn[nodeId + (fromSign?"+":"-")] = true;
				}
				else {
					indegm++ ;
					distinctRealOut[nodeId + (fromSign?"-":"+")] = true;
				}
				in.emplace_back(sign, nodeId, fromSign) ;
			}
			
			bool checkExistence(bool bSign, std::string toId, bool toSign){
				for(auto& i : out){
					if(i.baseSign == bSign and i.neighborSign == toSign){
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
			spp::sparse_hash_map<std::string, bool> distinctRealIn;
			spp::sparse_hash_map<std::string, bool> distinctRealOut;

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
