#include "FatPufferGraph.hpp"
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iterator>
#include <type_traits>


int main(int argc, char* argv[]){
	std::string gfa_file = argv[1] ;

	std::ifstream file(gfa_file) ;
	std::string ln ;

	pufg::Graph G ;

	while(std::getline(file,ln)){
		char firstC = ln[0] ;
		if(firstC != 'L') continue ;
		if(firstC == 'L'){
			stx::string_view lnview(ln);
			std::vector<stx::string_view> splited = util::split(lnview, '\t');

			std::string leftLinkId = splited[1].to_string() ;
			std::string leftLinkOre = splited[2].to_string() ;
			std::string rightLinkId = splited[3].to_string() ;
			std::string rightLinkOre = splited[4].to_string() ;

			bool fromSign = (leftLinkOre == "+") ? true : false ;
			bool toSign = (rightLinkOre == "+") ? true : false ;

			G.addEdge(leftLinkId, fromSign, rightLinkId, toSign) ;
		}
	}
	auto& Vertices = G.getVertices() ;

	std::cerr << "\n Number of vertices " << Vertices.size() << "\n" ;
	for(auto& node : Vertices){
		std::cout << " \n id " << node.second.getId() << " real out deg " << (int)node.second.getRealOutdeg() << "\n" ;
	}

	//std::cout << "\n out degree of 2 " << Vertices["2"].getRealIndeg() << "\n" ;

	return 0 ;
}
