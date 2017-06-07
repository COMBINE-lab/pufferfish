#include "FatPufferGraph.hpp"
#include "OurGFAReader.hpp"
#include "semiCompactedCompactor.hpp"

int main(int argc, char* argv[]){
	std::string gfa_infilename = argv[1] ;
	std::string gfa_outfilename = argv[2];

	PosFinder pf(argv[1], 31);
	pf.parseFile();
	
	SemiCompactedCompactor scc(pf.getNewSegments(), 
						  pf.getContigNameMap(),
						  pf.getPathStart(),
						  pf.getPathEnd(),
						  pf.getSemiCG());
	scc.compact();
	scc.updatePath(pf.getPaths());
	scc.writeFile(gfa_outfilename);

	return 0 ;
}
