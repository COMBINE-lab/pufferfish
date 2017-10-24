#include "FatPufferGraph.hpp"
#include "zstr/zstr.hpp"
#include <fstream>
#include <iostream>
//#include "OurGFAReader.hpp"
//#include "semiCompactedCompactor.hpp"
#include "CLI/CLI.hpp"
//#include "GFAConverter.hpp"
#include "Util.hpp"



void parseGFA(std::string& gfaFile, pufg::Graph& g){
  std::string ln;
  std::string tag, id, value;
  size_t contig_cnt{0};

  std::unique_ptr<zstr::ifstream> file;
  file.reset(new zstr::ifstream(gfaFile));

  std::cerr << "Start loading segments... \n";

  //we just need ids and not sequences
  //spp::sparse_hash_map<std::string, uint64_t> seq2id;

  uint64_t idCntr = 0;


  while (std::getline(*file, ln)) {
    char firstC = ln[0];
    //skipping everything other than path
		if (firstC != 'P')
			continue;
    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = util::split(lnview, '\t');
    if (splited.size() != 3) {
      std::cerr << "the path line should have 3 fields (tab separated) --- skipping!\n";
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

		//contigVec = convertOld2newPath(contigVec);

		// update graph and add the new set of segments to it
    for (size_t i = 0; i < contigVec.size()-1; i++) 
      g.addEdge(contigVec[i].first, contigVec[i].second, contigVec[i+1].first, contigVec[i+1].second);

		// set pathStart and pathEnd (do we need this)
		// setPathStart(id, contigVec[0]); (do we need this)
		// setDonTouchSide(contigVec); (do we need this)
	}
  std::cerr << "Done Reading gfa file\n";
 

}

int main(int argc, char* argv[]) {
  (void)argc;

  CLI::App app{"Calculate Edge density"};

  std::string gfaFile ;// = "/mnt/scratch7/pufferfish_data/gencode.v25.pc_transcripts_fixed.k27.pufferized.gfa";

  app.add_option("-g,--gfa", gfaFile, "path to the input GFA file")
    ->required();

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }



  pufg::Graph g ;
  parseGFA(gfaFile, g) ;

  //Go over all the vertices and 
  //calculate the statistics

  spp::sparse_hash_map<uint64_t, pufg::Node>& nodes = g.getVertices();
  double avgDensity{0} ;

  for(const auto& x : nodes){
    auto node = x.second ;
    auto inEdges = node.getPredecessors() ;
    auto outEdges = node.getSuccessors() ;

    size_t totalEdges = inEdges.size() + outEdges.size() ;
    avgDensity += double(totalEdges)/8.0 ; 
  }

  avgDensity = avgDensity / double(nodes.size()) ;

  std::cout << "Edge density for "<<gfaFile<<" "<<avgDensity << "\n" ;


}
