#include "FatPufferGraph.hpp"
//#include "OurGFAReader.hpp"
//#include "semiCompactedCompactor.hpp"
#include "CLI/CLI.hpp"
#include "GFAConverter.hpp"

struct PufferizeOpts {
  uint32_t k = 31;
  std::string gfaFile;
  std::string outFile;
};

int main(int argc, char* argv[]) {
  (void)argc;

  PufferizeOpts popts;
  CLI::App app{"pufferize : prepare a GFA file for pufferfish indexing."};
  app.add_option("-k,--klen", popts.k,
                 "length of the k-mer with which the compacted dBG was built",
                 static_cast<uint32_t>(31))
      ->required();
  app.add_option("-g,--gfa", popts.gfaFile, "path to the input GFA file")
      ->required();
  app.add_option("-o,--output", popts.outFile, "path to the output GFA file")
      ->required();

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }

  GFAConverter gc(popts.gfaFile.c_str(), popts.k);
  gc.parseFile();
  gc.buildGraph();
  gc.randomWalk();
  gc.writeFile(popts.outFile.c_str());
  /*	PosFinder pf(argv[1], k);
      pf.parseFile();

      SemiCompactedCompactor scc(pf.getNewSegments(),
                            pf.getContigNameMap(),
                            pf.getPathStart(),
                            pf.getPathEnd(),
                            pf.getSemiCG(),
                            k);
      scc.compact();
      scc.updatePath(pf.getPaths());
      scc.writeFile(gfa_outfilename);
  */

  return 0;
}
