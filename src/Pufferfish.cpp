//
// Pufferfish - An efficient compacted dBG index
// 
// Copyright (C) 2017 Rob Patro, Fatemeh Almodaresi, Hirak Sarkar
//
// This file is part of Pufferfish.
//
// RapMap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RapMap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RapMap.  If not, see <http://www.gnu.org/licenses/>.
//

#include <iostream>
#include <fstream>
#include <vector>
#include "CLI/CLI.hpp"
//#include <cereal/archives/json.hpp>

#include "Util.hpp"
#include "PufferfishConfig.hpp"
//#include "IndexHeader.hpp"

int pufferfishIndex(util::IndexOptions& indexOpts);//int argc, char* argv[]);
int pufferfishTest(int argc, char* argv[]);
//int rapMapMap(int argc, char* argv[]);
//int rapMapSAMap(int argc, char* argv[]);

void printUsage() {
    std::string versionString = pufferfish::version;
    std::cerr << "Pufferfish v"
              << versionString << '\n';
    std::cerr << "=====================================\n";
    auto usage =
        R"(
There are currently 2 Pufferfish subcommands
    index --- builds a Pufferfish index
    test --- tests k-mer lookup in the index

Run a corresponding command "pufferrish <cmd> -h" for
more information on each of the possible Pufferfish
commands.)";
    std::cerr << usage << '\n';
}


int testIndexApp(CLI::App* app) {
  std::cerr << "woot!\n";
  auto opts = app->parse_order();
  for (auto& o : opts) {
    std::cerr << "option " << o->results().front() << "\n";
  }
}

int main(int argc, char* argv[]) {

  CLI::App app{"App description"};
  auto indexApp = app.add_subcommand("index", "build the pufferfish index");
  auto testApp = app.add_subcommand("test", "test k-mer lookup in the index");

  util::IndexOptions indexOpt;
  indexApp->add_option("-k,--klen", indexOpt.k, "length of the k-mer with which the compacted dBG was built", static_cast<uint32_t>(31));
  indexApp->add_option("-g,--gfa", indexOpt.gfa_file, "path to the GFA file");
  indexApp->add_option("-o,--output", indexOpt.outdir, "directory where index is written");
  indexApp->add_option("-f,--fa", indexOpt.rfile, "path to the Fasta file with reads");

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  if (app.got_subcommand(indexApp)) {
    return pufferfishIndex(indexOpt);
    //return testIndexApp(indexApp);
  } else if (app.got_subcommand(testApp)) {
    std::cerr << "not yet implemented";
  }
  /*
  uint32_t k;
  std::string gfa_file;
  std::string cfile;
  std::string rfile;
  std::string outdir;
  popl::Switch helpOption("h", "help", "produce help message");
  popl::Value<uint32_t> kOpt("k", "klen", "length of the k-mer with which the compacted dBG was built", 31, &k);
  popl::Value<std::string> gfaOpt("g", "gfa", "path to the GFA file");
  popl::Value<std::string> outOpt("o", "output", "directory where index is written");
  popl::Value<std::string> readOpt("f", "fa", "path to the Fasta file with reads");

  popl::OptionParser op("Allowed options");
	op.add(helpOption)
    .add(kOpt)
    .add(gfaOpt)
    .add(readOpt)
    .add(outOpt);

  op.parse(argc, argv);
  if (helpOption.isSet()) {
    std::cout << op << '\n';
    std::exit(0);
  }

  gfa_file = gfaOpt.getValue();
  rfile = readOpt.getValue();
  outdir = outOpt.getValue();



    std::vector<char*> args;
    args.push_back(argv[0]);

    if (argc < 2) {
        printUsage();
        std::exit(0);
    }

    for (int i = 2; i < argc; ++i) {
        args.push_back(argv[i]);
    }

    if (std::string(argv[1]) == "-h" or
        std::string(argv[1]) == "--help") {
        printUsage();
        std::exit(0);
    }

    if (std::string(argv[1]) == "index") {
        return pufferfishIndex(argc - 1, args.data());
    } else if (std::string(argv[1]) == "test") {
        return pufferfishTest(argc - 1, args.data());
    } else {
        std::cerr << "the command " << argv[1]
                  << " is not yet implemented\n";
      return 1;
    }
  */
    return 0;
}
