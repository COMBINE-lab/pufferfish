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

#include "CLI/CLI.hpp"
#include <fstream>
#include <iostream>
#include <vector>
//#include <cereal/archives/json.hpp>

#include "PufferfishConfig.hpp"
#include "Util.hpp"
//#include "IndexHeader.hpp"

int pufferfishIndex(util::IndexOptions& indexOpts); // int argc, char* argv[]);
int pufferfishTest(util::TestOptions& testOpts);    // int argc, char* argv[]);
int pufferfishValidate(
    util::ValidateOptions& validateOpts); // int argc, char* argv[]);
int pufferfishTestLookup(
    util::ValidateOptions& lookupOpts); // int argc, char* argv[]);
int pufferfishAligner(util::AlignmentOpts& alignmentOpts) ;
// int rapMapMap(int argc, char* argv[]);
// int rapMapSAMap(int argc, char* argv[]);

void printUsage() {
  std::string versionString = pufferfish::version;
  std::cerr << "Pufferfish v" << versionString << '\n';
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

int main(int argc, char* argv[]) {

  CLI::App app{"Pufferfish : An efficient dBG index."};
  auto indexApp = app.add_subcommand("index", "build the pufferfish index");
  auto testApp = app.add_subcommand("test", "test k-mer lookup in the index");
  auto validateApp = app.add_subcommand(
      "validate", "test k-mer lookup for reference sequences");
  auto lookupApp = app.add_subcommand("lookup", "test k-mer lookup");
  auto alignmentApp = app.add_subcommand("align", "align paired end RNA-seq reads") ;

  util::IndexOptions indexOpt;
  indexApp
      ->add_option("-k,--klen", indexOpt.k,
                   "length of the k-mer with which the compacted dBG was built",
                   static_cast<uint32_t>(31))
      ->required();
  indexApp->add_option("-g,--gfa", indexOpt.gfa_file, "path to the GFA file")
      ->required();
  indexApp
      ->add_option("-o,--output", indexOpt.outdir,
                   "directory where index is written")
      ->required();
  indexApp
    ->add_flag("-s,--sparse", indexOpt.isSparse,
               "use the sparse pufferfish index (less space, but slower lookup)");
  indexApp
    ->add_option("-e,--extension", indexOpt.extensionSize,
                 "length of the extension to store in the sparse index",
                 static_cast<uint32_t>(4));

  util::TestOptions testOpt;

  util::ValidateOptions validateOpt;
  validateApp
      ->add_option("-i,--index", validateOpt.indexDir,
                   "directory where the pufferfish index is stored")
      ->required();
  validateApp
      ->add_option("-r,--ref", validateOpt.refFile,
                   "fasta file with reference sequences")
      ->required();

  util::ValidateOptions lookupOpt;
  lookupApp
      ->add_option("-i,--index", lookupOpt.indexDir,
                   "directory where the pufferfish index is stored")
      ->required();
  lookupApp
      ->add_option("-r,--ref", lookupOpt.refFile,
                   "fasta file with reference sequences")
      ->required();

  util::AlignmentOpts alignmentOpt ;
  alignmentApp
    ->add_option("-i,--index", alignmentOpt.indexDir,
                 "directory where the pufferfish index is stored")
    ->required() ;
  alignmentApp
    ->add_option(",--mate1", alignmentOpt.read1,
                 "path to left end of the read files")
    ->required() ;
  alignmentApp
    ->add_option(",--mate2", alignmentOpt.read2,
                 "path to right end of the read files")
    ->required() ;
  alignmentApp
    ->add_option("-p,--threads", alignmentOpt.numThreads,
                 "specfy number of threads") ;
  alignmentApp
    ->add_option("-o,--outdir", alignmentOpt.outname,
                 "output directory where the alignment results would get stored")
    ->required() ;

  alignmentApp
    ->add_flag(",--writeOrphans", alignmentOpt.writeOrphans,
                 "write Orphans flag");

  alignmentApp
    ->add_flag(",--noOutput", alignmentOpt.noOutput,
                 "run without writing sam");


  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }

  if (app.got_subcommand(indexApp)) {
    return pufferfishIndex(indexOpt);
  } else if (app.got_subcommand(testApp)) {
    return pufferfishTest(testOpt);
  } else if (app.got_subcommand(validateApp)) {
    return pufferfishValidate(validateOpt);
  } else if (app.got_subcommand(lookupApp)) {
    return pufferfishTestLookup(lookupOpt);
  } else if (app.got_subcommand(alignmentApp)) {
    return pufferfishAligner(alignmentOpt);
  }
  else {
    std::cerr << "I don't know the requested sub-command\n";
    return 1;
  }
  return 0;
}
