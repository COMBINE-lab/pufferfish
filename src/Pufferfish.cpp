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
#include "clipp.h"
#include <fstream>
#include <iostream>
#include <vector>
//#include <cereal/archives/json.hpp>

#include "PufferfishConfig.hpp"
#include "ProgOpts.hpp"
#include "Util.hpp"
//#include "IndexHeader.hpp"

int pufferfishIndex(IndexOptions& indexOpts); // int argc, char* argv[]);
int pufferfishTest(TestOptions& testOpts);    // int argc, char* argv[]);
int pufferfishValidate(
    ValidateOptions& validateOpts); // int argc, char* argv[]);
int pufferfishTestLookup(
    ValidateOptions& lookupOpts); // int argc, char* argv[]);
int pufferfishAligner(AlignmentOpts& alignmentOpts) ;
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
  using namespace clipp;
  using std::cout;
  enum class mode {help, index, validate, lookup, align};
  mode selected = mode::help;
  AlignmentOpts alignmentOpt ;
  IndexOptions indexOpt;
  //TestOptions testOpt;
  ValidateOptions validateOpt;
  ValidateOptions lookupOpt;

  auto indexMode = (
                    command("index").set(selected, mode::index),
                    (required("-o", "--output").call([]{cout << "parsing --output\n\n";}) & value("output_dir", indexOpt.outdir)) % "directory where index is written",
                    (required("-g", "--gfa").call([]{cout << "parsing --gfa\n\n";}) & value("gfa_file", indexOpt.gfa_file)) % "path to the GFA file",
                    (required("-r", "--ref").call([]{cout << "parsing --ref\n\n";}) & value("ref_file", indexOpt.rfile)) % "path to the reference fasta file",
                    (option("-k", "--klen") & value("kmer_length", indexOpt.k))  % "length of the k-mer with which the dBG was built (default = 31)",
                    (option("-s", "--sparse").set(indexOpt.isSparse, true)) % "use the sparse pufferfish index (less space, but slower lookup)",
                    (option("-e", "--extension") & value("extension_size", indexOpt.extensionSize)) % "length of the extension to store in the sparse index (default = 4)"
                    );

  /*
  auto testMode = (
                   command("test").set(selected, mode::test)
                   );
  */
  auto validateMode = (
                       command("validate").set(selected, mode::validate),
                       (required("-i", "--index") & value("index", validateOpt.indexDir)) % "directory where the pufferfish index is stored",
                       (required("-r", "--ref") & value("ref", validateOpt.refFile)) % "fasta file with reference sequences",
                       (required("-g", "--gfa") & value("gfa", validateOpt.gfaFileName)) % "GFA file name needed for edge table validation"
                       );
  auto lookupMode = (
                     command("lookup").set(selected, mode::lookup),
                     (required("-i", "--index") & value("index", lookupOpt.indexDir)) % "directory where the pufferfish index is stored",
                     (required("-r", "--ref") & value("ref", lookupOpt.refFile)) % "fasta file with reference sequences"
                     );

  std::string throwaway;
  auto alignMode = (
                    command("align").set(selected, mode::align),
                    (required("-i", "--index") & value("index", alignmentOpt.indexDir)) % "directory where the pufferfish index is stored",
                    (
                     (
                     ((required("--mate1") & value("mate 1", alignmentOpt.read1)) % "path to the left end of the read files"),
                     ((required("--mate2") & value("mate 2", alignmentOpt.read2)) % "path to the right end of the read files")
                     ) |
                     ((required("--read") & value("reads", throwaway)) % "path to single-end read files")
                    ),
                    (option("-m", "--just-mapping") & value("just map", alignmentOpt.justMap)) % "don't attempt alignment validation; just do mapping",
                    (option("-p", "--threads") & value("num threads", alignmentOpt.numThreads)) % "specify the number of threads (default=8)",
                    (required("-o", "--outdir") & value("output file", alignmentOpt.outname)) % "output file where the alignment results will be stored",
                    (option("--maxSpliceGap") & value("max splice gap", alignmentOpt.maxSpliceGap)) % "specify maximum splice gap that two uni-MEMs should have",
                    (option("--maxFragmentLength") & value("max frag length", alignmentOpt.maxFragmentLength)) % "specify the maximum distance between the last uni-MEM of the left and first uni-MEM of the right end of the read pairs",
                    (option("--writeOrphans").set(alignmentOpt.writeOrphans, true)) % "write Orphans flag",
                    (option("--noOutput").set(alignmentOpt.noOutput, true)) % "run without writing SAM file"
                    );

  auto cli = (
              (indexMode | validateMode | lookupMode | alignMode | command("help").set(selected,mode::help) ),
              option("-v", "--version").call([]{std::cout << "version 0.1.0\n\n";}).doc("show version"));

  decltype(parse(argc, argv, cli)) res;
  try {
    res = parse(argc, argv, cli);
  } catch (std::exception& e) {
    std::cout << "\n\nParsing command line failed with exception: " << e.what() << "\n";
    std::cout << "\n\n";
    std::cout << make_man_page(cli, "mantis");
    return 1;
  }


  if(res) {
    switch(selected) {
    case mode::index: pufferfishIndex(indexOpt);  break;
    case mode::validate: pufferfishValidate(validateOpt);  break;
    case mode::lookup: pufferfishTestLookup(lookupOpt); break;
    case mode::align: pufferfishAligner(alignmentOpt); break;
    case mode::help: std::cout << make_man_page(cli, "pufferfish"); break;
    }
  } else {
    debug::print(std::cerr, res);
    auto b = res.begin();
    auto e = res.end();
    std::cerr << "any blocked " << res.any_blocked() << "\n";
    std::cerr << "any conflict " << res.any_conflict() << "\n";
    std::cerr << "any bad repeat " << res.any_bad_repeat() << "\n";
    std::cerr << "any error " << res.any_error() << "\n";
    for( auto& m : res.missing() ) {
      std::cerr << "missing " << m.param()->label() << "\n";
    }
    if (std::distance(b,e) > 0) {
      if (b->arg() == "index") {
        std::cout << make_man_page(indexMode, "pufferfish");
        while (b != e) {
          std::cerr << b->arg() << "\n";
          ++b;
        }
      } else if (b->arg() == "validate") {
        std::cout << make_man_page(validateMode, "pufferfish");
      } else if (b->arg() == "lookup") {
        std::cout << make_man_page(lookupMode, "pufferfish");
      } else if (b->arg() == "align") {
        std::cout << make_man_page(alignMode, "pufferfish");
      } else {
        std::cout << "There is no command \"" << b->arg() << "\"\n";
        std::cout << usage_lines(cli, "pufferfish") << '\n';
        return 1;
      }
    } else {
      std::cout << usage_lines(cli, "pufferfish") << '\n';
      return 1;
    }
  }


  /*
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

  AlignmentOpts alignmentOpt ;
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
    ->add_flag("-m,--just-mapping", alignmentOpt.justMap,
               "don't attempt alignment validation; just do mapping");
  alignmentApp
    ->add_option("-p,--threads", alignmentOpt.numThreads,
                 "specfy number of threads") ;
  alignmentApp
    ->add_option("-o,--outdir", alignmentOpt.outname,
                 "output directory where the alignment results would get stored")
    ->required() ;

  alignmentApp
    ->add_option(",--maxSpliceGap", alignmentOpt.maxSpliceGap,
               "specify maximum splice gap that two uni-mems should have");
  alignmentApp
    ->add_option(",--maxFragmentLength", alignmentOpt.maxFragmentLength,
               "specify maximum distance between last uni-mem of the left end and first uni-mem of the right end of the read pairs");
  alignmentApp
    ->add_flag(",--writeOrphans", alignmentOpt.writeOrphans,
                 "write Orphans flag");

  alignmentApp
    ->add_flag(",--noOutput", alignmentOpt.noOutput,
                 "run without writing sam");


  CLI::App app{"Pufferfish : An efficient dBG index."};
  auto indexApp = app.add_subcommand("index", "build the pufferfish index");
  auto testApp = app.add_subcommand("test", "test k-mer lookup in the index");
  auto validateApp = app.add_subcommand(
      "validate", "test k-mer lookup for reference sequences");
  auto lookupApp = app.add_subcommand("lookup", "test k-mer lookup");
  auto alignmentApp = app.add_subcommand("align", "align paired end RNA-seq reads") ;

  IndexOptions indexOpt;
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
  IndexOptions indexOpt;

  TestOptions testOpt;

  ValidateOptions validateOpt;
  validateApp
      ->add_option("-i,--index", validateOpt.indexDir,
                   "directory where the pufferfish index is stored")
      ->required();
  validateApp
      ->add_option("-r,--ref", validateOpt.refFile,
                   "fasta file with reference sequences")
      ->required();
  validateApp
    ->add_option("-g,--gfa", validateOpt.gfaFileName,
                 "GFA file name needed for edge table validation") ;

  ValidateOptions lookupOpt;
  lookupApp
      ->add_option("-i,--index", lookupOpt.indexDir,
                   "directory where the pufferfish index is stored")
      ->required();
  lookupApp
      ->add_option("-r,--ref", lookupOpt.refFile,
                   "fasta file with reference sequences")
      ->required();


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
  */
  return 0;
}
