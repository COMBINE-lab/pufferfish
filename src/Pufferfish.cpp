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

#include "clipp.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <clocale>
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

int main(int argc, char* argv[]) {
  using namespace clipp;
  using std::cout;
  std::setlocale(LC_ALL, "en_US.UTF-8");

  enum class mode {help, index, validate, lookup, align};
  mode selected = mode::help;
  AlignmentOpts alignmentOpt ;
  IndexOptions indexOpt;
  //TestOptions testOpt;
  ValidateOptions validateOpt;
  ValidateOptions lookupOpt;

  auto indexMode = (
                    command("index").set(selected, mode::index),
                    (required("-o", "--output") & value("output_dir", indexOpt.outdir)) % "directory where index is written",
                    (required("-g", "--gfa") & value("gfa_file", indexOpt.gfa_file)) % "path to the GFA file",
                    (option("-r", "--ref") & value("ref_file", indexOpt.rfile)) % "path to the reference fasta file",
                    (option("-k", "--klen") & value("kmer_length", indexOpt.k))  % "length of the k-mer with which the dBG was built (default = 31)",
                    (option("-p", "--threads") & value("threads", indexOpt.p))  % "total number of threads to use for building MPHF (default = 16)",
                    (option("-l", "--build-edges").set(indexOpt.buildEdgeVec, true) % "build and record explicit edge table for the contaigs of the ccdBG (default = false)"),
                    (((option("-s", "--sparse").set(indexOpt.isSparse, true)) % "use the sparse pufferfish index (less space, but slower lookup)",
                     ((option("-e", "--extension") & value("extension_size", indexOpt.extensionSize)) % "length of the extension to store in the sparse index (default = 4)")) |
                     ((option("-x", "--lossy-rate").set(indexOpt.lossySampling, true)) & value("lossy_rate", indexOpt.lossy_rate) % "use the lossy sampling index with a sampling rate of x (less space and fast, but lower sensitivity)"))
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
  auto isValidRatio = [](const char* s) -> void {
    float r{0.0};
    std::string sv(s);
    try {
      r = std::stof(sv);
    } catch (std::exception& e) {
      std::string m = "Could not convert " + sv + " to a valid ratio\n";
      throw std::domain_error(m);
    }
    if (r <= 0 or r > 1) {
      std::string m = "The --scoreRatio you provided was " + sv + ", it must be in (0,1]\n";
      throw std::domain_error(m);
    }
  };

  auto alignMode = (
                    command("align").set(selected, mode::align),
                    (required("-i", "--index") & value("index", alignmentOpt.indexDir)) % "directory where the pufferfish index is stored",
                    (
                     (
                      ((required("--mate1", "-1") & value("mate 1", alignmentOpt.read1)) % "path to the left end of the read files"),
                      ((required("--mate2", "-2") & value("mate 2", alignmentOpt.read2)) % "path to the right end of the read files")
                     ) |
                     ((required("--read").set(alignmentOpt.singleEnd, true) & value("reads", alignmentOpt.unmatedReads)) % "path to single-end read files")
                    ),
                    (option("--scoreRatio") & value("score ratio", alignmentOpt.scoreRatio).call(isValidRatio)) % "mappings with a score < scoreRatio * OPT are discarded (default=1.0)",
                    (option("-p", "--threads") & value("num threads", alignmentOpt.numThreads)) % "specify the number of threads (default=8)",
                    (option("-m", "--just-mapping").set(alignmentOpt.justMap, true)) % "don't attempt alignment validation; just do mapping",
                    (
                      (required("--noOutput").set(alignmentOpt.noOutput, true)) % "run without writing SAM file"
                        |
                      (required("-o", "--outdir") & value("output file", alignmentOpt.outname)) % "output file where the alignment results will be stored"
                    ),
                    (option("--maxSpliceGap") & value("max splice gap", alignmentOpt.maxSpliceGap)) % "specify maximum splice gap that two uni-MEMs should have",
                    (option("--maxFragmentLength") & value("max frag length", alignmentOpt.maxFragmentLength)) % "specify the maximum distance between the last uni-MEM of the left and first uni-MEM of the right end of the read pairs",
                    (option("--noOrphans").set(alignmentOpt.noOrphan, true)) % "write Orphans flag",
                    (option("--noDiscordant").set(alignmentOpt.noDiscordant, true)) % "write Orphans flag",
		    (option("-z", "--compressedOutput").set(alignmentOpt.compressedOutput, true)) % "compress (gzip) the output file",
                    (
                    (option("-k", "--krakOut").set(alignmentOpt.krakOut, true)) % "write output in the format required for krakMap"
                    |
                    (option("-s", "--salmon").set(alignmentOpt.salmonOut, true)) % "write output in the format required for salmon"
                    ),
                    (option("--minScoreFraction") & value("minScoreFraction", alignmentOpt.minScoreFraction)) % "minScoreFraction"
                    );

  auto cli = (
              (indexMode | validateMode | lookupMode | alignMode | command("help").set(selected,mode::help) ),
              option("-v", "--version").call([]{std::cout << "version " << pufferfish::version << "\n"; std::exit(0);}).doc("show version"));

  decltype(parse(argc, argv, cli)) res;
  try {
    res = parse(argc, argv, cli);
  } catch (std::exception& e) {
    std::cout << "\n\nParsing command line failed with exception: " << e.what() << "\n";
    std::cout << "\n\n";
    std::cout << make_man_page(cli, pufferfish::progname);
    return 1;
  }


  if(res) {
    switch(selected) {
    case mode::index: pufferfishIndex(indexOpt);  break;
    case mode::validate: pufferfishValidate(validateOpt);  break;
    case mode::lookup: pufferfishTestLookup(lookupOpt); break;
    case mode::align: pufferfishAligner(alignmentOpt); break;
    case mode::help: std::cout << make_man_page(cli, pufferfish::progname); break;
    }
  } else {
    auto b = res.begin();
    auto e = res.end();
    if (std::distance(b,e) > 0) {
      if (b->arg() == "index") {
        std::cout << make_man_page(indexMode, pufferfish::progname);
      } else if (b->arg() == "validate") {
        std::cout << make_man_page(validateMode, pufferfish::progname);
      } else if (b->arg() == "lookup") {
        std::cout << make_man_page(lookupMode, pufferfish::progname);
      } else if (b->arg() == "align") {
        std::cout << make_man_page(alignMode, pufferfish::progname);
      } else {
        std::cout << "There is no command \"" << b->arg() << "\"\n";
        std::cout << usage_lines(cli, pufferfish::progname) << '\n';
        return 1;
      }
    } else {
      std::cout << usage_lines(cli, pufferfish::progname) << '\n';
      return 1;
    }
  }

  return 0;
}
