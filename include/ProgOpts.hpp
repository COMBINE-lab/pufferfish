#ifndef __PUFFERFISH_PROG_OPTS_HPP__
#define __PUFFERFISH_PROG_OPTS_HPP__
#include <cstdint>

class IndexOptions {
public:
  uint32_t k{31};
  std::string gfa_file;
  std::string cfile;
  std::string rfile;
  std::string outdir;
  bool isSparse{false};
  uint32_t extensionSize{4};
  uint32_t sampleSize{9};
};

class TestOptions {
public:
};

class ValidateOptions {
public:
  std::string indexDir;
  std::string refFile;
  std::string gfaFileName ;
};

class AlignmentOpts{
public:

  std::string indexDir;
  std::string read1;
  std::string read2;
  std::string unmatedReads;
  uint32_t numThreads{1};
  uint32_t maxNumHits{200};
  uint32_t maxSpliceGap{100};
  uint32_t maxFragmentLength{100000};
  double scoreRatio{0.5};
  std::string outname;
  double quasiCov{0.0};
  bool pairedEnd{false};
  bool noOutput{false};
  bool sensitive{false};
  bool strictCheck{false};
  bool fuzzy{false};
  bool consistentHits{false};
  bool quiet{false};
	bool writeOrphans{false} ;
  bool justMap{false};
};


#endif
