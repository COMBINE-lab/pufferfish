#ifndef __PUFFERFISH_PROG_OPTS_HPP__
#define __PUFFERFISH_PROG_OPTS_HPP__
#include <cstdint>

class IndexOptions {
public:
  uint32_t k{31};
  uint32_t p{16};
  std::string gfa_file;
  std::string cfile;
  std::string rfile;
  std::string outdir;
  bool isSparse{false};
  uint32_t extensionSize{4};
  uint32_t sampleSize{9};
  bool lossySampling{false};
  uint32_t lossy_rate{5};
  bool buildEdgeVec{false};
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
  bool singleEnd{false};
  uint32_t numThreads{1};
  uint32_t maxNumHits{200};
  uint32_t maxSpliceGap{100};
  uint32_t maxFragmentLength{1000};
  double scoreRatio{0.6};
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
  bool krakOut{false};
  bool salmonOut{false};
  bool noDiscordant{false};
  bool noOrphan{false};
  bool compressedOutput{false};
  bool verbose{false};
  bool validateMappings{false};
  bool strictFilter{false};
  int32_t gapOpenPenalty{5};
  int32_t gapExtendPenalty{3};
  int32_t matchScore{2};
  int32_t missMatchScore{-4};
  uint32_t refExtendLength{20};
  double minScoreFraction{0.65};
  bool fullAlignment{false};
  bool heuristicChaining{true};
  bool mergeMems{true};
  bool genomicReads{false};
  std::string genesNamesFile{""};
  bool filterGenomics{false};
  bool filterMicrobiom{false};
  bool primaryAlignment{false};
  bool listOfReads{false};
};


#endif
