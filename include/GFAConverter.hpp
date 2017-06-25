#ifndef OUR_GFA_READER_H
#define OUR_GFA_READER_H

#include "FatPufferGraph.hpp"
#include "Util.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/vector.hpp"
#include "sparsepp/spp.h"
#include "string_view.hpp"
#include "zstr/zstr.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

class GFAConverter {
private:
  std::unique_ptr<zstr::ifstream> file;
  size_t k;
  struct Contig {
    std::string seq;
    std::string id;
  };
  spp::sparse_hash_map<
      uint64_t,
      std::pair<std::string, std::vector<std::pair<uint64_t, bool>>>>
      new2seqAoldids;
  std::vector<uint64_t> ksizeContig;
  spp::sparse_hash_map<uint64_t, std::vector<std::pair<uint64_t, bool>>>
      old2newids;
  // path maps each transcript_id to a pair of <contig_id, orientation>
  // orientation : +/true main, -/false reverse
  spp::sparse_hash_map<std::string, std::vector<std::pair<uint64_t, bool>>>
      path;

  std::map<std::pair<uint64_t, bool>, bool, util::cmpByPair> pathStart;
  std::map<std::pair<uint64_t, bool>, bool, util::cmpByPair> pathEnd;

  pufg::Graph semiCG;

  void
  processContigSeq(uint64_t& contigId, std::string& contigSeq,
                   spp::sparse_hash_map<std::string, uint64_t>& seq2newid,
                   uint64_t& idCntr);
  bool is_start(uint64_t& nodeId);
  bool is_end(uint64_t& nodeId);
  void update_start(uint64_t& newId, bool newOri);
  void update_end(uint64_t& newId, bool newOri);
  // bool isCornerCase(pufg::Node& n, pufg::edgetuple, bool mergeIn);
  bool isCornerCase(pufg::Node& n, bool mergeIn);
  void mergeIn(pufg::Node& n);
  void mergeOut(pufg::Node& n);
  void eraseFromOldList(uint64_t nodeId);

public:
  GFAConverter(const char* gfaFileName, size_t input_k);
  void parseFile();
  void buildGraph();
  void randomWalk();
  void writeFile(const char* gfaFileName);
};

#endif
