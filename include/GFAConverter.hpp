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

enum class Side : uint8_t {
		LEFT = 0,
		RIGHT = 1,
		LEFT_RIGHT = 2
};

class GFAConverter {
private:
		std::string gfaFileName_;
  std::unique_ptr<zstr::ifstream> file;
  size_t k;
  spp::sparse_hash_map<uint64_t, std::string> id2seq;
  std::vector<uint64_t> ksizeContig;
  // path maps each transcript_id to a pair of <contig_id, orientation>
  // orientation : +/true main, -/false reverse

  std::map<std::pair<uint64_t, bool>, std::vector<std::string>, util::cmpByPair> start2path;
  spp::sparse_hash_map<uint64_t, Side> donTouch;
  spp::sparse_hash_map<std::string, std::pair<uint64_t, bool> > path2start;
  spp::sparse_hash_map<uint64_t, std::vector<std::pair<uint64_t, bool> > > old2newids;
  pufg::Graph semiCG;

  void
  processContigSeq(uint64_t& contigId, std::string& contigSeq,
                   spp::sparse_hash_map<std::string, uint64_t>& seq2id,
                   uint64_t& idCntr);
  std::vector<std::pair<uint64_t, bool> >& convertOld2newPath(std::vector<std::pair<uint64_t, bool> >& oContigVec);
  void setDonTouchSide(std::vector<std::pair<uint64_t, bool> >& contigVec); 
  void updateDonTouchSide(uint64_t& id, Side side);
  void mergeDonTouchSides(uint64_t& mergedId, bool mergedOri, uint64_t& mergingId, bool mergingOri);
  bool dont_touch(uint64_t& nodeId, Side side);
  void setPathStart(std::string& pathId, std::pair<uint64_t, bool>& startNode);
  void updatePathStart(uint64_t& oStartNodeId, bool oStartNodeOri, uint64_t& nStartNodeId, bool nStartNodeOri);
  bool isCornerCase(pufg::Node& n, bool mergeLeft);
  void mergeLeft(pufg::Node& n);
  void mergeRight(pufg::Node& n);
  void updateGraph(std::vector<std::pair<uint64_t, bool> >& contigVec);

public:
  GFAConverter(const char* gfaFileName, size_t input_k);
  void parseFile();
//  void buildGraph();
  void randomWalk();
  void reconstructPathAndWrite(const char* gfaOutFileName, const char* refFileName);
};

#endif
