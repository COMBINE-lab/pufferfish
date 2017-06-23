#ifndef OUR_GFA_READER_H
#define OUR_GFA_READER_H

#include "FatPufferGraph.hpp"
#include "Util.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/vector.hpp"
#include "sdsl/int_vector.hpp"
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

class PosFinder {
private:
  std::string filename_;
  std::unique_ptr<zstr::ifstream> file;
  size_t k;
  struct Contig {
    std::string seq;
    std::string id;
  };

  spp::sparse_hash_map<uint64_t, util::PackedContigInfo>
      contigid2seq; // map of contig_id to # of letters in contig (contig
                    // length)
  spp::sparse_hash_map<std::string, std::string> seq2contigid;
  // path maps each transcript_id to a pair of <contig_id, orientation>
  // orientation : +/true main, -/false reverse
  spp::sparse_hash_map<uint64_t, std::vector<std::pair<uint64_t, bool>>> path;
  spp::sparse_hash_map<uint64_t, uint32_t> refIDs;

  // spp::sparse_hash_map<uint64_t, std::string> refMap;
  std::vector<std::string> refMap;
  // maps each contig to a list of positions in different transcripts
  std::vector<std::pair<uint64_t, bool>> explode(const stx::string_view str,
                                                 const char& ch);

  /** Added for the new requirements to build semiCompacted DBG **/
  // spp::sparse_hash_map<std::string, bool> pathStart;
  // spp::sparse_hash_map<std::string, bool> pathEnd;
  std::map<std::pair<std::string, bool>, bool, util::cmpByPair> pathStart;
  std::map<std::pair<std::string, bool>, bool, util::cmpByPair> pathEnd;

  sdsl::int_vector<2> seqVec_;
  std::vector<std::pair<std::string, std::string>> newSegments;
  pufg::Graph semiCG;

  size_t fillContigInfoMap_();
  bool is_number(const std::string& s);
  void encodeSeq(sdsl::int_vector<2>& seqVec, size_t offset,
                 stx::string_view str);

  // Avoiding un-necessary stream creation + replacing strings with string view
  // is a bit > than a 2x win!
  // implementation from : https://marcoarena.wordpress.com/tag/string_view/
  std::vector<stx::string_view> split(stx::string_view str, char delims);

public:
  spp::sparse_hash_map<uint64_t, std::vector<util::Position>> contig2pos;
  PosFinder(const char* gfaFileName, size_t input_k);
  // spp::sparse_hash_map<uint64_t, std::string>& getContigNameMap();
  spp::sparse_hash_map<uint64_t, util::PackedContigInfo>& getContigNameMap();

  spp::sparse_hash_map<std::string, std::string>& getContigIDMap();
  // spp::sparse_hash_map<uint32_t, std::string>& getRefIDs();
  std::vector<std::string>& getRefIDs();
  std::map<std::pair<std::string, bool>, bool, util::cmpByPair>& getPathStart();
  std::map<std::pair<std::string, bool>, bool, util::cmpByPair>& getPathEnd();
  std::vector<std::pair<std::string, std::string>>& getNewSegments();
  sdsl::int_vector<2>& getContigSeqVec();
  // spp::sparse_hash_map<std::string, std::vector< std::pair<std::string, bool>
  // > >& getPaths() {return path;}
  // spp::sparse_hash_map<std::string, std::vector< std::pair<std::string, bool>
  // > >& getPaths() {return path;}
  pufg::Graph& getSemiCG();

  void parseFile();
  void mapContig2Pos();
  void clearContigTable();
  void serializeContigTable(const std::string& odir);
  void deserializeContigTable();
  // void writeFile(std::string fileName);
};

#endif
