#ifndef REF_SEQ_CONSTRUCTOR_HPP
#define REF_SEQ_CONSTRUCTOR_HPP


#include <sparsepp/spp.h>

#include "Util.hpp"
#include "CanonicalKmer.hpp"

enum Task {
           SUCCESS,
           FAILURE
};

enum SearchType {
  PREDECESSOR,
  SUCCESSOR
};

struct nextCompatibleStruct {
  util::ContigBlock cntg ;
  size_t tpos ;
  uint32_t cpos ;
  bool isCurContigFw ;

  nextCompatibleStruct(util::ContigBlock cntgIn, size_t tposIn, uint32_t cposIn, bool mFw) : cntg(cntgIn), tpos(tposIn), cpos(cposIn), isCurContigFw(mFw) {} 
} ;


template <typename PufferfishIndexT>
class RefSeqConstructor {

public:
  RefSeqConstructor(PufferfishIndexT* pfi, spp::sparse_hash_map<uint32_t, util::ContigBlock>* contigSeqCache);
  Task fillSeq(size_t tid,
                                               size_t tpos,
                                               bool isCurContigFw,
                                               util::ContigBlock& curContig,
                                               size_t startp,
                                               util::ContigBlock& endContig,
                                               size_t endp,
                                               bool isEndContigFw,
                                               uint32_t txpDist,
                                               std::string& seq);
  Task doBFS(size_t tid,
             size_t tpos,
             bool isCurContigFw,
             util::ContigBlock& curContig,
             size_t startp,
             util::ContigBlock& endContig,
             bool isEndContigFw,
             uint32_t threshold,
             std::string& seq);

  //search predecessors
  Task fillSeqLeft(size_t tid,
                   size_t tpos,
                   util::ContigBlock& curContig,
                   bool isCurContigFw,
                   uint32_t cstart,
                   std::string& seq);

  Task doRevBFS(size_t tid,
                size_t tpos,
                util::ContigBlock& curContig,
                bool isCurContigFw,
                uint32_t cstart,
                uint32_t threshold,
                std::string& seq);

private:
  PufferfishIndexT* pfi_ ;
  size_t k ;
  spp::sparse_hash_map<uint32_t, util::ContigBlock>* contigSeqCache_;



  size_t remainingLen(util::ContigBlock& contig, size_t startp, bool isCurContigFw, bool fromTheEnd);
  void append(std::string& seq, util::ContigBlock& contig, size_t startp, size_t endp, bool isCurContigFw);
  void appendByLen(std::string& seq, util::ContigBlock& contig, size_t startp, size_t len, bool isCurContigFw, bool appendSuffix);
  //TODO 
  //void prependByLen(std::string& seq, util::ContigBlock& contig, size_t startp, size_t len, bool isCurContigFw, bool appendSuffix);
  std::string getRemSeq(util::ContigBlock& contig, size_t len, bool isCurContigFw, bool appendSuffix);
  void cutoff(std::string& seq, size_t len);
  std::string rc(std::string str);
  char rev(const char& c);
  std::vector<nextCompatibleStruct> fetchSuccessors(util::ContigBlock& contig,
                                                 bool isCurContigFw,
                                                 size_t tid,
                                                 size_t tpos);
  std::vector<nextCompatibleStruct> fetchPredecessors(util::ContigBlock& contig,
                                                    bool isCurContigFw,
                                                    size_t tid,
                                                    size_t tpos);
};

#endif
