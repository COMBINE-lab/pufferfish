#ifndef GRAPH_TRAVERSER_HPP
#define GRAPH_TRAVERSER_HPP


#include "Util.hpp"
#include "CanonicalKmer.hpp"

enum Task {
           SUCCESS,
           FAILURE
};





template <typename PufferfishIndexT> class GraphTraverser {
public:
  GraphTraverser(PufferfishIndexT& pfi) : pfi_(pfi) { k = pfi_->k(); }

  Task doBFS(size_t tid, size_t tpos, bool moveFw, util::ContigBlock curContig, size_t startp, util::ContigBlock endContig, size_t endp, uint32_t threshold, std::string& seq) {}


private:
  std::vector<util::ContigBlock> fetchSuccessors(util::ContigBlock& contig,
                                                 bool moveFw,
                                                 size_t tid,
                                                 size_t tpos) {
    std::vector<util::ContigBlock> successors ;
    CanonicalKmer::k(k) ;

    auto& edges = pfi_->getEdge() ;
    util::Direction dir = moveFw?util::Direction::FORWARD:util::Direction::BACKWORD ;

    uint8_t edgeVec = edges[contig.contigIdx_] ;
    std::vector<util::extension> ext = util::getExts(edgeVec) ;

    if(!ext.empty()){
      CanonicalKmer kb = contig.kb ;
      CanonicalKmer ke = contig.ke ;
      CanonicalKmer kt ;

      for(auto& ed : ext){
        if(ed.dir == dir){
          (dir == util::Direction::FORWARD)?ke.shiftFw(ed.c):kb.shiftBw(ed.c) ;
          (dir == util::Direction::FORWARD)?kt.fromNum(ke.getCanonicalWord()):kt.fromNum(kb.getCanonicalWord()) ;

          auto& nextHit = pfi_->getRefPos(kt) ;
          util::ContigBlock ctgBlock = pfi_->getContigBlock(nextHit.contigIdx_) ;
          if(isCompatible(ctgBlock,tid,tpos,moveFw))
            successors.push_back(ctgBlock) ;

        }
      }
    }

    return successors;
  }

  bool isCompatible(util::ContigBlock& contig, size_t tid, size_t tpos, bool fw) {
    auto& pvec = pfi_->refList(contig.contigIdx_) ;
    bool comp{false} ;
    
    for(auto& posIt : pvec){
      if(posIt.transcript_id() == tid){
        comp = fw?(posIt.pos() > tpos):(posIt.pos() < tpos) ;
        if(comp){
          return true ;
        }
      }
    }

    return false;
  }


  PufferfishIndexT* pfi_ ;
  size_t k ;

};

#endif
