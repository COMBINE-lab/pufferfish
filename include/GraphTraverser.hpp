#ifndef GRAPH_TRAVERSER_HPP
#define GRAPH_TRAVERSER_HPP


#include "Util.hpp"

enum Task {
           SUCCESS,
           FAILURE
};





template <typename PufferfishIndexT> class GraphTraverser {
public:
  GraphTraverser(PufferfishIndexT& pfi) : pfi_(pfi) { k = pfi_->k(); }

  Task doBFS(size_t tid, size_t tpos, bool moveFw, util::ContigBlock curContig, size_t startp, util::ContigBlock endContig, size_t endp, uint32_t threshold, std::string& seq) {}


private:
  std::vector<util::ContigBlock> fetchSuccessors(util::ContigBlock& contig, bool moveFw) { return NULL;}

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
