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
  Task doBFS(size_t tid, size_t tpos, bool moveFw, util::ContigBlock& curContig, size_t startp, util::ContigBlock& endContig, size_t endp, uint32_t threshold, std::string& seq) {
    if (curContig.contigIdx_ == endContig.contigIdx_) {
      if (!endContig.isDummy()) {
        append(seq, curContig, startp, endp, moveFw);
        return Task::SUCCESS;
      }
      else { // if the end contig is a dummy
        
      }
    }

    if (!isCompatible(curContig, tid, tpos)) return Task::FAILURE;

    return Task::SUCCESS;
 }


private:
  PufferfishIndexT* pfi_ ;
  size_t k ;


  size_t distance(size_t startp, size_t endp, bool moveFw) {if(moveFw)return endp-startp; else return startp-endp;}
  void append(std::string& seq, util::ContigBlock& contig, size_t startp, size_t endp, bool moveFw) {
    if(moveFw)
      seq += contig.substrSeq(startp, endp-startp);
    else {
      seq += rc(contig.substrSeq(endp, startp-endp)); // we are always building the seq by moving forward in transcript, so we always append any substring that we construct
    }
  }
  void cutoff(std::string& seq, size_t len) {
    seq = seq.substr(0, seq.length()-len);
  }

  std::string& rc(std::string& str) {
    for (int i = 0; i < str.length()/2; i++) {
      char tmp = str[i];
      str[i] = rev(str[str.length()-1-i]);
      str[str.length()-1-i] = rev(tmp);
    }
    return str;
  }

  char rev(const char& c) {
    switch(c){
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    }
  }

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


};

#endif
