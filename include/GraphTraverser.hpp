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
  //NOTE: Since we are ALWAYS traversing the mems and hence the unmapped sequences in transcript in forward direction,
  //whenever we are moving backward in a contig that means that the contig is in reverse orientation wrt the transcript, so the sequence we fetch from the contig
  //should be reverse-complemented!!!
  //TODO: make sure about the fact above
  Task doBFS(size_t tid, size_t tpos, bool moveFw, util::ContigBlock& curContig, size_t startp, util::ContigBlock& endContig, size_t endp, uint32_t threshold, std::string& seq) {
    if (curContig.contigIdx_ == endContig.contigIdx_) {
        append(seq, curContig, startp, endp, moveFw);
        return Task::SUCCESS;
    }
    if (endContig.isDummy()) { // if the end of the path is open
      auto& dist = distance(startp, endp, moveFw);
      if (dist>threshold) {
        appendByLen(seq, curContig, startp, threshold, moveFw);
        return Task::SUCCESS;
      } else {
        appendByLen(seq, curContig, startp, dist, moveFw); // update seq
        threshold -= dist; // update threshold
        for (auto& c : fetchSuccessors(curContig, moveFw, tid, tpos)) {
          // act greedily and return with the first successfully constructed sequence.
          if (doBFS(tid, c.tpos, c.moveFw, c.cntg, c.cpos, endContig, endp, threshold, seq) == Task::SUCCESS)
            return Task::SUCCESS;
        }
        // If couldn't find any successful successors, then this path was a dead end. Revert your append and return with failure
        cutoff(seq, dist);
        return Task::FAILURE;
      }
    }

    // used for the last two conditions
    auto& remLen = remainingLen(curContig, startp, moveFw);
    // DON'T GET STUCK IN INFINITE LOOPS
    // if we have not reached the last contigId
    // and also end of the path is NOT open
    // then if the remaining of the contig from this position is less than threshold it should be counted as a failure
    // because we couldn't find the path between start and end that is shorter than some threshold
    if ( remLen < threshold) {
      return Task::FAILURE; // I'm in the middle of no where!! lost!!
    }

    // last but not least
    // fetch al the remaining string from the current contig and do BFS for its successors to either get to the end contig or exhaust the threshold
    appendByLen(seq, curContig, startp, remLen, moveFw);
    threshold -= remLen; // update threshold
    for (auto& c : fetchSuccessors(curContig, moveFw, tid, tpos)) {
      // act greedily and return with the first successfully constructed sequence.
      if (doBFS(tid, c.tpos, c.moveFw, c.cntg, c.cpos, endContig, endp, threshold, seq) == Task::SUCCESS)
        return Task::SUCCESS;
    }
    // If couldn't find any successful successors, then this path was a dead end. Revert your append and return with failure
    cutoff(seq, remLen);
    return Task::FAILURE;
 }


private:
  PufferfishIndexT* pfi_ ;
  size_t k ;

  struct nextCompatibleStruct{
    util::ContigBlock cntg ;
    uint32_t tpos ;
    uint32_t cpos ;
    bool moveFw ;

   nextCompatibleStruct(util::ContigBlock cntgIn, uint32_t tposIn, uint32_t cposIn, bool mFw) : cntg(cntgIn), tpos(tposIn), cpos(cposIn), moveFw(mFw) {} 
  } ;


  size_t distance(size_t startp, size_t endp, bool moveFw) {if(moveFw)return endp-startp; else return startp-endp;}
  size_t remainingLen(util::ContigBlock& contig, size_t startp, bool moveFw) {
    if (moveFw)
      return contig.contigLen_ - startp;
    else
      //TODO For backward walk, check the inclusion/exclusion of the current base and have a concensus on the assumption
      return startp+1;
  }
  void append(std::string& seq, util::ContigBlock& contig, size_t startp, size_t endp, bool moveFw) {
    if(moveFw)
      seq += contig.substrSeq(startp, endp-startp);
    else {
      // we are always building the seq by moving forward in transcript, so we always append any substring that we construct
      seq += rc(contig.substrSeq(endp, startp-endp)); 
    }
  }
  void appendByLen(std::string& seq, util::ContigBlock& contig, size_t startp, size_t len, bool moveFw) {
    if(moveFw)
      seq += contig.substrSeq(startp, len);
    else
      seq += rc(contig.substrSeq(startp-len, len));
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

  std::vector<nextCompatibleStruct> fetchSuccessors(util::ContigBlock& contig,
                                                 bool moveFw,
                                                 size_t tid,
                                                 size_t tpos) {
    std::vector<nextCompatibleStruct> successors ;
    CanonicalKmer::k(k) ;

    auto& edges = pfi_->getEdge() ;
    util::Direction dir = moveFw?util::Direction::FORWARD:util::Direction::BACKWORD ;

    uint8_t edgeVec = edges[contig.contigIdx_] ;
    std::vector<util::extension> ext = util::getExts(edgeVec) ;

    if(!ext.empty()){
      CanonicalKmer kb ;
      kb.fromStr(contig.seq.substr(0,k)) ;
      CanonicalKmer ke ;
      ke.fromStr(contig.seq.substr(contig.contigLen_-k-1,k)) ;
      CanonicalKmer kt ;

      for(auto& ed : ext){
        if(ed.dir == dir){
          (dir == util::Direction::FORWARD)?ke.shiftFw(ed.c):kb.shiftBw(ed.c) ;
          (dir == util::Direction::FORWARD)?kt.fromNum(ke.getCanonicalWord()):kt.fromNum(kb.getCanonicalWord()) ;

          auto& nextHit = pfi_->getRefPos(kt) ;
          auto& refs = nextHit.refRange ;
          std::string contigseq = pfi_->getSeqStr(nextHit.globalPos_,nextHit.contigLen_);
          util::ContigBlock cb({nextHit.contigIdx_, nextHit.globalPos_, nextHit.contigLen_, contigseq}) ;
          //TODO seq for cb
          for(auto& posIt: refs){
            if(posIt.transcript_id() == tid and tpos < posIt.position()){
              successors.emplace_back({cb, posIt.position(), posIt.orientation()?k-1:cb.contigLen_-k, posIt.orientation()}) ;
            }
          }
        }
      }
    }

    return successors;
  }

};

#endif
