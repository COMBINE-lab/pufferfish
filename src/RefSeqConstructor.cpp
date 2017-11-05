#include "RefSeqConstructor.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"

#include <sparsepp/spp.h>

#define verbose false

#define suffixIfFw true
#define prefixIfFw false

/* story 
NOTE: Since we are ALWAYS traversing the mems --and hence
the unmapped sequences-- in the forward direction wrt the transcript,
whenever we are moving backward in a contig that means that the contig
is in reverse orientation wrt the transcript, so the sequence we fetch from the contig
should be reverse-complemented!!!
EXCEPTION ::::: constructing the string before the first unimem :(((
TODO: make sure about the fact above startp : relative position in curContig
that toBeAligned string starts endp : relative position in endContig that toBeAligned string ends
txpDist : number of bases that should be fetched before we stop spreading the path through BFS
seq : will be appended through out the traversal
*/


template <typename PufferfishIndexT>
RefSeqConstructor<PufferfishIndexT>::RefSeqConstructor(PufferfishIndexT* pfi,
                                                       spp::sparse_hash_map<uint32_t,
                                                       util::ContigBlock>* contigSeqCache)
                                                       : pfi_(pfi), contigSeqCache_(contigSeqCache) { k = pfi_->k(); }
template <typename PufferfishIndexT>
Task RefSeqConstructor<PufferfishIndexT>::fillSeq(size_t tid,
                                             size_t tpos,
                                             bool isCurContigFw,
                                             util::ContigBlock& curContig,
                                             size_t startp,
                                             util::ContigBlock& endContig,
                                             size_t endp,
                                             bool isEndContigFw,
                                             uint32_t txpDist,
                                             std::string& seq) {
  if(verbose) std::cout << "GO time\n" << "called with curr contig: "
                        << curContig.contigIdx_ << " of length " << curContig.contigLen_
                        << " start " << startp << "\n"
                        << "end contig index "<< endContig.contigIdx_ << " ending at " << endp << "\n" ;



  if (curContig.contigIdx_ == endContig.contigIdx_) {
    append(seq, curContig, startp, endp, isCurContigFw);
    return Task::SUCCESS;
  }

  if (!endContig.isDummy()) {
    auto endRemLen = remainingLen(endContig, endp, isEndContigFw, prefixIfFw);
    if(verbose) std::cout << "END remaining length " << endRemLen << " txpDist:" << txpDist << " endp " <<endp << " in direction "<<isEndContigFw<<"\n" ;
    if (endRemLen >= txpDist) {
      appendByLen(seq, endContig, endp, txpDist, isEndContigFw, prefixIfFw);
      return Task::SUCCESS;
    } else if (endRemLen > 0) {
      appendByLen(seq, endContig, endp, endRemLen, isEndContigFw, prefixIfFw);
      txpDist -= endRemLen;
    }
  }

  if (verbose)
    std::cout << "\n\nWOOOOOT!! Got to bfs\n";
  return doBFS(tid, tpos, isCurContigFw, curContig, startp, endContig, isEndContigFw, txpDist, seq);

}

template <typename PufferfishIndexT>
Task RefSeqConstructor<PufferfishIndexT>::doBFS(size_t tid,
             size_t tpos,
             bool isCurContigFw,
             util::ContigBlock& curContig,
             size_t startp,
             util::ContigBlock& endContig,
             bool isEndContigFw,
             uint32_t txpDist,
             std::string& seq) {




    if (curContig.isDummy()) {
      //TODO should be taken care of (hard hard)
      return Task::FAILURE;
    }

    if (startp >= curContig.contigLen_)
      std::cerr << "ERROR!!! shouldn't happen ---> startp >= curContig.contigLen_ : " << startp << ">" << curContig.contigLen_ << "\n";
    // used in all the following terminal conditions
    auto remLen = remainingLen(curContig, startp, isCurContigFw, suffixIfFw);
    if(verbose) std::cout << "[doBFS] : remaining length " << remLen << " txpDist " << txpDist << " startp " <<startp << " in direction "<<isCurContigFw<<"\n" ;

    if (remLen >= txpDist) {
      if (endContig.isDummy()){
        appendByLen(seq, curContig, startp, txpDist, isCurContigFw, suffixIfFw);
        return Task::SUCCESS;
      }
      // DON'T GET STUCK IN INFINITE LOOPS
      // if we have not reached the last contigId
      // and also end of the path is NOT open
      // then if the remaining of the contig from this position is less than txpDist it should be counted as a failure
      // because we couldn't find the path between start and end that is shorter than some txpDist
      else {
        if (remLen > txpDist) {

          if (remLen-txpDist < endContig.contigLen_ &&
              getRemSeq(curContig, remLen-txpDist, isCurContigFw, suffixIfFw) == getRemSeq(endContig, remLen-txpDist, isEndContigFw, prefixIfFw)) {
            if (verbose) std::cout << "[doBFS] remSequences are equal. safe to append and exit\n"; 
            appendByLen(seq, curContig, startp, txpDist, isCurContigFw, suffixIfFw);
            return Task::SUCCESS;
          }
          else {
            if(verbose) std::cout << "[doBFS] returning failure\n";
            return Task::FAILURE;
          }
        }
        // terminal condition
        // called even when txpDist == 0
        // remLen == txpDist
        for (auto& c : fetchSuccessors(curContig, isCurContigFw, tid, tpos)) {
          if (c.cntg.contigLen_-(k-1) <= endContig.contigLen_ &&
              getRemSeq(c.cntg, c.cntg.contigLen_-(k-1), c.isCurContigFw, suffixIfFw) == getRemSeq(endContig, c.cntg.contigLen_-(k-1), isEndContigFw, prefixIfFw))
            return Task::SUCCESS;
        }
        return Task::FAILURE; // I'm in the middle of no where!! lost!!
      }
    }

    // if we are here, this means remLen < txpDist
    
    // If we didn't meet any terminal conditions, we need to dig deeper into the tree through BFS
    // The approach is the same for both valid and dummy end nodes
    if(verbose) std::cout << "[doBFS] remLen < txpDist : " << remLen << " < " << txpDist << "\n" ;

    appendByLen(seq, curContig, startp, remLen, isCurContigFw, suffixIfFw);
    
    txpDist -= remLen; // update txpDist
    for (auto& c : fetchSuccessors(curContig, isCurContigFw, tid, tpos)) {
      // act greedily and return with the first successfully constructed sequence.
      if (doBFS(tid, c.tpos, c.isCurContigFw, c.cntg, c.cpos, endContig, isEndContigFw, txpDist, seq) == Task::SUCCESS) {
        if(verbose) std::cout <<"[doBFS] got to success!!\n";
        return Task::SUCCESS;
      }
    }

    // If couldn't find any valid/compatible successors,
    //then this path was a dead end. Revert your append and return with failure
    cutoff(seq, remLen);
    return Task::FAILURE;
}


template <typename PufferfishIndexT>
      size_t RefSeqConstructor<PufferfishIndexT>::remainingLen(util::ContigBlock& contig, size_t startp, bool isCurContigFw, bool fromTheEnd) {
      if ( isCurContigFw == fromTheEnd)
      return contig.contigLen_ - startp - 1;
    else
      return startp ;
}


template <typename PufferfishIndexT>
void RefSeqConstructor<PufferfishIndexT>::append(std::string& seq,
                                                 util::ContigBlock& contig,
                                                 size_t startp, size_t endp,
                                                 bool isCurContigFw) {

  if(isCurContigFw) {
    if(verbose) std::cout << "\t[append] 1 " << seq << " clipping by pos: from "<<startp+1<< " to " << endp << " in a contig with len "<<contig.contigLen_<<"\n" ;
      seq += contig.substrSeq(startp+1, endp-startp-1);
  }
    else {
      if(verbose) std::cout << "\t[append] 2 rc " << seq << " clipping by pos: from "<<endp+1<< " to " << startp << " in a contig with len "<<contig.contigLen_<<"\n" ;
      // we are always building the seq by moving forward in transcript, so we always append (& never prepend) any substring that we construct
      seq += rc(contig.substrSeq(endp+1, startp-endp-1)); 
    }
}


template <typename PufferfishIndexT>
void RefSeqConstructor<PufferfishIndexT>::appendByLen(std::string& seq, util::ContigBlock& contig, size_t startp, size_t len, bool isCurContigFw, bool appendSuffix) {
  if (len == 0)
    return;
  if (isCurContigFw && appendSuffix) { // append suffix
    if(verbose) std::cout << "\t[appendByLen] 1 from " << startp+1 << " to " << startp+1+len << " total length " << contig.contigLen_ << "\n";
    seq += contig.substrSeq(startp+1, len);
  }
  else if (isCurContigFw && !appendSuffix) {// append prefix
    if(verbose) std::cout << "\t[appendByLen] 2 from " << startp-len << " to " << startp << " total length " << contig.contigLen_ << "\n";
    seq = contig.substrSeq(startp-len, len) + seq;
  }
  else if (!isCurContigFw && appendSuffix) {// append rc of the seq from the other end as a suffix
    if(verbose) std::cout << "\t[appendByLen] 3 rc from " << startp-len << " to " << startp << " total length " << contig.contigLen_ << "\n";
    seq += rc(contig.substrSeq(startp-len, len));
  }
  else if (!isCurContigFw && !appendSuffix) {// append rc of the seq as prefix
    if(verbose) std::cout << "\t[appendByLen] 4 rc from " << startp+1 << " to " << startp+1+len << " total length " << contig.contigLen_ << "\n";
    seq = rc(contig.substrSeq(startp+1, len)) + seq;
  }
}

template <typename PufferfishIndexT>
std::string RefSeqConstructor<PufferfishIndexT>::getRemSeq(util::ContigBlock& contig, size_t len, bool isCurContigFw, bool appendSuffix) {
  std::string seq = "";
  if (len == 0)
    return "";
  if (isCurContigFw && appendSuffix) {// append suffix
    if(verbose) std::cout << "\t[getRemSeq] 1 from " << contig.contigLen_-len << " to " << contig.contigLen_ << " total length " << contig.contigLen_ << "\n";
    seq += contig.substrSeq(contig.contigLen_-len, len);
  }
  else if (isCurContigFw && !appendSuffix) {// append prefix
    if(verbose) std::cout << "\t[getRemSeq] 2 from " << 0 << " to " << len << " total length " << contig.contigLen_ << "\n";
    seq = contig.substrSeq(0, len) + seq;
  }
  else if (!isCurContigFw && appendSuffix) {// append rc of the seq from the other end as a suffix
    if(verbose) std::cout << "\t[getRemSeq] 3 rc from " << 0 << " to " << len << " total length " << contig.contigLen_ << "\n";
    seq += rc(contig.substrSeq(0, len));
  }
  else if (!isCurContigFw && !appendSuffix) {// append rc of the seq as prefix
    if(verbose) std::cout << "\t[getRemSeq] 4 rc from " << contig.contigLen_-len << " to " << contig.contigLen_ << " total length " << contig.contigLen_ << "\n";
    seq = rc(contig.substrSeq(contig.contigLen_-len, len)) + seq;
  }
  return seq;
}

template <typename PufferfishIndexT>
void RefSeqConstructor<PufferfishIndexT>::cutoff(std::string& seq, size_t len) {
    seq = seq.substr(0, seq.length()-len);
}


//TODO couldn't make it reference bc of the section seq += rc(str)
// now eachtime calling this we are copying a string twice which is not good
template <typename PufferfishIndexT>
std::string RefSeqConstructor<PufferfishIndexT>::rc(std::string str) {
  if(verbose)std::cout << "\t[rc] of " << str << " is ";
    for (uint32_t i = 0; i < str.length()/2; i++) {
      char tmp = str[i];
      str[i] = rev(str[str.length()-1-i]);
      str[str.length()-1-i] = rev(tmp);
    }
    if (verbose) std::cout << str << "\n";
    return str;
}


template <typename PufferfishIndexT>
char RefSeqConstructor<PufferfishIndexT>::rev(const char& c) {
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
    return 'N';
}

//TODO the vector should be passed by reference as argument
template <typename PufferfishIndexT>
std::vector<nextCompatibleStruct> RefSeqConstructor<PufferfishIndexT>::fetchSuccessors(util::ContigBlock& contig,
                                                 bool isCurContigFw,
                                                 size_t tid,
                                                 size_t tpos) {

    if(verbose) std::cout << "\t[Fetch successors] \n" ;
    std::vector<nextCompatibleStruct> successors ;
    CanonicalKmer::k(k) ;

    auto& edges = pfi_->getEdge() ;
    util::Direction dir = isCurContigFw?util::Direction::FORWARD:util::Direction::BACKWORD ;

    uint8_t edgeVec = edges[contig.contigIdx_] ;
    std::vector<util::extension> ext = util::getExts(edgeVec) ;

    if(!ext.empty()){
      CanonicalKmer kb ;
      kb.fromStr(contig.seq.substr(0,k)) ;
      CanonicalKmer ke ;
      ke.fromStr(contig.seq.substr(contig.contigLen_-k,k)) ;
      CanonicalKmer kt ;

      for(auto& ed : ext){
        if(ed.dir == dir){
          (dir == util::Direction::FORWARD)?ke.shiftFw(ed.c):kb.shiftBw(ed.c) ;
          (dir == util::Direction::FORWARD)?kt.fromNum(ke.getCanonicalWord()):kt.fromNum(kb.getCanonicalWord()) ;

          auto nextHit = pfi_->getRefPos(kt) ;

          if(contigSeqCache_->find(nextHit.contigIdx_) == contigSeqCache_->end()){
            std::string contigseq = pfi_->getSeqStr(nextHit.globalPos_,nextHit.contigLen_);
            (*contigSeqCache_)[nextHit.contigIdx_] = {nextHit.contigIdx_, nextHit.globalPos_, nextHit.contigLen_, contigseq} ;
          }

          util::ContigBlock cb = (*contigSeqCache_)[nextHit.contigIdx_] ;
          //TODO seq for cb
          for(auto& posIt: nextHit.refRange){
            if(posIt.transcript_id() == tid and tpos < posIt.pos()){
              successors.emplace_back(cb, posIt.pos(), posIt.orientation()?k-1:cb.contigLen_-k, posIt.orientation()) ;
            }
          }
        }
      }
    }

    return successors;
}


template class RefSeqConstructor<PufferfishIndex>;
template class RefSeqConstructor<PufferfishSparseIndex>;
