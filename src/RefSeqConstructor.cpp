#include "RefSeqConstructor.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"

#include <sparsepp/spp.h>

#define verbose true

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
                                                       util::ContigBlock>& contigCache)
                                                       : pfi_(pfi), contigCache_(contigCache) { k = pfi_->k(); }
template <typename PufferfishIndexT>
Task RefSeqConstructor<PufferfishIndexT>::doBFS(size_t tid,
             size_t tpos,
             bool isCurContigFw,
             util::ContigBlock& curContig,
             size_t startp,
             util::ContigBlock& endContig,
             size_t endp,
                                                bool isEndContigFw,
             uint32_t txpDist,
             std::string& seq) {



  if(verbose) std::cout << "Doing BFS, GO time\n" << "called with curr contig: "
            << curContig.contigIdx_ << " of length " << curContig.contigLen_
            << " start " << startp << "\n"
            << "end contig index "<< endContig.contigIdx_ << " ending at " << endp << "\n" ;



    if (curContig.contigIdx_ == endContig.contigIdx_) {
        append(seq, curContig, startp, endp, isCurContigFw);
        return Task::SUCCESS;
    }

    if (curContig.isDummy()) {
      //TODO should be taken care of (hard hard)
      return Task::FAILURE;
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
        endp = 0;
      }
    }
    
    // used in all the following terminal conditions
    auto remLen = remainingLen(curContig, startp, isCurContigFw, suffixIfFw);
    if(verbose) std::cout << "remaining length " << remLen <<" startp " <<startp << " in direction "<<isCurContigFw<<"\n" ;

    if (remLen >= txpDist) {
      if (endContig.isDummy()){// || curContig.getRemSeq(startp, remLen, false) == endContig.getRemSeq(endp, remLen, true)) {// if the end of the path is open
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
          if (getSubstr(curContig, curContig.contigLen_ - remLen, remLen, isCurContigFw, suffixIfFw) == getSubstr(endContig, 0, remLen, isEndContigFw, prefixIfFw)) {
            appendByLen(seq, curContig, startp, txpDist, isCurContigFw, suffixIfFw);
            return Task::SUCCESS;
          }
          else
            return Task::FAILURE;
        }
        // remLen == txpDist
        for (auto& c : fetchSuccessors(curContig, isCurContigFw, tid, tpos)) {
          if (getSubstr(c.cntg, c.cntg.contigLen_ - remLen, remLen, c.isCurContigFw, suffixIfFw) == getSubstr(endContig, 0, remLen, isEndContigFw, prefixIfFw))
            return Task::SUCCESS;
        }
        return Task::FAILURE; // I'm in the middle of no where!! lost!!
      }
    }

    // if we are here, this means remLen < txpDist
    
    // If we didn't meet any terminal conditions, we need to dig deeper into the tree through BFS
    // The approach is the same for both valid and dummy end nodes
    if(verbose) std::cout << " Okay until here \n" ;

    appendByLen(seq, curContig, startp, remLen, isCurContigFw, suffixIfFw);
    
    if(verbose) std::cout << " never got printed in case of seg fault  \n" ;

    txpDist -= remLen; // update txpDist
    for (auto& c : fetchSuccessors(curContig, isCurContigFw, tid, tpos)) {
      // act greedily and return with the first successfully constructed sequence.
      if (doBFS(tid, c.tpos, c.isCurContigFw, c.cntg, c.cpos, endContig, endp, isEndContigFw, txpDist, seq) == Task::SUCCESS)
        return Task::SUCCESS;
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

  if(verbose) std::cout << "clipping by pos "<<startp<< " to " << endp << " from "<<contig.contigLen_<<"\n" ;
   if(isCurContigFw)
      seq += contig.substrSeq(startp+1, endp-startp-1);
    else {
      // we are always building the seq by moving forward in transcript, so we always append (& never prepend) any substring that we construct
      seq += rc(contig.substrSeq(endp+1, startp-endp-1)); 
    }
}


template <typename PufferfishIndexT>
void RefSeqConstructor<PufferfishIndexT>::appendByLen(std::string& seq, util::ContigBlock& contig, size_t startp, size_t len, bool isCurContigFw, bool appendSuffix) {
  if(verbose) std::cout << "clipping by len "<<len<< " from " << startp << " total length: "<<contig.contigLen_<<"\n" ;
  if (len == 0)
    return;
  if (isCurContigFw && appendSuffix) // append suffix
    seq += contig.substrSeq(startp+1, len);
  else if (isCurContigFw && !appendSuffix) // append prefix
    seq = contig.substrSeq(startp-len, len) + seq;
  else if (!isCurContigFw && appendSuffix) // append rc of the seq from the other end as a suffix
    seq += rc(contig.substrSeq(startp-len, len));
  else if (!isCurContigFw && !appendSuffix) // append rc of the seq as prefix
    seq = rc(contig.substrSeq(startp+1, len)) + seq;
}

template <typename PufferfishIndexT>
std::string RefSeqConstructor<PufferfishIndexT>::getSubstr(util::ContigBlock& contig, size_t startp, size_t len, bool isCurContigFw, bool appendSuffix) {
  std::string seq;
  if(verbose) std::cout << "clipping by len "<<len<< " from " << startp << " total length: "<<contig.contigLen_<<"\n" ;
  if (len == 0)
    return "";
  if (isCurContigFw && appendSuffix) // append suffix
    seq += contig.substrSeq(startp+1, len);
  else if (isCurContigFw && !appendSuffix) // append prefix
    seq = contig.substrSeq(startp-len, len) + seq;
  else if (!isCurContigFw && appendSuffix) // append rc of the seq from the other end as a suffix
    seq += rc(contig.substrSeq(startp-len, len));
  else if (!isCurContigFw && !appendSuffix) // append rc of the seq as prefix
    seq = rc(contig.substrSeq(startp+1, len)) + seq;
}

template <typename PufferfishIndexT>
void RefSeqConstructor<PufferfishIndexT>::cutoff(std::string& seq, size_t len) {
    seq = seq.substr(0, seq.length()-len);
}


//TODO couldn't make it reference bc of the section seq += rc(str)
// now eachtime calling this we are copying a string twice which is not good
template <typename PufferfishIndexT>
std::string RefSeqConstructor<PufferfishIndexT>::rc(std::string str) {
    for (uint32_t i = 0; i < str.length()/2; i++) {
      char tmp = str[i];
      str[i] = rev(str[str.length()-1-i]);
      str[str.length()-1-i] = rev(tmp);
    }
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

    if(verbose) std::cout << "Fetch successors \n" ;
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
          std::string contigseq = pfi_->getSeqStr(nextHit.globalPos_,nextHit.contigLen_);
          util::ContigBlock cb({nextHit.contigIdx_, nextHit.globalPos_, nextHit.contigLen_, contigseq}) ;
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
