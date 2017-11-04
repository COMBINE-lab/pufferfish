#include "RefSeqConstructor.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"

#include <sparsepp/spp.h>

#define verbose true

/* story 
NOTE: Since we are ALWAYS traversing the mems --and hence
the unmapped sequences-- in the forward direction wrt the transcript,
whenever we are moving backward in a contig that means that the contig
is in reverse orientation wrt the transcript, so the sequence we fetch from the contig
should be reverse-complemented!!!
EXCEPTION ::::: constructing the string before the first unimem :(((
TODO: make sure about the fact above startp : relative position in curContig
that toBeAligned string starts endp : relative position in endContig that toBeAligned string ends
threshold : number of bases that should be fetched before we stop spreading the path through BFS
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
             bool moveFw,
             util::ContigBlock& curContig,
             size_t startp,
             util::ContigBlock& endContig,
             size_t endp,
             uint32_t threshold,
             std::string& seq) {



  if(verbose) std::cout << "Doing BFS, GO time\n" << "called with curr contig: "
            << curContig.contigIdx_ << " of length " << curContig.contigLen_
            << " start " << startp << "\n"
            << "end contig index "<< endContig.contigIdx_ << " ending at " << endp << "\n" ;



    if (curContig.contigIdx_ == endContig.contigIdx_) {
        append(seq, curContig, startp, endp, moveFw);
        return Task::SUCCESS;
    }

    if (curContig.isDummy()) {
      //TODO should be taken care of (hard hard)
      return Task::FAILURE;
    }

    // used in all the following terminal conditions
    auto remLen = remainingLen(curContig, startp, moveFw);
    if(verbose) std::cout << "remaining length " << remLen <<" startp " <<startp << " in direction "<<moveFw<<"\n" ;

    if (remLen >= threshold) {
      if (endContig.isDummy()) {// if the end of the path is open
        appendByLen(seq, curContig, startp, threshold, moveFw);
        return Task::SUCCESS;
      }
      // DON'T GET STUCK IN INFINITE LOOPS
      // if we have not reached the last contigId
      // and also end of the path is NOT open
      // then if the remaining of the contig from this position is less than threshold it should be counted as a failure
      // because we couldn't find the path between start and end that is shorter than some threshold
      else
        return Task::FAILURE; // I'm in the middle of no where!! lost!!
    }

    // if we are here, this means threshold > remLen
    
    // If we didn't meet any terminal conditions, we need to dig deeper into the tree through BFS
    // The approach is the same for both valid and dummy end nodes
    if(verbose) std::cout << " Okay until here \n" ;

    appendByLen(seq, curContig, startp, remLen, moveFw);
    
    if(verbose) std::cout << " never got printed in case of seg fault  \n" ;

    threshold -= remLen; // update threshold
    for (auto& c : fetchSuccessors(curContig, moveFw, tid, tpos)) {
      // act greedily and return with the first successfully constructed sequence.
      if (doBFS(tid, c.tpos, c.moveFw, c.cntg, c.cpos, endContig, endp, threshold, seq) == Task::SUCCESS)
        return Task::SUCCESS;
    }

    // If couldn't find any valid/compatible successors,
    //then this path was a dead end. Revert your append and return with failure
    cutoff(seq, remLen);
    return Task::FAILURE;
}


template <typename PufferfishIndexT>
size_t RefSeqConstructor<PufferfishIndexT>::distance(size_t startp,
                                                     size_t endp, bool moveFw)
{if(moveFw)return endp-startp; else return startp-endp;}


template <typename PufferfishIndexT>
size_t RefSeqConstructor<PufferfishIndexT>::remainingLen(util::ContigBlock& contig, size_t startp, bool moveFw) {
    if (moveFw)
      return contig.contigLen_ - startp;
    else
      //TODO For backward walk, check the inclusion/exclusion of the current base and have a concensus on the assumption
      //this is a bad idea we might fall over the edge 
      //return startp+1;
      return startp ;
}


template <typename PufferfishIndexT>
void RefSeqConstructor<PufferfishIndexT>::append(std::string& seq,
                                                 util::ContigBlock& contig,
                                                 size_t startp, size_t endp,
                                                 bool moveFw) {

  if(verbose) std::cout << "clipping by pos "<<startp<< " to " << endp << " from "<<contig.contigLen_<<"\n" ;
   if(moveFw)
      seq += contig.substrSeq(startp, endp-startp);
    else {
      // we are always building the seq by moving forward in transcript, so we always append any substring that we construct
      seq += rc(contig.substrSeq(endp, startp-endp)); 
    }
}


template <typename PufferfishIndexT>
void RefSeqConstructor<PufferfishIndexT>::appendByLen(std::string& seq, util::ContigBlock& contig, size_t startp, size_t len, bool moveFw) {
  if(verbose) std::cout << "clipping by len "<<len<< " from " << startp << " total length: "<<contig.contigLen_<<"\n" ;
    if(moveFw)
      seq += contig.substrSeq(startp, len);
    else
      seq += rc(contig.substrSeq(startp-len, len));
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
                                                 bool moveFw,
                                                 size_t tid,
                                                 size_t tpos) {

    if(verbose) std::cout << "Fetch successors \n" ;
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
