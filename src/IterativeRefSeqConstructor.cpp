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
                                                       util::ContigBlock>* contigSeqCache)
                                                       : pfi_(pfi), contigSeqCache_(contigSeqCache) { k = pfi_->k(); }

/**
   tid : txp id
   tpos : projected txp start position of the sequence we want to fetch
   isCurContigFw : If the start contig is forward
   startp : position of the base before the region we want to fetch the sequence for in current contig (position on the contig)

 **/
template <typename PufferfishIndexT>
std::string RefSeqConstructor<PufferfishIndexT>::fillSeq(size_t tid,
                                             size_t tpos,
                                             bool isCurContigFw,
                                             util::ContigBlock& curContig,
                                             uint32_t startp,
                                             util::ContigBlock& endContig,
                                             uint32_t endp,
                                             bool isEndContigFw,
                                             uint32_t txpDist) {
 

  std::string seq = "";
  if (curContig.contigIdx_ == endContig.contigIdx_) {
    append(seq, curContig, startp, endp, isCurContigFw);
    return seq;
  }

  if (!endContig.isDummy()) {
    auto endRemLen = remainingLen(endContig, endp, isEndContigFw, prefixIfFw);
    //if(curContig.isDummy()) std::cerr << std::this_thread::get_id() << " "  << "END remaining length " << endRemLen << " txpDist:" << txpDist << " endp " <<endp << " in direction "<<isEndContigFw<<"\n" ;
    if (endRemLen >= txpDist) {
      auto pair = appendByLen(endContig, endp, txpDist, isEndContigFw, prefixIfFw);
      seq += pair.first;
      return seq;
    } else if (endRemLen > 0) {
      auto pair = appendByLen(endContig, endp, endRemLen, isEndContigFw, prefixIfFw);
      if (pair.second)
        seq += pair.first;
      else
        seq = pair.first + seq;
      txpDist -= endRemLen;
    }
  }

  if (verbose) std::cerr << std::this_thread::get_id() << " "  << "\n\nWOOOOOT!! Got to bfs\n";
  if (curContig.isDummy()) {
    return seq;
  //std::cerr << std::this_thread::get_id() << " "  << "left Dummy\n";
    //return doRevBFS(tid, tpos, isCurContigFw, curContig, startp, endContig, isEndContigFw, txpDist, seq);
  }
  std::string temp = doDFS(tid, tpos, isCurContigFw, curContig, startp, endContig, isEndContigFw, txpDist);
  seq += temp;
  return seq;
}

template <typename PufferfishIndexT>
std::string RefSeqConstructor<PufferfishIndexT>::doDFS(size_t tid,
             size_t tstart,
             bool isCurContigFw,
             util::ContigBlock& curContig,
             uint32_t cstart,
             util::ContigBlock& endContig,
             bool isEndContigFw,
             uint32_t txpDist) {

  std::stack<nextCompatibleStruct> cidStack;
  std::stack<strInfo> path;
  //uint32_t originalTxpDist = txpDist;

  if (verbose) std::cerr << "end contig index "<< endContig.contigIdx_ << " is dummy: " << endContig.isDummy() << "\ntxpDist: " << txpDist << "\n";


  cidStack.push({curContig.contigIdx_, tstart, cstart, isCurContigFw});
  bool success = false;
  while (!cidStack.empty()) {
    if (!path.empty() && path.top().pid == cidStack.size())
      path.pop();
    auto cInfo = cidStack.top();
    cidStack.pop();
    size_t pid = cidStack.size();
    util::ContigBlock& cb = (*contigSeqCache_)[cInfo.cid] ;
    if ( cInfo.cpos >= cb.contigLen_)
      std::cerr << std::this_thread::get_id() << " "  << "ERROR!!! shouldn't happen ---> startp >= curContig.contigLen_ : " << cInfo.cpos << ">" << cb.contigLen_ << "\n";
    // used in all the following terminal conditions
    auto remLen = remainingLen(cb, cInfo.cpos, cInfo.isCurContigFw, suffixIfFw);
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\n[doBFS]\n" << "called for txp " << tid << " with pos " << cInfo.tpos << " with curr contig: "
                          << cb.contigIdx_ << " of length " << cb.contigLen_ << " ori " << cInfo.isCurContigFw
                          << " start " << cInfo.cpos <<  " remLen:" << remLen << " txpDist: " << txpDist << "\n";
    if (remLen >= txpDist) {
      if (endContig.isDummy()){
        auto cpair = appendByLen(cb, cInfo.cpos, txpDist, cInfo.isCurContigFw, suffixIfFw);
        path.push({pid, cpair.first, cpair.second});
        if (verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[After append] " << cpair.first << "\n";
        success = true;
        break;
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
            auto cpair = appendByLen(cb, cInfo.cpos, txpDist, cInfo.isCurContigFw, suffixIfFw);
            path.push({pid, cpair.first, cpair.second});
            if (verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[After append2] " << cpair.first << "\n";
            success = true;
            break;
          }
          else {
            if(verbose) std::cerr << std::this_thread::get_id() << " "  << "[doBFS] returning failure\n";
            continue;
          }
        }
        // terminal condition
        // called even when txpDist == 0
        // remLen == txpDist
        bool breakFromOuter = false;
        for (auto& c : fetchSuccessors(cb, cInfo.isCurContigFw, tid, tstart, txpDist)) {
          util::ContigBlock& nextCB = (*contigSeqCache_)[c.cid];
          if (nextCB.contigLen_-(k-1) <= endContig.contigLen_ &&
              getRemSeq(nextCB, nextCB.contigLen_-(k-1), c.isCurContigFw, suffixIfFw) == getRemSeq(endContig, nextCB.contigLen_-(k-1), isEndContigFw, prefixIfFw)) {
            auto cpair = appendByLen(cb, cInfo.cpos, txpDist, cInfo.isCurContigFw, suffixIfFw);
            path.push({pid, cpair.first, cpair.second});
            breakFromOuter = true;
            break;
          }
        }
        if (breakFromOuter) {
          success = true;
          break;
        }
        continue; // I'm in the middle of no where!! lost!!
      }
    }

    // if we are here, this means remLen < txpDist

    // If we didn't meet any terminal conditions, we need to dig deeper into the tree through BFS
    // The approach is the same for both valid and dummy end nodes
    //if(verbose) std::cerr << std::this_thread::get_id() << " "  << "[doBFS] remLen < txpDist : " << remLen << " < " << txpDist << "\n" ;
    
    auto cpair = appendByLen(cb, cInfo.cpos, remLen, cInfo.isCurContigFw, suffixIfFw);
    path.push({pid, cpair.first, cpair.second});

    txpDist -= remLen; // update txpDist
    for (auto& c : fetchSuccessors(cb, cInfo.isCurContigFw, tid, cInfo.tpos, txpDist)) {
      // act greedily and return with the first successfully constructed sequence.
      if (verbose) std::cerr << std::this_thread::get_id() << " "  << "\nRecursive txp " << tid << " prev tpos " << tstart << " new tpos " << c.tpos << " txpDist: " << txpDist << "\n";
      cidStack.push(c);
    }
  }
  std::string seq = "";
  if (success) {
    while (!path.empty()) {
      auto s = path.top();
      path.pop();
      if (s.shouldAppend)
        seq += s.str;
      else
        seq = s.str + seq;
    }
  }
  return seq;
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
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[append] 1 " << seq << " clipping by pos: from "<<startp+1<< " to " << endp << " in a contig with len "<<contig.contigLen_
                          << " str len: " << contig.seq.length() << "\n" ;
      seq += contig.substrSeq(startp+1, endp-startp-1);
  }
    else {
      if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[append] 2 rc " << seq << " clipping by pos: from "<<endp+1<< " to " << startp << " in a contig with len "<<contig.contigLen_<<"\n" ;
      // we are always building the seq by moving forward in transcript, so we always append (& never prepend) any substring that we construct
      seq += rc(contig.substrSeq(endp+1, startp-endp-1)); 
    }
}


template <typename PufferfishIndexT>
std::pair<std::string, bool> RefSeqConstructor<PufferfishIndexT>::appendByLen(util::ContigBlock& contig, size_t startp, size_t len, bool isCurContigFw, bool appendSuffix) {
  if (len == 0)
    return std::make_pair("", true);
  if (isCurContigFw && appendSuffix) { // append suffix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[appendByLen] 1 from " << startp+1 << " to " << startp+1+len << " total length " << contig.contigLen_ << "\n";
    return std::make_pair(contig.substrSeq(startp+1, len), true);
  }
  else if (isCurContigFw && !appendSuffix) {// append prefix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[appendByLen] 2 from " << startp-len << " to " << startp << " total length " << contig.contigLen_ << "\n";
    return std::make_pair(contig.substrSeq(startp-len, len), false);
  }
  else if (!isCurContigFw && appendSuffix) {// append rc of the seq from the other end as a suffix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[appendByLen] 3 rc from " << startp-len << " to " << startp << " total length " << contig.contigLen_ << "\n";
    return std::make_pair(rc(contig.substrSeq(startp-len, len)), true);
  }
  else if (!isCurContigFw && !appendSuffix) {// append rc of the seq as prefix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[appendByLen] 4 rc from " << startp+1 << " to " << startp+1+len << " total length " << contig.contigLen_ << "\n";
    return std::make_pair(rc(contig.substrSeq(startp+1, len)), false);
  }
  return std::make_pair("", true);
}


//This will give warning b/c of empty function
/*
template <typename PufferfishIndexT>
void RefSeqConstructor<PufferfishIndexT>::prependByLen(std::string& seq, util::ContigBlock& contig, size_t startp, size_t len, bool isCurContigFw, bool appendSuffix){

}*/







template <typename PufferfishIndexT>
std::string RefSeqConstructor<PufferfishIndexT>::getRemSeq(util::ContigBlock& contig, size_t len, bool isCurContigFw, bool appendSuffix) {
  std::string seq = "";
  if (len == 0)
    return "";
  if (isCurContigFw && appendSuffix) {// append suffix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[getRemSeq] 1 from " << contig.contigLen_-len << " to " << contig.contigLen_ << " total length " << contig.contigLen_ << "\n";
    seq += contig.substrSeq(contig.contigLen_-len, len);
  }
  else if (isCurContigFw && !appendSuffix) {// append prefix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[getRemSeq] 2 from " << 0 << " to " << len << " total length " << contig.contigLen_ << "\n";
    seq = contig.substrSeq(0, len) + seq;
  }
  else if (!isCurContigFw && appendSuffix) {// append rc of the seq from the other end as a suffix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[getRemSeq] 3 rc from " << 0 << " to " << len << " total length " << contig.contigLen_ << "\n";
    seq += rc(contig.substrSeq(0, len));
  }
  else if (!isCurContigFw && !appendSuffix) {// append rc of the seq as prefix
    if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[getRemSeq] 4 rc from " << contig.contigLen_-len << " to " << contig.contigLen_ << " total length " << contig.contigLen_ << "\n";
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
  //if(verbose)std::cerr << std::this_thread::get_id() << " "  << "\t[rc] of " << str << " is ";
    for (uint32_t i = 0; i < str.length()/2; i++) {
      char tmp = str[i];
      str[i] = rev(str[str.length()-1-i]);
      str[str.length()-1-i] = rev(tmp);
    }
    //if (verbose) std::cerr << std::this_thread::get_id() << " "  << str << "\n";
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
                                                                                       size_t tpos,
                                                                                       size_t txpDist) {

  //if(verbose) std::cerr << std::this_thread::get_id() << " "  << "\t[Fetch successors] \n" ;
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
          nextCompatibleStruct theBest {cb.contigIdx_, tpos, 0, true};
          bool isBestValid = false;
          for(auto& posIt: nextHit.refRange){
            size_t succLastBaseTpos = posIt.pos() + cb.contigLen_ - 1;
            int overlap = tpos - posIt.pos() + 1;
            if (posIt.transcript_id() == tid and tpos < succLastBaseTpos and overlap >= 0) {
              //heuristic : TODO make sure it's right
              // As a successor, if the txp distance is not less than txpDist, then they are not successors in this txp
              // that 1 base in tpos + txpDist + 1 was kept as an overlap with end contig to have sanity check
              if (theBest.tpos < succLastBaseTpos) {
                isBestValid = true;
                theBest.tpos = succLastBaseTpos; // NOTE: careful
                theBest.isCurContigFw = posIt.orientation();
                theBest.cpos=theBest.isCurContigFw?(overlap-1):(cb.contigLen_-overlap);
                if (txpDist < succLastBaseTpos - tpos)
                  break;
              }
            }
          }
          if (isBestValid)
            successors.push_back(theBest) ;
        }
      }
    }

    if (verbose) std::cerr << " RETURNING SUCCESSOR!!\n";
    return successors;
}

template class RefSeqConstructor<PufferfishIndex>;
template class RefSeqConstructor<PufferfishSparseIndex>;
