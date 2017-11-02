//I am not sure what are the things I would need soon
//let's go with what we have
#include <iostream>
#include <mutex>
#include <vector>
#include <random>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <thread>
#include <tuple>
#include <sstream>
#include <fstream>
#include <iostream>
#include <tuple>
#include <memory>
#include <cstring>
#include <queue>

//we already have timers

//cereal include
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <sparsepp/spp.h>


#include "spdlog/spdlog.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/sinks/stdout_sinks.h"
#include "spdlog/sinks/ansicolor_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"


//jellyfish 2 include
#include "FastxParser.hpp"
#include "jellyfish/mer_dna.hpp"

//index header
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "SpinLock.hpp"
#include "MemCollector.hpp"
#include "SAMWriter.hpp"
#include "RefSeqConstructor.hpp"

#define START_CONTIG_ID ((uint32_t)-1) 
#define END_CONTIG_ID ((uint32_t)-2)

#define THRESHOLD 10

using paired_parser = fastx_parser::FastxParser<fastx_parser::ReadPair>;
using AlignmentOpts = util::AlignmentOpts ;
using HitCounters = util::HitCounters ;
using QuasiAlignment = util::QuasiAlignment ;
using MateStatus = util::MateStatus ;

using SpinLockT = std::mutex ;


void joinReadsAndFilter(spp::sparse_hash_map<size_t,std::vector<util::MemCluster>>& leftMemClusters,
                        spp::sparse_hash_map<size_t, std::vector<util::MemCluster>>& rightMemClusters,
                        std::vector<util::JointMems>& jointMemsList,
                        uint32_t maxFragmentLength,
                        uint32_t readLen,
                        bool verbose=false) {
  //orphan reads should be taken care of maybe with a flag!
  uint32_t maxCoverage{0};
  //std::cout << "txp count:" << leftMemClusters.size() << "\n";
  for (auto& leftClustItr : leftMemClusters) {
    // reference id
    size_t tid = leftClustItr.first;
    // left mem clusters
    auto& lClusts = leftClustItr.second;
    // right mem clusters for the same reference id
    auto& rClusts = rightMemClusters[tid];
    //if ((lClusts.size() > 5 || rClusts.size() > 5) && (lClusts.size()>0 && rClusts.size()>0))
    //std::cout << "\t" << tid << ": lClusts.size:" << lClusts.size() << " , rClusts.size:" << rClusts.size() << "\n";
    // Compare the left clusters to the right clusters to filter by positional constraints
    for (auto lclust =  lClusts.begin(); lclust != lClusts.end(); lclust++) {
      for (auto rclust =  rClusts.begin(); rclust != rClusts.end(); rclust++) {
        // if both the left and right clusters are oriented in the same direction, skip this pair
        // NOTE: This should be optional as some libraries could allow this.
        if (lclust->isFw == rclust->isFw) {
          continue;
        }

        // FILTER 1
        // filter read pairs based on the fragment length which is approximated by the distance between the left most start and right most hit end
        size_t fragmentLen = rclust->lastRefPos() + rclust->lastMemLen() - lclust->firstRefPos();
        if (lclust->firstRefPos() > rclust->firstRefPos()) {
          fragmentLen = lclust->lastRefPos() + lclust->lastMemLen() - rclust->firstRefPos();
        }

        if ( fragmentLen < maxFragmentLength) {
          // This will add a new potential mapping. Coverage of a mapping for read pairs is left->coverage + right->coverage
          // If we found a perfect coverage, we would only add those mappings that have the same perfect coverage
          if (maxCoverage < 2 * readLen || (lclust->coverage + rclust->coverage) == maxCoverage) {
            jointMemsList.emplace_back(tid, lclust, rclust, fragmentLen);
            if (verbose) {
              std::cout <<"\ntid:"<<tid<<"\n";
              std::cout <<"left:" << lclust->isFw << " size:" << lclust->mems.size() << " cov:" << lclust->coverage << "\n";
              for (size_t i = 0; i < lclust->mems.size(); i++){
                std::cout << "--- t" << lclust->mems[i].tpos << " r" << lclust->mems[i].memInfo->rpos << " cid:" << lclust->mems[i].memInfo->cid << " cpos: " << lclust->mems[i].memInfo->cpos << " len:" << lclust->mems[i].memInfo->memlen;
              }
              std::cout << "\nright:" << rclust->isFw << " size:" << rclust->mems.size() << " cov:" << rclust->coverage << "\n";
              for (size_t i = 0; i < rclust->mems.size(); i++){
                std::cout << "--- t" << rclust->mems[i].tpos << " r" << rclust->mems[i].memInfo->rpos << " cid:" << rclust->mems[i].memInfo->cid << " cpos: " << rclust->mems[i].memInfo->cpos << " len:" << rclust->mems[i].memInfo->memlen;
              }
            }
            uint32_t currCoverage =  jointMemsList.back().coverage();
            if (maxCoverage < currCoverage) {
              maxCoverage = currCoverage;
            }
          }
        } else {
          break;
        }
      }
    }
  }
  if (verbose) {
    std::cout << "\nBefore filter " << jointMemsList.size() << " maxCov:" << maxCoverage << "\n";
    std::cerr << "\n" << jointMemsList.size() << " maxCov:" << maxCoverage << "\n";
  }
  // FILTER 2
  // filter read pairs that don't have enough base coverage (i.e. their coverage is less than half of the maximum coverage for this read)
  double coverageRatio = 0.5;
  // if we've found a perfect match, we will erase any match that is not perfect
  if (maxCoverage == 2*readLen) {
    jointMemsList.erase(std::remove_if(jointMemsList.begin(), jointMemsList.end(),
                                        [&maxCoverage](util::JointMems& pairedReadMems) -> bool {
                                          return pairedReadMems.coverage() != maxCoverage;
                                        }),
                         jointMemsList.end());
  }
  else {// ow, we will go with the heuristic that just keep those mappings that have at least half of the maximum coverage
    //std::cout << "no!\n";
    jointMemsList.erase(std::remove_if(jointMemsList.begin(), jointMemsList.end(),
                                     [&maxCoverage, coverageRatio](util::JointMems& pairedReadMems) -> bool {
                                       return pairedReadMems.coverage() < coverageRatio*maxCoverage ;
                                     }),
                      jointMemsList.end());
  }
  if (verbose) {
    std::cout << "\nFinal stat: " << jointMemsList.size() << "\n";
    std::cerr << jointMemsList.size() << "\n";
  }
}


/*
template<typename PufferfishIndexT>
void updateCacheBlock(util::ContigCecheBlock& cacheBlock, uint32_t startClip, uint32_t endClip, PufferfishIndexT& pfi, uint32_t cid){

  uint32_t globalPos = std::numeric_limits<uint64_t>max() ;

  if(startClip >= cacheBlock.cpos and endClip < cacheBlock.cpos + cacheBlock.blockLen){
    tmp = cacheBlock.cseq.substr(startClip-cacheBlock.cpos, THRESHOLD + offset) ;
  }
  else{
    uint32_t newCpos, newLen ;
    //if doisjoint
    globalPos = pfi.getGlobalPos(firstMem.memInfo->cid) ;
    if(startClip < cacheBlock.cpos and endClip < cacheBlock.cpos){
      //prepend
      int gap = cacheBlock.cpos - endClip ;
      if(gap < THRESHOLD){
        tmp = pfi.getSeqStr(globalPos+startClip, cacheBlock.cpos - startClip) ;
        cacheBlock.blockLen += (cacheBlock.cpos - startClip + 1);
        cacheBlock.cseq = tmp + cacheBlock.cseq ;
      }else{
        //get new
        tmp = pfi.getSeqStr(globalPos+startClip, endClip - startClip +1) ;
        cacheBlock.cpos = startClip ;
        cacheBlock.cseq = tmp ;
        cacheBlock.blockLen = endClip - startClip + 1 ;

      }
    }else if((startClip >= cacheBlock.cpos + cacheBlock.blockLen) and (endClip > cacheBlock.cpos + cacheBlock.blockLen)){
      //append 
      int gap = startClip - cacheBlock.cpos ;
      if(gap < THRESHOLD){
        tmp = pfi.getSeqStr(globalPos+cacheBlock.cpos+cacheBlock.blockLen , cacheBlock.cpos - startClip) ;
        cacheBlock.cseq = tmp + cacheBlock.cseq ;
      }else{
        //get new
        tmp = pfi.getSeqStr(globalPos+startClip, endClip - startClip +1) ;
        cacheBlock.cpos = startClip ;
        cacheBlock.cseq = tmp ;
        cacheBlock.blockLen = endClip - startClip + 1 ;

      }

    }
  }

}*/


//this is super heuristic
//not a complete BFS
template <typename PufferfishIndexT>
void populatePaths(util::MemInfo& smem,
                   util::MemInfo& emem,
                   std::string& path,
                   PufferfishIndexT& pfi,
                   uint32_t tid,
                   std::map<uint32_t, util::ContigCecheBlock>& contigSeqCache,
                   size_t readGapDist){

  //In graph searching
  //std::cerr << "Graph searching: go time\n" ;

  CanonicalKmer::k(pfi.k()) ;

  auto s = *smem.memInfo ;
  auto e = *emem.memInfo ;

  auto stpos = smem.tpos ;

  auto sCid = s.cid ;
  auto eCid = e.cid ;

  auto& edges = pfi.getEdge() ;

  std::map<uint32_t,bool> visited ;
  uint32_t numOfNodesVisited{0} ;


  //Initialization
  std::queue<uint32_t> queue ;
  queue.push(sCid) ;
  path = "" ;

  //decide direction
  //util::Direction toGo;

  int remainingLen = readGapDist ;
  bool start{true} ;

  //at worst we would visit readGapDist number of nodes 
  while(numOfNodesVisited < readGapDist and !queue.empty()){
    auto cid = queue.front() ;
    numOfNodesVisited += 1 ;
    visited[cid] = true ;
    queue.pop() ;
    auto contigLen = pfi.getContigLen(cid) ;

    //traverse this node
    int toClip{0} ;
    if(!start){
      toClip = (contigLen < remainingLen)?(contigLen):(remainingLen) ;
      //std::cerr << "middle clipping " << toClip << " from cid "<<cid << "\n" ;
      path += pfi.getSeqStr(pfi.getGlobalPos(cid),toClip) ;
    }
    else{
      toClip = std::min(static_cast<int>(contigLen - (s.cpos+s.memlen)),remainingLen) ;
      if(toClip>0){
        //std::cerr << "start clipping " << toClip << " from cid "<<cid << " glob pos start "<< pfi.getGlobalPos(cid)<< " cpos: " << s.cpos << " memlen: " << s.memlen << "\n" ;
        path += pfi.getSeqStr(pfi.getGlobalPos(cid)+s.cpos+s.memlen,toClip) ;
      }
      start = false ;
    }
    //read finished 
    remainingLen -= static_cast<int> (toClip) ;
    if(!remainingLen)
      return ;

    //explore edges, only one edge will be chosen, the extension is greedy
    uint8_t edgeVec = edges[cid] ;
    std::vector<util::extension> ext = util::getExts(edgeVec) ;

    if(!ext.empty()){
      //go over the ext
      CanonicalKmer kb = pfi.getStartKmer(cid) ;
      CanonicalKmer ke = pfi.getEndKmer(cid) ;
      for(auto& ed : ext){
        auto kbtmp = kb ;
        auto ketmp = ke ;
        char c = ed.c ;
        CanonicalKmer kt ;

        if(ed.dir == util::Direction::FORWARD){
          ke.shiftFw(c) ;
          kt.fromNum(ke.getCanonicalWord()) ;
        }else{
          kb.shiftBw(c) ;
          kt.fromNum(kb.getCanonicalWord()) ;
        }


        auto nextHit = pfi.getRefPos(kt) ;
        uint32_t nextCid{0} ;
        if(!nextHit.empty()){
           nextCid = nextHit.contigID() ;
        }else{
          std::cerr << "should not happen\n" ;
        }

        //1. transcript consistensy
        //2. position consistensy

        if(visited.find(nextCid) == visited.end()){
          for(auto& posIt : nextHit.refRange){
            //check for transcript consistensy
            //auto refPosOri = nextHit.decodeHit(posIt) ;
            if(posIt.transcript_id() == tid){
              //check for position consistensy
              if((posIt.pos() > stpos and posIt.orientation()) or (posIt.pos() < stpos and !posIt.orientation())){
                if(eCid != nextCid){
                  queue.push(nextCid) ;
                  stpos = posIt.pos() ;
                  goto StopEdge ;
                }
                //found an edge 
                break ;

              }
            }
          }
        
        }
        kb = kbtmp ;
        ke = ketmp ;
      }// end of edge list

    }else{
      //the previous node is the end node
      return ;
    }
   
   StopEdge: ;

  }

}

//NOTE: not doing mem chains 
//I did not take the reference since I
//want to sort them according to finishing time
//without disrupting the order
/*
std::vector<int> dynamicMemChain(std::vector<util::MemInfo> mems){
  std::vector<int> finalChain ;

  return finalChain ;
}



std::vector<int> greedyMemChain(std::vector<util::MemInfo>& mems, bool fw){
  //start with the starting chain and go forward or
  //backward to get any overlappping chain
  //can be shitty
  std::vector<int> finalChain ;
  finalChain.push_back(0) ;
  //little heuristic
  if(mems.size() == 1){
    return finalChain ;
  }

  auto prevBlock = *mems[0].memInfo ;
  for(size_t i = 1; i < mems.size(); ++i){
    auto thisBlock = *mems[i].memInfo ;
    auto prevEnd = (fw)?(prevBlock.rpos + prevBlock.memlen):(prevBlock.rpos) ;
    auto thisBegin = (fw)?(thisBlock.rpos):(thisBlock.rpos + thisBlock.memlen) ;
    bool overlap = (fw)?(thisBegin < prevEnd):(thisBegin > prevEnd) ;
    if(!overlap){
      prevBlock = *mems[i].memInfo ;
      finalChain.push_back(i) ; 
    }
  }
  return finalChain ;
  }*/




/*
template<typename PufferfishIndexT>

void goOverClust2(PufferfishIndexT& pfi,
                 std::vector<util::MemCluster>::iterator clust,
                 fastx_parser::ReadSeq& read,
                 std::map<uint32_t,util::ContigCecheBlock>& contigSeqCache,
                 uint32_t tid,
                 bool  verbose){

  std::string readSeq = read.seq ;
  std::string readName = read.name ;

  size_t clustSize = clust->mems.size() ;
  bool rdFwd = clust->isFw ;

  for(int it = 0 ; it < clustSize ; ++it){
    auto& thismem = clust->mems[it] ;
    auto& nextmem = (it < clustSize-1)?clust->mems[it+1]:clust->mems[it] ;
    //check if there is an insert
    //no read gap in corresponding tpos
    auto readGapDist = rdFwd?(nextmem.memInfo->rpos - (thismem.memInfo->rpos + thismem.memInfo->memlen)):(thismem.memInfo->rpos - (nextmem.memInfo->rpos+nextmem.memInfo))

  }



  }*/


template <typename PufferfishIndexT>
void goOverClust(PufferfishIndexT& pfi,
                 std::vector<util::MemCluster>::iterator clust,
                 fastx_parser::ReadSeq& read,
                 std::map<uint32_t,util::ContigCecheBlock>& contigSeqCache,
                 uint32_t tid,
                 bool  verbose){

  //size_t readLen = readSeq.length() ;

  std::string readSeq = read.seq ;
  std::string readName = read.name ;
  auto readLen = readSeq.length() ;

  size_t overhangLeft ;
  size_t overhangRight ;
  size_t clustSize = clust->mems.size() ;

  size_t startIndex = 0;
  size_t endIndex = clustSize -1 ;

  if(!clust->isFw) std::swap(startIndex,endIndex) ;
  clust->isVisited = true ;

  if(clustSize > 1){
    bool overlap = (clust->isFw)?((clust->mems[0].memInfo->rpos + clust->mems[0].memInfo->memlen) > clust->mems[clustSize-1].memInfo->rpos):((clust->mems[clustSize-1].memInfo->rpos + clust->mems[clustSize-1].memInfo->memlen) > clust->mems[0].memInfo->rpos);
    if(overlap)
      return ;
  }


  //if(!clust->isFw) std::reverse(clust->mems.begin(),clust->mems.end()) ;
  //NOTE: Taking care of the part between beginning
  //of the read block and 0
  /*for(size_t it = 0 ; it < clustSize ; ++it){
    util::MemInfo *smem{nullptr} ;
    util::MemInfo *emem{nullptr} ;
    int readGapDist ;
    int contigGapDist ;
    util::Direction togo ;

    

    if(it == 0){
      emem = &clust->mems[it] ;
      togo = (emem->memInfo->cIsFw)?util::Direction::FORWARD:util::Direction::BACKWORD ;
      readGapDist = emem->memInfo->rpos ;
    } else if(it > 0 and it < clustSize-1){
      smem = &clust->mems[it] ;
      emem = &clust->mems[it+1] ;
      togo = (smem->memInfo->cIsFw)?util::Direction::FORWARD:util::Direction::BACKWORD ;
    } else{
      smem = &clust->mems[it] ;
      togo = (smem->memInfo->cIsFw)?util::Direction::FORWARD:util::Direction::BACKWORD ;
    }
    //left part of read
    if(!smem){
      readGapDist = emem->memInfo->rpos ;
    }

    if(it == 0){
      smem = &clust->mems[it] ;
      auto readGapDist = smem->memInfo->rpos ;
      uint32_t remainingContigLen{0} ;

      util::Direction togo = (smem->memInfo->cIsFw)?util::Direction::FORWARD:util::Direction::BACKWORD ;
      if(togo == util::Direction::FORWARD){
        remainingContigLen = pfi.getContigLen(smem->memInfo->cid) - (pfi.getContigLen(smem->memInfo->cid) + smem->memInfo->memlen) ;
      }else{
        remainingContigLen = smem->memInfo->cpos ;
      }
      if(remainingContigLen < readGapDist){
        //clip contig 
        readGapDist -= remainingContigLen ;
        //call graph 
      }
    }
    if(it > 0 and it < clustSize-1){
      smem = &clust->mems[it] ;
      emem = &clust->mems[it+1] ;
      util::Direction togo = (smem->memInfo->cIsFw)?util::Direction::FORWARD:util::Direction::BACKWORD ;
    }

  }*/


  overhangLeft = clust->mems[startIndex].memInfo->rpos ;

  if(overhangLeft > 0){
    auto sInfo = *clust->mems[startIndex].memInfo ;
    bool torc = clust->mems[startIndex].memInfo->cIsFw ;
    std::string tmp(readSeq.substr(0,overhangLeft));
    std::string rdtmp = (torc) ? tmp : util::reverseComplement(tmp)  ;

    auto cstart = (static_cast<int>(sInfo.cpos - overhangLeft) > 0)?(sInfo.cpos - overhangLeft):0 ;
    size_t toClip = std::min(static_cast<uint32_t>(overhangLeft), pfi.getContigLen(sInfo.cid)) ;
    //std::cerr << "left clipping "<<toClip<<" from "<<cstart<<"\n" ;
    std::string contigtmp = pfi.getSeqStr(pfi.getGlobalPos(sInfo.cid) + cstart, toClip, sInfo.cIsFw) ;
    clust->alignableStrings.push_back({rdtmp, contigtmp}) ;
  }

  overhangRight = static_cast<int>(readLen - (clust->mems[endIndex].memInfo->rpos + clust->mems[endIndex].memInfo->memlen)) ;

  if(overhangRight > 0){
    auto eInfo = *clust->mems[endIndex].memInfo ;
    bool torc = eInfo.cIsFw ;

    std::string tmp(readSeq.substr(overhangRight));
    std::string rdtmp = (torc) ? tmp : util::reverseComplement(tmp)  ;

    auto cstart = eInfo.cpos + eInfo.memlen ;
    size_t toClip = std::min(pfi.getContigLen(eInfo.cid) - cstart, static_cast<uint32_t>(overhangRight)) ;

    //std::cerr << "right clipping "<<toClip<<" from "<<cstart<<"\n" ;
    std::string contigtmp = pfi.getSeqStr(pfi.getGlobalPos(eInfo.cid) + cstart, toClip, eInfo.cIsFw) ;
    clust->alignableStrings.push_back({rdtmp, contigtmp}) ;
  }

  if(verbose){
    std::cerr << readName << "\n" ;
    for(auto& m: clust->mems){
      std::cerr << "cpos "<<m.memInfo->cpos << "\n" ;
    }
  }

  //if start and end has overlap then return NOTE: heuristic
  if(clustSize > 1){
    bool overlap = (clust->isFw)?((clust->mems[0].memInfo->rpos + clust->mems[0].memInfo->memlen) > clust->mems[clustSize-1].memInfo->rpos):((clust->mems[clustSize-1].memInfo->rpos + clust->mems[clustSize-1].memInfo->memlen) > clust->mems[0].memInfo->rpos);
    if(overlap)
      return ;
  }

  for(size_t it=0 ; it < clustSize -1 ; ++it){
    //search if contig sequence is in the cacheBlock and to clipon the same contig explore paths between unimems from graphs 
    auto startMem = clust->mems[it] ;
    auto endMem = clust->mems[it+1] ;
    auto sInfo = *startMem.memInfo ;
    auto eInfo = *endMem.memInfo ;

    if(!clust->isFw) std::swap(sInfo,eInfo) ;

    std::string readgap = "" ;
    std::string path = "";
    
    auto readGapDist = static_cast<int>(eInfo.rpos - (sInfo.rpos + sInfo.memlen)) ;

    if(readGapDist >= 0){
      //clip the appropriate read sequence
      readgap = readSeq.substr(sInfo.rpos + sInfo.memlen, readGapDist) ;
      if(startMem.memInfo->cid  != endMem.memInfo->cid){
        //going into graph search
        //std::cerr << "read name  "<< readName<<"\n" ;
        //std::cerr << "start cpos "<<startMem.memInfo->cpos<<"\n" ;
        populatePaths(startMem, endMem, path, pfi, tid, contigSeqCache, readGapDist) ;
      } else if(readGapDist > 0){

        //it's on the same contig just insert the sequence from the contig
        //std::cerr << "c"<<sInfo.cid <<"\t cpos:" << sInfo.cpos <<"\t rpos:" << sInfo.rpos << "\n" ;
        auto cstart = sInfo.cpos + sInfo.memlen;
        auto toClip = readGapDist ;
        //std::cerr << "middle clipping "<<toClip<<" from "<<cstart<<"\n" ;
        path = pfi.getSeqStr(pfi.getGlobalPos(sInfo.cid)+cstart, toClip) ;
      }
    }

    clust->alignableStrings.push_back({readgap, path}) ;

  }

}

//template <typename ReadPairT ,typename PufferfishIndexT>
std::string extractReadSeq(const std::string readSeq, uint32_t rstart, uint32_t rend, bool isFw) {
  //DEBUG: tracking substr error
  //std::cerr << "checking boundaries " << readSeq.length() << " " << rstart << " " << rend << "\n" ;
  std::string subseq = readSeq.substr(rstart, rend-rstart);
  if (isFw)
    return subseq;
  return util::reverseComplement(subseq); //reverse-complement the substring
}

template <typename PufferfishIndexT>
void createSeqPairs(PufferfishIndexT* pfi,
                    std::vector<util::MemCluster>::iterator clust,
                    const std::string& readSeq,
                    RefSeqConstructor<PufferfishIndexT>& refSeqConstructor,
                    uint32_t tid,
                    bool verbose) {

  (void)verbose;
  for(size_t it=0 ; it < clust->mems.size() -1 ; ++it) {
    // while mems overlap, continue
    if (clust->mems[it+1].tpos < clust->mems[it].tpos + clust->mems[it].memInfo->memlen){
      //std::cerr<<"Overlapping mems skipping\n" ;
      continue;
    }

    uint32_t rstart = (clust->mems[it].memInfo->rpos + clust->mems[it].memInfo->memlen);
    uint32_t rend = clust->mems[it+1].memInfo->rpos;
    if (!clust->isFw) {
      rstart = (clust->mems[it+1].memInfo->rpos + clust->mems[it+1].memInfo->memlen);
      rend = clust->mems[it].memInfo->rpos;
    }
    //Insertion in read
    if (clust->mems[it+1].tpos == clust->mems[it].tpos + clust->mems[it].memInfo->memlen) {
        if (rend-rstart > 0) {
          std::string tmp = extractReadSeq(readSeq, rstart, rend, clust->isFw) ;
          clust->alignableStrings.push_back(std::make_pair(extractReadSeq(readSeq, rstart, rend, clust->isFw), ""));
        }
    }
    else {// we have a gap on transcript and need to do graph traversal
      //std::cerr << "Need to do graph traversal \n" ; 
      std::string refSeq;
      bool firstContigDirWRTref = clust->mems[it].memInfo->cIsFw == clust->isFw; // I am changing the logic 
      bool secondContigDirWRTref = clust->mems[it+1].memInfo->cIsFw == clust->isFw;
      util::ContigBlock scb = {
                               clust->mems[it].memInfo->cid,
                               clust->mems[it].memInfo->cGlobalPos,
                               clust->mems[it].memInfo->clen,
                               pfi->getSeqStr(clust->mems[it].memInfo->cGlobalPos, clust->mems[it].memInfo->clen)};
      util::ContigBlock ecb = {
                               clust->mems[it+1].memInfo->cid,
                               clust->mems[it+1].memInfo->cGlobalPos,
                               clust->mems[it+1].memInfo->clen,
                               pfi->getSeqStr(clust->mems[it+1].memInfo->cGlobalPos, clust->mems[it+1].memInfo->clen)};
      //TODO -1 needed or not? This is the problem
      //We should check if cpos does not beome negative, since we are doing -1
      
      uint32_t cstart = firstContigDirWRTref?(clust->mems[it].memInfo->cpos + clust->mems[it].memInfo->memlen):(clust->mems[it].memInfo->cpos > 0?clust->mems[it].memInfo->cpos-1:0);
      //TODO also we have to check if the next does
      //not start from 0
      uint32_t cend = secondContigDirWRTref?(clust->mems[it+1].memInfo->cpos > 0 ? clust->mems[it+1].memInfo->cpos-1:0):(clust->mems[it+1].memInfo->cpos + clust->mems[it+1].memInfo->memlen);

      //NOTE: if they the mem cstart and cend are back to back and on same contig
      //do we need to do graph ? probably not
      if(scb.contigIdx_ != ecb.contigIdx_ or cstart != cend)
        refSeqConstructor.doBFS(tid, clust->mems[it].tpos, firstContigDirWRTref, scb, cstart, ecb, cend, rend-rstart+THRESHOLD, refSeq);
      if (rend-rstart>0) {
        clust->alignableStrings.push_back(std::make_pair(extractReadSeq(readSeq, rstart, rend, clust->isFw), refSeq));
      }
    }
  }

  //TODO take care of
  //gap before first unimem
  //gap after last unimem

  clust->isVisited = true;
}

template <typename ReadPairT ,typename PufferfishIndexT>
void traverseGraph(ReadPairT& rpair,
                   util::JointMems& hit,
                   PufferfishIndexT& pfi,
                   RefSeqConstructor<PufferfishIndexT>& refSeqConstructor,
                   bool verbose){

  size_t tid = hit.tid ;
  auto readLen = rpair.first.seq.length() ;


  if(!hit.leftClust->isVisited && hit.leftClust->coverage < readLen)
    createSeqPairs(&pfi, hit.leftClust, rpair.first.seq, refSeqConstructor, tid, verbose);
    //goOverClust(pfi, hit.leftClust, rpair.first, contigSeqCache, tid, verbose) ;
  if(!hit.rightClust->isVisited && hit.rightClust->coverage < readLen)
    createSeqPairs(&pfi, hit.rightClust, rpair.second.seq, refSeqConstructor, tid, verbose);
    //goOverClust(pfi, hit.rightClust, rpair.second, contigSeqCache, tid, verbose) ;
}

template <typename PufferfishIndexT>
void processReadsPair(paired_parser* parser,
                     PufferfishIndexT& pfi,
                     SpinLockT* iomutex,
                     std::shared_ptr<spdlog::logger> outQueue,
                     HitCounters& hctr,
                     AlignmentOpts* mopts){
  MemCollector<PufferfishIndexT> memCollector(&pfi) ;


  spp::sparse_hash_map<uint32_t, util::ContigBlock> contigSeqCache ;
  RefSeqConstructor<PufferfishIndexT> refSeqConstructor(&pfi, contigSeqCache);

  //std::cerr << "\n In process reads pair\n" ;
    //TODO create a memory layout to store
    //strings then will allocate alignment to them
    //accordingly
  std::vector<std::string> refBlocks ;

  auto logger = spdlog::get("stderrLog") ;
  fmt::MemoryWriter sstream ;
  //size_t batchSize{2500} ;
  size_t readLen{0} ;

  spp::sparse_hash_map<size_t, std::vector<util::MemCluster>> leftHits ;
  spp::sparse_hash_map<size_t, std::vector<util::MemCluster>> rightHits ;
  std::vector<util::JointMems> jointHits ;
  PairedAlignmentFormatter<PufferfishIndexT*> formatter(&pfi);

  
  auto rg = parser->getReadGroup() ;
  while(parser->refill(rg)){
    for(auto& rpair : rg){
      readLen = rpair.first.seq.length() ;
      //std::cout << readLen << "\n";
      //std::cout << rpair.first.name << "\n";
      bool verbose = false;//rpair.first.name == "read25634671/ENST00000396859;mate1:392-491;mate2:535-633";

      ++hctr.numReads ;

      jointHits.clear() ;
      leftHits.clear() ;
      rightHits.clear() ;
      memCollector.clear();

      //help me to debug, will deprecate later
      //std::cerr << "\n first seq in pair " << rpair.first.seq << "\n" ;
      //std::cerr << "\n second seq in pair " << rpair.second.seq << "\n" ;

      //std::cerr << "\n going inside hit collector \n" ;
      //readLen = rpair.first.seq.length() ;
      bool lh = memCollector(rpair.first.seq,
                             leftHits,
                             mopts->maxSpliceGap,
                             MateStatus::PAIRED_END_LEFT,
                             verbose
                             /*
                             mopts->consistentHits,
                             refBlocks*/) ;
      bool rh = memCollector(rpair.second.seq,
                             rightHits,
                             mopts->maxSpliceGap,
                             MateStatus::PAIRED_END_RIGHT,
                             verbose
                             /*,
                             mopts->consistentHits,
                             refBlocks*/) ;
      //do intersection on the basis of
      //performance, or going towards selective alignment
      //otherwise orphan
      if(verbose){
        for(auto& l : leftHits){
          auto& lclust = l.second ;
          for(auto& clust : lclust)
            for(auto& m : clust.mems){
              std::cerr << "before join "<<m.memInfo->cid << " cpos "<< m.memInfo->cpos<< "\n" ;
            }
        }
        for(auto& l : rightHits){
          auto& lclust = l.second ;
          for(auto& clust : lclust)
            for(auto& m : clust.mems){
              std::cerr << "before join "<<m.memInfo->cid << " cpos "<< m.memInfo->cpos<< "\n" ;
            }
        }
      }

      if(lh && rh){
        joinReadsAndFilter(leftHits, rightHits, jointHits, mopts->maxFragmentLength, readLen, verbose) ;
      }
      else{
        //ignore orphans for now
      }
      //jointHits is a vector
      //this can be used for BFS

      //void traverseGraph(std::string& leftReadSeq, std::string& rightReadSeq, util::JointMems& hit, PufferfishIndexT& pfi,   std::map<uint32_t, std::string>& contigSeqCache){
      for(auto& h: jointHits){
        for(auto& m : h.leftClust->mems){
          auto cSeq = pfi.getSeqStr(pfi.getGlobalPos(m.memInfo->cid)+m.memInfo->cpos, m.memInfo->memlen) ;
          auto tmp = rpair.first.seq.substr(m.memInfo->rpos, m.memInfo->memlen) ;
          auto rseq = m.memInfo->cIsFw?tmp:util::reverseComplement(tmp);
          if(cSeq != rseq){
            std::cerr <<"rpos:"<<m.memInfo->rpos<<"\t"<<"c "<<static_cast<int>(m.memInfo->cid)<<"\t"<<"cpos:"<<m.memInfo->cpos<<"\t"<<m.memInfo->memlen<<" cfw: "<<int(m.memInfo->cIsFw)<<" isFw "<<h.leftClust->isFw<< " tpos: "<<m.tpos<<"\n" ;
            std::cerr<<rpair.first.name << "\n" ;
            std::cerr << cSeq << "\n" << rseq <<"\n" << tmp << "\n";
            std::exit(1) ;
          }
        }
        for(auto& m : h.rightClust->mems){
          auto cSeq = pfi.getSeqStr(pfi.getGlobalPos(m.memInfo->cid)+m.memInfo->cpos, m.memInfo->memlen) ;
          auto tmp = rpair.second.seq.substr(m.memInfo->rpos, m.memInfo->memlen) ;
          auto rseq = m.memInfo->cIsFw?tmp:util::reverseComplement(tmp);
          if(cSeq != rseq){
            std::cerr <<"rpos:"<<m.memInfo->rpos<<"\t"<<"c "<<static_cast<int>(m.memInfo->cid)<<"\t"<<"cpos:"<<m.memInfo->cpos<<"\t"<<m.memInfo->memlen<<" cfw: "<<int(m.memInfo->cIsFw)<<" isFw "<<h.rightClust->isFw<<" tpos: "<<m.tpos<<"\n" ;
            std::cerr<<rpair.second.name << "\n" ;
            std::cerr << cSeq << "\n" << rseq <<"\n" << tmp << "\n" ;
            std::exit(1) ;
          }
        }

      }

      bool doTraverse = true;
      if (doTraverse) {
        //TODO Have to make it per thread 
        //have to make write access thread safe
        //std::cout << "traversing graph \n" ;

        int hitNum{0} ;

        if(jointHits.empty() or jointHits.front().coverage() < 2*readLen){
          for(auto& hit : jointHits){
            traverseGraph(rpair, hit, pfi, refSeqConstructor, verbose) ;
            hitNum++ ;

          }
        }
      }



      //extractSuitableAligningPairs(joinHits);
      //TODO Write them to a sam file
      hctr.totHits += jointHits.size();
      hctr.peHits += jointHits.size();
      hctr.numMapped += !jointHits.empty() ? 1 : 0;

      //std::cerr << "\n Number of total joint hits" << jointHits.size() << "\n" ;
      //TODO When you get to this, you should be done aligning!!
      //fill the QuasiAlignment list
      std::vector<QuasiAlignment> jointAlignments;
      for (auto& jointHit : jointHits) {
        jointAlignments.emplace_back(jointHit.tid,           // reference id
                                     jointHit.leftClust->getTrFirstHitPos(),     // reference pos
                                     jointHit.leftClust->isFw ,     // fwd direction
                                     readLen, // read length
                                     jointHit.fragmentLen,       // fragment length
                                     true);         // properly paired
         // Fill in the mate info		         // Fill in the mate info
        auto& qaln = jointAlignments.back();
         qaln.mateLen = readLen;
         qaln.matePos = jointHit.rightClust->getTrFirstHitPos();
         qaln.mateIsFwd = jointHit.rightClust->isFw;
         qaln.mateStatus = MateStatus::PAIRED_END_PAIRED;
      }

      if(jointHits.size() > 0 and !mopts->noOutput){
        writeAlignmentsToStream(rpair, formatter, jointAlignments, sstream, mopts->writeOrphans) ;
      }




      //TODO write them on cmd
      if (hctr.numReads > hctr.lastPrint + 1000000) {
        hctr.lastPrint.store(hctr.numReads.load());
        if (!mopts->quiet and iomutex->try_lock()) {
          if (hctr.numReads > 0) {
            std::cerr << "\r\r";
          }
          std::cerr << "saw " << hctr.numReads << " reads : "
                    << "pe / read = " << hctr.peHits / static_cast<float>(hctr.numReads)
                    << " : se / read = " << hctr.seHits / static_cast<float>(hctr.numReads) << ' ';
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
          std::cerr << ": true hit \% = "
                    << (100.0 * (hctr.trueHits / static_cast<float>(hctr.numReads)));
#endif // __DEBUG__
          iomutex->unlock();
        }
      }
    } // for all reads in this job



     //TODO Dump Output
    // DUMP OUTPUT
    if (!mopts->noOutput) {
      std::string outStr(sstream.str());
      //std::cerr << "\n OutStream size "<< outStr.size() << "\n" ;
      // Get rid of last newline
      if (!outStr.empty()) {
        outStr.pop_back();
        outQueue->info(std::move(outStr));
      }
      sstream.clear();
      /*
        iomutex->lock();
        outStream << sstream.str();
        iomutex->unlock();
        sstream.clear();
      */
    }

  } // processed all reads
}

template <typename PufferfishIndexT>
bool spawnProcessReadsthreads(
                              uint32_t nthread,
                              paired_parser* parser,
                              PufferfishIndexT& pfi,
                              SpinLockT& iomutex,
                              std::shared_ptr<spdlog::logger> outQueue,
                              HitCounters& hctr,
                              AlignmentOpts* mopts){

  std::vector<std::thread> threads ;

  for(size_t i = 0; i < nthread ; ++i){

    threads.emplace_back(processReadsPair<PufferfishIndexT>,
                         parser,
                         std::ref(pfi),
                         &iomutex,
                         outQueue,
                         std::ref(hctr),
                         mopts);
  }
  for(auto& t : threads) { t.join(); }

  return true ;
}


void printAlignmentSummary(HitCounters& hctrs, std::shared_ptr<spdlog::logger> consoleLog) {
  consoleLog->info("Done mapping reads.");
  consoleLog->info("\n\n");
  consoleLog->info("=====");
  consoleLog->info("Observed {} reads", hctrs.numReads);
  consoleLog->info("Mapping rate : {:03.2f}%", (100.0 * static_cast<float>(hctrs.numMapped)) / hctrs.numReads);
  consoleLog->info("Average # hits per read : {}", hctrs.totHits / static_cast<float>(hctrs.numReads));
  consoleLog->info("=====");
}

template <typename PufferfishIndexT>
bool alignReads(
              PufferfishIndexT& pfi,
              std::shared_ptr<spdlog::logger> consoleLog,
              AlignmentOpts* mopts) {



  std::streambuf* outBuf ;
  std::ofstream outFile ;
  //bool haveOutputFile{false} ;
  if(mopts->outname == ""){

    outBuf = std::cout.rdbuf() ;
  }else{
    outFile.open(mopts->outname) ;
    outBuf = outFile.rdbuf() ;
    //haveOutputFile = true ;
  }

  //out stream to the buffer
  //it can be std::cout or a file

  std::ostream outStream(outBuf) ;

  //this is my async queue
  //a power of 2
  size_t queueSize{268435456};
  spdlog::set_async_mode(queueSize);

  auto outputSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(outStream) ;
  std::shared_ptr<spdlog::logger> outLog = std::make_shared<spdlog::logger>("puffer::outLog",outputSink) ;
  outLog->set_pattern("%v");

  //mopts->noOutput = true;

  uint32_t nthread = mopts->numThreads ;
  std::unique_ptr<paired_parser> pairParserPtr{nullptr} ;

  //write the SAMHeader
  //If nothing gets printed by this
  //time we are in trouble
  writeSAMHeader(pfi, outLog) ;


  size_t chunkSize{10000} ;
  SpinLockT iomutex ;
  {
    ScopedTimer timer(!mopts->quiet) ;
    HitCounters hctrs ;
    consoleLog->info("mapping reads ... \n\n\n") ;
    std::vector<std::string> read1Vec = util::tokenize(mopts->read1, ',') ;
    std::vector<std::string> read2Vec = util::tokenize(mopts->read2, ',') ;

    if(read1Vec.size() != read2Vec.size()){
      consoleLog->error("The number of provided files for"
                        "-1 and -2 are not same!") ;
      std::exit(1) ;
    }

    uint32_t nprod = (read1Vec.size() > 1) ? 2 : 1;
    pairParserPtr.reset(new paired_parser(read1Vec, read2Vec, nthread, nprod, chunkSize));
    pairParserPtr->start();

    spawnProcessReadsthreads(nthread, pairParserPtr.get(), pfi, iomutex,
                             outLog, hctrs, mopts) ;

    pairParserPtr->stop();
	consoleLog->info("flushing output queue.");
  printAlignmentSummary(hctrs, consoleLog);
	outLog->flush();
  }


  return true ;


}

int pufferfishAligner(AlignmentOpts& alnargs){

  //auto outputSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(outStream) ;
  //std::shared_ptr<spdlog::logger> outLog = std::make_shared<spdlog::logger>("puffer::outLog",outputSink) ;

//auto rawConsoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
//auto rawConsoleSink = spdlog::sinks::stderr_sink_mt::instance() ;
 // auto consoleSink = std::make_shared<spdlog::sinks::ansicolor_sink>(rawConsoleSink);
    auto consoleLog = spdlog::stderr_color_mt("console");


    bool success{false} ;

  auto indexDir = alnargs.indexDir ;

  std::string indexType;
  {
    std::ifstream infoStream(indexDir + "/info.json");
    cereal::JSONInputArchive infoArchive(infoStream);
    infoArchive(cereal::make_nvp("sampling_type", indexType));
    std::cerr << "Index type = " << indexType << '\n';
    infoStream.close();
  }

  if(indexType == "dense"){
    PufferfishIndex pfi(indexDir) ;
    success = alignReads(pfi, consoleLog, &alnargs) ;

  }else if(indexType == "sparse"){
    PufferfishSparseIndex pfi(indexDir) ;
    success = alignReads(pfi, consoleLog, &alnargs) ;
  }

  if (!success) {
    consoleLog->warn("Problem mapping.");
  }
  return 0 ;
}
