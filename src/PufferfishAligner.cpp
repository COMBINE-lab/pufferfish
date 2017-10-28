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
#include <sparsepp/spp.h>

#define START_CONTIG_ID ((uint32_t)-1) 
#define END_CONTIG_ID ((uint32_t)-2)

#define THRESHOLD 10

using paired_parser = fastx_parser::FastxParser<fastx_parser::ReadPair>;
using AlignmentOpts = util::AlignmentOpts ;
using HitCounters = util::HitCounters ;
using QuasiAlignment = util::QuasiAlignment ;
using MateStatus = util::MateStatus ;

using SpinLockT = std::mutex ;


void joinReadsAndFilter(spp::sparse_hash_map<size_t,
                        std::vector<util::MemCluster>>& leftMemClusters,
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

        // We will use the left and right iterators to point to the clusters that are
        // "actually" the leftmost and rightmost (rather than those that come from read 1 and read 2).
        decltype(lClusts.begin()) left = lclust;
        decltype(rClusts.begin()) right = rclust;
        // If read 1 comes after read 2, then:
        // left = read 2
        // right = read 1
        // otherwise, vice versa.
        if (lclust->firstRefPos() > rclust->firstRefPos()) {
          left = rclust;
          right = lclust;
        }

        // FILTER 1
        // filter read pairs based on the fragment length which is approximated by the distance between the left most start and right most hit end
        size_t fragmentLen = right->lastRefPos() + right->lastMemLen() - left->firstRefPos();
        if ( fragmentLen < maxFragmentLength) {
          // This will add a new potential mapping. Coverage of a mapping for read pairs is left->coverage + right->coverage
          // If we found a perfect coverage, we would only add those mappings that have the same perfect coverage
          if (maxCoverage < 2 * readLen || (left->coverage + right->coverage) == maxCoverage) {
            jointMemsList.emplace_back(tid, left, right, fragmentLen);
            if (verbose) {
              std::cout <<"tid:"<<tid<<"\n";
              std::cout <<"left\n";
              std::cout <<"leftsize = " << left->mems.size() << "\n";
              for (size_t i = 0; i < left->mems.size(); i++){
                std::cout << "--- t" << left->mems[i].tpos << " r" << left->mems[i].memInfo->rpos << " cid:" << left->mems[i].memInfo->cid << " ori:" << left->isFw;
              }
              std::cout << "\nright\n";
              std::cout <<"rightsize = " << left->mems.size() << "\n";
              for (size_t i = 0; i < right->mems.size(); i++){
                std::cout << "--- t" << right->mems[i].tpos << " r" << right->mems[i].memInfo->rpos << " cid:" << right->mems[i].memInfo->cid << " ori:" << right->isFw;
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
  // FILTER 2
  // filter read pairs that don't have enough base coverage (i.e. their coverage is less than half of the maximum coverage for this read)
  double coverageRatio = 0.5;
  // if we've found a perfect match, we will erase any match that is not perfect
  if (maxCoverage == 2*readLen)
    jointMemsList.erase(std::remove_if(jointMemsList.begin(), jointMemsList.end(),
                                        [&maxCoverage](util::JointMems& pairedReadMems) -> bool {
                                          return pairedReadMems.coverage() != maxCoverage;
                                        }),
                         jointMemsList.end());
  else // ow, we will go with the heuristic that just keep those mappings that have at least half of the maximum coverage
    jointMemsList.erase(std::remove_if(jointMemsList.begin(), jointMemsList.end(),
                                     [&maxCoverage, coverageRatio](util::JointMems& pairedReadMems) -> bool {
                                       return pairedReadMems.coverage() < coverageRatio*maxCoverage ;
                                     }),
                      jointMemsList.end());
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
void populatePaths(util::MemInfo& smem, util::MemInfo& emem, std::string& path, PufferfishIndexT& pfi, uint32_t tid, std::map<uint32_t, util::ContigCecheBlock>& contigSeqCache, size_t readGapDist){

  CanonicalKmer::k(pfi.k()) ;

  auto s = smem.memInfo ;
  auto e = emem.memInfo ;

  auto stpos = smem.tpos ;
  //auto etpos = emem.tpos ;
  auto sCid = s->cid ;

  auto eCid = e->cid ;

  //bool sCFwd = s->cIsFw ;

  //auto& seq = pfi.getSeq() ;
  auto& edges = pfi.getEdge() ;

  std::map<uint32_t,bool> visited ;
  //auto dist = std::abs(s->rpos - e->rpos) ;
  uint32_t numOfNodesVisited{0} ;


  std::queue<uint32_t> queue ;
  //bool first{false} ;
  queue.push(sCid) ;
  path = "" ;

  auto remainingLen = readGapDist ;

  //at worst we would visit readGapDist number of nodes 
  while(numOfNodesVisited < readGapDist and !queue.empty()){
    auto cid = queue.front() ;
    numOfNodesVisited += 1 ;
    visited[cid] = true ;
    queue.pop() ;
    auto contigLen = pfi.getContigLen(cid) ;

    //traverse this node 
    auto toClip = (contigLen <= remainingLen)?(contigLen):(remainingLen) ;
    remainingLen -= toClip ;
    path += pfi.getSeqStr(pfi.getGlobalPos(cid),toClip) ;
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
        auto nextCid = nextHit.contigID() ;

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



template <typename PufferfishIndexT>
void goOverClust(PufferfishIndexT& pfi, std::vector<util::MemCluster>::iterator clust, std::vector<std::pair<std::string,std::string>>& paths, std::string& readSeq, std::map<uint32_t, util::ContigCecheBlock>& contigSeqCache, uint32_t tid){

  //size_t readLen = readSeq.length() ;
  size_t overhangLeft ;
  //size_t overhangRight ;
  size_t clustSize = clust->mems.size() ;


  overhangLeft = (clust->isFw)?(clust->mems[0].memInfo->rpos):(clust->mems[clustSize-1].memInfo->rpos) ;
  //overhangRight = (clust->isFw)?(readLen - (clust->mems[clustSize-1].memInfo->rpos + clust->mems[clustSize-1].memInfo->memlen)):(readLen - (clust->mems[0].memInfo->rpos + clust->mems[0].memInfo->memlen)) ;

  //if start and end has overlap then return
  if(clustSize > 1){
    bool overlap = (clust->isFw)?((clust->mems[0].memInfo->rpos + clust->mems[0].memInfo->memlen) > clust->mems[clustSize-1].memInfo->rpos):((clust->mems[clustSize-1].memInfo->rpos + clust->mems[clustSize-1].memInfo->memlen) > clust->mems[0].memInfo->rpos);
    if(overlap)
      return ;
  }

  for(size_t i = 0; i < clust->mems.size() - 1; i++){

      //search if contig sequence is in the cacheBlock and to clipon the same contig explore paths between unimems from graphs 
    auto startMem = clust->mems[i] ;
    auto endMem = clust->mems[i+1] ;
    std::string readgap = "" ;
    std::string path = "";
    
    auto readGapDist = (clust->isFw)?(endMem.memInfo->rpos - (startMem.memInfo->rpos + startMem.memInfo->memlen)):(startMem.memInfo->rpos - (endMem.memInfo->rpos + endMem.memInfo->memlen)) ;

    if(readGapDist >= 0){
      //clip the appropriate read sequence
      if(clust->isFw){
        readgap = readSeq.substr(startMem.memInfo->rpos + startMem.memInfo->memlen, readGapDist) ;
      }
      else{
        std::string tmp(readSeq.substr(endMem.memInfo->rpos + endMem.memInfo->memlen, readGapDist)) ;
        readgap = util::reverseComplement(tmp) ;
      }

      if(startMem.memInfo->cid  != endMem.memInfo->cid){
        populatePaths(startMem, endMem, path, pfi, tid, contigSeqCache, readGapDist) ;
      }else if(readGapDist > 0){
        //it's on the same contig just insert the sequence from the contig
        auto clipStart = (clust->isFw)?(startMem.memInfo->cpos + startMem.memInfo->memlen):(endMem.memInfo->cpos + endMem.memInfo->memlen) ;
        std::string tmp = pfi.getSeqStr(pfi.getGlobalPos(startMem.memInfo->cid)+clipStart, readGapDist+THRESHOLD) ;
        path = tmp ;
      }

    }
    paths.push_back({readgap, path}) ;
  }

  //TODO overhangRight include in path 

}

template <typename PufferfishIndexT>
void traverseGraph(std::string& leftReadSeq,
                   std::string& rightReadSeq,
                   util::JointMems& hit,
                   PufferfishIndexT& pfi,
                   std::vector<std::pair<std::string,std::string>>& leftClustPaths,
                   std::vector<std::pair<std::string,std::string>>& rightClustPaths,
                   std::map<uint32_t, util::ContigCecheBlock>& contigSeqCache){

  //for all memes in left memcluster ;
  auto tid = hit.tid ;


  //void goOverClust(PufferfishIndexT& pfi, std::vector<util::MemCluster>::iterator clust, std::vector<std::pair<std::string,std::string>>& paths, std::string& readSeq, std::map<uint32_t, std::string>& contigSeqCache, uint32_t tid){

  goOverClust(pfi, hit.leftClust, leftClustPaths, leftReadSeq, contigSeqCache, tid) ;
  goOverClust(pfi, hit.rightClust, rightClustPaths, rightReadSeq, contigSeqCache, tid) ;

}

template <typename PufferfishIndexT>
void processReadsPair(paired_parser* parser,
                     PufferfishIndexT& pfi,
                     SpinLockT* iomutex,
                     std::shared_ptr<spdlog::logger> outQueue,
                     HitCounters& hctr,
                     AlignmentOpts* mopts){
  MemCollector<PufferfishIndexT> memCollector(&pfi) ;

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
      //std::cout << rpair.first.name << "\n";
      bool verbose = false;// rpair.first.name == "read17739350/ENST00000421512;mate1:506-605;mate2:614-712";
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

      if(lh && rh){
        joinReadsAndFilter(leftHits, rightHits, jointHits, mopts->maxFragmentLength, readLen, verbose) ;
      }
      else{
        //ignore orphans for now
      }
      //jointHits is a vector
      //this can be used for BFS

      //void traverseGraph(std::string& leftReadSeq, std::string& rightReadSeq, util::JointMems& hit, PufferfishIndexT& pfi,   std::map<uint32_t, std::string>& contigSeqCache){
      bool doTraverse = false;
      if (doTraverse) {
        //TODO Have to make it per thread 
        //have to make write access thread safe
        //std::cout << "traversing graph \n" ;
        std::map<uint32_t, util::ContigCecheBlock> contigSeqCache ;
        using PathType = std::vector<std::pair<std::string,std::string>>  ;
        std::map<uint32_t, std::pair<PathType,PathType>> alnPaths ;
        for(auto& hit : jointHits){
          PathType lpath ;
          PathType rpath ;
          traverseGraph(rpair.first.seq, rpair.second.seq, hit, pfi, lpath, rpath, contigSeqCache) ;
          alnPaths[hit.tid] = {lpath,rpath} ;
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
