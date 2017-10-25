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

using paired_parser = fastx_parser::FastxParser<fastx_parser::ReadPair>;
using AlignmentOpts = util::AlignmentOpts ;
using HitCounters = util::HitCounters ;
using QuasiAlignment = util::QuasiAlignment ;
using MateStatus = util::MateStatus ;

using SpinLockT = std::mutex ;


inline void joinReverseOrientationMems(size_t tid,
                                       std::vector<util::MemCluster>& fwClusters,
                                       std::vector<util::MemCluster>& rcClusters,
                                       std::vector<util::JointMems>& jointMemsList,
                                       uint32_t maxFragmentLength,
                                       uint32_t& maxCoverage, // pass it by reference since it should be updated by both pairs <fw, rc> and <rc, fw>
                                       bool verbose) {
  bool isLeftFw = true;
  for (auto lclust =  fwClusters.begin(); lclust != fwClusters.end(); lclust++) {
    for (auto rclust =  rcClusters.begin(); rclust != rcClusters.end(); rclust++) {
      decltype(fwClusters.begin()) left = lclust;
      decltype(rcClusters.begin()) right = rclust;
      isLeftFw = true;
      if (lclust->mems[0].tpos > rclust->mems[0].tpos) {
        std::cerr<<"isLeftFw?: " << isLeftFw << " , and this happened\n";
        //std::swap(lclust, rclust);
        left = rclust;
        right = lclust;
        isLeftFw = !isLeftFw;
      }
      // FILTER 1
      // filter read pairs based on the fragment length which is approximated by the distance between the left most start and right most hit end
      size_t fragmentLen = right->mems.back().tpos + right->mems.back().memlen - left->mems.front().tpos;
      if ( fragmentLen < maxFragmentLength) {
        if (verbose) {
        std::cout << "jointMemsList start\n";
        }
        jointMemsList.emplace_back(tid, isLeftFw, *left, *right, fragmentLen);
        if (verbose) {
        std::cout << "jointMemsList end.\n";
        }
        if (verbose) {
          const auto& last = jointMemsList.back();
          std::cout << isLeftFw << "\n";
          std::cout <<"left\n";
          std::cout <<"leftsize = " << left->mems.size() << ", last leftMems size = " << last.leftMems.mems.size() << "\n";
          for (size_t i = 0; i < left->mems.size(); i++){
            std::cout << " t" << left->mems[i].tpos << " t2:" << last.leftMems.mems[i].tpos << " l" << left->mems[i].memlen << "\n";
          }
          std::cout << "right\n";
          for (auto& mem : last.rightMems.mems) {
            std::cout <<  " t" << mem.tpos << " r" << mem.rpos << " l" << mem.memlen << "\n";
          }
        }
        uint32_t currCoverage =  jointMemsList.back().coverage;
        if (maxCoverage < currCoverage) {
          maxCoverage = currCoverage;
        }
      }
      else {
        break;
      }
    }
  }
}

void joinReadsAndFilter(spp::sparse_hash_map<size_t,
                        util::TrClusters>& leftMemClusters,
                        spp::sparse_hash_map<size_t, util::TrClusters>& rightMemClusters,
                        std::vector<util::JointMems>& jointMemsList,
                        uint32_t maxFragmentLength) {
  //orphan reads should be taken care of maybe with a flag!
  uint32_t maxCoverage{0};
  for (auto& leftClustItr : leftMemClusters) {
    size_t tid = leftClustItr.first;
    auto& lClusts = leftClustItr.second;
    auto& rClusts = rightMemClusters[tid];
    joinReverseOrientationMems(tid, lClusts.fwClusters, rClusts.rcClusters, jointMemsList, maxFragmentLength, maxCoverage, tid == 6);
    joinReverseOrientationMems(tid, rClusts.fwClusters, lClusts.rcClusters, jointMemsList, maxFragmentLength, maxCoverage, tid == 6);
  }
  //std::cerr << "mc:" << maxCoverage << "\n";
  // FILTER 2
  // filter read pairs that don't have enough base coverage (i.e. their coverage is less than half of the maximum coverage for that transcript)
  jointMemsList.erase(std::remove_if(jointMemsList.begin(), jointMemsList.end(),
                                     [&maxCoverage](util::JointMems& pairedReadMems) -> bool {
                                       return pairedReadMems.coverage < 0.5*maxCoverage ;
                                     }),
                      jointMemsList.end());
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


  spp::sparse_hash_map<size_t, util::TrClusters> leftHits ;
  spp::sparse_hash_map<size_t, util::TrClusters> rightHits ;
  std::vector<util::JointMems> jointHits ;
  PairedAlignmentFormatter<PufferfishIndexT*> formatter(&pfi);


  auto rg = parser->getReadGroup() ;

  while(parser->refill(rg)){
    for(auto& rpair : rg){
      readLen = rpair.first.seq.length() ;
      ++hctr.numReads ;

      jointHits.clear() ;
      leftHits.clear() ;
      rightHits.clear() ;


      //help me to debug, will deprecate later
      //std::cerr << "\n first seq in pair " << rpair.first.seq << "\n" ;
      //std::cerr << "\n second seq in pair " << rpair.second.seq << "\n" ;

      //std::cerr << "\n going inside hit collector \n" ;
      //readLen = rpair.first.seq.length() ;

      bool lh = memCollector(rpair.first.seq,
                             leftHits,
                             mopts->maxSpliceGap
                             /*,
                             MateStatus::PAIRED_END_LEFT,
                             mopts->consistentHits,
                             refBlocks*/) ;

      bool rh = memCollector(rpair.second.seq,
                             rightHits,
                             mopts->maxSpliceGap
                             /*,
                             MateStatus::PAIRED_END_RIGHT,
                             mopts->consistentHits,
                             refBlocks*/) ;
      //do intersection on the basis of
      //performance, or going towards selective alignment
      //otherwise orphan

      if(lh && rh){
        joinReadsAndFilter(leftHits, rightHits, jointHits, mopts->maxFragmentLength) ;
      }
      else{
        //ignore orphans for now
      }


      //TODO Write them to a sam file
      hctr.totHits += jointHits.size() ;

      //std::cerr << "\n Number of total joint hits" << jointHits.size() << "\n" ;
      //TODO When you get to this, you should be done aligning!!
      //fill the QuasiAlignment list
      std::vector<QuasiAlignment> jointAlignments;
      for (auto& jointHit : jointHits) {
        if (jointHit.tid == 6) {
          std::cout << "ltpos:" << jointHit.leftMems.getTrFirstHitPos() << " isFw: " <<  jointHit.isLeftFw <<
            " fragmentLen:" << jointHit.fragmentLen << " rtpos:" << jointHit.rightMems.getTrFirstHitPos() << "\n";
        }
        jointAlignments.emplace_back(jointHit.tid,           // reference id
                                     jointHit.leftMems.getTrFirstHitPos(),     // reference pos
                                     jointHit.isLeftFw ,     // fwd direction
                                     readLen, // read length
                                     jointHit.fragmentLen,       // fragment length
                                     true);         // properly paired
         // Fill in the mate info		         // Fill in the mate info
        auto& qaln = jointAlignments.back();
         qaln.mateLen = readLen;
         qaln.matePos = jointHit.rightMems.getTrFirstHitPos();
         qaln.mateIsFwd = !jointHit.isLeftFw;
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
  std::cerr << "\n In spawn threads\n" ;

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


  consoleLog->info("Done mapping reads.");
  consoleLog->info("In total saw {} reads.", hctrs.numReads);
  consoleLog->info("Final # hits per read = {}", hctrs.totHits / static_cast<float>(hctrs.numReads));
	consoleLog->info("flushing output queue.");
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
