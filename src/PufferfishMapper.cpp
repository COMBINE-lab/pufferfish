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
#include "HitCollector.hpp"
#include "SAMWriter.hpp"

using paired_parser = fastx_parser::FastxParser<fastx_parser::ReadPair>;
using MappingOpts = util::MappingOpts ;
using HitCounters = util::HitCounters ;
using QuasiAlignment = util::QuasiAlignment ;
using MateStatus = util::MateStatus ;

using SpinLockT = std::mutex ;


void mergeHits(std::vector<QuasiAlignment>& leftHits, std::vector<QuasiAlignment>& rightHits, std::vector<QuasiAlignment>& jointHits){
  /*
  //sort leftHits
  if(leftHits.size() > 1)
    std::sort(leftHits.begin(),leftHits.end(),
              [](QuasiAlignment& q1, QuasiAlignment& q2) -> bool{
                return q1.tid < q2.tid ;
              });

  //sort rightHits
  if(rightHits.size() > 1)
    std::sort(rightHits.begin(),rightHits.end(),
              [](QuasiAlignment& q1, QuasiAlignment& q2) -> bool{
                return q1.tid < q2.tid ;
              });
  */
  /**
   * based on potential standard implementation :http://en.cppreference.com/w/cpp/algorithm/set_intersection
   **/

  // We will not handle chimetric reads
  int32_t signedZero{0};
  auto rqIt = rightHits.begin();
  auto rqEnd = rightHits.end();
  auto lqIt = leftHits.begin();
  auto lqEnd = leftHits.end();
  while (lqIt != lqEnd and rqIt != rqEnd) {
    auto lid = lqIt->tid ;
    auto rid = rqIt->tid;
    if (lid < rid) {
      ++lqIt;
    } else {
    //for(auto rqIt = rightHits.begin(); rqIt != rightHits.end(); ++rqIt) {
      //for(auto& rq : rightHits){
      //for now let's ignore chimeric reads
      if(!(rid < lid)){
        int32_t startRead1 = std::max(lqIt->pos, signedZero);
        int32_t startRead2 = std::max(rqIt->pos, signedZero);
        bool read1First{(startRead1 < startRead2)};
        int32_t fragStartPos = read1First ? startRead1 : startRead2;
        int32_t fragEndPos = read1First ?
          (startRead2 + rqIt->readLen) : (startRead1 + lqIt->readLen);
        uint32_t fragLen = fragEndPos - fragStartPos;
        jointHits.emplace_back(lid,           // reference id
                               lqIt->pos,     // reference pos
                               lqIt->fwd,     // fwd direction
                               lqIt->readLen, // read length
                               fragLen,       // fragment length
                               true);         // properly paired

        // Fill in the mate info
        auto& qaln = jointHits.back();
        qaln.mateLen = rqIt->readLen;
        qaln.matePos = rqIt->pos;
        qaln.mateIsFwd = rqIt->fwd;
        qaln.mateStatus = MateStatus::PAIRED_END_PAIRED;
        ++lqIt;
      }
      ++rqIt;
    } // end else
  } // end while

    //if there is no jointHits then
    //add orphan reads to the list
    if(jointHits.size() == 0){
        if(leftHits.size() + rightHits.size() > 0){
	        //std::cerr<<"orphans here\n";
            jointHits.insert(jointHits.end(),
                    std::make_move_iterator(leftHits.begin()),
                    std::make_move_iterator(leftHits.end()));
            jointHits.insert(jointHits.end(),
                    std::make_move_iterator(rightHits.begin()),
                    std::make_move_iterator(rightHits.end()));

        }
    }

}

template <typename PufferfishIndexT>
void processReadsPair(paired_parser* parser,
                     PufferfishIndexT& pfi,
                     SpinLockT* iomutex,
                     std::shared_ptr<spdlog::logger> outQueue,
                     HitCounters& hctr,
                     MappingOpts* mopts){
  HitCollector<PufferfishIndexT> hitCollector(&pfi) ;

  //std::cerr << "\n In process reads pair\n" ;
    //TODO create a memory layout to store
    //strings then will allocate alignment to them
    //accordingly
  std::vector<std::string> refBlocks ;

  auto logger = spdlog::get("stderrLog") ;
  fmt::MemoryWriter sstream ;
  //size_t batchSize{2500} ;
  size_t readLen{0} ;


  std::vector<QuasiAlignment> leftHits ;
  std::vector<QuasiAlignment> rightHits ;
  std::vector<QuasiAlignment> jointHits ;
  PairedAlignmentFormatter<PufferfishIndexT*> formatter(&pfi);


  auto rg = parser->getReadGroup() ;

  while(parser->refill(rg)){
    for(auto& rpair : rg){
      //readLen = rpair.first.seq.length() ;
      ++hctr.numReads ;

      jointHits.clear() ;
      leftHits.clear() ;
      rightHits.clear() ;


      //help me to debug, will deprecate later
      //std::cerr << "\n first seq in pair " << rpair.first.seq << "\n" ;
      //std::cerr << "\n second seq in pair " << rpair.second.seq << "\n" ;

      //std::cerr << "\n going inside hit collector \n" ;
      readLen = rpair.first.seq.length() ;

      bool lh = hitCollector(rpair.first.seq,
                             leftHits,
                             MateStatus::PAIRED_END_LEFT,
                             mopts->consistentHits,
                             refBlocks) ;

      bool rh = hitCollector(rpair.second.seq,
                             rightHits,
                             MateStatus::PAIRED_END_RIGHT,
                             mopts->consistentHits,
                             refBlocks) ;
      //do intersection on the basis of
      //performance, or going towards selective alignment
      //otherwise orphan

      if(lh && rh){
        mergeHits(leftHits, rightHits, jointHits) ;
      }
      else{
        //ignore orphans for now
      }


      //TODO Write them to a sam file
      hctr.totHits += jointHits.size() ;

      //std::cerr << "\n Number of total joint hits" << jointHits.size() << "\n" ;
      if(jointHits.size() > 0 and !mopts->noOutput){
        writeAlignmentsToStream(rpair, formatter, jointHits, sstream, mopts->writeOrphans) ;
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
                              MappingOpts* mopts){

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
bool mapReads(
              PufferfishIndexT& pfi,
              std::shared_ptr<spdlog::logger> consoleLog,
              MappingOpts* mopts) {



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

int pufferfishMapper(MappingOpts& mapargs){

  //auto outputSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(outStream) ;
  //std::shared_ptr<spdlog::logger> outLog = std::make_shared<spdlog::logger>("puffer::outLog",outputSink) ;

//auto rawConsoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
//auto rawConsoleSink = spdlog::sinks::stderr_sink_mt::instance() ;
 // auto consoleSink = std::make_shared<spdlog::sinks::ansicolor_sink>(rawConsoleSink);
    auto consoleLog = spdlog::stderr_color_mt("console");


    bool success{false} ;

  auto indexDir = mapargs.indexDir ;

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
    success = mapReads(pfi, consoleLog, &mapargs) ;

  }else if(indexType == "sparse"){
    PufferfishSparseIndex pfi(indexDir) ;
    success = mapReads(pfi, consoleLog, &mapargs) ;
  }

  if (!success) {
    consoleLog->warn("Problem mapping.");
  }
  return 0 ;
}
