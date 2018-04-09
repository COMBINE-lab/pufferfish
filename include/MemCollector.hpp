#ifndef HIT_COLLECTOR_HPP
#define HIT_COLLECTOR_HPP

#include "CommonTypes.hpp"
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "Util.hpp"
#include "edlib.h"
#include "jellyfish/mer_dna.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <sparsepp/spp.h>
// using spp:sparse_hash_map;

// polute the namespace --- put this in the functions that need it.
namespace kmers = combinelib::kmers;

template <typename PufferfishIndexT> class MemCollector {
  enum class ExpansionTerminationType : uint8_t { MISMATCH = 0, CONTIG_END, READ_END };  

public:
  MemCollector(PufferfishIndexT* pfi) : pfi_(pfi) { k = pfi_->k(); }

  void recoverGaps(spp::sparse_hash_map<pufferfish::common_types::ReferenceID, 
                                        std::vector<util::MemCluster>>& memClustersMap,
                                        std::vector<util::UniMemInfo>& memCollection,
                                        size_t rlen, 
                                        bool verbose=false) {
    if (verbose)
      std::cerr << "\n[RECOVERGAPS]\n";
    uint64_t thresh = k/2;
    // a sample for <readPos,memLen> if read is mapped in rc against the reference
    // and the list is sorted based on reference pos
    //read poses 20,47  9,41  4,35  0,34
    //trns poses 100,47 117,41
    for (auto& kv : memClustersMap) {
      //auto& refId = kv.first;
      auto& memClusters = kv.second;
      for (auto& memCluster : memClusters) {
        // If we have reasonable coverage
        // at least 4 kmers here
        if (memCluster.coverage >= k + 3) { // 3 is a complete heuristic! 
         
          auto rfirstMem = memCluster.mems[0].memInfo;
          auto rlastMem = memCluster.mems.back().memInfo;
          size_t tfirstpos = memCluster.mems[0].tpos;
          size_t tlastpos = memCluster.mems.back().tpos;
            
          if (memCluster.isFw) {
            uint32_t firstGap = rfirstMem->rpos;
            uint32_t lastGap = rlen - (rlastMem->rpos + rlastMem->memlen);

            // we also check that the read is not an overhang
            // meaning that it hangs over the transcript. 
            // In that case, we should add as many nucleotides as left in transcript's end
            if (firstGap > 0 && firstGap <= thresh) {
              size_t tpos = tfirstpos - firstGap;
              int potentialtpos = tfirstpos - firstGap;
              if (potentialtpos < 0) {
                firstGap = firstGap + potentialtpos;
                tpos = 0;
              }
              if (verbose) {

                std::cerr << "\nbefore: " << memCluster.isFw << " " << memCluster.coverage <<"\n ";
                for (auto mem : memCluster.mems) {              
                        std::cerr << " r:" << mem.memInfo->rpos 
                        << " t:" << mem.tpos
                        << " len:" << mem.memInfo->memlen << " -- ";
                }
              }
              memCollection.emplace_back(rfirstMem->cid, rfirstMem->cIsFw,
                                    0, firstGap-1, rfirstMem->cpos,
                                    rfirstMem->cGlobalPos, rfirstMem->clen);
              memCluster.addMem(std::prev(memCollection.end()),
                                    tpos , 0);
            }
            if (lastGap > 0 && lastGap <= thresh) {
              // FIXME -> what if it passes the transcript length?
              // we need transcript len for this!!
              size_t tpos = tlastpos +  rlastMem->memlen + 1;
              if (verbose) {

                std::cerr << "\nbefore: " << memCluster.isFw << " " << memCluster.coverage <<"\n ";
                for (auto mem : memCluster.mems) {              
                        std::cerr << " r:" << mem.memInfo->rpos 
                        << " t:" << mem.tpos
                        << " len:" << mem.memInfo->memlen << " -- ";
                }
              }
              memCollection.emplace_back(rlastMem->cid, rlastMem->cIsFw,
                                    rlastMem->rpos+rlastMem->memlen+1, lastGap-1, rlastMem->cpos,
                                    rlastMem->cGlobalPos, rlastMem->clen);
              memCluster.addMem(std::prev(memCollection.end()),
                                    tpos, memCluster.mems.size());  
              if (verbose) {
                std::cerr << "\nafter: " << memCluster.coverage << "\n ";
                    for (auto mem : memCluster.mems) {
                  
                            std::cerr << " r:" << mem.memInfo->rpos 
                            << " t:" << mem.tpos
                            << " len:" << mem.memInfo->memlen << " -- ";
                }
                std::cerr << "\n";
              }
            }  
            std::vector<util::MemInfo>::size_type siz = memCluster.mems.size();
            for (size_t i = 0; i < siz - 1; i++) {
              auto& prev = memCluster.mems[i];
              auto& next = memCluster.mems[i+1];
              uint64_t middleGap = next.memInfo->rpos - (prev.memInfo->rpos+prev.memInfo->memlen);
              //uint64_t refGap = next.tpos - (prev.tpos+prev.memInfo->memlen);
              if ( middleGap > 2 && middleGap <= thresh+1) {
                // first of all, since it's READ coverage, you should see what the gap is in read
                // also, to make it closer to edit distance, 
                // we can consider 2 + extra characters in either read or reference as the edit distance
              if (verbose) {

                std::cerr << "\nbefore: " << memCluster.isFw << " " << memCluster.coverage <<"\n ";
                for (auto mem : memCluster.mems) {              
                        std::cerr << " r:" << mem.memInfo->rpos 
                        << " t:" << mem.tpos
                        << " len:" << mem.memInfo->memlen << " -- ";
                }
              }
                uint64_t readPos = prev.memInfo->rpos + prev.memInfo->memlen + 1;
                memCollection.emplace_back(prev.memInfo->cid, prev.memInfo->cIsFw,
                                    readPos, middleGap-2, prev.memInfo->cpos,
                                    prev.memInfo->cGlobalPos, prev.memInfo->clen);
                memCluster.addMem(std::prev(memCollection.end()),
                                    prev.tpos+prev.memInfo->memlen + 1, i+1);
              if (verbose) {
                std::cerr << "\nafter: " << memCluster.coverage << "\n ";
                    for (auto mem : memCluster.mems) {
                  
                            std::cerr << " r:" << mem.memInfo->rpos 
                            << " t:" << mem.tpos
                            << " len:" << mem.memInfo->memlen << " -- ";
                }
                std::cerr << "\n";
              }
                siz++;
                i++; 
                
              }
            }
          } else {
            uint32_t firstGap = rlen - (rfirstMem->rpos + rfirstMem->memlen);
            uint32_t lastGap = rlastMem->rpos;

            // also check that the read is not an overhang
            // meaning that it hangs over the transcript. 
            // In that case, we should add as many nucleotides as left in transcript's end
            if (firstGap > 0 && firstGap <= thresh) {
              size_t tpos = tfirstpos - firstGap;
              int potentialtpos = tfirstpos - firstGap;
              if (potentialtpos < 0) {
                firstGap = firstGap + potentialtpos;
                tpos = 0;
              }
              if (verbose) {

                std::cerr << "\nbefore: " << memCluster.isFw << " " << memCluster.coverage <<"\n ";
                for (auto mem : memCluster.mems) {              
                        std::cerr << " r:" << mem.memInfo->rpos 
                        << " t:" << mem.tpos
                        << " len:" << mem.memInfo->memlen << " -- ";
                }
              }

              memCollection.emplace_back(rfirstMem->cid, rfirstMem->cIsFw,
                                    rfirstMem->rpos+rfirstMem->memlen+1, 
                                    firstGap-1, rfirstMem->cpos,
                                    rfirstMem->cGlobalPos, rfirstMem->clen);
              memCluster.addMem(std::prev(memCollection.end()),
                                    tpos, 0);  
              if (verbose) {
                std::cerr << "\nafter: " << memCluster.coverage << "\n ";
                    for (auto mem : memCluster.mems) {
                  
                            std::cerr << " r:" << mem.memInfo->rpos 
                            << " t:" << mem.tpos
                            << " len:" << mem.memInfo->memlen << " -- ";
                }
                std::cerr << "\n";
              }
            }
            if (lastGap > 0 && lastGap <= thresh) {
              // FIXME -> what if it passes the transcript length?
              // we need transcript len for this!!
              if (verbose) {

                std::cerr << "\nbefore: " << memCluster.isFw << " " << memCluster.coverage <<"\n ";
                for (auto mem : memCluster.mems) {              
                        std::cerr << " r:" << mem.memInfo->rpos 
                        << " t:" << mem.tpos
                        << " len:" << mem.memInfo->memlen << " -- ";
                }
              }
              size_t tpos = tlastpos +  rlastMem->memlen + 1;
              memCollection.emplace_back(rlastMem->cid, rlastMem->cIsFw,
                                    0, lastGap-1, rlastMem->cpos,
                                    rlastMem->cGlobalPos, rlastMem->clen);
              memCluster.addMem(std::prev(memCollection.end()),
                                    tpos, memCluster.mems.size());  
              if (verbose) {
                std::cerr << "\nafter: " << memCluster.coverage << "\n ";
                    for (auto mem : memCluster.mems) {
                  
                            std::cerr << " r:" << mem.memInfo->rpos 
                            << " t:" << mem.tpos
                            << " len:" << mem.memInfo->memlen << " -- ";
                }
                std::cerr << "\n";
              }

              
            }  

            std::vector<util::MemInfo>::size_type siz = memCluster.mems.size();
            for (size_t i = siz-1; i > 0; i--) {
              auto& prev = memCluster.mems[i];
              auto& next = memCluster.mems[i-1];
              uint64_t middleGap = next.memInfo->rpos - (prev.memInfo->rpos+prev.memInfo->memlen);
              //uint64_t refGap = prev.tpos - (next.tpos+next.memInfo->memlen);
              if ( middleGap > 2 && middleGap <= thresh+1) {
              if (verbose) {

                std::cerr << "\nbefore: " << memCluster.isFw << " " << memCluster.coverage <<"\n ";
                for (auto mem : memCluster.mems) {              
                        std::cerr << " r:" << mem.memInfo->rpos 
                        << " t:" << mem.tpos
                        << " len:" << mem.memInfo->memlen << " -- ";
                }
              }
                uint64_t readPos = prev.memInfo->rpos + prev.memInfo->memlen + 1;
                memCollection.emplace_back(prev.memInfo->cid, prev.memInfo->cIsFw,
                                    readPos, middleGap-2, prev.memInfo->cpos,
                                    prev.memInfo->cGlobalPos, prev.memInfo->clen);
                memCluster.addMem(std::prev(memCollection.end()),
                                  prev.tpos - (middleGap-2) - 1, i);
               
              if (verbose) {
                std::cerr << "\nafter: " << memCluster.coverage << "\n ";
                    for (auto mem : memCluster.mems) {
                  
                            std::cerr << " r:" << mem.memInfo->rpos 
                            << " t:" << mem.tpos
                            << " len:" << mem.memInfo->memlen << " -- ";
                }
                std::cerr << "\n";
              }
                //size++;
                //i++;
              }
            }
          }
        }
      }
    }
    if (verbose) {
      std::cerr << "\nfinal memCollection size: " << memCollection.size() << "\n";
    }
    if (verbose)
      std::cerr << "\n[END OF RECOVERGAPS]\n";

  }

  bool clusterMems(std::vector<std::pair<int, util::ProjectedHits>>& hits,
                   spp::sparse_hash_map<pufferfish::common_types::ReferenceID, std::vector<util::MemCluster>>& memClusters,
                   uint32_t maxSpliceGap, std::vector<util::UniMemInfo>& memCollection, bool verbose = false) {
        if (verbose)
      std::cerr << "\n[CLUSTERMEMS]\n";

    using namespace pufferfish::common_types;
    namespace kmers = combinelib::kmers;
    //(void)verbose;

    if (hits.empty()) {
      return false;
    }

    // Map from (reference id, orientation) pair to a cluster of MEMs.
    std::map<std::pair<ReferenceID, bool>, std::vector<util::MemInfo>>
        trMemMap;
    // here we guarantee that even if later we fill up
    // every gap between the hits and before the first and after the last hit
    // still the memCollection size doesn't change and hence the pointers are valid
    memCollection.reserve(200*2*hits.size() + 1); 
    if (verbose)
      std::cerr<<"\nreserved memCollection size: " << 200*2*hits.size()+1 << "\n";
    for (auto& hit : core::range<decltype(hits.begin())>(hits.begin(), hits.end())) {
      auto& readPos = hit.first;
      auto& projHits = hit.second;
      // NOTE: here we rely on internal members of the ProjectedHit (i.e., member variables ending in "_").
      // Maybe we want to change the interface (make these members public or provide accessors)?
      auto& refs = projHits.refRange;
      if (refs.size() < 200) {
        memCollection.emplace_back(projHits.contigIdx_, projHits.contigOrientation_,
                                   readPos, projHits.k_, projHits.contigPos_,
                                   projHits.globalPos_-projHits.contigPos_, projHits.contigLen_);
        auto memItr = std::prev(memCollection.end());
        for (auto& posIt : refs) {
          auto refPosOri = projHits.decodeHit(posIt);
          trMemMap[std::make_pair(posIt.transcript_id(), refPosOri.isFW)]
            .emplace_back(memItr, refPosOri.pos);
        }
      }
    }

    for (auto& trMem : core::range<decltype(trMemMap.begin())>(trMemMap.begin(), trMemMap.end())) {
      auto& trOri = trMem.first;
      auto& tid = trOri.first;
      auto& isFw = trOri.second;
      auto& memList = trMem.second;
      // sort memList according to mem reference positions
      if (isFw)
        std::sort(memList.begin(), memList.end(),
                [](util::MemInfo& q1, util::MemInfo& q2) -> bool {
                  return q1.tpos == q2.tpos? q1.memInfo->rpos < q2.memInfo->rpos:q1.tpos < q2.tpos; // sort based on tpos
                  });
      else
        std::sort(memList.begin(), memList.end(),
                  [](util::MemInfo& q1, util::MemInfo& q2) -> bool {
                    return q1.tpos == q2.tpos? q1.memInfo->rpos > q2.memInfo->rpos:q1.tpos < q2.tpos; // sort based on tpos
                  });

      std::vector<util::MemCluster> currMemClusters;
      // cluster MEMs so that all the MEMs in one cluster are concordant.
      for (auto& hit : core::range<decltype(memList.begin())>(memList.begin(), memList.end())) {
        //bool foundAtLeastOneCluster = false;
        //bool gapIsSmall = false;
        bool addNewCluster = currMemClusters.size() == 0;
        bool foundAtLeastOne = false;
        //bool foundPrev = false;
        for (auto prevClus = currMemClusters.rbegin();
             prevClus != currMemClusters.rend(); prevClus++) {
          if (hit.tpos - prevClus->getTrLastHitPos() < maxSpliceGap) {
            // if the distance between last mem and the new one is NOT longer than maxSpliceGap
            //gapIsSmall = true;
            if (
                (hit.tpos >= (prevClus->getTrLastHitPos() + prevClus->getTrLastMemLen()) &&
                 (
                  (isFw && hit.memInfo->rpos >= (prevClus->getReadLastHitPos() + prevClus->getTrLastMemLen())) ||
                  (!isFw && (hit.memInfo->rpos + hit.memInfo->memlen) <= prevClus->getReadLastHitPos())
                  )) ||
                // then hit.tpos < (prevClus->getTrLastHitPos() + prevClus->getTrLastMemLen() and
                //(hit.tpos < (prevClus->getTrLastHitPos() + prevClus->getTrLastMemLen())) and
                 (isFw && (prevClus->getTrLastHitPos() + prevClus->getTrLastMemLen() - hit.tpos) == (prevClus->getReadLastHitPos() + prevClus->getTrLastMemLen()-hit.memInfo->rpos)) ||
                 (!isFw && (prevClus->getTrLastHitPos() + prevClus->getTrLastMemLen() - hit.tpos) == (hit.memInfo->rpos + hit.memInfo->memlen - prevClus->getReadLastHitPos()))
                ) {
              // NOTE: Adds a new mem to the list of cluster mems and updates the coverage
              prevClus->addMem(hit.memInfo, hit.tpos);
              //if(verbose){std::cerr << "len after putting in cluster: "<<hit.memInfo->memlen<<"\n";}
              foundAtLeastOne = true;
              //foundPrev = true;
              //addNewCluster = false;
              //break;
            }
          } else {
            if (!foundAtLeastOne)
              addNewCluster = true;
            break;
          }
        }

        if (addNewCluster) {
          //auto prevLastIndex = static_cast<int32_t>(currMemClusters.size()) - 1;
          // Create the new clusters on the end of the currMemClusters vector
          currMemClusters.emplace_back(isFw);
          auto& newClus = currMemClusters.back();
          /*if ((prevLastIndex > 0) and gapIsSmall) {
            auto& lastClus = currMemClusters[prevLastIndex];
            // add all previous compatable mems before this last one that was crossed
            for (auto mem = lastClus.mems.begin(); mem != lastClus.mems.end() && mem->memInfo->rpos < hit.memInfo->rpos; mem++) {
              newClus.mems.emplace_back(mem->memInfo, mem->tpos);
            }
            }*/
          // NOTE: Adds a new mem to the list of cluster mems and updates the coverage
          newClus.addMem(hit.memInfo, hit.tpos);
          //if(verbose){std::cerr << "len after putting in cluster: "<<hit.memInfo->memlen<<"\n";}
        }
      }
      /*for (auto& clus : currMemClusters) {
        if (clus.mems.size() > 69) {
          std::cerr << "mem size: " << clus.mems.size() << "\n";
          for (auto& mem : clus.mems) {
            std::cerr << "t" << mem.tpos << " r" << mem.memInfo->rpos << " cid" << mem.memInfo->cid << " -- ";
          }
          std::cerr << "\n";
        }
        }*/
      if (verbose) {
        std::cerr << "\ntid" << tid << " " << pfi_->refName(tid) << " , isFw:" << isFw << " cluster size:" << currMemClusters.size() << "\n";
        for (auto& clus : currMemClusters) {
          std::cerr << "mem size: " << clus.mems.size() << "\n";
          for (auto& mem : clus.mems) {
            std::cerr << "t" << mem.tpos << " r" << mem.memInfo->rpos << " cid" << mem.memInfo->cid << " -- ";
          }
          std::cerr << "\n";
        }
      }
      // This is kind of inefficient (copying the currMemClusters while probably we can build it on the fly)
      memClusters[tid].insert(memClusters[tid].end(), std::make_move_iterator(currMemClusters.begin()),
                              std::make_move_iterator(currMemClusters.end()));
    }
        if (verbose)
      std::cerr << "\n[END OF CLUSTERMEMS]\n";

    return true;
  }

  size_t expandHitEfficient(util::ProjectedHits& hit,
                            pufferfish::CanonicalKmerIterator& kit,
			    ExpansionTerminationType& et,
                            bool verbose=false) {

    if(verbose){
      std::cerr <<"\nBefore cpos " << hit.contigPos_ << "\n" ;
    }

    auto& allContigs = pfi_->getSeq();
    // startPos points to the next kmer in contig (which can be the left or
    // right based on the orientation of match)
    size_t cStartPos =
        hit.globalPos_ - hit.contigPos_   ; // next kmer in the read
    size_t cEndPos = cStartPos + hit.contigLen_;
    size_t cCurrPos = hit.globalPos_; // start from next character if fw match
    if (hit.contigOrientation_) { // if match is fw, go to the next k-mer in the
                                  // contig
      cCurrPos += k;
    }

    int currReadStart = kit->second + 1;
    auto readSeqView = kit.seq();
    auto readSeqLen = readSeqView.size();
    auto readSeqStart = currReadStart;
    auto readSeqOffset = currReadStart + k - 1;
    int fastNextReadCode{0};
    bool stillMatch = true;
    bool foundTermCondition = false;

    while (stillMatch and
           (cCurrPos < cEndPos) and
           (cCurrPos > cStartPos) and
           readSeqOffset < readSeqLen) { // over kmers in contig

      if (hit.contigOrientation_) { // if fw match, compare read last base with
                                    // contig first base and move fw in the
                                    // contig
        auto baseCnt = k < cEndPos - cCurrPos ? k : cEndPos - cCurrPos;
        uint64_t fk = allContigs.get_int(2 * (cCurrPos), 2 * baseCnt);
        cCurrPos += baseCnt;
        for (size_t i = 0; i < baseCnt && readSeqOffset < readSeqLen; i++) {
          // be dirty and peek into the underlying read
          fastNextReadCode = kmers::codeForChar(readSeqView[readSeqOffset]);
          int contigCode = (fk >> (2 * i)) & 0x3;
          if (fastNextReadCode != contigCode) {
            stillMatch = false;
	    et = ExpansionTerminationType::MISMATCH;
	    foundTermCondition = true;
            break;
          }
          hit.k_++;
          readSeqOffset++;
          readSeqStart++;
        }
      } else { // if rc match, compare read last base with contig last base and
               // move backward in the contig
        auto baseCnt = k < cCurrPos - cStartPos ? k : cCurrPos - cStartPos;
        uint64_t fk = allContigs.get_int(2 * (cCurrPos - baseCnt), 2 * baseCnt);
        cCurrPos -= baseCnt;
        for (int i = baseCnt - 1; i >= 0 && readSeqOffset < readSeqLen; i--) {
          // be dirty and peek into the underlying read
          fastNextReadCode = kmers::codeForChar(kmers::complement(readSeqView[readSeqOffset]));
          int contigCode = (fk >> (2 * i)) & 0x3;
          if (fastNextReadCode != contigCode) {
            stillMatch = false;
	    et = ExpansionTerminationType::MISMATCH;
	    foundTermCondition = true;
            break;
          }
          hit.k_++;
          readSeqOffset++;
          readSeqStart++;
        }
      }
    }

    if (!foundTermCondition) {
	    et = (cCurrPos >= cEndPos or cCurrPos <= cStartPos) ? 
		    ExpansionTerminationType::CONTIG_END : ExpansionTerminationType::READ_END;
    }

    if (!hit.contigOrientation_) {
      // if (verbose)
      //   std::cerr << hit.k_ << " prev contig pos:" << hit.contigPos_ << "\n";
      hit.contigPos_ -= (hit.k_ - k);
      hit.globalPos_ -= (hit.k_ - k);
      /* if (verbose){
        std::cerr << "after updating: " << hit.contigPos_ << "\n";
      } */
      
    }
    //std::cerr << "after updating coverage: " << hit.k_ << "\n";
    kit.jumpTo(readSeqStart);
    return currReadStart;
  }

  bool operator()(std::string& read,
                  spp::sparse_hash_map<size_t, std::vector<util::MemCluster>>& memClusters,
                  uint32_t maxSpliceGap,
                  util::MateStatus mateStatus,
                  util::QueryCache& qc,
                  bool verbose=false) {
    // currently unused:
    // uint32_t readLen = static_cast<uint32_t>(read.length()) ;
   /*  if (verbose) {
      std::cerr << (mateStatus == util::MateStatus::PAIRED_END_RIGHT) << "\n";
    } */

    util::ProjectedHits phits;
    std::vector<std::pair<int, util::ProjectedHits>> rawHits;

    CanonicalKmer::k(k);
    pufferfish::CanonicalKmerIterator kit_end;
    pufferfish::CanonicalKmerIterator kit1(read);

    /**
     *  Testing heuristic.  If we just succesfully matched a k-mer, and extended it to a uni-MEM, then
     *  the nucleotide after the end of the match is very likely to be a sequencing error (or a variant).
     *  In this case, the next "k" k-mers will still likely overlap that error and not match.  So, instead
     *  of looking at each of them (or skipping all of them, which might decrease sensitivity), we will 
     *  skip (by a Mohsen number ;P) until we are at least k bases away from the end of the valid match.
     *  Then, so as to maintain high sensitivity, we will start skipping by only 1 k-mer at a time.
     **/              

    // Start off pretending we are at least k bases away from the last hit
    uint32_t skip{1};
    uint32_t altSkip{1};
    int32_t signedK = static_cast<int32_t>(k);
    int32_t basesSinceLastHit{signedK};
    ExpansionTerminationType et {ExpansionTerminationType::MISMATCH};

    while (kit1 != kit_end) {
      auto phits = pfi_->getRefPos(kit1->first, qc);
      skip = (basesSinceLastHit >= signedK) ? 1 : altSkip;
      if (!phits.empty()) {
        // kit1 gets updated inside expandHitEfficient function
        // stamping the reasPos
        size_t readPosOld = kit1->second ;
        /* if(verbose){
          std::cerr<< "Index "<< phits.contigID() << " ContigLen "<<phits.contigLen_<< " GlobalPos " << phits.globalPos_ << " ore " << phits.contigOrientation_ << " ref size " << phits.refRange.size() <<"\n\n\n" ;
          std::cerr<<kit1->first.to_str() << "\n" ;
          for(auto& posIt : phits.refRange){
            auto refPosOri = phits.decodeHit(posIt);
            std::cerr << posIt.transcript_id() << "\t" <<  refPosOri.isFW << "\t" << refPosOri.pos << "\n" ;
          } 
        } */
        expandHitEfficient(phits, kit1, et, verbose);
        //if(verbose) std::cerr<<"len after expansion: "<<phits.k_<<"\n" ;
        
        rawHits.push_back(std::make_pair(readPosOld, phits));
        basesSinceLastHit = 1;
        skip = (et == ExpansionTerminationType::MISMATCH) ? altSkip : 1;
        kit1 += (skip-1);
       //} else {
       //  ++kit1;
       //}
       //++pos;
      } else {
       // ++pos;
       basesSinceLastHit += skip;
       kit1 += skip;
       //++kit1;
      }
    }

    // if this is the right end of a paired-end read, use memCollectionRight,
    // otherwise (left end or single end) use memCollectionLeft.
    auto* memCollection = (mateStatus == util::MateStatus::PAIRED_END_RIGHT) ?
      &memCollectionRight : &memCollectionLeft;
    if (rawHits.size() > 0) {
      clusterMems(rawHits, memClusters, maxSpliceGap, *memCollection, verbose);
      if (verbose) {
        std::cerr << "lets see what we have\n";
        for (auto kv : memClusters) {
          
        std::cerr <<"tid:" <<  kv.first << "\n";
          for (auto clust : kv.second) {
            auto lclust = &clust;
          
          std::cerr << lclust->isFw << " size:" << lclust->mems.size() << " cov:" << lclust->coverage << "\n";
          for (size_t i = 0; i < lclust->mems.size(); i++){
              std::cerr << "--- t" << lclust->mems[i].tpos << " r"
              << lclust->mems[i].memInfo->rpos << " cid:"
              << lclust->mems[i].memInfo->cid << " cpos: "
              << lclust->mems[i].memInfo->cpos << " len:"
              << lclust->mems[i].memInfo->memlen << " fw:"
              << lclust->mems[i].memInfo->cIsFw << "\n";
          }
  }
        }
        std::cerr << "\n and now after recover gaps\n\n";
      
      }
      //recoverGaps(memClusters, *memCollection, read.length(), verbose );
      
      if (verbose) {
        for (auto kv : memClusters) {
          
        std::cerr <<"tid:" <<  kv.first << "\n";
          for (auto clust : kv.second) {
            auto lclust = &clust;
          
          std::cerr << lclust->isFw << " size:" << lclust->mems.size() << " cov:" << lclust->coverage << "\n";
          for (size_t i = 0; i < lclust->mems.size(); i++){
              std::cerr << "--- t" << lclust->mems[i].tpos << " r"
              << lclust->mems[i].memInfo->rpos << " cid:"
              << lclust->mems[i].memInfo->cid << " cpos: "
              << lclust->mems[i].memInfo->cpos << " len:"
              << lclust->mems[i].memInfo->memlen << " fw:"
              << lclust->mems[i].memInfo->cIsFw << "\n";
          }
  }
        }
      
      }
      
      return true;
    }
    return false;
  }

  void clear() {
    memCollectionLeft.clear();
    memCollectionRight.clear();
  }

private:
  PufferfishIndexT* pfi_;
  size_t k;
  //AlignerEngine ae_;
  std::vector<util::UniMemInfo> memCollectionLeft;
  std::vector<util::UniMemInfo> memCollectionRight;
};
#endif
