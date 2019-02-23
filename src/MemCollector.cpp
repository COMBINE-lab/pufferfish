#include "MemCollector.hpp"

// using spp:sparse_hash_map;

// polute the namespace --- put this in the functions that need it.
namespace kmers = combinelib::kmers;

template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::recoverGaps(spp::sparse_hash_map<pufferfish::common_types::ReferenceID, 
                                        std::vector<util::MemCluster>>& memClustersMap,
                                        std::vector<util::UniMemInfo>& memCollection,
                                        size_t rlen, 
                                        bool verbose) {
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
                tpos, (size_t) 0);
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
                tpos, (size_t) 0);
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


template <typename PufferfishIndexT>
size_t MemCollector<PufferfishIndexT>::expandHitEfficient(util::ProjectedHits& hit, 
                       pufferfish::CanonicalKmerIterator& kit, 
                       ExpansionTerminationType& et) {

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
     readSeqOffset < readSeqLen) { 

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
  hit.contigPos_ -= (hit.k_ - k);
  hit.globalPos_ -= (hit.k_ - k);
  }
  kit.jumpTo(readSeqStart);
  return currReadStart;
}


template <typename PufferfishIndexT>
bool MemCollector<PufferfishIndexT>::operator()(std::string &read,
                  //spp::sparse_hash_map<size_t, std::vector<util::MemCluster>>& memClusters,
                  uint32_t maxSpliceGap,
                  util::MateStatus mateStatus,
                  util::QueryCache& qc,
                  bool verbose) {

  // currently unused:
  // uint32_t readLen = static_cast<uint32_t>(read.length()) ;
 /*  if (verbose) {
  std::cerr << (mateStatus == util::MateStatus::PAIRED_END_RIGHT) << "\n";
  } */
  (void) maxSpliceGap;
  util::ProjectedHits phits;
  std::vector<std::pair<int, util::ProjectedHits>> rawHits;

  CanonicalKmer::k(k);
  pufferfish::CanonicalKmerIterator kit_end;
  pufferfish::CanonicalKmerIterator kit1(read);
  if (verbose) {
  std::cerr << "ORIGINAL READ:\n";
  std::cerr << read << "\n";
  }

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
    // NOTE: expandHitEfficient advances kit1 by *at least* 1 base
    size_t readPosOld = kit1->second ;
    expandHitEfficient(phits, kit1, et);

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
  isSingleEnd = mateStatus == util::MateStatus::SINGLE_END;
  auto* memCollection = (mateStatus == util::MateStatus::PAIRED_END_RIGHT) ?
  &memCollectionRight : &memCollectionLeft;
  if (rawHits.size() > 0) {
  if (mateStatus == util::MateStatus::PAIRED_END_LEFT)
    mc.fillMemCollection(rawHits, trMemMap, *memCollection, util::ReadEnd::LEFT, verbose);
  else
    mc.fillMemCollection(rawHits, trMemMap, *memCollection, util::ReadEnd::RIGHT, verbose);

  //mc.findOptChain(rawHits, memClusters, maxSpliceGap, *memCollection, verbose);
  //mc.clusterMems(rawHits, memClusters, maxSpliceGap, *memCollection, verbose);
  if (verbose) {
    std::cerr << "lets see what we have\n";
    for (auto kv : trMemMap) {
    
    std::cerr <<"tid:" <<  kv.first.first << "\n";
    for (auto mem : kv.second) {
      //auto lclust = &clust;

      //std::cerr << lclust->isFw << " size:" << lclust->mems.size() << " cov:" << lclust->coverage << "\n";
      //for (size_t i = 0; i < lclust->mems.size(); i++) {
      std::cerr << "--- t" << mem.tpos << " r"
            << mem.memInfo->rpos << " cid:"
            << mem.memInfo->cid << " cpos: "
            << mem.memInfo->cpos << " len:"
            << mem.memInfo->memlen << " re:"
            << mem.memInfo->readEnd << " fw:"
            << mem.isFw << "\n";
      std::cerr << read.substr(mem.memInfo->rpos,mem.memInfo->memlen) << "\n";
      //}
    }
    }
  }
  /*recoverGaps(memClusters, *memCollection, read.length(), verbose );*/
  
  return true;
  }
  return false;
}

  

template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::findBestChain(std::vector<util::JointMems> &jointHits,
                       std::vector<util::MemCluster> &all,
                       uint32_t maxSpliceGap, uint32_t maxFragmentLength, bool verbose) {
      mc.findOptChain(trMemMap, jointHits, all, maxSpliceGap, maxFragmentLength, verbose);
}

template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::findOptChainAllowingOneJumpBetweenTheReadEnds(std::vector<util::JointMems> &jointHits,
                                                       std::vector<util::MemCluster> &all,
                                                       uint32_t maxSpliceGap, uint32_t maxFragmentLength,
                                                       uint32_t readLenLeft,
                                                       uint32_t readLenRight,
                                                       bool verbose) {
  mc.findOptChainAllowingOneJumpBetweenTheReadEnds(trMemMap, jointHits, all, maxSpliceGap, maxFragmentLength, pfi_,
              readLenLeft, readLenRight, verbose);
}

  

template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::clear() {
  memCollectionLeft.clear();
  memCollectionRight.clear();
  trMemMap.clear();
}

template class MemCollector<PufferfishIndex>;

template class MemCollector<PufferfishSparseIndex>;

template class MemCollector<PufferfishLossyIndex>;
