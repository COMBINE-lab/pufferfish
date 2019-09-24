#include "MemCollector.hpp"

// using spp:sparse_hash_map;

// polute the namespace --- put this in the functions that need it.
namespace kmers = combinelib::kmers;

template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::configureMemClusterer(uint32_t max) {
  mc.setMaxAllowedRefsPerHit(max);
}

template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::setHitFilterPolicy(pufferfish::util::HitFilterPolicy hfp) {
  mc.setHitFilterPolicy(hfp);
}

template <typename PufferfishIndexT>
pufferfish::util::HitFilterPolicy MemCollector<PufferfishIndexT>::getHitFilterPolicy() const {
  return mc.getHitFilterPolicy();
}


template <typename PufferfishIndexT>
size_t MemCollector<PufferfishIndexT>::expandHitEfficient(pufferfish::util::ProjectedHits& hit,
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
      uint64_t fk = allContigs.get_int(2*cCurrPos, 2*baseCnt);
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
      uint64_t fk = allContigs.get_int(2*(cCurrPos - baseCnt), 2*baseCnt);
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
void MemCollector<PufferfishIndexT>::setConsensusFraction(double cf) {
  mc.setConsensusFraction(cf);
}

template <typename PufferfishIndexT>
double MemCollector<PufferfishIndexT>::getConsensusFraction() const {
  return mc.getConsensusFraction();
}

template <typename PufferfishIndexT>
bool MemCollector<PufferfishIndexT>::operator()(std::string &read,
                  pufferfish::util::QueryCache& qc,
                  bool isLeft,
                  bool verbose) {

  // currently unused:
  // uint32_t readLen = static_cast<uint32_t>(read.length()) ;
  pufferfish::util::ProjectedHits phits;
  auto& rawHits = isLeft ? left_rawHits : right_rawHits;

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
  uint32_t altSkip{5};
  int32_t signedK = static_cast<int32_t>(k);
  int32_t basesSinceLastHit{signedK};
  ExpansionTerminationType et {ExpansionTerminationType::MISMATCH};

  while (kit1 != kit_end) {
    auto phits = pfi_->getRefPos(kit1->first, qc);
    skip = (basesSinceLastHit >= signedK) ? 1 : altSkip;
    if (!phits.empty()) {
      // kit1 gets updated inside expandHitEfficient function
      // stamping the readPos
      // NOTE: expandHitEfficient advances kit1 by *at least* 1 base
      size_t readPosOld = kit1->second;
      expandHitEfficient(phits, kit1, et);
			if (verbose){
			  std::cerr<<"after expansion\n";
				std::cerr<<"readPosOld:"<<readPosOld<<" kmer:"<< kit1->first.to_str() <<"\n";
			}
      rawHits.push_back(std::make_pair(readPosOld, phits));
    
      basesSinceLastHit = 1;
      skip = (et == ExpansionTerminationType::MISMATCH) ? altSkip : 1;
      kit1 += (skip-1);
    } else {
     basesSinceLastHit += skip;
     kit1 += skip;
    }
  }

  // To consider references this end maps to for allowing hits on the other end
  /*auto& ref_ids = isLeft ? left_refs : right_refs;
  for (auto &hit : core::range<decltype(rawHits.begin())>(rawHits.begin(), rawHits.end())) {
    auto &projHits = hit.second;
    auto &refs = projHits.refRange;
    if (refs.size() < mc.getMaxAllowedRefsPerHit()) {
      for (auto &posIt : refs) {
        auto refid = posIt.transcript_id();
        ref_ids[refid] = true;
      }
    }
  }*/
  return rawHits.size() != 0;
}

template <typename PufferfishIndexT>
bool MemCollector<PufferfishIndexT>::findChains(std::string &read,
                                                pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& memClusters,
                                                //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>>& memClusters,
                  uint32_t maxSpliceGap,
                  pufferfish::util::MateStatus mateStatus,
                  bool hChain,
                  bool isLeft,
                  bool verbose) {
  (void) maxSpliceGap;
  uint64_t firstDecoyIndex = pfi_->firstDecoyIndex();
  auto& rawHits = isLeft ? left_rawHits : right_rawHits;
  // if this is the right end of a paired-end read, use memCollectionRight,
  // otherwise (left end or single end) use memCollectionLeft.
  isSingleEnd = mateStatus == pufferfish::util::MateStatus::SINGLE_END;
  auto* memCollection = (mateStatus == pufferfish::util::MateStatus::PAIRED_END_RIGHT) ?
  &memCollectionRight : &memCollectionLeft;
  if (rawHits.size() > 0) {
    auto& other_end_refs = isLeft ? right_refs : left_refs;
    trMemMap.clear();
    mc.findOptChain(rawHits, memClusters, maxSpliceGap, *memCollection, read.length(), other_end_refs, hChain, trMemMap, firstDecoyIndex, verbose);
    /*
    if (verbose) {
      std::cerr << "lets see what we have\n";
      for (auto kv : trMemMap) {
        std::cerr <<"tid:" <<  kv.first.first << "\n";
        for (auto mem : kv.second) {
          std::cerr << "--- t" << mem.tpos << " r"
              << mem.rpos << " cid:"
              << mem.memInfo->cid << " cpos: "
              << mem.memInfo->cpos << " len:"
              << mem.extendedlen << " re:"
              << mem.memInfo->readEnd << " fw:"
              << mem.isFw << "\n";
          std::cerr << read.substr(mem.rpos,mem.extendedlen) << "\n";
        }
      }
    }
    */
    return true;
  }
  return false;
}


template <typename PufferfishIndexT>
void MemCollector<PufferfishIndexT>::clear() {
  memCollectionLeft.clear();
  memCollectionRight.clear();
  trMemMap.clear();
  left_refs.clear();
  right_refs.clear();
  left_rawHits.clear();
  right_rawHits.clear();
}

template class MemCollector<PufferfishIndex>;

template class MemCollector<PufferfishSparseIndex>;

template class MemCollector<PufferfishLossyIndex>;
