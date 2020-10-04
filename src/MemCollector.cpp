#include "MemCollector.hpp"
#include <cmath>
#include <limits>

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
void MemCollector<PufferfishIndexT>::setAltSkip(uint32_t as) {
  altSkip = as;
}

  
enum class LastSkipType : uint8_t { NO_HIT=0, SKIP_READ=1, SKIP_UNI=2 };

// Idea is to move the logic of the search into here.
// We should also consider the optimizations that can 
// be done here (like having small checks expected to)
// be on the same contig bypass a hash lookup.
struct SkipContext {

  SkipContext(std::string& read, int32_t k_in) : 
    kit1(read), kit_tmp(read), read_len(static_cast<int32_t>(read.length())),
    expected_cid(invalid_cid), safe_skip(1), k(k_in), 
    last_skip_type(LastSkipType::NO_HIT)
  { }
  
  inline bool is_exhausted() {
    return kit1 == kit_end;
  }

  inline CanonicalKmer& curr_kmer() {
    return kit1->first;
  }

  template <typename PufferfishIndexT>
  inline bool query_kmer(PufferfishIndexT* pfi, pufferfish::util::QueryCache& qc) {
    phits = pfi->getRefPos(kit1->first, qc);
    return !phits.empty();
  }

  inline bool hit_is_unexpected() {
    return (expected_cid != invalid_cid and phits.contigIdx_ != expected_cid);
  }

  inline int32_t read_pos() { return kit1->second; }

  inline pufferfish::util::ProjectedHits proj_hits() { return phits; }

  inline void clear_expectation() { expected_cid = invalid_cid; }

  inline int32_t compute_safe_skip() {
    // if we are calling this, then we *did* get 
    // an unexpected hit.  To avoid iterating to that
    // already queried hit again, the target pos
    // will be the position before that hit
    read_target_pos = kit1->second-1;
    read_prev_pos = kit_tmp->second;
    int32_t dist = read_target_pos - read_prev_pos;
    
    //if (dist < 0) { std::cerr << "dist = " << dist << " should not happen!\n"; }
    
    safe_skip = 1;
    // if the distance is large enough
    if (dist > 2) {
      safe_skip = (read_target_pos - read_prev_pos) / 2;
      safe_skip = static_cast<uint32_t>(safe_skip < 1 ? 1 : safe_skip);
    }
    return dist;
  }

  inline void rewind() {
    // start from the next position on the read and 
    // walk until we hit the read end or the position 
    // that we wanted to skip to
    if (kit_tmp->second >= kit1->second) {
      std::cerr << "rewinding to a greater position; shouldn't happen!\n";
    }
    kit_swap = kit1;
    kit1 = kit_tmp;
    kit_tmp = kit_swap;
  }

  inline void fast_forward() {
    // jump back to the place we were 
    // at when we called rewind()

    // NOTE: normally we would assume that 
    // kit1->second >= kit_tmp->second
    // but, this be violated if we skipped farther than we asked in the iterator
    // because there were 'N's in the read.
    // in that case case, ignore the fact that fast forward 
    // could be a rewind.  We always want to skip back to 
    // where we were before.  However, in this case we 
    // ASSUME that no hit *after* kit_tmp has been added
    // to the set of collected hits.
    kit1 = kit_tmp;
  }

  /**
   * Advance by the safe skip amount calculated in 
   * `compute_safe_skip`.  
   * Returns true if we are still before or at the expected
   * end point, and false otherwise.
   */
  inline bool advance_safe() {
    int32_t remaining_len = read_target_pos - kit1->second;
    safe_skip = (remaining_len <= safe_skip) ? remaining_len : safe_skip;
    safe_skip = (safe_skip < 1) ? 1 : safe_skip;
    kit1 += safe_skip;
    return kit1->second <= read_target_pos;
  }
        
  inline void advance_from_hit() {
      int32_t skip = 1;
      // the offset of the hit on the read
      int32_t read_offset = kit1->second;
      // the skip that would take us to the last base of the read
      int32_t read_skip = (read_len - read_offset + k - 1);
      // any valid skip should always be at least 1 base
      read_skip = (read_skip < 1) ? 1 : read_skip;
      
      size_t cStartPos = phits.globalPos_ - phits.contigPos_; 
      size_t cEndPos = cStartPos + phits.contigLen_;
      size_t cCurrPos = phits.globalPos_; 
      int32_t ctg_skip = 1;
      // fw ori
      if (phits.contigOrientation_) {
        ctg_skip = static_cast<int64_t>(cEndPos) - (static_cast<int64_t>(cCurrPos + k));
        if (ctg_skip < 0) { 
            std::cerr << "cStartPos =  " << cStartPos << "\n";
            std::cerr << "cEndPos = " << cEndPos << "\n";
            std::cerr << "cCurrPos = " << cCurrPos << "\n";
            std::cerr << "contig skip = " << ctg_skip << " < 0, should not happen!\n";
        }
      } else { // rc ori
        ctg_skip = static_cast<int32_t>(phits.contigPos_);
      }
      // we're already at the end of the contig
      bool at_contig_end = (ctg_skip == 0);
      // if we're at the end of the contig
      // we'll set the contig skip to be one 
      // (i.e. look for the next k-mer), but we won't
      // set an expectation on what that contig should be
      if (at_contig_end) {
        ctg_skip = 1;
        expected_cid = invalid_cid;
      } else { // otherwise, we expect the next contig to be the same
        expected_cid = phits.contigIdx_;
      }

      // check here if we are going to attempt to skip to the 
      // end of the read, and if this skip is non trivial (>1).
      // if that doesn't produce a hit, we'll want to "rewind"
      // to the current position and move forward by alt.
      // nontrivial_last_skip = (skip > ulast_skip) and ulast_skip > 1;
      
      // remember where we are coming from
      kit_tmp = kit1;
      
      // remember the reason we are skipping
      last_skip_type = (read_skip <= ctg_skip) ? LastSkipType::SKIP_READ : LastSkipType::SKIP_UNI;
      
      // skip will be min of read and contig
      skip = (last_skip_type == LastSkipType::SKIP_READ) ? read_skip : ctg_skip;

      // skip must be at least 1
      skip = skip < 1 ? 1 : skip;
      
      // onward
      kit1 += skip;
  }

  inline void advance_from_miss() {
      int32_t skip = 1;
      read_target_pos = kit1->second;

      // distance from backup position 
      int32_t dist = read_target_pos - kit_tmp->second;

      switch (last_skip_type) {
        // we could have not yet seen a hit
        // should move 1 at a time
        case LastSkipType::NO_HIT : {
          kit1 += 2;
          return;
        }
        break;

        // we could have seen a hit, and tried to jump to the end of the read
        // and the previous search was either that hit, or a miss
        case LastSkipType::SKIP_READ  :
        case LastSkipType::SKIP_UNI : {
          // if distance from backup is 1, that means we already tried
          // last position, so we should just move to the next position 
          if (dist <= 2) {
            kit1 += 1;
            return;
          } else {
            // otherwise move the backup position toward us
            // and move the current point to the backup
            skip = dist / 2;
            skip = (skip < 2) ? 2 : skip;
            kit_tmp += skip;
            kit1 = kit_tmp;
            return;
          }
        }
        break;
        
        // we could have seen a hit, and tried to jump to the end of the contig
        // and the previous search was either that hit, or a miss
        //case LastSkipType::SKIP_UNI : {
        //
        //
        // }
      }
      // we could have previously not seen a hit either
      // that doesn't matter since it is recursively one of the above cases.
  }

  pufferfish::CanonicalKmerIterator kit1;
  pufferfish::CanonicalKmerIterator kit_tmp;
  pufferfish::CanonicalKmerIterator kit_end;
  pufferfish::CanonicalKmerIterator kit_swap;
  int32_t read_len;
  int32_t read_target_pos;
  int32_t read_current_pos;
  int32_t read_prev_pos;
  int32_t safe_skip;
  int32_t k;
  static constexpr uint32_t invalid_cid{std::numeric_limits<uint32_t>::max()};
  uint32_t expected_cid;//  = invalid_cid;
  
  LastSkipType last_skip_type{LastSkipType::NO_HIT};

  pufferfish::util::ProjectedHits phits;
};

template <typename PufferfishIndexT>
bool MemCollector<PufferfishIndexT>::get_raw_hits_sketch(std::string &read,
                  pufferfish::util::QueryCache& qc,
                  bool isLeft,
                  bool verbose) {
  (void) verbose;
  pufferfish::util::ProjectedHits phits;
  auto& raw_hits = isLeft ? left_rawHits : right_rawHits;

  CanonicalKmer::k(k);
  int32_t k = static_cast<int32_t>(CanonicalKmer::k());
  SkipContext skip_ctx(read, k);
  
  // while it is possible to search further
  while (!skip_ctx.is_exhausted()) {
    
    // if we had a hit
    if (skip_ctx.query_kmer(pfi_, qc)) {

      // record this hit
      bool added_orig_hit = false;
      int32_t read_pos = skip_ctx.read_pos();
      auto proj_hits = skip_ctx.proj_hits();
      
      // if the hit was not inline with what we were 
      // expecting. 
      if (skip_ctx.hit_is_unexpected()){
        int32_t dist_to_target = skip_ctx.compute_safe_skip();

        // special case is the prev hit was the 
        // one right before this, in this case
        // no need to do this work, we already 
        // have both hits.  In that case, just
        // proceed as usual.
        if (dist_to_target > 0) {
            // start from the next position on the read and
            // walk until we hit the read end or the position
            // that we wanted to skip to
            skip_ctx.rewind();
            while (skip_ctx.advance_safe()) {
              if (skip_ctx.query_kmer(pfi_, qc)) {
                raw_hits.push_back(
                    std::make_pair(skip_ctx.read_pos(), skip_ctx.proj_hits()));
              }
            }

            // if we got here, then either we are done, or we 
            // reached our target position so we don't
            // have an expectation of what should come next.
            skip_ctx.clear_expectation();
            skip_ctx.fast_forward();
            // so that the next skip can be computed properly.
            skip_ctx.phits = proj_hits;
          }
      } 
      // push this hit and advance
      raw_hits.push_back(std::make_pair(read_pos, proj_hits));
      skip_ctx.advance_from_hit();
    } else {
      // if we got here, we looked for a match and didn't find one
      skip_ctx.advance_from_miss();
    }
  }
  
  return raw_hits.size() != 0;
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
  uint32_t homoPolymerSkip{1};
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

      if (kit1->first.is_homopolymer()) {
        kit1 += homoPolymerSkip;
        continue;
      }

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
void MemCollector<PufferfishIndexT>::setChainSubOptThresh(double pre_merge_chain_sub_thresh) {
  pre_merge_chain_sub_thresh_ = pre_merge_chain_sub_thresh;
  inv_pre_merge_chain_sub_thresh_ = std::nexttoward(1.0 / pre_merge_chain_sub_thresh, std::numeric_limits<long double>::infinity());
  mc.set_chain_sub_opt_thresh_(pre_merge_chain_sub_thresh_, inv_pre_merge_chain_sub_thresh_);
}


template <typename PufferfishIndexT>
double MemCollector<PufferfishIndexT>::chainSubOptThresh() const {
  return pre_merge_chain_sub_thresh_;
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
    mc.findOptChain(rawHits, memClusters, maxSpliceGap, *memCollection, read.length(), 
                    other_end_refs, hChain, trMemMap, firstDecoyIndex, verbose);
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
