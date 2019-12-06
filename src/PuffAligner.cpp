#include "nonstd/string_view.hpp"
#include "PuffAligner.hpp"
#include "Util.hpp"

std::string extractReadSeq(const std::string& readSeq, uint32_t rstart, uint32_t rend, bool isFw) {
    std::string subseq = readSeq.substr(rstart, rend - rstart);
    if (isFw)
        return subseq;
    return pufferfish::util::reverseComplement(subseq); //reverse-complement the substring
}

// decode CIGAR, just like : https://github.com/lh3/ksw2/blob/master/cli.c#L134
std::string cigar2str(const ksw_extz_t *ez) {
    std::string cigar;
    if (ez->n_cigar > 0) {
        cigar.resize(ez->n_cigar * 2);
        for (int i = 0; i < ez->n_cigar; ++i) {
            cigar += (std::to_string(ez->cigar[i] >> 4) + "MID"[ez->cigar[i] & 0xf]);
        }
    }
    return cigar;
}

int32_t addCigar(pufferfish::util::CIGARGenerator &cigarGen, const ksw_extz_t *ez, bool start_gap) {
  uint32_t gap_size = 0;
  for (int i = 0; i < ez->n_cigar; ++i) {
    char cigar_type = "MID"[ez->cigar[i] & 0xf];
    cigarGen.add_item(ez->cigar[i] >> 4, cigar_type);
    if (start_gap and (cigar_type == 'M' or cigar_type == 'D'))
      gap_size += ez->cigar[i] >> 4;
  }
  return gap_size;
}

int32_t PuffAligner::align_ungapped(const char* const query, const char* const target, int len) {
  int score = 0;
  bool computeCIGAR = !(aligner.config().flag & KSW_EZ_SCORE_ONLY);
  auto& cigarGen = cigarGen_;
  for (uint32_t i = 0; i < len; i++) {
    if (computeCIGAR) cigarGen.add_item(1, 'M');
    score += query[i] == target[i] ? mopts.matchScore : mopts.missMatchPenalty;
  }
  return score;
}

/*std::string getRefSeq(compact::vector<uint64_t, 2> &refseq, uint64_t refAccPos, size_t tpos, uint32_t memlen) {
    if (memlen == 0) return "";
    std::string tseq; tseq.reserve(memlen);
    uint64_t bucket_offset = (refAccPos + tpos) * 2;
    auto len_on_vector = memlen * 2;
    for (uint32_t w = 0; w <= len_on_vector / 64; w++) {
        uint32_t len = std::min((uint32_t) 64, len_on_vector - w * 64);
        if (len == 0) continue;
        uint64_t word = refseq.get_int(bucket_offset, len);
        for (uint32_t i = 0; i < len; i += 2) {
            uint8_t next_bits = ((word >> i) & 0x03);
            char next = 'A';
            if (next_bits == 1)
                next = 'C';
            else if (next_bits == 2)
                next = 'G';
            else if (next_bits == 3)
                next = 'T';
            tseq += next;
        }
        bucket_offset += len;
    }
    return tseq;
}*/

bool fillRefSeqBuffer(compact::vector<uint64_t, 2> &refseq, uint64_t refAccPos, size_t tpos, uint32_t memlen, std::string& refBuffer_) {
  refBuffer_.clear();
  if (memlen == 0) return false;
  uint64_t bucket_offset = (refAccPos + tpos) * 2;
  auto len_on_vector = memlen * 2;
  int32_t toFetch = len_on_vector;
  while (toFetch > 0) {
    uint32_t len = (toFetch >= 64) ? 64 : toFetch;
    toFetch -= len;
    uint64_t word = refseq.get_int(bucket_offset, len);
    for (uint32_t i = 0; i < len; i += 2) {
      uint8_t next_bits = ((word >> i) & 0x03);
      refBuffer_ += "ACGT"[next_bits];
    }
    bucket_offset += len;
  }
  return true;
}

// NOTE: This fills in refBuffer_ with the reference sequence from tpos to tpos+memlen in *reverse* order.
// refBuffer will contain the reverse of the reference substring, *NOT* the reverse-complement.
bool fillRefSeqBufferReverse(compact::vector<uint64_t, 2> &refseq, uint64_t refAccPos, size_t tpos, uint32_t memlen, std::string& refBuffer_) {
  refBuffer_.resize(memlen, 'N');
  if (memlen == 0) return false;
  uint64_t bucket_offset = (refAccPos + tpos) * 2;
  auto len_on_vector = memlen * 2;
  int32_t toFetch = len_on_vector;
  int32_t offset = static_cast<int32_t>(memlen)-1;
  while (toFetch > 0) {
    uint32_t len = (toFetch >= 64) ? 64 : toFetch;
    toFetch -= len;
    uint64_t word = refseq.get_int(bucket_offset, len);
    for (uint32_t i = 0; i < len; i += 2) {
      uint8_t next_bits = ((word >> i) & 0x03);
      refBuffer_[offset] = "ACGT"[next_bits];
      --offset;
    }
    bucket_offset += len;
  }
  return true;
}

/**
 *  Align the read `original_read`, whose mems consist of `mems` against the index and return the result
 *  in `arOut`.  How the alignment is computed (i.e. full vs between-mem and CIGAR vs. score only) depends
 *  on the parameters of how this PuffAligner object was constructed.
 **/
bool PuffAligner::alignRead(std::string& read, std::string& read_rc, const std::vector<pufferfish::util::MemInfo> &mems, bool perfectChain,
                            bool isFw, size_t tid, AlnCacheMap &alnCache, HitCounters &hctr, AlignmentResult& arOut, bool isLeft) {
  int32_t alignmentScore{std::numeric_limits<decltype(arOut.score)>::min()};
  if (mems.empty()) {
    arOut.score = alignmentScore;
    return false;
  }

  /*auto logger_ = spdlog::get("console");
  spdlog::set_level(spdlog::level::debug); // Set global log level to debug
  logger_->set_pattern("%v");*/

  int32_t refExtLength = static_cast<int32_t>(mopts.refExtendLength);
  bool firstMem = true;
  int32_t lastHitEnd_read = -1;
  int32_t currHitStart_read = 0;
  int64_t lastHitEnd_ref = -1;
  int64_t currHitStart_ref = 0;
  int32_t alignment{0};

  uint32_t openGapLen{0};

  // will we be computing CIGAR strings for this alignment
  bool computeCIGAR = !(aligner.config().flag & KSW_EZ_SCORE_ONLY);
  // the cigar generator we will use
  auto& cigarGen = cigarGen_;
  cigarGen.clear();

  std::string cigar = "";

  // where this reference starts, and its length.
  int64_t refAccPos = tid > 0 ? refAccumLengths[tid - 1] : 0;
  int64_t refTotalLength = refAccumLengths[tid] - refAccPos;

  auto& frontMem = mems.front();
  auto rpos = frontMem.rpos;
  auto memlen = frontMem.extendedlen;
  auto readLen = read.length();
  auto tpos = frontMem.tpos;

  // do full alignment if we are in that mode, or if the
  // current read was recovered via orphan recovery.
  // @mohsen & @fataltes : we need a better signal than memlen == 1
  // to designate this was a recovered orphan.
  bool recoveredOrphan = memlen == 1;
  bool doFullAlignment = mopts.fullAlignment or recoveredOrphan;
  currHitStart_read = isFw ? rpos : readLen - (rpos + memlen);

  if (currHitStart_read < 0 or currHitStart_read >= (int32_t) readLen) {
    std::cerr << "[ERROR in PuffAligner::alignRead :] currHitStart_read is invalid; this should not happen!\n";
    return false;
  }

  currHitStart_ref = tpos;

  // the length of the reference string we will use to check the
  // alignment cache (or to compute a full alignment).
  int32_t keyLen = 0;

  std::string tseq;

  uint64_t hashKey{0};
  bool didHash{false};
  uint32_t refStart, readStart{0};
  //We want the maximum length of the buffer that if we include those number of indels we can still achieve a good quality alignment
  //For the read length of readLen, the following formula will create a buffer of size 21
  uint32_t buff = (readLen*mopts.matchScore*(1-mopts.minScoreFraction)-mopts.gapOpenPenalty)/mopts.gapExtendPenalty;

  int32_t signedRefStartPos = currHitStart_ref - currHitStart_read;
  int32_t signedRefEndPos = signedRefStartPos + read.length();
  bool invalidStart = (signedRefStartPos < 0);
  bool invalidEnd = (signedRefEndPos > refTotalLength);

  bool allowOverhangSoftclip = mopts.allowOverhangSoftclip;

  if (mopts.mimicBT2Strict and (invalidStart or invalidEnd)) {
    arOut.score = std::numeric_limits<decltype(arOut.score)>::min();
    arOut.cigar = "";
    arOut.openGapLen = 0;
    return false;
  }

  ksw_reset_extz(&ez);
  // If we are only aligning between MEMs
  if (!doFullAlignment) {
    refStart = (currHitStart_ref >= currHitStart_read) ? currHitStart_ref - currHitStart_read : 0;
    //keyLen = (refStart + readLen < refTotalLength) ? readLen : refTotalLength - refStart;
    keyLen = (refStart + readLen + buff < refTotalLength) ? readLen + buff : refTotalLength - refStart;
  } else { // we are aligning from the start of the read
    // If the first hit starts further in the reference than in the
    // read, then we align from the beginning of the read and (ref_start - read_start) on
    // the reference.
    if (currHitStart_ref >= currHitStart_read) {
      refStart = currHitStart_ref - currHitStart_read;
      readStart = 0;
    } else if (currHitStart_ref < currHitStart_read) {
      // If the first hit starts further in the read than in the reference, then
      // the read overhangs the beginning of the reference and we start aligning
      // from the beginning of the reference and from position (read_start - ref_start)
      // in the read.
      readStart = currHitStart_read - currHitStart_ref;
      refStart = 0;
    }
    keyLen = (refStart + readLen + buff < refTotalLength) ? readLen + buff : refTotalLength - refStart;
  }

  if (perfectChain) {
    arOut.score = alignmentScore = readLen * mopts.matchScore;
    if (computeCIGAR) { cigarGen.add_item(readLen, 'M'); }
    hctr.skippedAlignments_byCov += 1;
    SPDLOG_DEBUG(logger_,"[[");
    //SPDLOG_DEBUG(logger_,"read sequence ({}) : {}", (isFw ? "FW" : "RC"), readView);
    SPDLOG_DEBUG(logger_,"ref  sequence      : {}", (doFullAlignment ? tseq : refSeqBuffer_));
    SPDLOG_DEBUG(logger_,"perfect chain!\n]]\n");
    arOut.cigar = cigar;
    arOut.openGapLen = openGapLen;
    return true;
  }

  fillRefSeqBuffer(allRefSeq, refAccPos, refStart, keyLen, refSeqBuffer_);
  int32_t originalRefSeqLen = static_cast<int32_t>(refSeqBuffer_.length());
  // If we're not using fullAlignment, we'll need the full reference sequence later
  // so copy it into tseq.
  if (!doFullAlignment) { tseq = refSeqBuffer_; }

  if (!alnCache.empty() and isMultimapping_) {
    // hash the reference string
    MetroHash64::Hash(reinterpret_cast<uint8_t *>(const_cast<char*>(refSeqBuffer_.data())), keyLen, reinterpret_cast<uint8_t *>(&hashKey), 0);
    didHash = true;
    // see if we have this hash
    auto hit = alnCache.find(hashKey);
    // if so, we know the alignment score
    if (hit != alnCache.end() ) {
      hctr.skippedAlignments_byCache += 1;
      arOut.score = hit->second.score;
      if (computeCIGAR) { arOut.cigar = hit->second.cigar; }
      arOut.openGapLen = hit->second.openGapLen;
      return true;
    }
  }

  // @mohsen & @fataltes --- we should figure out how to
  // avoid computing the rc of a read if we've already done it.
  if (!isFw and read_rc.empty()) { read_rc = pufferfish::util::reverseComplement(read); }
  nonstd::string_view readView = (isFw) ? read : read_rc;

  int32_t minLengthGapRequired = ( 2 * (mopts.gapOpenPenalty + mopts.gapExtendPenalty) + mopts.matchScore ) / (mopts.matchScore - mopts.missMatchPenalty);
  int32_t minAcceptedScore = scoreStatus_.getCutoff(read.length());
  auto alignable = [&minAcceptedScore, &readLen] (uint32_t alignedLen, uint32_t matchScore, int32_t alignmentScore) -> bool {
    return minAcceptedScore <= alignmentScore + matchScore*(readLen-alignedLen);
  };

  if (doFullAlignment) {
    // if we allow softclipping of overhanging bases, then we can cut off the part of the read
    // before the start of the reference
    decltype(readStart) readOffset = allowOverhangSoftclip ? readStart : 0;
    if (computeCIGAR)
      cigarGen.begin_softClipped_len = readOffset;
    nonstd::string_view readSeq = readView.substr(readOffset);
    hctr.aligner_calls_count += 1;
    ksw_reset_extz(&ez);
    aligner(readSeq.data(), readSeq.length(), refSeqBuffer_.data(), refSeqBuffer_.length(), &ez,
            allowOverhangSoftclip, ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
    // if we allow softclipping of overhaning bases, then we only care about the best
    // score to the end of the query or the end of the reference.  Otherwise, we care about
    // the best score all the way until the end of the query.
    alignmentScore = allowOverhangSoftclip ? std::max(ez.mqe, ez.mte) : ez.mqe;
    SPDLOG_DEBUG(logger_,"readSeq ({}) : {}\nrefSeq  : {}\nscore   : {}\nreadStart : {}", (isFw ? "FW" : "RC"), readSeq, refSeqBuffer_, alignmentScore, readStart);
    SPDLOG_DEBUG(logger_,"currHitStart_read : {}, currHitStart_ref : {}\nmqe : {}, mte : {}\n", currHitStart_read, currHitStart_ref, ez.mqe, ez.mte);
    // can make start pos negative, but sam writer deals with this right now
    if (computeCIGAR) {
      addCigar(cigarGen, &ez, false);
      if (allowOverhangSoftclip and ez.mte > ez.mqe) {
        cigarGen.end_softClipped_len = readSeq.length() - ez.max_q - 1;
      }
    }
    SPDLOG_DEBUG(logger_,"The current cigar is: {}", cigarGen.get_cigar());
  } else {
    // ROB ATTEMPT ---
    alignmentScore = 0;

    //std::stringstream ss;
    SPDLOG_DEBUG(logger_,"[[");
    SPDLOG_DEBUG(logger_,"read sequence ({}) : {}", (isFw ? "FW" : "RC"), readView);
    SPDLOG_DEBUG(logger_,"ref  sequence      : {}\nrefID : {}", tseq, tid);

    // If the first mem does not start at the beginning of the
    // read, then there is a gap to align.
    int32_t firstMemStart_read = (isFw) ? rpos : readLen - (rpos + memlen);
    if (firstMemStart_read > 0) {
      // align the part before the first mem

      // the gap is of length firstMemStart_read, so grab that much (plus buffer) before the
      // first start position on the reference.
      int32_t firstMemStart_ref = tpos;
      int32_t readStartPosOnRef = tpos - firstMemStart_read;
      int32_t dataDepBuff = std::min(refExtLength, 5*firstMemStart_read);

      //int32_t refWindowStart = (readStartPosOnRef - refExtLength) > 0 ? (readStartPosOnRef - refExtLength) : 0;
      int32_t refWindowStart = (readStartPosOnRef - dataDepBuff) > 0 ? (readStartPosOnRef - dataDepBuff) : 0;
      int32_t refWindowLength = tpos - refWindowStart;
      fillRefSeqBufferReverse(allRefSeq, refAccPos, refWindowStart, refWindowLength, refSeqBuffer_);

      if (refSeqBuffer_.length() > 0) {
        auto readWindow = readView.substr(0, firstMemStart_read).to_string();
        std::reverse(readWindow.begin(), readWindow.end());
        SPDLOG_DEBUG(logger_,"PRE:\nreadStartPosOnRef : {}\nrefWindowStart : {}", readStartPosOnRef, refWindowStart);
        SPDLOG_DEBUG(logger_,"refWindowLength : {}\nread : [{}]\nref : [{}]", refWindowLength, readWindow, refSeqBuffer_);
        ksw_reset_extz(&ez);
        hctr.aligner_calls_count += 1;
        aligner(readWindow.data(), readWindow.length(), refSeqBuffer_.data(), refSeqBuffer_.length(), &ez,
                allowOverhangSoftclip, ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
        alignmentScore += allowOverhangSoftclip ? std::max(ez.mqe, ez.mte) : ez.mqe;
        if (ez.mte > ez.mqe and computeCIGAR and allowOverhangSoftclip) {
            cigarGen.begin_softClipped_len = readWindow.length() - ez.max_q - 1;
        }
        openGapLen = computeCIGAR ? addCigar(cigarGen, &ez, true) : (ez.mqe_t + 1);

        SPDLOG_DEBUG(logger_,"The current cigar is: {}", cigarGen.get_cigar());
        SPDLOG_DEBUG(logger_,"score : {}", std::max(ez.mqe, ez.mte));
      } else {
        // do any special soft-clipping penalty here if we want
        if (computeCIGAR and allowOverhangSoftclip)
          cigarGen.begin_softClipped_len = firstMemStart_read;
        else if (computeCIGAR)
          cigarGen.add_item(firstMemStart_read, 'I');
        alignmentScore += allowOverhangSoftclip ? 0 :
          (-1 * mopts.gapOpenPenalty + -1 * mopts.gapExtendPenalty * firstMemStart_read);
        openGapLen = firstMemStart_read;
      }
    }

    int32_t prevMemEnd_read = firstMemStart_read - 1; //isFw ? rpos : readLen - (rpos+memlen);
    int32_t prevMemEnd_ref = tpos;
    if (!alignable(prevMemEnd_read + 1, mopts.matchScore, alignmentScore)) {
      hctr.not_alignable_skips+=1;
      ez.stopped = 1;
    }
    SPDLOG_DEBUG(logger_,"\t Aligning through MEM chain : ");
    // for the second through the last mem
    for(auto it = mems.begin(); it != mems.end() and !ez.stopped; ++it) {
      auto& mem = *it;
      rpos = mem.rpos;
      memlen = mem.extendedlen;
      tpos = mem.tpos;

      // first score the mem match
      int32_t score = mopts.matchScore * memlen;

      int32_t currMemStart_ref = tpos;
      int32_t currMemStart_read = isFw ? rpos : readLen - (rpos + memlen);
      int32_t gapRef = currMemStart_ref - prevMemEnd_ref - 1; // both inclusive
      int32_t gapRead = currMemStart_read - prevMemEnd_read - 1; // both inclusive

      if ((gapRef <= 0 or gapRead <= 0) and gapRef != gapRead) {
        int32_t gapDiff = std::abs(gapRef - gapRead);
        // subtract off extra matches
        score += mopts.matchScore * std::min(gapRead, gapRef);
        if (computeCIGAR) cigarGen.remove_match(std::abs(std::min(gapRead, gapRef)));
        // Add Insertions to accomodate for the extra gap in the query
        score += -1 * mopts.gapOpenPenalty + -1 * mopts.gapExtendPenalty * gapDiff;
        if (computeCIGAR) cigarGen.add_item(gapDiff, gapRead < gapRef ? 'D' : 'I');

        SPDLOG_DEBUG(logger_,"\t GAP NOT THE SAME:\n\t gapRef : {}, gapRead : {}", gapRef, gapRead);
        SPDLOG_DEBUG(logger_,"\t totalScore (MEM + gapDiff) : {}", score);
      } else if (gapRead > 0 and gapRef > 0) {
        SPDLOG_DEBUG(logger_,"\t\t overlaps : \n\t\t gapRef : {}, gapRead : {}", gapRef, gapRead);

        auto readWindow = readView.substr(prevMemEnd_read + 1, gapRead);
        const char* refSeq1 = tseq.data() + (prevMemEnd_ref) - refStart + 1;

        SPDLOG_DEBUG(logger_,"\t\t aligning\n\t\t [{}]\n\t\t\ [{}]", readWindow, nonstd::string_view(refSeq1, gapRef));
        if (prevMemEnd_ref - refStart + 1 + gapRef >= tseq.size()) {
          SPDLOG_DEBUG(logger_,"\t\t tseq was not long enough; need to fetch more!");
        }
        if (readWindow.length() <= minLengthGapRequired and gapRead == gapRef ) {
          score += align_ungapped(readWindow.data(), refSeq1, gapRead);
        } else {
          hctr.aligner_calls_count += 1;
          hctr.between_aligner_calls_count += 1;
          ksw_reset_extz(&ez);
          score += aligner(readWindow.data(), readWindow.length(), refSeq1, gapRef, &ez,
                    ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::GLOBAL>());
          if (computeCIGAR) { addCigar(cigarGen, &ez, false); }
        }
      } else if ( it > mems.begin() and ((currMemStart_read <= prevMemEnd_read) or (currMemStart_ref <= prevMemEnd_ref)) ) {
        SPDLOG_DEBUG(logger_,"]]\n");
        std::cerr << "[ERROR in PuffAligner, between-MEM alignment] : Improperly compacted MEM chain.  Should not happen!\n";
        std::cerr << "gapRef : " << gapRef << ", gapRead : " << gapRead << ", memlen : " << memlen << "\n";
        std::cerr << "prevMemEnd_read : " << prevMemEnd_read << ", currMemStart_read : " << currMemStart_read << "\n";
        std::cerr << "prevMemEnd_ref  : " << prevMemEnd_ref <<  ", currMemStart_ref  : " << currMemStart_ref << "\n";
        std::exit(1);
      }
      if (computeCIGAR) cigarGen.add_item(memlen, 'M');

      SPDLOG_DEBUG(logger_,"\t The mem matched with the score: {}, CIGAR: {}M", mopts.matchScore * memlen, memlen);
      SPDLOG_DEBUG(logger_,"\t MEM (rpos : {}, memlen : {}, tpos : {})", rpos, memlen, tpos);
      SPDLOG_DEBUG(logger_,"\t gapRef : {}, gapRead : {}", gapRef, gapRead);
      //auto printView = readView.substr(currMemStart_read, memlen);
      //auto refView = nonstd::string_view(tseq.c_str() + tpos - refStart, memlen);

      SPDLOG_DEBUG(logger_,"\t read [{}], pos : {}, len : {}, ori : {}", printView, currMemStart_read, memlen, (isFw ? "FW" : "RC"));
      SPDLOG_DEBUG(logger_,"\t ref  [{}], pos : {}, len : {}", refView, currMemStart_ref, memlen);
      /*if (printView.length() != refView.length()) {
        SPDLOG_DEBUG(logger_,"\t readView length != refView length; should not happen!");
        std::exit(1);
      }*/

      prevMemEnd_read = currMemStart_read + memlen - 1;
      prevMemEnd_ref = tpos + memlen - 1;
      alignmentScore += score;

      if (!alignable(prevMemEnd_read+1, mopts.matchScore, alignmentScore)) {
        hctr.not_alignable_skips+=1;
        ez.stopped = 1;
      }
    }
    // If we got to the end, and there is a read gap left, then align that as well
    SPDLOG_DEBUG(logger_,"prevMemEnd_read : {}, readLen : {}", prevMemEnd_read, readLen);
    bool gapAtEnd = (prevMemEnd_read + 1) <= (static_cast<int32_t>(readLen) - 1);
    if (gapAtEnd) {
      int32_t gapRead = (readLen - 1) - (prevMemEnd_read + 1) + 1;
      int32_t refTailStart = prevMemEnd_ref + 1;
      int32_t dataDepBuff = std::min(refExtLength, 5*gapRead);
      int32_t refTailEnd = refTailStart + gapRead + dataDepBuff;
      if (refTailEnd >= refTotalLength) {refTailEnd = refTotalLength - 1;}
      int32_t refLen = (refTailEnd > refTailStart) ? refTailEnd - refTailStart + 1 : 0;
      auto readWindow = readView.substr(prevMemEnd_read + 1);
      fillRefSeqBuffer(allRefSeq, refAccPos, refTailStart, refLen, refSeqBuffer_);

      SPDLOG_DEBUG(logger_,"POST:");
      SPDLOG_DEBUG(logger_,"read : [{}]", readWindow);
      SPDLOG_DEBUG(logger_,"ref  : [{}]", refSeqBuffer_);
      SPDLOG_DEBUG(logger_,"gapRead : {}, refLen : {}, refBuffer_.size() : {}, refTotalLength : {}", gapRead, refLen, refSeqBuffer_.size(), refTotalLength);

      if (refLen > 0) {
        ksw_reset_extz(&ez);
        hctr.aligner_calls_count += 1;
        aligner(readWindow.data(), readWindow.length(), refSeqBuffer_.data(), refLen, &ez,
                allowOverhangSoftclip, ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
        int32_t alnCost = allowOverhangSoftclip ? std::max(ez.mqe, ez.mte) : ez.mqe;
        int32_t insertCost = (-1 * mopts.gapOpenPenalty + -1 * mopts.gapExtendPenalty * readWindow.length());
        alignmentScore += std::max(alnCost, insertCost);
        SPDLOG_DEBUG(logger_,"POST score : {}", std::max(ez.mqe, ez.mte));
        if (computeCIGAR) {
          SPDLOG_DEBUG(logger_,"The current cigar is: {}", cigarGen.get_cigar());
          if (alnCost >= insertCost) {
            addCigar(cigarGen, &ez, false);
            if (ez.mte > ez.mqe and allowOverhangSoftclip) 
              cigarGen.end_softClipped_len = readWindow.length() - ez.max_q - 1;
          } else {
            cigarGen.add_item(readWindow.length(), 'I');
          }
          SPDLOG_DEBUG(logger_,"The current cigar is: {}", cigarGen.get_cigar());
        }
      } else {
        // do any special soft-clipping penalty here if we want
        if (allowOverhangSoftclip and computeCIGAR)
          cigarGen.end_softClipped_len = readWindow.length();
        else if (computeCIGAR)
          cigarGen.add_item(readWindow.length(), 'I');
        alignmentScore += allowOverhangSoftclip ? 0 :
          (-1 * mopts.gapOpenPenalty + -1 * mopts.gapExtendPenalty * readWindow.length());
      }
    }

    SPDLOG_DEBUG(logger_,"score : {}\n]]\n", alignmentScore);
  }

  if (computeCIGAR) {
     cigar = cigarGen.get_cigar();
     cigarGen.clear();
     SPDLOG_DEBUG(logger_,"score: {}\tcigar : {}\n", alignmentScore, cigar);
  }
  if (isMultimapping_) { // don't bother to fill up a cache unless this is a multi-mapping read
    if (!didHash) {
      // We want the alignment cache to be on the hash of the full underlying reference sequence.
      // If we are using fullAlignment, this is in refSeqBuffer_, but if we are using between-mem alignment
      // then refSeqBuffer_ could have been used to store shorter portions of the reference during
      // the alignment procedure.  In that case, get the original reference sequence from tseq, which
      // was copied from the full reference sequence in the beginning of the function.
      char* ptr = const_cast<char*>((doFullAlignment ? refSeqBuffer_.data() : tseq.data()));
      MetroHash64::Hash(reinterpret_cast<uint8_t *>(ptr), keyLen, reinterpret_cast<uint8_t *>(&hashKey), 0);
    }
    AlignmentResult aln;
    aln.score = alignmentScore;
    if (computeCIGAR) { aln.cigar = cigar; }
    aln.openGapLen = openGapLen;
    alnCache[hashKey] = aln;
  }
  arOut.score = alignmentScore;
  arOut.cigar = cigar;
  arOut.openGapLen = openGapLen;
  return true;
}


/**
 *  Align read_left and read_right, filling the relevant alignment information into the output joinHit structure.
 *  The behavior of alignment (whether the alignment is done only between MEMs or over the full read length, and
 *  if CIGAR strings are computed or just scores, is controlled by the configuration that has been passed to this
 *  PuffAligner object).
 **/
int32_t PuffAligner::calculateAlignments(std::string& read_left, std::string& read_right, pufferfish::util::JointMems& jointHit,
                                         HitCounters& hctr, bool isMultimapping, bool verbose) {
    isMultimapping_ = isMultimapping;
    auto tid = jointHit.tid;
    double optFrac{mopts.minScoreFraction};
    bool computeCIGAR = !(aligner.config().flag & KSW_EZ_SCORE_ONLY);
    auto threshold = [&, optFrac] (uint64_t len) -> double {
        return (!mopts.matchScore)?(-0.6+-0.6*len):optFrac*mopts.matchScore*len;
    };
    constexpr const auto invalidScore = std::numeric_limits<decltype(ar_left.score)>::min();

    // If this read is in an orphaned mapping
    if (jointHit.isOrphan()) {
        hctr.totalAlignmentAttempts += 1;

        // If this mapping was an orphan, then this is the orphaned read
        bool isLeft = jointHit.isLeftAvailable();
        std::string& read_orphan = isLeft ? read_left : read_right;
        std::string& rc_orphan = isLeft ? read_left_rc_ : read_right_rc_;
        auto& ar_orphan = isLeft ? ar_left : ar_right;
        auto& orphan_aln_cache = isLeft ? alnCacheLeft : alnCacheRight;

        ar_orphan.score = invalidScore;
        alignRead(read_orphan, rc_orphan, jointHit.orphanClust()->mems,
                  jointHit.orphanClust()->perfectChain,
                  jointHit.orphanClust()->isFw, tid, orphan_aln_cache, hctr, ar_orphan, true);
        jointHit.alignmentScore =
          ar_orphan.score > threshold(read_orphan.length())  ? ar_orphan.score : invalidScore;
        jointHit.orphanClust()->cigar = (computeCIGAR) ? ar_orphan.cigar : "";
        jointHit.orphanClust()->openGapLen = ar_orphan.openGapLen;
        //jointHit.orphanClust()->coverage = jointHit.alignmentScore;
        if (jointHit.alignmentScore < 0 and verbose) {
          std::cerr << read_orphan.length() << " " << threshold(read_orphan.length()) << " " << ar_left.score << "\n";
        }
        return jointHit.alignmentScore;
    } else {
        hctr.totalAlignmentAttempts += 2;
        ar_left.score = ar_right.score = invalidScore;
        if (verbose) { std::cerr << "left\n"; }
        alignRead(read_left, read_left_rc_, jointHit.leftClust->mems, jointHit.leftClust->perfectChain,
                                            jointHit.leftClust->isFw, tid, alnCacheLeft, hctr, ar_left, true);

        if (verbose) { std::cerr << "right\n"; }
        alignRead(read_right, read_right_rc_, jointHit.rightClust->mems, jointHit.rightClust->perfectChain,
                                             jointHit.rightClust->isFw, tid, alnCacheRight, hctr, ar_right, false);

        jointHit.alignmentScore = ar_left.score > threshold(read_left.length()) ? ar_left.score : invalidScore;
        jointHit.mateAlignmentScore = ar_right.score > threshold(read_right.length()) ? ar_right.score : invalidScore;
/*
        jointHit.alignmentScore = (score_left == invalidScore or score_right == invalidScore)?
                                  invalidScore : score_left + score_right;
*/
//      jointHit.leftClust->coverage = score_left;
        jointHit.leftClust->openGapLen = ar_left.openGapLen;
//      jointHit.rightClust->coverage = score_right;
        jointHit.rightClust->openGapLen = ar_right.openGapLen;
        if (computeCIGAR) {
          jointHit.leftClust->cigar = ar_left.cigar;
          jointHit.rightClust->cigar = ar_right.cigar;
        } else {
          jointHit.leftClust->cigar = "";
          jointHit.rightClust->cigar = "";
        }
        return (jointHit.alignmentScore == invalidScore or jointHit.mateAlignmentScore == invalidScore) ?
        invalidScore:jointHit.alignmentScore+jointHit.mateAlignmentScore;
    }
}



/**
 *  Align read, filling the relevant alignment information into the output joinHit structure.
 *  The behavior of alignment (whether the alignment is done only between MEMs or over the full read length, and
 *  if CIGAR strings are computed or just scores, is controlled by the configuration that has been passed to this
 *  PuffAligner object).
 **/
int32_t PuffAligner::calculateAlignments(std::string& read, pufferfish::util::JointMems& jointHit, HitCounters& hctr, bool isMultimapping, bool verbose) {
  isMultimapping_ = isMultimapping;
    auto tid = jointHit.tid;
    double optFrac{mopts.minScoreFraction};
    bool computeCIGAR = !(aligner.config().flag & KSW_EZ_SCORE_ONLY);
    auto threshold = [&, optFrac] (uint64_t len) -> double {
        return (!mopts.matchScore)?(-0.6+-0.6*len):optFrac*mopts.matchScore*len;
    };
    constexpr const auto invalidScore = std::numeric_limits<decltype(ar_left.score)>::min();

    hctr.totalAlignmentAttempts += 1;
    ar_left.score = invalidScore;
    const auto& oc = jointHit.orphanClust();
    alignRead(read, read_left_rc_, oc->mems, oc->perfectChain, oc->isFw, tid, alnCacheLeft, hctr, ar_left, true);
    jointHit.alignmentScore =
      ar_left.score > threshold(read.length())  ? ar_left.score : invalidScore;
    jointHit.orphanClust()->cigar = (computeCIGAR) ? ar_left.cigar : "";
    jointHit.orphanClust()->openGapLen = ar_left.openGapLen;
    //jointHit.orphanClust()->coverage = jointHit.alignmentScore;
    if (jointHit.alignmentScore < 0 and verbose) {
      std::cerr << read.length() << " " << threshold(read.length()) << " " << ar_left.score << "\n";
    }
    return jointHit.alignmentScore;
}

bool PuffAligner::recoverSingleOrphan(std::string& read_left, std::string& read_right, pufferfish::util::MemCluster& clust, std::vector<pufferfish::util::MemCluster> &recoveredMemClusters, uint32_t tid, bool anchorIsLeft, bool verbose) {
  int32_t anchorLen = anchorIsLeft ? read_left.length() : read_right.length();
  auto tpos = clust.mems[0].tpos;
  auto anchorStart = clust.mems[0].isFw ? clust.mems[0].rpos : anchorLen - (clust.mems[0].rpos + clust.mems[0].extendedlen);
  uint32_t anchorPos = tpos >= anchorStart ? tpos - anchorStart : 0;

  bool recovered_fwd;
  uint32_t recovered_pos=-1;

  auto* r1 = read_left.data();
  auto* r2 = read_right.data();
  auto l1 = static_cast<int32_t>(read_left.length());
  auto l2 = static_cast<int32_t>(read_right.length());
  const char* rptr{nullptr};
  bool anchorFwd{clust.isFw};
  int32_t startPos = -1, maxDist = -1, otherLen = -1, rlen = -1;
  std::string* otherReadPtr{nullptr};
  const char* otherRead{nullptr};
  char* otherReadRC{nullptr};
  std::string* otherRCSpace{nullptr};
  char* r1rc = nullptr;
  char* r2rc = nullptr;

  std::unique_ptr<char[]> windowSeq{nullptr};
  int32_t windowLength = -1;

  int32_t maxDistRight = l2 / 4;
  int32_t maxDistLeft = l1 / 4;
  constexpr const int32_t signedZero{0};

  if (anchorIsLeft) {
    anchorLen = l1;
    otherLen = l2;
    maxDist = maxDistRight;
    otherReadPtr = &read_right;
    otherRCSpace = &rc2_;
    otherRead = r2;
    otherReadRC = r2rc;
    /* from rapmap
    anchorLen = l1;
    otherLen = l2;
    maxDist = maxDistRight;
    lpos = anchorPos;
    rpos = -1;
    lfwd = anchorFwd;
    rfwd = !lfwd;
    otherReadPtr = &rightRead;
    otherRCSpace = &rc2;
    otherRead = r2;
    otherReadRC = r2rc;
    leftChainStatus = anchorHit.chainStatus.getLeft();
    */
  } else {
    anchorLen = l2;
    otherLen = l1;
    maxDist = maxDistLeft;
    otherReadPtr = &read_left;
    otherRCSpace = &rc1_;
    otherRead = r1;
    otherReadRC = r1rc;
  }

  uint64_t refAccPos = tid > 0 ? refAccumLengths[tid - 1] : 0;
  uint64_t refLength = refAccumLengths[tid] - refAccPos;

  if (anchorFwd) {
    if (!otherReadRC){
      pufferfish::util::reverseRead(*otherReadPtr, *otherRCSpace);
      otherReadRC = const_cast<char*>(otherRCSpace->data());
    }
    rptr = otherReadRC;
    rlen = otherLen;
    startPos = std::max(signedZero, static_cast<int32_t>(anchorPos));
    windowLength = std::min(500, static_cast<int32_t>(refLength - startPos));
  } else {
    rptr = otherRead;
    rlen = otherLen;
    int32_t endPos = std::min(static_cast<int32_t>(refLength), static_cast<int32_t>(anchorPos) + anchorLen);
    startPos = std::max(signedZero,  endPos - 500);
    windowLength = std::min(500, endPos);
  }

  if (verbose) { std::cerr<< anchorPos<< "\n"; }
  fillRefSeqBuffer(allRefSeq, refAccPos, startPos, windowLength, refSeqBuffer_);
  /*windowSeq.reset(new char[tseq.length() + 1]);
  strcpy(windowSeq.get(), tseq.c_str());
  */

  // Note -- we use score only mode to find approx end position in rapmap, can we
  // do the same here?
  EdlibAlignResult result = edlibAlign(rptr, rlen, refSeqBuffer_.data(), windowLength,
                                       edlibNewAlignConfig(maxDist, EDLIB_MODE_HW, EDLIB_TASK_LOC));

  if (result.editDistance > -1) {
    recovered_fwd = !anchorFwd;
    recovered_pos = startPos + result.startLocations[0];
    recoveredMemClusters.push_back(pufferfish::util::MemCluster(recovered_fwd, rlen));
    auto it = recoveredMemClusters.begin() + recoveredMemClusters.size() - 1;
    if (verbose) {
      std::cerr<< anchorIsLeft << " " << anchorFwd << " " <<  anchorPos << " " << startPos + result.startLocations[0] <<" " << result.editDistance << "\n";
    }
    orphanRecoveryMemCollection.push_back(pufferfish::util::UniMemInfo());
    auto memItr = orphanRecoveryMemCollection.begin() + orphanRecoveryMemCollection.size() - 1;
    if (verbose) { std::cerr<<recovered_fwd << " " << orphanRecoveryMemCollection.size()<<"\n"; }
    it->addMem(memItr, recovered_pos, 1, recovered_fwd ? 1 : rlen-1, recovered_fwd);

    //delete windowSeq;
    edlibFreeAlignResult(result);
    return true;
  } else {
    //delete windowSeq;
    edlibFreeAlignResult(result);
    return false;
  }
}
