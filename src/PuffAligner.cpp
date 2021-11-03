#include "nonstd/string_view.hpp"
#include "PuffAligner.hpp"
#include "Util.hpp"
#include "libdivide/libdivide.h"

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

int32_t addCigar(pufferfish::util::CIGARGenerator &cigarGen, ksw_extz_t ez, bool beginGap) {
    if (!beginGap) {
        int32_t insertionDeletion = 0;
        for (int i = 0; i < ez.n_cigar; ++i) {
            char cigar_type = 'N';
            if ((ez.cigar[i] & 0xf) == 0) cigar_type = 'M';
            if ((ez.cigar[i] & 0xf) == 1) {
                cigar_type = 'I';
                insertionDeletion += ez.cigar[i] >> 4;
            }
            if ((ez.cigar[i] & 0xf) == 2) {
                cigar_type = 'D';
                insertionDeletion -= ez.cigar[i] >> 4;
            }
            cigarGen.add_item(ez.cigar[i] >> 4, cigar_type);
        }
        return insertionDeletion;
    } else {
        uint32_t gapSize = 0;
        for (int i = ez.n_cigar - 1; i >= 0; --i) {
            char cigar_type = 'N';
            if ((ez.cigar[i] & 0xf) == 0) {
                cigar_type = 'M';
                gapSize += ez.cigar[i] >> 4;
            }
            if ((ez.cigar[i] & 0xf) == 1) { cigar_type = 'I'; }
            if ((ez.cigar[i] & 0xf) == 2) {
                cigar_type = 'D';
                gapSize += ez.cigar[i] >> 4;
            }
            cigarGen.add_item(ez.cigar[i] >> 4, cigar_type);
        }
        return gapSize;
    }
}

std::string getRefSeq(compact::vector<uint64_t, 2> &refseq, uint64_t refAccPos, size_t tpos, uint32_t memlen) {
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
}

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


// Ungapped Alignment of two sequences with the same length 
int32_t PuffAligner::align_ungapped(const char* const query, const char* const target, int32_t len) {
  int32_t score = 0;
  // bool computeCIGAR = !(aligner.config().flag & KSW_EZ_SCORE_ONLY);
  auto& cigarGen = cigarGen_;
  for (int32_t i = 0; i < len; i++) {
    if (mopts.computeCIGAR) cigarGen.add_item(1, 'M');
    if (query[i] != 'N' and target[i] != 'N')
      score += query[i] == target[i] ? mopts.matchScore : mopts.mismatchPenalty;
  }
  return score;
}

// FOR DEBUGGING PURPOSES. Builds CIGAR string from the expanded CIGAR (contains =/X instead of M)
std::string collapse_cigar(std::string exp_cigar)
{
  if(exp_cigar.size() < 1) return "";
  std::string collapsed = "";
  char last_char = exp_cigar[0];
  int cnt = 1;
  for(unsigned int i = 1; i < exp_cigar.size(); i++)
  {
    if(exp_cigar[i] == last_char)
    {
      cnt++;
    }
    else
    {
      collapsed += std::to_string(cnt) + last_char;
      last_char = exp_cigar[i];
      cnt = 1;
    }
  }
  collapsed += std::to_string(cnt) + last_char;
  return collapsed;
}

// FOR DEBUGGING PURPOSES. Builds the CIGAR string as well as the expanded CIGAR string (contains =/X and needs to be collapsed with collapse_ciga() function)
std::string ksw2_to_cigar(const char *qseq, const char *tseq, ksw_extz_t *ez, std::string &exp_cigar)
{
  int qpos = 0;
  int tpos = 0;
  std::string cigar = "";
  exp_cigar = "";
  if (ez->n_cigar > 0) {
    for (int i = 0; i < ez->n_cigar; ++i)
    {
      if((ez->cigar[i]&0xf) == 1)
      {
        qpos += ez->cigar[i]>>4;
        exp_cigar += std::string(ez->cigar[i]>>4, 'I');
      }
      else if ((ez->cigar[i]&0xf) == 2)
      {
        tpos += ez->cigar[i]>>4;
        exp_cigar += std::string(ez->cigar[i]>>4, 'D');
      }
      else
      {
        for(unsigned int j = 0; j < ez->cigar[i]>>4; j++)
        {
          if(tseq[tpos] == qseq[qpos])
            exp_cigar += '=';
          else
            exp_cigar += 'X';
          tpos++;
          qpos++;
        }
      }
      cigar += std::to_string(ez->cigar[i]>>4) + "MID"[ez->cigar[i]&0xf];
    }
  }
  return cigar;
}

/**
 *  Align the read `original_read`, whose mems consist of `mems` against the index and return the result
 *  in `arOut`.  How the alignment is computed (i.e. full vs between-mem and CIGAR vs. score only) depends
 *  on the parameters of how this PuffAligner object was constructed.
 **/
bool PuffAligner::alignRead(std::string& read, std::string& read_rc, const std::vector<pufferfish::util::MemInfo> &mems, 
                            uint64_t queryChainHash, bool perfectChain, bool isFw, size_t tid, AlnCacheMap &alnCache, 
                            HitCounters &hctr, AlignmentResult& arOut, std::shared_ptr<spdlog::logger>& logger_, bool /*verbose*/) {

  int32_t alignmentScore{std::numeric_limits<decltype(arOut.score)>::min()};
  if (mems.empty()) {
    arOut.score = alignmentScore;
    return false;
  }

  int32_t refExtLength = static_cast<int32_t>(mopts.refExtendLength);
  // bool firstMem = true;
  // int32_t lastHitEnd_read = -1;
  int32_t currHitStart_read = 0;
  // int64_t lastHitEnd_ref = -1;
  int64_t currHitStart_ref = 0;
  // int32_t alignment{0};
  uint32_t openGapLen{0};

  // the cigar generator we will use
  auto& cigarGen = cigarGen_;
  cigarGen.clear();

  ksw_reset_extz(&ez);

  // where this reference starts, and its length.
  int64_t refAccPos = tid > 0 ? refAccumLengths[tid - 1] : 0;
  int64_t refTotalLength = refAccumLengths[tid] - refAccPos;

  const auto& frontMem = mems.front();
  auto rpos = frontMem.rpos;
  auto memlen = frontMem.extendedlen;
  auto readLen = read.length();
  auto tpos = frontMem.tpos;

  if (perfectChain) {
    logger_->debug("\t\tperfect chain!");
    hctr.skippedAlignments_byCov += 1;
    arOut.score = alignmentScore = readLen * mopts.matchScore;
    arOut.cigar = (mopts.computeCIGAR ? std::to_string(readLen) + 'M' : "");
    arOut.openGapLen = openGapLen;
    return true;
  }


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

  int32_t signedRefStartPos = currHitStart_ref - currHitStart_read;
  int32_t signedRefEndPos = signedRefStartPos + read.length();
  bool invalidStart = (signedRefStartPos < 0);
  bool invalidEnd = (signedRefEndPos > refTotalLength);

  uint32_t maxSoftclipLenGeneral = mopts.maxSoftclipFractionGeneral * readLen;
  uint32_t maxSoftclipLenOverhang = mopts.maxSoftclipFractionOverhang * readLen;
  
  if (mopts.end2end) {
    maxSoftclipLenGeneral = 0;
    maxSoftclipLenOverhang = 0;
  }
  uint32_t maxSoftclipLen = std::max(maxSoftclipLenGeneral, maxSoftclipLenOverhang);

  int32_t remainedSoftClipLenGeneral = maxSoftclipLenGeneral;
  int32_t remainedSoftClipLenOverhang = maxSoftclipLenOverhang;

  if (mopts.mimicBT2Strict and (invalidStart or invalidEnd)) {
    arOut.score = std::numeric_limits<decltype(arOut.score)>::min();
    arOut.cigar = "";
    arOut.openGapLen = 0;
    return false;
  }

  auto& bandwidth = aligner.config().bandwidth;
  auto& aligner_config = aligner.config();
  libdivide::divider<int32_t> gapExtDivisor(static_cast<int32_t>(mopts.gapExtendPenalty));
  const int32_t minAcceptedScore = scoreStatus_.getCutoff(read.length() - maxSoftclipLen); //mopts.minScoreFraction * mopts.matchScore * (readLen - maxSoftclipLen);
  logger_->debug("\t\tNOTE mems.size(): {}, isRev: {}, minAcceptedScore: {}", mems.size(), !isFw, minAcceptedScore);
  // compute the maximum gap length that would be allowed given the length of read aligned so far and the current 
  // alignment score.
  auto maxAllowedGaps = [&minAcceptedScore, &readLen, &gapExtDivisor, this] (uint32_t alignedLen, int32_t alignmentScore) -> int {
    int maxAllowedGaps = (alignmentScore + static_cast<int32_t>(mopts.matchScore) * static_cast<int32_t>(readLen - alignedLen) - minAcceptedScore - static_cast<int32_t>(mopts.gapOpenPenalty)) / gapExtDivisor;
    return std::max(maxAllowedGaps + 1, 1);                                                                       
  };
  auto alignable = [&minAcceptedScore, &readLen] (uint32_t alignedLen, uint32_t matchScore, int32_t alignmentScore) -> bool {
    return minAcceptedScore <= alignmentScore + static_cast<int32_t>(matchScore*(readLen-alignedLen));
  };

  const auto& lastMem = mems.back();
  int32_t last_rpos = isFw ? lastMem.rpos : readLen - (lastMem.rpos + lastMem.extendedlen);
  int32_t leftEndRead = readLen - (last_rpos + lastMem.extendedlen);
  int32_t dataDepBuff = maxAllowedGaps(0, 0) + 1; // std::min(refExtLength, 5 * leftEndRead);
  bool overhangingEnd = static_cast<int64_t>(lastMem.tpos + lastMem.extendedlen + leftEndRead + dataDepBuff) > refTotalLength;
  bool overhangingStart = currHitStart_ref < currHitStart_read;
  uint32_t remLen = 0;
  // If we are only aligning between MEMs
  if (!doFullAlignment) {
    refStart = (currHitStart_ref >= currHitStart_read) ? currHitStart_ref - currHitStart_read : 0;
    overhangingStart = static_cast<int32_t>(refStart) < refExtLength;
    int32_t refStartExtLength = (refExtLength < static_cast<int32_t>(refStart)) ? refExtLength : refStart;
    refStart = static_cast<int32_t>(refStart) > refExtLength ? static_cast<int32_t>(refStart) - refExtLength : 0;
    remLen = (currHitStart_ref >= currHitStart_read) ? readLen : readLen - (currHitStart_read - currHitStart_ref);
    keyLen = (static_cast<int64_t>(refStart + refStartExtLength + remLen + refExtLength) < refTotalLength) ? refStartExtLength + remLen + refExtLength : refTotalLength - refStart;
  } else { // we are aligning from the start of the read
    // If the first hit starts further in the reference than in the
    // read, then we align from the beginning of the read and (ref_start - read_start) on
    // the reference.
    if (currHitStart_ref > currHitStart_read) {
      refStart = currHitStart_ref - currHitStart_read;
      readStart = 0;
      remLen = readLen;
    } else if (currHitStart_ref < currHitStart_read) {
      // If the first his starts further in the read than in the reference, than
      // the read overhangs the beginning of the reference and we start aligning
      // from the beginning of the reference and from position (read_start - ref_start)
      // in the read.
      readStart = currHitStart_read - currHitStart_ref;
      refStart = 0;
      remLen = readLen - readStart;
    } else {
      // If the first hit starts at the same position in the reference and the read
      // ... what is this case?
      readStart = 0;
      refStart = currHitStart_ref - currHitStart_read;
      remLen = readLen;
    }
    keyLen = (static_cast<int64_t>(refStart + remLen + refExtLength) < refTotalLength) ? remLen + refExtLength : refTotalLength - refStart;
  }

  fillRefSeqBuffer(allRefSeq, refAccPos, refStart, keyLen, refSeqBuffer_);
  // int32_t originalRefSeqLen = static_cast<int32_t>(refSeqBuffer_.length());
  // If we're not using fullAlignment, we'll need the full reference sequence later
  // so copy it into tseq.
  if (!doFullAlignment) { tseq = refSeqBuffer_; }

  // Mohsen says that we avoid using the cache only when both ends are overhanging (including extensions)
  bool useAlnCache = mopts.useAlignmentCache and isMultimapping_ and !(overhangingEnd and overhangingStart);

  if (useAlnCache and !alnCache.empty()) {
    logger_->debug("\t\tNOTE checking cache");
    // hash the reference string
    MetroHash64::Hash(reinterpret_cast<uint8_t *>(const_cast<char*>(refSeqBuffer_.data())), keyLen, reinterpret_cast<uint8_t *>(&hashKey), 0);
    hashKey ^= queryChainHash;
    didHash = true;
    // see if we have this hash
    auto hit = alnCache.find(hashKey);
    // if so, we know the alignment score
    if (hit != alnCache.end() and (isFw == hit->second.isFw) ) {//}and refStart + readLen + refExtLength < refTotalLength) {
      logger_->debug("\t\tNOTE cache hit... done!");
      hctr.skippedAlignments_byCache += 1;
      arOut.score = hit->second.score;
      if (mopts.computeCIGAR) { arOut.cigar = hit->second.cigar; }
      if (mopts.allowSoftclip) {
        arOut.softclip_start = hit->second.softclip_start;
        arOut.softclip_end = hit->second.softclip_end;
      }
      arOut.openGapLen = hit->second.openGapLen;
      return true;
    }
    logger_->debug("\t\tNOTE cache miss... continue");
  }

  // @mohsen & @fataltes --- we should figure out how to
  // avoid computing the rc of a read if we've already done it.
  if (!isFw and read_rc.empty()) { read_rc = pufferfish::util::reverseComplement(read); }
  nonstd::string_view readView = (isFw) ? read : read_rc;

  logger_->debug("\t\t+++ NEXT MEM readPos: {}, refPos: {}, memLen: {}", (isFw ? rpos : readLen - (rpos + memlen)), tpos, memlen);

  // NOTE: We don't worry about soft clipping within the `doFullAlignment`
  // branch, since these two options are currently incompatible.
  if (doFullAlignment) {
    // if we allow softclipping of overhanging bases, then we can cut off the
    // part of the read before the start of the reference
    // Note: a more accurate way is to perform an extensention to the left as well
    decltype(readStart) readOffset = (mopts.allowSoftclip && (int32_t)(readStart) <= remainedSoftClipLenOverhang) ? readStart : 0;
    remainedSoftClipLenOverhang -= readOffset;
    remainedSoftClipLenGeneral -= readOffset;
    arOut.softclip_start = readOffset;

    if (mopts.computeCIGAR && readOffset > 0) {
      cigarGen.add_item(readOffset, 'S');
    }
    
    nonstd::string_view readSeq = readView.substr(readOffset);
    auto cutoff = minAcceptedScore - mopts.matchScore * read.length();
    aligner(readSeq.data(), readSeq.length(), refSeqBuffer_.data(),
            refSeqBuffer_.length(), &ez, cutoff, remainedSoftClipLenGeneral,
            ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
    // if we allow softclipping of overhaning bases, then we only care about
    // the best score to the end of the query or the end of the reference.
    // Otherwise, we care about the best score all the way until the end of
    // the query.
    // modify this to select between mqe and max like MEMs
    int32_t numSoftClipped = readSeq.length() - (ez.max_q + 1);

    if (remainedSoftClipLenGeneral < numSoftClipped || ez.mqe + aligner_config.end_bonus > ez.max)
    {
      alignmentScore = ez.mqe;
      numSoftClipped = 0;
    }
    else
    {
      alignmentScore = ez.max;
    }
    
    arOut.softclip_end = numSoftClipped;

    logger_->debug("readSeq : {}\nrefSeq  : {}\nscore   : {}\nreadStart : {}", readSeq, refSeqBuffer_, alignmentScore, readStart);
    logger_->debug("currHitStart_read : {}, currHitStart_ref : {}\nmqe : {}, mte : {}\n", currHitStart_read, currHitStart_ref, ez.mqe, ez.mte);

    if (mopts.computeCIGAR) {
      openGapLen = addCigar(cigarGen, ez, false);
    } else {
      // can make start pos negative, but sam writer deals with this right now
      openGapLen = currHitStart_read;
      //((currHitStart_ref - currHitStart_read) < 0) ? 0 :
      //(currHitStart_read);
    }
  } else {
    alignmentScore = 0;

    logger_->debug("\t\tread sequence ({}) : {}", (isFw ? "FW" : "RC"), readView);
    logger_->debug("\t\tref  sequence      : {}", tseq);

    // If the first mem does not start at the beginning of the
    // read, then there is a gap to align.
    int32_t firstMemStart_read = (isFw) ? rpos : readLen - (rpos + memlen);
    if (firstMemStart_read > 0) {
      logger_->debug("\t\t### EXTEND LEFT");
      // align the part before the first mem

      // the gap is of length firstMemStart_read, so grab that much (plus
      // buffer) before the first start position on the reference. int32_t
      // firstMemStart_ref = tpos;
      int32_t readStartPosOnRef = tpos - firstMemStart_read;
      // int32_t dataDepBuff = std::min(refExtLength, 5*firstMemStart_read);

      int32_t maxGaps = maxAllowedGaps(0, 0) + 1;

      int32_t refWindowStart = (readStartPosOnRef - maxGaps) > 0
                                    ? (readStartPosOnRef - maxGaps)
                                    : 0;
      int32_t refWindowLength = tpos - refWindowStart;
      fillRefSeqBufferReverse(allRefSeq, refAccPos, refWindowStart,
                              refWindowLength, refSeqBuffer_);

      int32_t numSoftClipped = 0;

      if (refSeqBuffer_.length() > 0) {
        logger_->debug("\t\t\tCASE 1: some reference bases left to align");
        auto readWindow = readView.substr(0, firstMemStart_read).to_string();
        std::reverse(readWindow.begin(), readWindow.end());
        logger_->debug("\t\t\treadStartPosOnRef: {}\trefWindowStart: {}\trefWindowLength: {}", readStartPosOnRef, refWindowStart, refWindowLength);
        logger_->debug("\t\t\tread: [{}]", readWindow);
        logger_->debug("\t\t\tref:  [{}]", refSeqBuffer_);
        
        bandwidth = maxGaps;
        logger_->debug("\t\t\tksw2_parameters: bandwidth={}, end_bonus={}, zdrop={}", bandwidth, aligner_config.end_bonus, aligner_config.dropoff);
        auto cutoff = minAcceptedScore - mopts.matchScore * read.length();
        aligner(readWindow.data(), readWindow.length(), refSeqBuffer_.data(),
                refSeqBuffer_.length(), &ez, cutoff, remainedSoftClipLenGeneral,
                ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
        logger_->debug("\t\t\tksw2_results:");
        logger_->debug("\t\t\t\tmax={}, max_q={}, max_t={}", ez.max, ez.max_q, ez.max_t);
        logger_->debug("\t\t\t\tmqe={}, mqe_t={}", ez.mqe, ez.mqe_t);
        logger_->debug("\t\t\t\tmte={}, mte_q={}", ez.mte, ez.mte_q);
        logger_->debug("\t\t\t\tstopped={}, reach_end={}, zdropped={}, score={}", ez.stopped == 1, ez.reach_end == 1, ez.zdropped == 1, ez.score);
        // FOR DEBUGGING PURPOSES ONLY
        // if(mopts.computeCIGAR) {
        //   std::string tmp_cg_exp;
        //   logger_->debug("\t\t\t\talnCigar: {}", ksw2_to_cigar(readWindow.c_str(), refSeqBuffer_.c_str(), &ez, tmp_cg_exp));
        //   logger_->debug("\t\t\t\talnCigar: {}", collapse_cigar(tmp_cg_exp));
        // } 
        if (ez.stopped) hctr.skippedAlignments_notAlignable += 1;
        if (ez.mqe != KSW_NEG_INF and ez.stopped>0) ez.stopped = 0;

        // NOTE: pre soft-clip code for adjusting the alignment score.
        // alignmentScore += allowOverhangSoftclip ? std::max(ez.mqe, ez.mte) : ez.mqe;

        // NOTE: If we are not writing down the CIGAR, we assume *for the
        // purpose of computing positions* that the prefix of the read we are
        // aligning here is optimally aligned in a gapless fashion. This,
        // obviously, may not be the case.  However, if the aligned prefix
        // does contain indels, and we don't write down the appropriate CIGAR
        // string to designate this fact, then downstream tools will not be
        // able to properly adjust (so that e.g. the MEMs in the chain are
        // properly positioned). This approximation means that the position we
        // report may not be the one corresponding to the optimal alignment
        // score that we actually calculate.  However, if we are not writing
        // out the actual CIGAR string, this is the best we can do. (addresses
        // https://github.com/COMBINE-lab/salmon/issues/475).

        decltype(alignmentScore) part_score = ez.max;
        numSoftClipped = readWindow.length() - (ez.max_q + 1);
        if(remainedSoftClipLenGeneral < numSoftClipped || ez.mqe + aligner_config.end_bonus > ez.max)
        {
          part_score = ez.mqe;
          openGapLen = ez.mqe_t + 1;
          numSoftClipped = 0;
        }
        else
        {
          part_score = ez.max;
          openGapLen = ez.max_t + 1;
          cigarGen.add_item(numSoftClipped, 'S');
        }
        remainedSoftClipLenGeneral -= numSoftClipped;
        remainedSoftClipLenOverhang -= numSoftClipped;
        logger_->debug("\t\t\t\tremainedSoftClipLen={}->{};", remainedSoftClipLenGeneral + numSoftClipped, remainedSoftClipLenGeneral);

        alignmentScore += part_score;
        addCigar(cigarGen, ez, true);
      } else {
        logger_->debug("\t\t\tCASE 2: no reference bases left (start)");
        // do any special soft clipping penalty here if we want
        numSoftClipped = 
            mopts.allowSoftclip
            ? firstMemStart_read
            : 0;
        alignmentScore +=
            mopts.allowSoftclip && remainedSoftClipLenOverhang >= numSoftClipped
                ? 0
                : (-1 * mopts.gapOpenPenalty +
                    -1 * mopts.gapExtendPenalty * firstMemStart_read);
        openGapLen = firstMemStart_read;

        char cigar_char;
        if (mopts.allowSoftclip && remainedSoftClipLenOverhang >= numSoftClipped) {
          cigar_char = 'S';
          openGapLen = 0;
        } else {
          cigar_char = 'I';
          numSoftClipped = 0;
        }

        if (mopts.computeCIGAR) {
          cigarGen.add_item(firstMemStart_read, cigar_char);
        }
      }
      arOut.softclip_start = numSoftClipped;
      logger_->debug("\t\t\tscore_sofar: {}", alignmentScore);
      logger_->debug("\t\t\tcigar_sofar: {}", cigarGen.get_cigar());
    }

    int32_t prevMemEnd_read = firstMemStart_read - 1; //isFw ? rpos : readLen - (rpos + memlen);
    int32_t prevMemEnd_ref = tpos - 1;
    
    if (!alignable(prevMemEnd_read + 1, mopts.matchScore, alignmentScore)) {
      hctr.skippedAlignments_notAlignable += 1;
      ez.stopped = 1;
      logger_->debug("ez stopped");
    }

    logger_->debug("\t\t### ADD MEM SCORE");
    alignmentScore += mopts.matchScore * memlen;
    cigarGen.add_item(memlen, 'M');
    prevMemEnd_read = firstMemStart_read + memlen - 1;
    prevMemEnd_ref = tpos + memlen - 1;

    // for the second through the last mem
    for(auto it = mems.begin() + 1; it != mems.end() and !ez.stopped; ++it) {
      logger_->debug("\t\t\tscore_sofar: {}", alignmentScore);
      logger_->debug("\t\t\tcigar_sofar: {}", cigarGen.get_cigar());
      logger_->debug("\t\t### FILL GAP");
      auto& mem = *it;
      rpos = mem.rpos;
      memlen = mem.extendedlen;
      tpos = mem.tpos;
      logger_->debug("\t\t+++ NEXT MEM rpos: {}, tpos: {}, memLen: {}", (isFw ? rpos : readLen - (rpos + memlen)), tpos, memlen);

      int32_t score = 0;

      int32_t currMemStart_ref = tpos;
      int32_t currMemStart_read = isFw ? rpos : readLen - (rpos + memlen);
      int32_t gapRef =
          currMemStart_ref - prevMemEnd_ref - 1; // both inclusive
      int32_t gapRead =
          currMemStart_read - prevMemEnd_read - 1; // both inclusive
      
      logger_->debug("\t\t\tgapRef: {}, gapRead: {}", gapRef, gapRead);

      if ((gapRef <= 0 or gapRead <= 0) and gapRef != gapRead) {
        logger_->debug("\t\t\tCASE 1: overlapping MEMs");
        int32_t gapDiff = std::abs(gapRef - gapRead);
        score += (-1 * mopts.gapOpenPenalty +
                  -1 * mopts.gapExtendPenalty * gapDiff);
        // subtract off extra matches
        score += mopts.matchScore * std::min(gapRead, gapRef);
        cigarGen.add_item(std::min(gapRead, gapRef), 'M');
        if(gapRef > gapRead)
        {
          cigarGen.add_item(gapDiff, 'D');
        }
        else // gapRef < gapRead
        {
          cigarGen.add_item(gapDiff, 'I');
        }

        logger_->debug("\t\t\t\tscore_adjustment: {}", score);
      } else if (gapRead > 0 and gapRef > 0) {
        logger_->debug("\t\t\tCASE 2: gaps to align");

        auto readWindow = readView.substr(prevMemEnd_read + 1, gapRead);
        const char* refSeq1 = tseq.data() + (prevMemEnd_ref)-refStart + 1;

        logger_->debug("\t\t\t\t[{}]\n\t\t\t\t[{}]", readWindow, nonstd::string_view(refSeq1, gapRef));
        if (prevMemEnd_ref - refStart + 1 + gapRef >= tseq.size()) {
          logger_->debug("\t\t WARNING: tseq was not long enough; need to fetch more!");
        }

        if (readWindow.length() <= minLengthGapRequired and gapRead == gapRef) {
          logger_->debug("\t\t\t\taligning with align_ungapped (readWindow.length()={} < minLengthGapRequired={})", readWindow.length(), minLengthGapRequired);
          score += align_ungapped(readWindow.data(), refSeq1, gapRead);
        } else {
          logger_->debug("\t\t\t\taligning with ksw2");
          bandwidth = maxAllowedGaps(prevMemEnd_read + 1, alignmentScore) + 1;
          score += aligner(
            readWindow.data(), readWindow.length(), refSeq1, gapRef, &ez,
            ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::GLOBAL>());
            
          if (mopts.computeCIGAR) {
            addCigar(cigarGen, ez, false);
          }
        }

        logger_->debug("\t\t\t\tscore_adjustment: {}", score);
      } else if (it > mems.begin() and
                  ((currMemStart_read <= prevMemEnd_read) or
                  (currMemStart_ref <= prevMemEnd_ref))) {
        logger_->debug("\t\t\tCASE 3: Improperly compacted MEM chain. Should not happen!");
        std::cerr << "[ERROR in PuffAligner, between-MEM alignment]: "
                  << "Improperly compacted MEM chain. Should not happen!\n";
        std::cerr << "gapRef : " << gapRef << ", gapRead : " << gapRead
                  << ", memlen : " << memlen << "\n";
        std::cerr << "prevMemEnd_read : " << prevMemEnd_read
                  << ", currMemStart_read : " << currMemStart_read << "\n";
        std::cerr << "prevMemEnd_ref  : " << prevMemEnd_ref
                  << ", currMemStart_ref  : " << currMemStart_ref << "\n";
        std::exit(1);
      }

      alignmentScore += score;
      logger_->debug("\t\t\tscore_sofar: {}", alignmentScore);
      logger_->debug("\t\t\tcigar_sofar: {}", cigarGen.get_cigar());

      logger_->debug("\t\t### ADD MEM SCORE");
      auto printView = readView.substr(currMemStart_read, memlen);
      auto refView = nonstd::string_view(tseq.c_str() + tpos - refStart, memlen);
      // logger_->debug("\t\t\tread [{}], pos: {}, len: {}, ori: {}", printView, currMemStart_read, memlen, (isFw ? "FW" : "RC"));
      // logger_->debug("\t\t\tref  [{}], pos: {}, len: {}", refView, currMemStart_ref, memlen);
      if (printView.length() != refView.length()) {
        logger_->debug("\t\t\treadView length != refView length; should not happen!");
        std::cerr << "[ERROR in PuffAligner, between-MEM alignment]: "
                  << "readView length != refView length; should not happen!" << "\t";
        std::exit(1);
      }
      // finally score the mem match
      alignmentScore += mopts.matchScore * memlen;
      cigarGen.add_item(memlen, 'M');

      prevMemEnd_read = currMemStart_read + memlen - 1;
      prevMemEnd_ref = tpos + memlen - 1;

      if (!alignable(prevMemEnd_read + 1, mopts.matchScore, alignmentScore)) {
        hctr.skippedAlignments_notAlignable += 1;
        ez.stopped = 1;
        logger_->debug("ez stopped");
      }
    }
    
    logger_->debug("\t\t\tscore_sofar: {}", alignmentScore);
    logger_->debug("\t\t\tcigar_sofar: {}", cigarGen.get_cigar());

    // If we got to the end, and there is a read gap left, then align that as
    // well
    logger_->debug("\t\tprevMemEnd_read : {}, readLen : {}", prevMemEnd_read, readLen);
    bool gapAtEnd = (prevMemEnd_read + 1) <= (static_cast<int32_t>(readLen) - 1);
    if (gapAtEnd and !ez.stopped) {
      logger_->debug("\t\t### EXTEND RIGHT");
      int32_t gapRead = (readLen - 1) - (prevMemEnd_read + 1) + 1;
      int32_t refTailStart = prevMemEnd_ref + 1;
      bandwidth = maxAllowedGaps(prevMemEnd_read + 1, alignmentScore) + 1;
      int32_t dataDepBuff = bandwidth; // std::min(refExtLength, 5 * gapRead);
      // refTailEnd is not included in the target sequence
      int32_t refTailEnd = refTailStart + gapRead + dataDepBuff;
      //overhangingEnd =  (refTailStart + gapRead + dataDepBuff) > refTotalLength;
      if (refTailEnd > refTotalLength) {
        refTailEnd = refTotalLength;
      }
      int32_t refLen =
          (refTailEnd > refTailStart) ? refTailEnd - refTailStart: 0;
      auto readWindow = readView.substr(prevMemEnd_read + 1).to_string();
      fillRefSeqBuffer(allRefSeq, refAccPos, refTailStart, refLen, refSeqBuffer_);

      int32_t numSoftClipped = 0; 

      logger_->debug("\t\t\tgapRead : {}, refLen : {}, refBuffer_.size() : {}, refTotalLength : {}", gapRead, refLen, refSeqBuffer_.size(), refTotalLength);
      if (refLen > 0) {
        logger_->debug("\t\t\tCASE 1: some reference bases left to align");
        logger_->debug("\t\t\tread: [{}]", readWindow);
        logger_->debug("\t\t\tref:  [{}]", refSeqBuffer_);

        logger_->debug("\t\t\tksw2_parameters: bandwidth={}, end_bonus={}, zdrop={}", bandwidth, aligner_config.end_bonus, aligner_config.dropoff);
        auto cutoff = minAcceptedScore - alignmentScore - mopts.matchScore * readWindow.length();
        aligner(readWindow.data(), readWindow.length(), refSeqBuffer_.data(),
                refLen, &ez, cutoff, remainedSoftClipLenGeneral,
                ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
        logger_->debug("\t\t\tksw2_results:");
        logger_->debug("\t\t\t\tmax={}, max_q={}, max_t={}", ez.max, ez.max_q, ez.max_t);
        logger_->debug("\t\t\t\tmqe={}, mqe_t={}", ez.mqe, ez.mqe_t);
        logger_->debug("\t\t\t\tmte={}, mte_q={}", ez.mte, ez.mte_q);
        logger_->debug("\t\t\t\tstopped={}, reach_end={}, zdropped={}, score={}", ez.stopped == 1, ez.reach_end == 1, ez.zdropped == 1, ez.score);
        // FOR DEBUGGING PURPOSES ONLY
        // if(mopts.computeCIGAR) {
        //   std::string tmp_cg_exp;
        //   logger_->debug("\t\t\t\talnCigar: {}", ksw2_to_cigar(readWindow.c_str(), refSeqBuffer_.c_str(), &ez, tmp_cg_exp));
        //   logger_->debug("\t\t\t\talnCigar: {}", collapse_cigar(tmp_cg_exp));
        // } 
        if (ez.stopped) hctr.skippedAlignments_notAlignable += 1;
        if (ez.mqe != KSW_NEG_INF) ez.stopped = 0;

        decltype(alignmentScore) part_score = ez.max;
        
        numSoftClipped = readWindow.length() - (ez.max_q + 1);
        addCigar(cigarGen, ez, false);

        if(remainedSoftClipLenGeneral < numSoftClipped || ez.mqe + aligner_config.end_bonus > ez.max)
        {
          part_score = ez.mqe;
          numSoftClipped = 0;
        }
        else
        {
          part_score = ez.max;
          cigarGen.add_item(numSoftClipped, 'S');
        }
        remainedSoftClipLenGeneral -= numSoftClipped;
        remainedSoftClipLenOverhang -= numSoftClipped;
        logger_->debug("\t\t\t\tremainedSoftClipLen={}->{};", remainedSoftClipLenGeneral + numSoftClipped, remainedSoftClipLenGeneral);

        alignmentScore += part_score;

        // NOTE: pre soft-clip code for adjusting the alignment score.
        // int32_t alnCost = allowOverhangSoftclip ? std::max(ez.mqe, ez.mte)
        // : ez.mqe; int32_t delCost = (-1 * mopts.gapOpenPenalty + -1 *
        // mopts.gapExtendPenalty * readWindow.length()); alignmentScore +=
        // std::max(alnCost, delCost);
      } else {
        logger_->debug("\t\t\tCASE 2: no reference bases left (end)");
        //overhangingEnd = true;
        // do any special soft clipping penalty here if we want
        numSoftClipped = 
            mopts.allowSoftclip
            ? readWindow.length()
            : 0;
        alignmentScore +=
            mopts.allowSoftclip && remainedSoftClipLenOverhang >= numSoftClipped
                ? 0
                : (-1 * mopts.gapOpenPenalty +
                    -1 * mopts.gapExtendPenalty * readWindow.length());

        char cigar_char;
        if (mopts.allowSoftclip && remainedSoftClipLenOverhang >= numSoftClipped) {
          cigar_char = 'S';
        } else {
          cigar_char = 'I';
          numSoftClipped = 0;
        }

        if (mopts.computeCIGAR) {
          cigarGen.add_item(readWindow.length(), cigar_char);
        }
      }
      arOut.softclip_end = numSoftClipped;
      logger_->debug("\t\t\tscore_sofar: {}", alignmentScore);
      logger_->debug("\t\t\tcigar_sofar: {}", cigarGen.get_cigar());
    }
  }

  std::string cigar = "";
  bool cigar_fixed{false};
  if (mopts.computeCIGAR) { cigar = cigarGen.get_cigar(); }
  if (cigar_fixed) { hctr.cigar_fixed_count++; }
  if (useAlnCache) { // don't bother to fill up a cache unless this is a multi-mapping read
    //mopts.useAlignmentCache and isMultimapping_ and !perfectChain and !overhangingEnd) { // don't bother to fill up a cache unless this is a multi-mapping read
    if (!didHash) {
      // We want the alignment cache to be on the hash of the full underlying reference sequence.
      // If we are using fullAlignment, this is in refSeqBuffer_, but if we are using between-mem alignment
      // then refSeqBuffer_ could have been used to store shorter portions of the reference during
      // the alignment procedure.  In that case, get the original reference sequence from tseq, which
      // was copied from the full reference sequence in the beginning of the function.
      char* ptr = const_cast<char*>((doFullAlignment ? refSeqBuffer_.data() : tseq.data()));
      MetroHash64::Hash(reinterpret_cast<uint8_t *>(ptr), keyLen, reinterpret_cast<uint8_t *>(&hashKey), 0);
      hashKey ^= queryChainHash;
    }
    AlignmentResult aln;
    aln.isFw = isFw;
    aln.score = alignmentScore;
    if (mopts.computeCIGAR) { aln.cigar = cigar; }
    if (mopts.allowSoftclip) {
      aln.softclip_start = arOut.softclip_start;
      aln.softclip_end = arOut.softclip_end;
    }
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
                                         HitCounters& hctr, bool isMultimapping, std::shared_ptr<spdlog::logger>& logger_, bool verbose) {
  isMultimapping_ = isMultimapping;
    auto tid = jointHit.tid;
    double optFrac{mopts.minScoreFraction};
    auto threshold = [&, optFrac] (uint64_t len) -> int32_t {
         return static_cast<int32_t>(std::floor((!mopts.matchScore)?(-0.6+-0.6*len):optFrac*mopts.matchScore*len));
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
        ar_orphan.softclip_start = 0;
        ar_orphan.softclip_end = 0;
        alignRead(read_orphan, rc_orphan, jointHit.orphanClust()->mems, jointHit.orphanClust()->queryChainHash,
                  jointHit.orphanClust()->perfectChain,
                  jointHit.orphanClust()->isFw, tid, orphan_aln_cache, hctr, ar_orphan, logger_, verbose);
        auto orphan_total_softclip_len = ar_orphan.softclip_start + ar_orphan.softclip_end;
        jointHit.alignmentScore =
          ar_orphan.score > threshold(read_orphan.length() - orphan_total_softclip_len)  ? ar_orphan.score : invalidScore;
        jointHit.orphanClust()->cigar = (mopts.computeCIGAR) ? ar_orphan.cigar : "";
        jointHit.orphanClust()->openGapLen = ar_orphan.openGapLen;
        jointHit.orphanClust()->softClipStart = ar_orphan.softclip_start;
//        jointHit.orphanClust()->coverage = jointHit.alignmentScore;
        if (jointHit.alignmentScore < 0 and verbose) {
          std::cerr << read_orphan.length() << " " << threshold(read_orphan.length() - orphan_total_softclip_len) << " " << ar_left.score << "\n";
        }
        return jointHit.alignmentScore;
    } else {
        hctr.totalAlignmentAttempts += 2;
        ar_left.score = ar_right.score = invalidScore;
        ar_left.softclip_start = ar_right.softclip_start = 0;
        ar_left.softclip_end = ar_right.softclip_end = 0;
        if (verbose) { std::cerr << "left\n"; }
        alignRead(read_left, read_left_rc_, jointHit.leftClust->mems, jointHit.leftClust->queryChainHash, jointHit.leftClust->perfectChain,
                                            jointHit.leftClust->isFw, tid, alnCacheLeft, hctr, ar_left, logger_, verbose);
        if (verbose) { std::cerr << "right\n"; }
        alignRead(read_right, read_right_rc_, jointHit.rightClust->mems, jointHit.rightClust->queryChainHash, jointHit.rightClust->perfectChain,
                                             jointHit.rightClust->isFw, tid, alnCacheRight, hctr, ar_right, logger_, verbose);

        auto left_total_softclip_len = ar_left.softclip_start + ar_left.softclip_end;
        jointHit.alignmentScore = ar_left.score > threshold(read_left.length() - left_total_softclip_len) ? ar_left.score : invalidScore;
        auto right_total_softclip_len = ar_right.softclip_start + ar_right.softclip_end;
        jointHit.mateAlignmentScore = ar_right.score > threshold(read_right.length() - right_total_softclip_len) ? ar_right.score : invalidScore;
/*
        jointHit.alignmentScore = (score_left == invalidScore or score_right == invalidScore)?
                                  invalidScore : score_left + score_right;
*/
//        jointHit.leftClust->coverage = score_left;
        jointHit.leftClust->openGapLen = ar_left.openGapLen;
        jointHit.leftClust->softClipStart = ar_left.softclip_start;

//        jointHit.rightClust->coverage = score_right;
        jointHit.rightClust->openGapLen = ar_right.openGapLen;
        jointHit.rightClust->softClipStart = ar_right.softclip_start;

        if (mopts.computeCIGAR) {
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
int32_t PuffAligner::calculateAlignments(std::string& read, pufferfish::util::JointMems& jointHit, HitCounters& hctr, bool isMultimapping, std::shared_ptr<spdlog::logger>& logger_, bool verbose) {
  isMultimapping_ = isMultimapping;
    auto tid = jointHit.tid;
    double optFrac{mopts.minScoreFraction};
    // bool computeCIGAR = !(aligner.config().flag & KSW_EZ_SCORE_ONLY);
    // bool approximateCIGAR = mopts.alignmentMode == pufferfish::util::PuffAlignmentMode::APPROXIMATE_CIGAR;
    auto threshold = [&, optFrac] (uint64_t len) -> int32_t {
        return static_cast<int32_t>(std::floor((!mopts.matchScore)?(-0.6+-0.6*len):optFrac*mopts.matchScore*len));
    };
    constexpr const auto invalidScore = std::numeric_limits<decltype(ar_left.score)>::min();

    hctr.totalAlignmentAttempts += 1;
    ar_left.score = invalidScore;
    ar_left.softclip_start = 0;
    ar_left.softclip_end = 0;
    const auto& oc = jointHit.orphanClust();
    alignRead(read, read_left_rc_, oc->mems, oc->queryChainHash, oc->perfectChain, oc->isFw, tid, alnCacheLeft, hctr, ar_left, logger_, verbose);
    auto total_softclip_len = ar_left.softclip_start + ar_left.softclip_end;
    jointHit.alignmentScore =
      ar_left.score > threshold(read.length() - total_softclip_len)  ? ar_left.score : invalidScore;
    jointHit.orphanClust()->cigar = (mopts.computeCIGAR) ? ar_left.cigar : "";
    jointHit.orphanClust()->openGapLen = ar_left.openGapLen;
    //        jointHit.orphanClust()->coverage = jointHit.alignmentScore;
    if (jointHit.alignmentScore < 0 and verbose) {
      std::cerr << read.length() << " " << threshold(read.length() - total_softclip_len) << " " << ar_left.score << "\n";
    }
    return jointHit.alignmentScore;
}

bool PuffAligner::recoverSingleOrphan(std::string& read_left, std::string& read_right, pufferfish::util::MemCluster& clust, std::vector<pufferfish::util::MemCluster> &recoveredMemClusters, uint32_t tid, bool anchorIsLeft, bool verbose) {
  int32_t anchorLen = anchorIsLeft ? read_left.length() : read_right.length();
  /*auto tpos = clust.mems[0].tpos;
  auto anchorStart = clust.mems[0].isFw ? clust.mems[0].rpos : anchorLen - (clust.mems[0].rpos + clust.mems[0].extendedlen);
  uint32_t anchorPos = tpos >= anchorStart ? tpos - anchorStart : 0;*/
  uint32_t anchorPos = clust.approxReadStartPos();

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

  int32_t maxDistRight = std::max(1, l2 / 4);
  int32_t maxDistLeft = std::max(1, l1 / 4);
  constexpr const int32_t signedZero{0};

  bool noDovetail = mopts.noDovetail;

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
    windowLength = std::min(static_cast<int32_t>(mopts.maxFragmentLength), static_cast<int32_t>(refLength - startPos));
    if (!noDovetail) {
      startPos = std::max(signedZero, static_cast<int32_t>(anchorPos - mopts.maxFragmentLength));
      windowLength = std::min(2*static_cast<int32_t>(mopts.maxFragmentLength), static_cast<int32_t>(refLength - startPos));
    }
  } else {
    rptr = otherRead;
    rlen = otherLen;
    startPos = std::max(signedZero, static_cast<int32_t>(anchorPos + anchorLen - mopts.maxFragmentLength));
    windowLength = std::min(static_cast<int32_t>(mopts.maxFragmentLength),  static_cast<int32_t>(anchorPos + anchorLen));
    if (!noDovetail) {
      windowLength = std::min(2*static_cast<int32_t>(mopts.maxFragmentLength), static_cast<int32_t>(refLength - startPos));
    }
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
    if (noDovetail and (recovered_fwd and static_cast<int32_t>(recovered_pos + rlen) > static_cast<int32_t>(anchorPos + anchorLen))) {
        edlibFreeAlignResult(result);
        return false;
    }
    recoveredMemClusters.push_back(pufferfish::util::MemCluster(recovered_fwd, rlen));
    auto it = recoveredMemClusters.begin() + recoveredMemClusters.size() - 1;
    if (verbose) {
      std::cerr<< anchorIsLeft << " " << anchorFwd << " " <<  anchorPos << " " << startPos + result.startLocations[0] <<" " << result.editDistance << "\n";
    }
    orphanRecoveryMemCollection.push_back(pufferfish::util::UniMemInfo());
    auto memItr = orphanRecoveryMemCollection.begin() + orphanRecoveryMemCollection.size() - 1;
    if (verbose) { std::cerr<<recovered_fwd << " " << orphanRecoveryMemCollection.size()<<"\n"; }
    it->addMem(memItr, recovered_pos, 1, recovered_fwd ? 0 : rlen-1, recovered_fwd);

    //delete windowSeq;
    edlibFreeAlignResult(result);
    return true;
  } else {
    //delete windowSeq;
    edlibFreeAlignResult(result);
    return false;
  }
}
