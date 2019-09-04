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

/**
 *  Align the read `original_read`, whose mems consist of `mems` against the index and return the result
 *  in `arOut`.  How the alignment is computed (i.e. full vs between-mem and CIGAR vs. score only) depends
 *  on the parameters of how this PuffAligner object was constructed.
 **/
bool PuffAligner::alignRead(std::string& read, std::string& read_rc, const std::vector<pufferfish::util::MemInfo> &mems, bool perfectChain,
                            bool isFw, size_t tid, AlnCacheMap &alnCache, HitCounters &hctr, AlignmentResult& arOut, bool verbose) {

  int32_t alignmentScore{std::numeric_limits<decltype(arOut.score)>::min()};
  if (mems.empty()) {
    arOut.score = alignmentScore;
    return false;
  }

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
  memset(&ez, 0, sizeof(ksw_extz_t));

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
    std::cerr << "[ERROR in PuffAligner::alignRead :] currHitStart_read is invalid; this hould not happen!\n";
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
  uint32_t buff{20};

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

  // If we are only aligning between MEMs
  if (!doFullAlignment) {
    refStart = (currHitStart_ref >= currHitStart_read) ? currHitStart_ref - currHitStart_read : 0;
    //keyLen = (refStart + readLen < refTotalLength) ? readLen : refTotalLength - refStart;
    keyLen = (refStart + readLen + buff < refTotalLength) ? readLen + buff : refTotalLength - refStart;
  } else { // we are aligning from the start of the read
    // If the first hit starts further in the reference than in the
    // read, then we align from the beginning of the read and (ref_start - read_start) on
    // the reference.
    if (currHitStart_ref > currHitStart_read) {
      refStart = currHitStart_ref - currHitStart_read;
      readStart = 0;
    } else if (currHitStart_ref < currHitStart_read) {
      // If the first his starts further in the read than in the reference, than
      // the read overhangs the beginning of the reference and we start aligning
      // from the beginning of the reference and from position (read_start - ref_start)
      // in the read.
      readStart = currHitStart_read - currHitStart_ref;
      refStart = 0;
    } else {
      // If the first hit starts at the same position in the reference and the read
      // ... what is this case?
      readStart = 0;
      refStart = currHitStart_ref - currHitStart_read;
    }
    keyLen = (refStart + readLen + buff < refTotalLength) ? readLen + buff : refTotalLength - refStart;
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
    if (hit != alnCache.end() ) {//}and refStart + readLen + refExtLength < refTotalLength) {
      hctr.skippedAlignments_byCache += 1;
      arOut.score /*= alignmentScore*/ = hit->second.score;
      if (computeCIGAR) { arOut.cigar = hit->second.cigar; }
      arOut.openGapLen = hit->second.openGapLen;
      return true;
    }
  }

  //auto logger_ = spdlog::get("console");
  //spdlog::set_level(spdlog::level::debug); // Set global log level to debug
  //logger_->set_pattern("%v");

  // @mohsen & @fataltes --- we should figure out how to
  // avoid computing the rc of a read if we've already done it.
  if (!isFw and read_rc.empty()) { read_rc = pufferfish::util::reverseComplement(read); }
  nonstd::string_view readView = (isFw) ? read : read_rc;
  /*
    std::string read = original_read;
    if (!isfw) {
    read = pufferfish::util::reversecomplement(read);
    }
  */

  if (perfectChain) {
    arOut.score /*= alignment*/ = alignmentScore = readLen * mopts.matchScore;
    if (computeCIGAR) { cigarGen.add_item(readLen, 'M'); }
    hctr.skippedAlignments_byCov += 1;
    SPDLOG_DEBUG(logger_,"[[");
    SPDLOG_DEBUG(logger_,"read sequence ({}) : {}", (isFw ? "FW" : "RC"), readView);
    SPDLOG_DEBUG(logger_,"ref  sequence      : {}", (doFullAlignment ? tseq : refSeqBuffer_));
    SPDLOG_DEBUG(logger_,"perfect chain!\n]]\n");
  } else if (doFullAlignment) {
    // if we allow softclipping of overhanging bases, then we can cut off the part of the read
    // before the start of the reference
    decltype(readStart) readOffset = allowOverhangSoftclip ? readStart : 0;
    nonstd::string_view readSeq = readView.substr(readOffset);
    aligner(readSeq.data(), readSeq.length(), refSeqBuffer_.data(), refSeqBuffer_.length(), &ez,
            ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
    // if we allow softclipping of overhaning bases, then we only care about the best
    // score to the end of the query or the end of the reference.  Otherwise, we care about
    // the best score all the way until the end of the query.
    alignmentScore = allowOverhangSoftclip ? std::max(ez.mqe, ez.mte) : ez.mqe;

    SPDLOG_DEBUG(logger_,"readSeq : {}\nrefSeq  : {}\nscore   : {}\nreadStart : {}", readSeq, refSeqBuffer_, alignmentScore, readStart);
    SPDLOG_DEBUG(logger_,"currHitStart_read : {}, currHitStart_ref : {}\nmqe : {}, mte : {}\n", currHitStart_read, currHitStart_ref, ez.mqe, ez.mte);

    if (computeCIGAR) { openGapLen = addCigar(cigarGen, ez, false); }
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

      int32_t refWindowStart = (readStartPosOnRef - refExtLength) > 0 ? (readStartPosOnRef - refExtLength) : 0;
      int32_t refWindowLength = tpos - refWindowStart;
      fillRefSeqBufferReverse(allRefSeq, refAccPos, refWindowStart, refWindowLength, refSeqBuffer_);

      if (refSeqBuffer_.length() > 0) {
        auto readWindow = readView.substr(0, firstMemStart_read).to_string();
        std::reverse(readWindow.begin(), readWindow.end());
        SPDLOG_DEBUG(logger_,"PRE:\nreadStartPosOnRef : {}\nrefWindowStart : {}", readStartPosOnRef, refWindowStart);
        SPDLOG_DEBUG(logger_,"refWindowLength : {}\nread : [{}]\nref : [{}]", refWindowLength, readWindow, refSeqBuffer_);
        ksw_reset_extz(&ez);
        aligner(readWindow.data(), readWindow.length(), refSeqBuffer_.data(), refSeqBuffer_.length(), &ez,
                ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
        alignmentScore += allowOverhangSoftclip ? std::max(ez.mqe, ez.mte) : ez.mqe;
        openGapLen = addCigar(cigarGen, ez, true);
        SPDLOG_DEBUG(logger_,"score : {}", std::max(ez.mqe, ez.mte));
      } else {
        // do any special soft clipping penalty here if we want
        alignmentScore += allowOverhangSoftclip ? 0 :
          (-1 * mopts.gapOpenPenalty + -1 * mopts.gapExtendPenalty * firstMemStart_read);
      }
    }

    int32_t prevMemEnd_read = isFw ? rpos : readLen - (rpos+memlen);
    int32_t prevMemEnd_ref = tpos;

    SPDLOG_DEBUG(logger_,"\t Aligning through MEM chain : ");
    // for the second through the last mem
    for(auto it = mems.begin(); it != mems.end(); ++it) {
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
        score += (-1 * mopts.gapOpenPenalty + -1 * mopts.gapExtendPenalty * gapDiff);
        if (gapRead < 0) { // subtract off extra matches
          score += mopts.matchScore * gapRead;
        }
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

        score += aligner(readWindow.data(), readWindow.length(), refSeq1, gapRef, &ez,
                        ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::GLOBAL>());
        addCigar(cigarGen, ez, false);
      } else if ( it > mems.begin() and ((currMemStart_read <= prevMemEnd_read) or (currMemStart_ref <= prevMemEnd_ref)) ){
        SPDLOG_DEBUG(logger_,"]]\n");
        std::cerr << "[ERROR in PuffAligner, between-MEM alignment] : Improperly compacted MEM chain.  Should not happen!\n";
        std::cerr << "gapRef : " << gapRef << ", gapRead : " << gapRead << ", memlen : " << memlen << "\n";
        std::cerr << "prevMemEnd_read : " << prevMemEnd_read << ", currMemStart_read : " << currMemStart_read << "\n";
        std::cerr << "prevMemEnd_ref  : " << prevMemEnd_ref <<  ", currMemStart_ref  : " << currMemStart_ref << "\n";
        std::exit(1);
      }

      SPDLOG_DEBUG(logger_,"\t MEM (rpos : {}, memlen : {}, tpos : {})", rpos, memlen, tpos);
      SPDLOG_DEBUG(logger_,"\t gapRef : {}, gapRead : {}", gapRef, gapRead);
      auto printView = readView.substr(currMemStart_read, memlen);
      auto refView = nonstd::string_view(tseq.c_str() + tpos - refStart, memlen);

      SPDLOG_DEBUG(logger_,"\t read [{}], pos : {}, len : {}, ori : {}", printView, currMemStart_read, memlen, (isFw ? "FW" : "RC"));
      SPDLOG_DEBUG(logger_,"\t ref  [{}], pos : {}, len : {}", refView, currMemStart_ref, memlen);
      if (printView.length() != refView.length()) {
        SPDLOG_DEBUG(logger_,"\t readView length != refView length; should not happen!");
        std::exit(1);
      }

      prevMemEnd_read = currMemStart_read + memlen - 1;
      prevMemEnd_ref = tpos + memlen - 1;
      alignmentScore += score;
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
        aligner(readWindow.data(), readWindow.length(), refSeqBuffer_.data(), refLen, &ez,
                ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
        int32_t alnCost = allowOverhangSoftclip ? std::max(ez.mqe, ez.mte) : ez.mqe;
        int32_t delCost = (-1 * mopts.gapOpenPenalty + -1 * mopts.gapExtendPenalty * readWindow.length());
        alignmentScore += std::max(alnCost, delCost);
        SPDLOG_DEBUG(logger_,"POST score : {}", std::max(ez.mqe, ez.mte));
      } else {
        // do any special soft-clipping penalty here if we want
        // do any special soft clipping penalty here if we want
        alignmentScore += allowOverhangSoftclip ? 0 :
          (-1 * mopts.gapOpenPenalty + -1 * mopts.gapExtendPenalty * readWindow.length());
      }
    }

    SPDLOG_DEBUG(logger_,"score : {}\n]]\n", alignmentScore);

    /*
    alignmentScore = 0;
    for (auto&& mem : mems) {
      rpos = mem.rpos;
      memlen = mem.extendedlen;
      tpos = mem.tpos;

      // To work around seg fault from mems
      // TODO: What?  why should we need a special case here?
      //if (!firstMem and tpos + memlen - refStart > readLen) {
      //  break;
      //}

      int score = 0;
      currHitStart_read = isFw ? rpos : readLen - (rpos + memlen);
      currHitStart_ref = tpos;

      if (verbose) { std::cerr << "currHitStart_read:" << currHitStart_read << "\tcurrHitStart_ref:" << currHitStart_ref << "\n"; }
      if (verbose) { std::cerr << "lastHitStart_read:" << lastHitEnd_read << "\tlastHitStart_ref:" << lastHitEnd_ref << "\n"; }

      if (firstMem) {
        lastHitEnd_ref =
          currHitStart_ref - currHitStart_read >= 0 ? currHitStart_ref - currHitStart_read - 1 : -1;
      } else {
        if (currHitStart_read <= lastHitEnd_read - (int32_t) k + 1) {
          if (verbose) { std::cerr << "The mems on read and reference are not consistent(1)\t" << read << "\n"; }
          break;
        }
        if (currHitStart_ref <= lastHitEnd_ref - (int32_t) k + 1) {
          if (verbose) { std::cerr << "The mems on read and reference are not consistent(2)\t" << read << "\n"; }
          break;
        }
        if (currHitStart_ref > lastHitEnd_ref + (int32_t) readLen) {
          if (verbose) { std::cerr << "The mems on read and reference are not consistent(3)\t" << read << "\n"; }
          break;
        }
        if (currHitStart_read <= lastHitEnd_read and
            lastHitEnd_read - currHitStart_read != lastHitEnd_ref - currHitStart_ref) {
          if (verbose) { std::cerr << "to align reads with indels\t" << read << "\n"; }
          break;
        }
        if (currHitStart_ref <= lastHitEnd_ref and
            lastHitEnd_read - currHitStart_read != lastHitEnd_ref - currHitStart_ref) {
          if (verbose) { std::cerr << "to align reads with indels\t" << read << "\n"; }
          break;
        }
      }

      
      // // Performing full alignment to validate the scores
      // if (verbose and firstMem) {
      // uint32_t refStart, readStart;
      // if (currHitStart_ref > currHitStart_read) {
      // refStart = currHitStart_ref - currHitStart_read;
      // readStart = 0;
      // } else if (currHitStart_ref < currHitStart_read) {
      // readStart = currHitStart_read - currHitStart_ref;
      // refStart = 0;
      // } else {
      // readStart = currHitStart_read;
      // refStart = currHitStart_ref;
      // }
      // auto refLength = (refStart + readLen + buff < refTotalLength) ? readLen + buff : refTotalLength -
      // refStart;

      // tseq = getRefSeq(allRefSeq, refAccPos, refStart, refLength);
      // auto readSeq = readStart > 0 ? extractReadSeq(read, readStart, readLen, 1) : read;

      // aligner(readSeq.c_str(), readSeq.length(), tseq.c_str(), tseq.length(), &ez,
      // ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
      // alignment = ez.mqe; //std::max(ez.mqe, ez.mte);

      // {
      // std::cerr << "Original read seq:\t" << original_read << "\n";
      // std::cerr << "Total alignment with the score\t" << alignment << "\t from position\t"
      // << readStart << "\t on the read:\n" << readSeq << "\n";
      // std::cerr << "At\t" << refStart << "\t for length of\t" << refLength << "\t on reference:\n"
      // << tseq << "\n";
      // for (int i = 0; i < ez.n_cigar; ++i) // print CIGAR
      // std::cerr<< (ez.cigar[i]>>4) << "MID"[ez.cigar[i] & 0xf];
      // std::cerr << "\n";
      // }
      // }
      

      // If there is a gap in the read or reference or both
      bool gapAligned = false;
      if (currHitStart_read > (lastHitEnd_read + 1) or currHitStart_ref > (lastHitEnd_ref + 1)) {
        gapAligned = true;
        auto readGapLength = currHitStart_read - lastHitEnd_read - 1;
        auto refGapLength = currHitStart_ref - lastHitEnd_ref - 1;

        if (verbose) { std::cerr << "readGapLength " << readGapLength << "\trefGapLength " << refGapLength << "\n"; }

        // get the next segment of the read sequence we wish to align
        auto readSeq = readView.substr(lastHitEnd_read + 1, readGapLength).to_string();

        if (verbose) {
          for (size_t x = 0; x < 50; ++x) {
            std::cerr << "+";
          }
          std::cerr << " start gap alignment\n";
          std::cerr << "read from " << lastHitEnd_read + 1 << " for length " << readGapLength
                    << " of the read:\t " << readSeq << "\n";
          std::cerr << "currHitStart_read: " << currHitStart_read << " lastHitEnd_read: "
                    << lastHitEnd_read << " currHitStart_ref: " << currHitStart_ref << " lastHitEnd_ref: "
                    << lastHitEnd_ref << "\n";
          std::cerr << mems[0].rpos << "\t" << readLen << "\t" << rpos << "\t" << memlen << "\n";
        }
          
        if (firstMem) {
          ///// no need because of line 512: lastHitEnd_ref = currHitStart_ref > currHitStart_read ? currHitStart_ref - currHitStart_read -1: -1;
          // Not extending the reference beyond its beginning
          auto refStartPos = lastHitEnd_ref > refExtLength ? lastHitEnd_ref + 1 - refExtLength : 0;
          refGapLength = currHitStart_ref - refStartPos;

          // SOFT CLIPPING for reads mapping with a left hanger at the beginning of a transcript
          if (refGapLength > 0) {
            // We reverse the strings because of ksw forces alignment from the beginning of the sequences
            std::reverse(readSeq.begin(), readSeq.end());
            fillRefSeqBufferReverse(allRefSeq, refAccPos, refStartPos, refGapLength, refSeqBuffer_);
            ksw_reset_extz(&ez);
            aligner(readSeq.data(), readSeq.length(), refSeqBuffer_.data(), refSeqBuffer_.length(), &ez,
                    ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
            // TODO : @fataltes & @mohsen --- I changed from the line below to use the max.
            // is there a reason we should not be using max? Please check this.
            //score = ez.mqe;
            score = std::max(ez.mqe, ez.mte);
            openGapLen = addCigar(cigarGen, ez, true);

            if (verbose) {
              std::cerr << "It's a gap at the beginning, so the reverse read from "
                        << lastHitEnd_read + 1 << " is\t" << readSeq << "\n";
              std::cerr << "The reverse reference from " << refStartPos << " with the length "
                        << refGapLength << " is\t" << tseq << "\n";
              std::cerr << "mqe: " << ez.mqe << "\t mte: " << ez.mte << "\n";
              for (int i = 0; i < ez.n_cigar; ++i) // print CIGAR
                std::cerr<< (ez.cigar[i]>>4) << "MID"[ez.cigar[i] & 0xf];
              std::cerr << "\n" << score << "\n";
            }

            int32_t insertionDeletion = 0; //TODO should be fixed for cigars
            (void) insertionDeletion;
          } else {
            score = 0;
            int32_t startHang = readGapLength - refGapLength;
            alignmentScore +=
              (-1) * mopts->gapOpenPenalty + (-1) * mopts->gapExtendPenalty * (startHang);
            cigarGen.add_item(startHang, 'I');
          }
        } else {
          //Considering the case where the gap in either the read or reference is empty
          if (readGapLength <= 0) {
            if (verbose) {
              std::cerr << "read gap lengths was of length 0 for the ref gap length of \t"
                        << refGapLength << "\n";
            }
            score = (-1) * mopts->gapOpenPenalty +
              (-1) * mopts->gapExtendPenalty * (refGapLength + std::abs(readGapLength));
            cigarGen.add_item(refGapLength, 'D');
          } else if (refGapLength <= 0) {
            if (verbose) {
              std::cerr << "ref gap lengths was of length 0 for the read gap length of \t"
                        << readGapLength << "\n";
            }
            score = (-1) * mopts->gapOpenPenalty +
              (-1) * mopts->gapExtendPenalty * (readGapLength + std::abs(refGapLength));
            cigarGen.add_item(readGapLength, 'I');
          } else {
            if (lastHitEnd_ref + 1 - refStart < 0 or
                lastHitEnd_ref + 1 - refStart >= (originalRefSeqLen + buff)) {
              std::cerr << "Should not happen: lastHitEnd_ref is " << lastHitEnd_ref
                        << "and refStart is " << refStart << ", but refSeq length is "
                        << originalRefSeqLen << "\n";
            }
            //char *refSeq1 = refSeq.get() + lastHitEnd_ref + 1 - refStart;
            const char* refSeq1 = tseq.data() + (lastHitEnd_ref + 1 - refStart);
            score = aligner(readSeq.c_str(), readSeq.length(), refSeq1, refGapLength, &ez,
                            ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::GLOBAL>());
            addCigar(cigarGen, ez, false);
            if (verbose) {
              tseq = getRefSeq(allRefSeq, refAccPos, lastHitEnd_ref + 1, refGapLength);
              std::stringstream out;

              out << "[[\n";
              out << "read from " << lastHitEnd_read + 1 << " with the gap length of "
                        << readGapLength << " is \t" << readSeq << "\n";
              out << "reference from " << lastHitEnd_ref + 1 << " with the length "
                        << refGapLength << " is \t" << nonstd::string_view(refSeq1, refGapLength) << "\n";
              out << "score is " << score << "\n]]";
              std::cerr << out.str() << "\n";
              //for (int i = 0; i < ez.n_cigar; ++i) { // print CIGAR 
              //  std::cerr<< (ez.cigar[i]>>4) << "MID"[ez.cigar[i] & 0xf];
              //}
              std::cout << "\n";
             }
          }
        }
        alignmentScore += score;
        if (verbose) {
          std::cerr << isFw << " at " << currHitStart_ref << " on " << tid << " for " << refGapLength
                    << "\n";
          std::cerr << "curr score2:\t" << score << "\n";
          for (size_t x = 0; x < 50; ++x)
            std::cerr << "+";
          std::cerr << " end gap alignments\n";
        }
      }

      int32_t extension_ref = currHitStart_ref + memlen - lastHitEnd_ref;
      int32_t extension_read = currHitStart_read + memlen - lastHitEnd_read;

      auto longest_match_extension = std::max(extension_ref, extension_read);
      //currHitStart_ref + memlen - lastHitEnd_ref > currHitStart_read + memlen - lastHitEnd_read ?
      //currHitStart_ref + memlen - lastHitEnd_ref : currHitStart_read + memlen - lastHitEnd_read;

      // @mohsen & @fataltes --- why don't we always get mopts->matchScore * memlen value for the current mem?
      auto num_match = (currHitStart_ref <= lastHitEnd_ref) ? (longest_match_extension - 1) : memlen;
      score = mopts->matchScore * num_match;

      if (lastHitEnd_ref - currHitStart_ref == lastHitEnd_read - currHitStart_read or firstMem) {
        if (verbose) {
          auto readSeq = extractReadSeq(read, currHitStart_read, currHitStart_read + memlen, 1);
          auto tseq1 = getRefSeq(allRefSeq, refAccPos, currHitStart_ref, memlen);
          std::cerr << "read from " << currHitStart_read << "\t with the sequence \n" << readSeq << "\n";
          std::cerr << "orientation of\t" << isFw << " at " << currHitStart_ref << " on reference " << tid
                    << " for the length of " << memlen << "\t with sequence:\n" << tseq1 << "\n";
          std::cerr << "with the score1:\t" << score << "\n";
        }
      } else {
        int32_t gapLength = std::abs(
                                     lastHitEnd_ref - currHitStart_ref - lastHitEnd_read + currHitStart_read);

        int32_t penalty = mopts->gapOpenPenalty + mopts->gapExtendPenalty * gapLength;
        if (!gapAligned) {
          score -= penalty;
          if (currHitStart_ref - lastHitEnd_ref > currHitStart_read - lastHitEnd_read) {
            cigarGen.add_item(gapLength, 'D');
          } else {
            cigarGen.add_item(gapLength, 'I');
          }
        }
        if (verbose) {
          if (!gapAligned) { std::cerr << "penalty\t" << penalty << "\n"; }
          auto readSeq = extractReadSeq(read, currHitStart_read, currHitStart_read + memlen, 1);
          auto tseq1 = getRefSeq(allRefSeq, refAccPos, currHitStart_ref, memlen);
          std::cerr << "read from " << currHitStart_read << "\t with the sequence \n" << readSeq << "\n";
          std::cerr << "orientation of\t" << isFw << " at " << currHitStart_ref << " on reference " << tid
                    << " for the length of " << memlen << "\t with sequence:\n" << tseq1 << "\n";
          std::cerr << "with the score1:\t" << score << "\n";
        }
      }
      alignmentScore += score;
      cigarGen.add_item(num_match, 'M');

      lastHitEnd_read = currHitStart_read + memlen - 1;
      lastHitEnd_ref = tpos + memlen - 1;
      if (lastHitEnd_ref - refStart + 1 > originalRefSeqLen + 1 + buff) {
        std::cerr << "Should not happen: lastHitEnd_ref is " << lastHitEnd_ref << " and refStart is "
                  << refStart << ", but refSeq length is " << originalRefSeqLen << "\n";
      }
      firstMem = false;
    }

    // Try and align any remaining sequence on the read
    if (lastHitEnd_read < static_cast<int32_t>(readLen) - 1) {
      auto readGapLength = readLen - 1 - lastHitEnd_read;
      auto refGapLength = lastHitEnd_ref + 1 + readGapLength + refExtLength < refTotalLength ?
                                                                              readGapLength + refExtLength : refTotalLength - lastHitEnd_ref - 1;

      if (lastHitEnd_ref + readGapLength > refTotalLength) {
        refGapLength = refTotalLength - (lastHitEnd_ref + 1);
        // @mohsen & @fataltes --- I commented out the below way of
        // dealing with this in favor of what is above.  Why would
        // we not want to align to what remains of the ref in this case
        // to see if the best alignment score is reasonable?

        //AlignmentResult aln;
        //arOut.score = -1;
        //if (computeCIGAR) { arOut.cigar = cigar; }
        //arOut.openGapLen = openGapLen;
        //return false;
      }

      // SOFT CLIPPING for reads mapping with a right hanger at the beginning of a transcript
      if (refGapLength != 0) {
        auto readSeq = readView.substr(lastHitEnd_read + 1, readGapLength).to_string();
        auto refStartPos = lastHitEnd_ref + 1;
        fillRefSeqBuffer(allRefSeq, refAccPos, refStartPos, refGapLength, refSeqBuffer_);
        aligner(readSeq.data(), readSeq.length(), refSeqBuffer_.data(), refGapLength, &ez,
                ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
        // TODO : @fataltes & @mohsen --- I changed from the line below to use the max.
        // is there a reason we should not be using max? Please check this.
        //auto score = ez.mqe;
        auto score = std::max(ez.mqe, ez.mte);

        alignmentScore += score;
        auto insertionDeletion = addCigar(cigarGen, ez, false);

        if (verbose) {
          std::cerr << "read from " << lastHitEnd_read + 1 << "\t\t\t" << readSeq << "\n";
          std::cerr << "at " << lastHitEnd_ref << " for " << refGapLength << "\t" << tseq << "\n";
          for (int i = 0; i < ez.n_cigar; ++i) // print CIGAR
            std::cerr<< (ez.cigar[i]>>4) << "MID"[ez.cigar[i] & 0xf];
          std::cerr << "\ncurr score3:\t" << score << "\n";
          std::cerr << ez.mqe << "\t " << ez.mte << "\n";
        }

        if (readGapLength > refGapLength + insertionDeletion) {
          int32_t endHang = readGapLength - refGapLength - insertionDeletion;
          alignmentScore += (-1) * mopts->gapOpenPenalty + (-1) * mopts->gapExtendPenalty * (endHang);
          cigarGen.add_item(endHang, 'I');
        }
      } else {
        int32_t endHang = readGapLength;
        alignmentScore += (-1) * mopts->gapOpenPenalty + (-1) * mopts->gapExtendPenalty * (endHang);
        cigarGen.add_item(endHang, 'I');
      }
    }
    */
  }
  bool cigar_fixed{false};
  if (computeCIGAR) { cigar = cigarGen.get_cigar(readLen, cigar_fixed); }
  if (cigar_fixed) { hctr.cigar_fixed_count++; }
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
                  jointHit.orphanClust()->isFw, tid, orphan_aln_cache, hctr, ar_orphan, verbose);
        jointHit.alignmentScore =
          ar_orphan.score > threshold(read_orphan.length())  ? ar_orphan.score : invalidScore;
        jointHit.orphanClust()->cigar = (computeCIGAR) ? ar_orphan.cigar : "";
        jointHit.orphanClust()->openGapLen = ar_orphan.openGapLen;
//        jointHit.orphanClust()->coverage = jointHit.alignmentScore;
        if (jointHit.alignmentScore < 0 and verbose) {
          std::cerr << read_orphan.length() << " " << threshold(read_orphan.length()) << " " << ar_left.score << "\n";
        }
        return jointHit.alignmentScore;
    } else {
        hctr.totalAlignmentAttempts += 2;
        ar_left.score = ar_right.score = invalidScore;
        if (verbose) { std::cerr << "left\n"; }
        alignRead(read_left, read_left_rc_, jointHit.leftClust->mems, jointHit.leftClust->perfectChain,
                                            jointHit.leftClust->isFw, tid, alnCacheLeft, hctr, ar_left, verbose);
        if (verbose) { std::cerr << "right\n"; }
        alignRead(read_right, read_right_rc_, jointHit.rightClust->mems, jointHit.rightClust->perfectChain,
                                             jointHit.rightClust->isFw, tid, alnCacheRight, hctr, ar_right, verbose);

        jointHit.alignmentScore = ar_left.score > threshold(read_left.length()) ? ar_left.score : invalidScore;
        jointHit.mateAlignmentScore = ar_right.score > threshold(read_right.length()) ? ar_right.score : invalidScore;
/*
        jointHit.alignmentScore = (score_left == invalidScore or score_right == invalidScore)?
                                  invalidScore : score_left + score_right;
*/
//        jointHit.leftClust->coverage = score_left;
        jointHit.leftClust->openGapLen = ar_left.openGapLen;
//        jointHit.rightClust->coverage = score_right;
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
    alignRead(read, read_left_rc_, oc->mems, oc->perfectChain, oc->isFw, tid, alnCacheLeft, hctr, ar_left, verbose);
    jointHit.alignmentScore =
      ar_left.score > threshold(read.length())  ? ar_left.score : invalidScore;
    jointHit.orphanClust()->cigar = (computeCIGAR) ? ar_left.cigar : "";
    jointHit.orphanClust()->openGapLen = ar_left.openGapLen;
    //        jointHit.orphanClust()->coverage = jointHit.alignmentScore;
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
