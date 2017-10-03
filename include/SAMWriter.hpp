#ifndef SAM_WRITER_HPP
#define SAM_WRITER_HPP
#include "PairedAlignmentFormatter.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "Util.hpp"

inline void getSamFlags(const util::QuasiAlignment& qaln, uint16_t& flags) {
  /*
    constexpr uint16_t pairedInSeq = 0x1;
    constexpr uint16_t mappedInProperPair = 0x2;
    constexpr uint16_t unmapped = 0x4;
    constexpr uint16_t mateUnmapped = 0x8;
  */
  constexpr uint16_t isRC = 0x10;
  /*
  constexpr uint16_t mateIsRC = 0x20;
  constexpr uint16_t isRead1 = 0x40;
  constexpr uint16_t isRead2 = 0x80;
  constexpr uint16_t isSecondaryAlignment = 0x100;
  constexpr uint16_t failedQC = 0x200;
  constexpr uint16_t isPCRDup = 0x400;
  constexpr uint16_t supplementaryAln = 0x800;
*/

  flags = 0;
  // Not paired in sequencing
  // flags1 = (peInput) ? pairedInSeq : 0;
  // flags |= properlyAligned;
  // we don't output unmapped yet
  // flags |= unmapped
  // flags |= mateUnmapped
  flags |= (qaln.fwd) ? 0 : isRC;
  // Mate flag meaningless
  // flags1 |= (qaln.mateIsFwd) ? 0 : mateIsRC;
  // flags |= isRead1;
  // flags2 |= isRead2;
}

// get the sam flags for the quasialignment qaln.
// peinput is true if the read is paired in *sequencing*; false otherwise
// the sam flags for mate 1 are written into flags1 and for mate2 into flags2
inline void getSamFlags(const util::QuasiAlignment& qaln, bool peInput,
                        uint16_t& flags1, uint16_t& flags2) {
  constexpr uint16_t pairedInSeq = 0x1;
  constexpr uint16_t properlyAligned = 0x2;
  constexpr uint16_t unmapped = 0x4;
  constexpr uint16_t mateUnmapped = 0x8;
  constexpr uint16_t isRC = 0x10;
  constexpr uint16_t mateIsRC = 0x20;
  constexpr uint16_t isRead1 = 0x40;
  constexpr uint16_t isRead2 = 0x80;
  // constexpr uint16_t isSecondaryAlignment = 0x100;
  // constexpr uint16_t failedQC = 0x200;
  // onstexpr uint16_t isPCRDup = 0x400;
  // constexpr uint16_t supplementaryAln = 0x800;

  flags1 = flags2 = 0;
  flags1 = (peInput) ? pairedInSeq : 0;
  flags1 |= (qaln.isPaired) ? properlyAligned : 0;
  flags2 = flags1;
  // we don't output unmapped yet
  bool read1Unaligned = qaln.mateStatus == util::MateStatus::PAIRED_END_RIGHT;
  bool read2Unaligned = qaln.mateStatus == util::MateStatus::PAIRED_END_LEFT;
  // If read 1 is unaligned, flags1 gets "unmapped" and flags2 gets "mate
  // unmapped"
  flags1 |= (read1Unaligned) ? unmapped : 0;
  flags2 |= (read1Unaligned) ? mateUnmapped : 0;
  // If read 2 is unaligned, flags2 gets "unmapped" and flags1 gets "mate
  // unmapped"
  flags2 |= (read2Unaligned) ? unmapped : 0;
  flags1 |= (read2Unaligned) ? mateUnmapped : 0;

  flags1 |= (qaln.fwd) ? 0 : isRC;
  flags1 |= (qaln.mateIsFwd) ? 0 : mateIsRC;
  flags2 |= (qaln.mateIsFwd) ? 0 : isRC;
  flags2 |= (qaln.fwd) ? 0 : mateIsRC;
  flags1 |= isRead1;
  flags2 |= isRead2;
}

template <typename IndexT>
inline void writeSAMHeader(IndexT& pfi, std::shared_ptr<spdlog::logger> out) {
  fmt::MemoryWriter hd;
  hd.write("@HD\tVN:1.0\tSO:unknown\n");

  auto& txpNames = pfi.getRefNames();
  // auto& txpLens = rmi.txpLens;

  auto numRef = txpNames.size();
  // for now go with constant length
  // TODO read reference information
  // while reading the index
  int txpLen = 1000;
  for (size_t i = 0; i < numRef; ++i) {
    hd.write("@SQ\tSN:{}\tLN:{:d}\n", txpNames[i], txpLen);
  }
  // Eventually output a @PG line
  // some other version nnumber for now,
  // will think about it later
  std::string version = "0.1.0";
  hd.write("@PG\tID:rapmap\tPN:rapmap\tVN:{}\n", version);
  std::string headerStr(hd.str());
  // Don't include the last '\n', since the logger will do it for us.
  headerStr.pop_back();
  out->info(headerStr);
}

// Declarations for functions dealing with SAM formatting and output
//
inline void adjustOverhang(int32_t& pos, uint32_t readLen, uint32_t txpLen,
                           util::FixedWriter& cigarStr) {
  cigarStr.clear();
  int32_t readLenS = static_cast<int32_t>(readLen);
  int32_t txpLenS = static_cast<int32_t>(txpLen);
  if (pos + readLenS < 0) {
    cigarStr.write("{}S", readLen);
    pos = 0;
  } else if (pos < 0) {
    int32_t matchLen = readLenS + pos;
    int32_t clipLen = readLenS - matchLen;
    cigarStr.write("{}S{}M", clipLen, matchLen);
    // Now adjust the mapping position
    pos = 0;
  } else if (pos > txpLenS) {
    cigarStr.write("{}S", readLen);
  } else if (pos + readLenS > txpLenS) {
    int32_t matchLen = txpLenS - pos;
    int32_t clipLen = readLenS - matchLen;
    cigarStr.write("{}M{}S", matchLen, clipLen);
  } else {
    cigarStr.write("{}M", readLen);
  }
}

inline void adjustOverhang(util::QuasiAlignment& qa, uint32_t txpLen,
                           util::FixedWriter& cigarStr1,
                           util::FixedWriter& cigarStr2) {
  if (qa.isPaired) { // both mapped
    adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr1);
    adjustOverhang(qa.matePos, qa.mateLen, txpLen, cigarStr2);
  } else if (qa.mateStatus == util::MateStatus::PAIRED_END_LEFT) {
    // left read mapped
    adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr1);
    // right read un-mapped will just be read length * S
    cigarStr2.clear();
    cigarStr2.write("{}S", qa.mateLen);
  } else if (qa.mateStatus == util::MateStatus::PAIRED_END_RIGHT) {
    // right read mapped
    adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr2);
    // left read un-mapped will just be read length * S
    cigarStr1.clear();
    cigarStr1.write("{}S", qa.readLen);
  }
}

template <typename ReadPairT, typename IndexT>
inline uint32_t writeAlignmentsToStream(
    ReadPairT& r, PairedAlignmentFormatter<IndexT>& formatter,
    std::vector<util::QuasiAlignment>& jointHits, fmt::MemoryWriter& sstream) {

  auto& read1Temp = formatter.read1Temp;
  auto& read2Temp = formatter.read2Temp;
  // auto& qual1Temp = formatter.qual1Temp;
  // auto& qual2Temp = formatter.qual2Temp;
  auto& cigarStr1 = formatter.cigarStr1;
  auto& cigarStr2 = formatter.cigarStr2;

  cigarStr1.clear();
  cigarStr2.clear();
  cigarStr1.write("{}M", r.first.seq.length());
  cigarStr2.write("{}M", r.second.seq.length());
  //std::cerr << cigarStr1.c_str() << "\n";
  uint16_t flags1, flags2;

  auto& readName = r.first.name;
  // If the read name contains multiple space-separated parts,
  // print only the first
  size_t splitPos = readName.find(' ');
  if (splitPos < readName.length()) {
    readName[splitPos] = '\0';
  } else {
    splitPos = readName.length();
  }

  if (splitPos > 2 and readName[splitPos - 2] == '/') {
    readName[splitPos - 2] = '\0';
  }

  auto& mateName = r.second.name;
  // If the read name contains multiple space-separated parts,
  // print only the first
  splitPos = mateName.find(' ');
  if (splitPos < mateName.length()) {
    mateName[splitPos] = '\0';
  } else {
    splitPos = mateName.length();
  }

  // trim /2 from the pe read
  if (splitPos > 2 and mateName[splitPos - 2] == '/') {
    mateName[splitPos - 2] = '\0';
  }

  std::string numHitFlag = fmt::format("NH:i:{}", jointHits.size());
  uint32_t alnCtr{0};
  // uint32_t trueHitCtr{0};
  // util::QuasiAlignment* firstTrueHit{nullptr};
  bool haveRev1{false};
  bool haveRev2{false};
  // bool* haveRev = nullptr;

  size_t i{0};
  for (auto& qa : jointHits) {
    ++i;
    auto& refName = formatter.index->refName(qa.tid);
    // === SAM
    if (qa.isPaired) {
      getSamFlags(qa, true, flags1, flags2);
      if (alnCtr != 0) {
        flags1 |= 0x100;
        flags2 |= 0x100;
      }

      /** NOTE : WHY IS txpLen 100 here --- we should store and read this from the index. **/
      // uint32_t txpLen = 100;
      // rapmap::utils::adjustOverhang(qa, txpLens[qa.tid], cigarStr1,
      // cigarStr2);
      // Reverse complement the read and reverse
      // the quality string if we need to
      std::string* readSeq1 = &(r.first.seq);
      // std::string* qstr1 = &(r.first.qual);
      if (!qa.fwd) {
        if (!haveRev1) {
          util::reverseRead(*readSeq1, read1Temp);
          haveRev1 = true;
        }
        readSeq1 = &(read1Temp);
        // qstr1 = &(qual1Temp);
      }

      std::string* readSeq2 = &(r.second.seq);
      // std::string* qstr2 = &(r.second.qual);
      if (!qa.mateIsFwd) {
        if (!haveRev2) {
          util::reverseRead(*readSeq2, read2Temp);
          haveRev2 = true;
        }
        readSeq2 = &(read2Temp);
        // qstr2 = &(qual2Temp);
      }
      // If the fragment overhangs the right end of the reference
      // adjust fragLen (overhanging the left end is already handled).
      int32_t read1Pos = qa.pos;
      int32_t read2Pos = qa.matePos;
      const bool read1First{read1Pos < read2Pos};

      // TODO : We don't have access to the txp len yet
      // const int32_t minPos = read1First ? read1Pos : read2Pos;
      // if (minPos + qa.fragLen > txpLen) { qa.fragLen = txpLen - minPos; }

      // get the fragment length as a signed int
      const int32_t fragLen = static_cast<int32_t>(qa.fragLen);

      sstream << readName.c_str() << '\t'                    // QNAME
              << flags1 << '\t'                              // FLAGS
              << refName << '\t'                             // RNAME
              << qa.pos + 1 << '\t'                          // POS (1-based)
              << 1 << '\t'                                   // MAPQ
              << cigarStr1.c_str() << '\t'                   // CIGAR
              << '=' << '\t'                                 // RNEXT
              << qa.matePos + 1 << '\t'                      // PNEXT
              << ((read1First) ? fragLen : -fragLen) << '\t' // TLEN
              << *readSeq1 << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\n';

      sstream << mateName.c_str() << '\t'                    // QNAME
              << flags2 << '\t'                              // FLAGS
              << refName << '\t'                             // RNAME
              << qa.matePos + 1 << '\t'                      // POS (1-based)
              << 1 << '\t'                                   // MAPQ
              << cigarStr2.c_str() << '\t'                   // CIGAR
              << '=' << '\t'                                 // RNEXT
              << qa.pos + 1 << '\t'                          // PNEXT
              << ((read1First) ? -fragLen : fragLen) << '\t' // TLEN
              << *readSeq2 << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\n';

    } else {
    }
    ++alnCtr;
  }

  return 0;
}

#endif
