#ifndef SAM_WRITER_HPP
#define SAM_WRITER_HPP

#include "PairedAlignmentFormatter.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishConfig.hpp"
#include "PufferfishSparseIndex.hpp"
#include "Util.hpp"
#include "BinWriter.hpp"
#include "parallel_hashmap/phmap.h"
#include "nonstd/string_view.hpp"

typedef uint16_t rLenType;
typedef uint32_t refLenType;


inline void getSamFlags(const pufferfish::util::QuasiAlignment& qaln, uint16_t& flags) {
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
inline void getSamFlags(const pufferfish::util::QuasiAlignment& qaln, bool peInput,
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
  bool read1Unaligned = qaln.mateStatus == pufferfish::util::MateStatus::PAIRED_END_RIGHT;
  bool read2Unaligned = qaln.mateStatus == pufferfish::util::MateStatus::PAIRED_END_LEFT;
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
inline void writeKrakOutHeader(IndexT& pfi, std::shared_ptr<spdlog::logger> out, pufferfish::AlignmentOpts* mopts) {
  BinWriter bw(100000);
  bw << !mopts->singleEnd; // isPaired (bool)
  auto& txpNames = pfi.getFullRefNames();
  // TODO: Determine which ref lengths we should use here
  auto& txpLens = pfi.getFullRefLengths();
  auto numRef = txpNames.size();
  bw << static_cast<uint64_t>(numRef); // refCount (size_t)
  std::cerr << "is paired: " << !mopts->singleEnd << "\n";
  std::cerr << "numRef: " << numRef << "\n";
            
  for (size_t i = 0; i < numRef; ++i) {
    bw << txpNames[i] << txpLens[i]; //txpName (string) , txpLength (size_t)
  }
  out->info("{}",bw);
}


template <typename IndexT>
inline void writeSAMHeader(IndexT& pfi, std::shared_ptr<spdlog::logger> out) {
  fmt::MemoryWriter hd;
  hd.write("@HD\tVN:1.0\tSO:unknown\n");

  auto& txpNames = pfi.getFullRefNames();
  auto& txpLens = pfi.getFullRefLengthsComplete();

  auto numRef = txpNames.size();
  for (size_t i = 0; i < numRef; ++i) {
    hd.write("@SQ\tSN:{}\tLN:{:d}\n", txpNames[i], txpLens[i]);
  }
  // Eventually output a @PG line
  hd.write("@PG\tID:pufferfish\tPN:pufferfish\tVN:{}\n", pufferfish::version);
  std::string headerStr(hd.str());
  // Don't include the last '\n', since the logger will do it for us.
  headerStr.pop_back();
  out->info(headerStr);
}


template <typename IndexT>
inline void writeSAMHeader(IndexT& pfi, std::ostream& outStream) {
  fmt::MemoryWriter hd;
  hd.write("@HD\tVN:1.0\tSO:unknown\n");

  auto& txpNames = pfi.getFullRefNames();
  auto& txpLens = pfi.getFullRefLengthsComplete();

  auto numRef = txpNames.size();
  for (size_t i = 0; i < numRef; ++i) {
    hd.write("@SQ\tSN:{}\tLN:{:d}\n", txpNames[i], txpLens[i]);
  }
  // Eventually output a @PG line
  hd.write("@PG\tID:pufferfish\tPN:pufferfish\tVN:{}\n", pufferfish::version);
  outStream << hd.str();
}

template <typename IndexT>
inline void writeSAMHeader(IndexT& pfi, std::shared_ptr<spdlog::logger> out,
        bool filterGenomics,
        phmap::flat_hash_set<std::string> gene_names,
        phmap::flat_hash_set<std::string> rrna_names) {
  fmt::MemoryWriter hd;
  hd.write("@HD\tVN:1.0\tSO:unknown\n");

  auto& txpNames = pfi.getFullRefNames();
  auto& txpLens = pfi.getFullRefLengthsComplete();

  auto numRef = txpNames.size();

  // for now go with constant length
  // TODO read reference information
  // while reading the index
  for (size_t i = 0; i < numRef; ++i) {
    if (filterGenomics and
    (gene_names.find(txpNames[i]) != gene_names.end() or rrna_names.find(txpNames[i]) != rrna_names.end()))
      continue; 
    hd.write("@SQ\tSN:{}\tLN:{:d}\n", txpNames[i], txpLens[i]);
  }
  // Eventually output a @PG line
  // some other version number for now,
  // will think about it later
  std::string version = "1.0.0";
  hd.write("@PG\tID:pufferfish\tPN:pufferfish\tVN:{}\n", pufferfish::version);
  std::string headerStr(hd.str());
  // Don't include the last '\n', since the logger will do it for us.
  headerStr.pop_back();
  out->info(headerStr);
}

// Declarations for functions dealing with SAM formatting and output
//
inline void adjustOverhang(int32_t& pos, uint32_t readLen, uint32_t txpLen,
                           pufferfish::util::FixedWriter& cigarStr) {
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

inline void adjustOverhang(pufferfish::util::QuasiAlignment& qa, uint32_t txpLen,
                           pufferfish::util::FixedWriter& cigarStr1,
                           pufferfish::util::FixedWriter& cigarStr2) {
  if (qa.isPaired) { // both mapped
    adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr1);
    adjustOverhang(qa.matePos, qa.mateLen, txpLen, cigarStr2);
  } else if (qa.mateStatus == pufferfish::util::MateStatus::PAIRED_END_LEFT) {
    // left read mapped
    adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr1);
    // right read un-mapped will just be read length * S
    cigarStr2.clear();
    cigarStr2.write("{}S", qa.mateLen);
  } else if (qa.mateStatus == pufferfish::util::MateStatus::PAIRED_END_RIGHT) {
    // right read mapped
    adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr2);
    // left read un-mapped will just be read length * S
    cigarStr1.clear();
    cigarStr1.write("{}S", qa.readLen);
  }
}

// dump paired end read
template <typename ReadT , typename IndexT >
inline uint32_t writeAlignmentsToKrakenDump(ReadT& r,
                                   PairedAlignmentFormatter<IndexT>& formatter,
                                   std::vector<pufferfish::util::JointMems>& validJointHits,
                                   BinWriter& bstream,
                                   bool justMap,
                                   bool wrtIntervals=true) {
  if (validJointHits.empty()) return 0;
  auto& readName = r.first.name;
  //std::cout << readName << " || ";
  // If the read name contains multiple space-separated parts,
  // print only the first
  size_t nameLen = readName.length();
  size_t splitPos = readName.find(' ');
  if (splitPos < readName.length()) {
    readName[splitPos] = '\0';
    nameLen = splitPos;
  } else {
    splitPos = readName.length();
  } 
  
  if (splitPos > 2 and readName[splitPos - 2] == '/') {
    readName[splitPos - 2] = '\0';
    nameLen = splitPos - 2;
  }
  //std::cout << strlen(readName.c_str() << "\n";
  bstream   << stx::string_view(readName.data(), nameLen) 
            << static_cast<uint32_t>(validJointHits.size())
            << static_cast<rLenType>(r.first.seq.length()) 
            << static_cast<rLenType>(r.second.seq.length());
  for (auto& qa : validJointHits) {

      /* auto& refName = formatter.index->refName(qa.tid);
      uint32_t refLength = formatter.index->refLength(qa.tid);
       */
    auto& clustLeft = qa.leftClust;
    auto& clustRight = qa.rightClust;
    size_t leftNumOfIntervals = 0;
    size_t rightNumOfIntervals = 0;
    refLenType leftRefPos = 0;
    refLenType rightRefPos = 0;
    if (qa.isLeftAvailable()) {
      uint64_t pos = clustLeft->getTrFirstHitPos() < 0 ? 0 : clustLeft->getTrFirstHitPos();
      leftNumOfIntervals = clustLeft->mems.size();
      leftRefPos = pos | (static_cast<refLenType>(clustLeft->isFw) << (sizeof(refLenType)*8-1));
    }
    if (qa.isRightAvailable()) {
      uint64_t pos = clustRight->getTrFirstHitPos() < 0 ? 0 : clustRight->getTrFirstHitPos();
      rightNumOfIntervals = clustRight->mems.size();
      rightRefPos = pos | (static_cast<refLenType>(clustRight->isFw) << (sizeof(refLenType)*8-1));
    }
    bstream << static_cast<uint32_t>(formatter.index->getRefId(qa.tid));
      /*std::cerr << "r " << readName << " puffid: " << qa.tid << " lcnt:" << leftNumOfIntervals
                << " rcnt:" << rightNumOfIntervals;*/
//	if (wrtIntervals) {
		bstream << static_cast<rLenType>(leftNumOfIntervals) 
            	<< static_cast<rLenType>(rightNumOfIntervals);
//	}

      if (qa.isLeftAvailable()) {
          //std::cerr << " lrefpos: " << leftRefPos;
        if (justMap) {
          bstream << qa.alignmentScore;
        } else {
          bstream << clustLeft->coverage;
        }
        bstream << static_cast<refLenType>(leftRefPos);
    	if (wrtIntervals) {
        	for (auto& mem: clustLeft->mems) {
         	   	bstream << static_cast<rLenType>(mem.rpos) 
                  		<< static_cast<rLenType>(mem.extendedlen);
        	}
      	}
	  }
      if (qa.isRightAvailable()) {
          //std::cerr << " rrefpos: " << rightRefPos;
          if (justMap) {
            bstream << qa.mateAlignmentScore;
          } else {
            bstream << clustRight->coverage;
          }
        bstream << static_cast<refLenType>(rightRefPos);
    	if (wrtIntervals) {
        	for (auto& mem: clustRight->mems) {
          		bstream << static_cast<rLenType>(mem.rpos)
                  		<< static_cast<rLenType>(mem.extendedlen);
        	}
      	}
      }
      //std::cerr << "\n";
  }
  return 0;
}
// dump single end read
template <typename ReadT , typename IndexT >
inline uint32_t writeAlignmentsToKrakenDump(ReadT& r,
                                   PairedAlignmentFormatter<IndexT>& formatter,
                                   std::vector<std::pair<uint32_t, std::vector<pufferfish::util::MemCluster>::iterator>>& validHits,
                                   BinWriter& binStream,
                                   bool wrtIntervals=true) {
    if (validHits.empty()) return 0;

  auto& readName = r.name;
  //std::cerr << readName  << " -- ";
  size_t nameLen = readName.length();

  // If the read name contains multiple space-separated parts,
  // print only the first
  size_t splitPos = readName.find(' ');
  if (splitPos < readName.length()) {
    readName[splitPos] = '\0';
    nameLen = splitPos;
  } else {
    splitPos = readName.length();
  }

  if (splitPos > 2 and readName[splitPos - 2] == '/') {
    readName[splitPos - 2] = '\0';
    nameLen = splitPos - 2;
  }
    //std::cerr << readName << " " << nameLen << "\n";
  binStream << stx::string_view(readName.data(), nameLen) 
            << static_cast<uint32_t>(validHits.size()) 
            << static_cast<rLenType>(r.seq.length());
  //std::cerr << readName << "\t" << validHits.size() << "\t" << r.seq.length() << "\n";
  for (auto& qa : validHits) {
    auto& clust = qa.second;
    binStream << static_cast<uint32_t>(formatter.index->getRefId(qa.first));
//    if (wrtIntervals) {
      binStream << static_cast<rLenType>(clust->mems.size());
//    }
    //std::cerr << qa.first << "\t" << clust->mems.size() << "\n";
    binStream << clust->coverage;

    binStream << static_cast<refLenType>(clust->getTrFirstHitPos() | (static_cast<refLenType>(clust->isFw) << (sizeof(refLenType)*8-1)));
    if (wrtIntervals) {
      for (auto& mem: clust->mems){
        binStream << static_cast<rLenType>(mem.rpos) 
                  << static_cast<rLenType>(mem.extendedlen);
      }
    }
  }
  return 0;

}

inline uint32_t writeUnalignedPairToStream(fastx_parser::ReadPair& r,
                                          fmt::MemoryWriter& sstream) {
        constexpr uint16_t flags1 = 0x1 | 0x4 | 0x8 | 0x40;
        constexpr uint16_t flags2 = 0x1 | 0x4 | 0x8 | 0x80;

        auto processReadName = [](const std::string& name) -> fmt::StringRef {
                                 nonstd::string_view readNameView(name);
                                 // If the read name contains multiple space-separated parts,
                                 // print only the first
                                 size_t splitPos = readNameView.find(' ');
                                 if (splitPos < readNameView.length()) {
                                   readNameView.remove_suffix(readNameView.length() - splitPos);
                                 } else {
                                   splitPos = readNameView.length();
                                 }

                                 // trim /1 from the pe read
                                 if (splitPos > 2 and readNameView[splitPos - 2] == '/') {
                                   readNameView.remove_suffix(2);
                                   //readName[splitPos - 2] = '\0';
                                 }
                                 return fmt::StringRef(readNameView.data(), readNameView.size());
                               };

        auto readNameView = processReadName(r.first.name);
        auto mateNameView = processReadName(r.second.name);
        std::string* readSeq1 = &(r.first.seq);
        std::string* readSeq2 = &(r.second.seq);

        sstream << readNameView << '\t' // QNAME
                << flags1 << '\t'       // FLAGS
                << "*\t"                // RNAME
                << "0\t"                // POS (1-based)
                << "255\t"              // MAPQ
                << "*\t"                // CIGAR
                << "*\t"                // RNEXT
                << "0\t"                // PNEXT
                << "0\t"                // TLEN
                << *readSeq1 << '\t'    // SEQ
                << "*\t"                // QUAL
                << "NH:i:0\t"
                << "HI:i:0\t"
                << "AS:i:0\n";

        sstream << mateNameView << '\t' // QNAME
                << flags2 << '\t'       // FLAGS
                << "*\t"                // RNAME
                << "0\t"                // POS (1-based)
                << "255\t"              // MAPQ
                << "*\t"                // CIGAR
                << "*\t"                // RNEXT
                << "0\t"                // PNEXT
                << "0\t"                // TLEN
                << *readSeq2 << '\t'    // SEQ
                << "*\t"                // QUAL
                << "NH:i:0\t"
                << "HI:i:0\t"
                << "AS:i:0\n";
        return 0;
      }

inline uint32_t writeUnalignedSingleToStream(fastx_parser::ReadSeq& r,
                                            fmt::MemoryWriter& sstream) {
        constexpr uint16_t flags = 0x4;

        nonstd::string_view readNameViewSV(r.name);
        // If the read name contains multiple space-separated parts, print
        // only the first
        size_t splitPos = readNameViewSV.find(' ');
        if (splitPos < readNameViewSV.length()) {
          readNameViewSV.remove_suffix(readNameViewSV.size() - splitPos);
        }
        fmt::StringRef readNameView(readNameViewSV.data(), readNameViewSV.size());
        std::string* readSeq = &(r.seq);

        sstream << readNameView << '\t' // QNAME
                << flags << '\t'    // FLAGS
                << "*\t"            // RNAME
                << 0 << '\t'        // POS (1-based)
                << 255 << '\t'      // MAPQ
                << "*\t"            // CIGAR
                << "*\t"            // MATE NAME
                << "0\t"            // MATE POS
                << "0\t"            // TLEN
                << *readSeq << '\t' // SEQ
                << "*\t"            // QSTR
                << "NH:i:0\t"
                << "HI:i:0\t"
                << "AS:i:0\n";
        return 0;
      }

template <typename ReadT, typename IndexT>
inline uint32_t writeAlignmentsToStreamSingle(
    ReadT& r, PairedAlignmentFormatter<IndexT>& formatter,
    std::vector<pufferfish::util::QuasiAlignment>& jointHits, fmt::MemoryWriter& sstream, bool writeOrphans, 
    bool tidsAlreadyDecoded = false) {
  (void) writeOrphans;


  auto& readTemp = formatter.read1Temp;
  auto& cigarStr = formatter.cigarStr1;

  uint16_t flags;

  nonstd::string_view readNameViewSV(r.name);

  // If the read name contains multiple space-separated parts, print
  // only the first
  size_t splitPos = readNameViewSV.find(' ');
  if (splitPos < readNameViewSV.length()) {
    readNameViewSV.remove_suffix(readNameViewSV.size() - splitPos);
  }
  fmt::StringRef readNameView(readNameViewSV.data(), readNameViewSV.size());

  std::string numHitFlag = fmt::format("NH:i:{}", jointHits.size());
  uint32_t alnCtr{0};
  bool haveRev{false};
  size_t i{0};

  auto* fullRefNames = tidsAlreadyDecoded ? &formatter.index->getFullRefNames() : nullptr;
  auto* fullRefLengths = tidsAlreadyDecoded ? &formatter.index->getFullRefLengths() : nullptr;
  for (auto& qa : jointHits) {
    ++i;
    auto& refName = tidsAlreadyDecoded ? (*fullRefNames)[qa.tid] : formatter.index->refName(qa.tid);
    uint32_t txpLen = tidsAlreadyDecoded ? (*fullRefLengths)[qa.tid] : formatter.index->refLength(qa.tid);
    // === SAM
      getSamFlags(qa, flags);
      if (alnCtr != 0) {
        flags |= 0x100;
      }

      // Reverse complement the read and reverse
      // the quality string if we need to
      std::string* readSeq = &(r.seq);
      // std::string* qstr1 = &(r.first.qual);
      if (!qa.fwd) {
        if (!haveRev) {
          pufferfish::util::reverseRead(*readSeq, readTemp);
          haveRev = true;
        }
        readSeq = &(readTemp);
      }

      adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr);

      sstream << readNameView << '\t' // QNAME
              << flags << '\t' // FLAGS
              << refName << '\t' // RNAME
              << qa.pos + 1 << '\t' // POS (1-based)
              << 255 << '\t' // MAPQ
              << (qa.cigar.empty() ? cigarStr.c_str() : qa.cigar) << '\t' // CIGAR
              << '*' << '\t' // MATE NAME
              << 0 << '\t' // MATE POS
              << qa.fragLen << '\t' // TLEN
              << *readSeq << '\t' // SEQ
              << "*\t" // QSTR
              << numHitFlag << '\t'
              << "HI:i:" << i << '\t'
              << "AS:i:" << qa.score << '\n';
    ++alnCtr;
  }
  return 0;
}

template <typename ReadPairT, typename IndexT>
inline uint32_t writeAlignmentsToStream(
    ReadPairT& r, PairedAlignmentFormatter<IndexT>& formatter,
    std::vector<pufferfish::util::QuasiAlignment>& jointHits, fmt::MemoryWriter& sstream,
    bool writeOrphans,
    bool tidsAlreadyDecoded = false) {

  auto& read1Temp = formatter.read1Temp;
  auto& read2Temp = formatter.read2Temp;
  auto& cigarStr1 = formatter.cigarStr1;
  auto& cigarStr2 = formatter.cigarStr2;

  uint16_t flags1, flags2;

  auto processReadName = [](const std::string& name) -> fmt::StringRef {
                           nonstd::string_view readNameView(name);
                           // If the read name contains multiple space-separated parts,
                           // print only the first
                           size_t splitPos = readNameView.find(' ');
                           if (splitPos < readNameView.length()) {
                             readNameView.remove_suffix(readNameView.length() - splitPos);
                           } else {
                             splitPos = readNameView.length();
                           }

                           // trim /1 from the pe read
                           if (splitPos > 2 and readNameView[splitPos - 2] == '/') {
                             readNameView.remove_suffix(2);
                             //readName[splitPos - 2] = '\0';
                           }
                           return fmt::StringRef(readNameView.data(), readNameView.size());
                         };

  auto readNameView = processReadName(r.first.name);
  auto mateNameView = processReadName(r.second.name);

  cigarStr1.clear();
  cigarStr2.clear();
  cigarStr1.write("{}M", r.first.seq.length());
  cigarStr2.write("{}M", r.second.seq.length());

  std::string numHitFlag = fmt::format("NH:i:{}", jointHits.size());
  uint32_t alnCtr{0};
  // uint32_t trueHitCtr{0};
  // pufferfish::util::QuasiAlignment* firstTrueHit{nullptr};
  bool haveRev1{false};
  bool haveRev2{false};
  bool* haveRev = nullptr;
  size_t i{0};

  auto* fullRefNames = tidsAlreadyDecoded ? &formatter.index->getFullRefNames() : nullptr;
  auto* fullRefLengths = tidsAlreadyDecoded ? &formatter.index->getFullRefLengths() : nullptr;

  for (auto& qa : jointHits) {
    ++i;
    auto& refName = tidsAlreadyDecoded ? (*fullRefNames)[qa.tid] : formatter.index->refName(qa.tid);
    uint32_t txpLen = tidsAlreadyDecoded ? (*fullRefLengths)[qa.tid] : formatter.index->refLength(qa.tid);
    // === SAM
    if (qa.isPaired) {
      getSamFlags(qa, true, flags1, flags2);
      if (alnCtr != 0) {
        flags1 |= 0x100;
        flags2 |= 0x100;
      }

      adjustOverhang(qa, txpLen, cigarStr1, cigarStr2);
      // Reverse complement the read and reverse
      // the quality string if we need to
      std::string* readSeq1 = &(r.first.seq);
      // std::string* qstr1 = &(r.first.qual);
      if (!qa.fwd) {
        if (!haveRev1) {
          pufferfish::util::reverseRead(*readSeq1, read1Temp);
          haveRev1 = true;
        }
        readSeq1 = &(read1Temp);
        // qstr1 = &(qual1Temp);
      }

      std::string* readSeq2 = &(r.second.seq);
      // std::string* qstr2 = &(r.second.qual);
      if (!qa.mateIsFwd) {
        if (!haveRev2) {
          pufferfish::util::reverseRead(*readSeq2, read2Temp);
          haveRev2 = true;
        }
        readSeq2 = &(read2Temp);
        // qstr2 = &(qual2Temp);
      }


      // If the fragment overhangs the right end of the transcript
      // adjust fragLen (overhanging the left end is already handled).
      int32_t read1Pos = qa.pos;
      int32_t read2Pos = qa.matePos;
      const bool read1First{read1Pos < read2Pos};
      const int32_t minPos = read1First ? read1Pos : read2Pos;
      if ((minPos + static_cast<int32_t>(qa.fragLen)) > static_cast<int32_t>(txpLen)) { qa.fragLen = txpLen - minPos; }
      // get the fragment length as a signed int
      const int32_t fragLen = static_cast<int32_t>(qa.fragLen);

      std::stringstream ss;

      sstream << readNameView << '\t'                    // QNAME
              //<< qa.numHits << '\t'
              << flags1 << '\t'                              // FLAGS
              << refName << '\t'                             // RNAME
              << qa.pos + 1 << '\t'                          // POS (1-based)
              << 1 << '\t'                                   // MAPQ
              << (qa.cigar.empty() ? cigarStr1.c_str() : qa.cigar) << '\t'                   // CIGAR
              //<< "100M" << '\t'
              //<< qa.cigar << '\t'                   // CIGAR
              << '=' << '\t'                                 // RNEXT
              << qa.matePos + 1 << '\t'                      // PNEXT
              << ((read1First) ? fragLen : -fragLen) << '\t' // TLEN
              << *readSeq1 << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\t'
              << "HI:i:" << i << '\t'
              << "AS:i:" << qa.score << '\n';

      sstream << mateNameView << '\t'                    // QNAME
              //<< qa.numHits << '\t'
              << flags2 << '\t'                              // FLAGS
              << refName << '\t'                             // RNAME
              << qa.matePos + 1 << '\t'                      // POS (1-based)
              << 1 << '\t'                                   // MAPQ
              << (qa.mateCigar.empty() ? cigarStr2.c_str() : qa.mateCigar) << '\t'                   // CIGAR
              //<< "100M" << '\t'
              //<< qa.mateCigar << '\t'                   // CIGAR
              << '=' << '\t'                                 // RNEXT
              << qa.pos + 1 << '\t'                          // PNEXT
              << ((read1First) ? -fragLen : fragLen) << '\t' // TLEN
              << *readSeq2 << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\t'
              << "HI:i:" << i << '\t'
              << "AS:i:" << qa.mateScore << '\n';
    } else if(writeOrphans) {
		//added orphan support
	  //std::cerr<<"orphans here";

      getSamFlags(qa, true, flags1, flags2);
      if (alnCtr != 0) {
        flags1 |= 0x100;
        flags2 |= 0x100;
      }
      /** NOTE : WHY IS txpLen 100 here --- we should store and read this from the index. **/
      adjustOverhang(qa, txpLen, cigarStr1, cigarStr2);
      // Reverse complement the read and reverse
      // the quality string if we need to

      std::string* readSeq{nullptr} ;
      std::string* unalignedSeq{nullptr} ;

      uint32_t flags, unalignedFlags ;

      fmt::StringRef* alignedName{nullptr};
      fmt::StringRef* unalignedName{nullptr};
      std::string* readTemp{nullptr};

      auto* cigarStr = &formatter.cigarStr1;
      cigarStr->clear();
      cigarStr->write("{}M", r.first.seq.length());

      //logic for assigning orphans
      if(qa.mateStatus == pufferfish::util::MateStatus::PAIRED_END_LEFT){ //left read
        alignedName = &readNameView;
        unalignedName = &mateNameView;

        readSeq = &(r.first.seq) ;
        unalignedSeq = &(r.second.seq) ;

        flags = flags1 ;
        unalignedFlags = flags2 ;

        haveRev = &haveRev1 ;
        readTemp = &read1Temp ;
      } else {
        alignedName = &mateNameView;
        unalignedName = &readNameView;

        cigarStr = &formatter.cigarStr2;
        cigarStr->clear();
        cigarStr->write("{}M", r.second.seq.length());

        readSeq = &(r.second.seq) ;
        unalignedSeq = &(r.first.seq) ;

        flags = flags2 ;
        unalignedFlags = flags1 ;

        haveRev = &haveRev2 ;
        readTemp = &read2Temp ;

      }


      // std::string* qstr1 = &(r.first.qual);
      if (!qa.fwd) {
        if (!*haveRev) {
          pufferfish::util::reverseRead(*readSeq, *readTemp);
          *haveRev = true;
        }
        readSeq = readTemp;
      }

      // If the fragment overhangs the right end of the reference
      // adjust fragLen (overhanging the left end is already handled).

      sstream << *(alignedName) << '\t'                    // QNAME
              << flags << '\t'                              // FLAGS
              << refName << '\t'                             // RNAME
              << qa.pos + 1 << '\t'                          // POS (1-based)
              << 1 << '\t'                                   // MAPQ
              << (qa.cigar.empty() ? (*cigarStr).c_str() : qa.cigar) << '\t'                   // CIGAR
              //<< "100M" << '\t' 
              << '=' << '\t'                                 // RNEXT
              << /* qa.matePos */ qa.pos + 1 << '\t'                      // PNEXT
              << 0 << '\t'                                     // TLEN
              << *readSeq << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\t'
              << "HI:i:" << i << '\t'
              << "AS:i:" << qa.score << '\n';


      sstream << *(unalignedName) << '\t'                    // QNAME
              << unalignedFlags << '\t'                              // FLAGS
              << refName << '\t'                             // RNAME
              << qa.pos + 1 << '\t'                      // POS (1-based)
              << 0 << '\t'                                   // MAPQ
              << "*\t"                   // CIGAR
              << '=' << '\t'                                 // RNEXT
              << qa.pos + 1 << '\t'                          // PNEXT
              << 0 << '\t'                                   // TLEN
              << *unalignedSeq << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\t'
              << "HI:i:" << i << '\t'
              << "AS:i:" << qa.mateScore << '\n';

    }
    ++alnCtr;
  }

  return 0;
}

#endif
