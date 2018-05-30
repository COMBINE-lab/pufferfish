#ifndef SAM_WRITER_HPP
#define SAM_WRITER_HPP

#include "PairedAlignmentFormatter.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "Util.hpp"
#include "BinWriter.hpp"

typedef uint16_t rLenType;
typedef uint32_t refLenType;


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
inline void writeKrakOutHeader(IndexT& pfi, std::shared_ptr<spdlog::logger> out, AlignmentOpts* mopts) {
  BinWriter bw(100000);
  bw << !mopts->singleEnd; // isPaired (bool)
  auto& txpNames = pfi.getRefNames();
  auto& txpLens = pfi.getRefLengths();
  auto numRef = txpNames.size();
  bw << numRef; // refCount (size_t)
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

  auto& txpNames = pfi.getRefNames();
  auto& txpLens = pfi.getRefLengths();

  auto numRef = txpNames.size();
  // for now go with constant length
  // TODO read reference information
  // while reading the index
  for (size_t i = 0; i < numRef; ++i) {
    hd.write("@SQ\tSN:{}\tLN:{:d}\n", txpNames[i], txpLens[i]);
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

// dump paired end read
template <typename ReadT/* , typename IndexT */>
inline uint32_t writeAlignmentsToKrakenDump(ReadT& r,
                                   //PairedAlignmentFormatter<IndexT>& formatter,
                                   std::vector<util::JointMems>& validJointHits,
                                   BinWriter& bstream,
                                   bool wrtIntervals=true) {

    if (validJointHits.size() == 0) return 0;
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
      leftNumOfIntervals = clustLeft->mems.size();
      leftRefPos = clustLeft->getTrFirstHitPos() | (static_cast<refLenType>(clustLeft->isFw) << (sizeof(refLenType)*8-1));
    }
    if (qa.isRightAvailable()) {
      rightNumOfIntervals = clustRight->mems.size();
      rightRefPos = clustRight->getTrFirstHitPos() | (static_cast<refLenType>(clustRight->isFw) << (sizeof(refLenType)*8-1));
    }
    bstream << static_cast<uint32_t>(qa.tid);
      /*std::cerr << "r " << readName << " puffid: " << qa.tid << " lcnt:" << leftNumOfIntervals
                << " rcnt:" << rightNumOfIntervals;*/
//	if (wrtIntervals) {
		bstream << static_cast<rLenType>(leftNumOfIntervals) 
            	<< static_cast<rLenType>(rightNumOfIntervals);
//	}

      if (qa.isLeftAvailable()) {
          //std::cerr << " lrefpos: " << leftRefPos;
        bstream << static_cast<refLenType>(leftRefPos);
    	if (wrtIntervals) {
        	for (auto& mem: clustLeft->mems) {
         	   	bstream << static_cast<rLenType>(mem.memInfo->rpos) 
                  		<< static_cast<rLenType>(mem.memInfo->memlen);
        	}
      	}
	  }
      if (qa.isRightAvailable()) {
          //std::cerr << " rrefpos: " << rightRefPos;
        bstream << static_cast<refLenType>(rightRefPos);
    	if (wrtIntervals) {
        	for (auto& mem: clustRight->mems) {
          		bstream << static_cast<rLenType>(mem.memInfo->rpos) 
                  		<< static_cast<rLenType>(mem.memInfo->memlen);
        	}
      	}
      }
      //std::cerr << "\n";
  }
  return 0;
}

// dump single end read
template <typename ReadT/* , typename IndexT */>
inline uint32_t writeAlignmentsToKrakenDump(ReadT& r,
                                   //PairedAlignmentFormatter<IndexT>& formatter,
                                   std::vector<std::pair<uint32_t, std::vector<util::MemCluster>::iterator>>& validHits,
                                   BinWriter& binStream,
                                   bool wrtIntervals=true) {
    if (validHits.size() == 0) return 0;
  auto& readName = r.name;
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

  binStream << stx::string_view(readName.data(), nameLen) 
            << static_cast<uint32_t>(validHits.size()) 
            << static_cast<rLenType>(r.seq.length());
  
  for (auto& qa : validHits) {
    auto& clust = qa.second;
    binStream << static_cast<uint32_t>(qa.first);
//    if (wrtIntervals) {
      binStream << static_cast<rLenType>(clust->mems.size());
//    }
    binStream << static_cast<refLenType>(clust->getTrFirstHitPos() | (static_cast<refLenType>(clust->isFw) << (sizeof(refLenType)*8-1)));
    if (wrtIntervals) {
      for (auto& mem: clust->mems){
        binStream << static_cast<rLenType>(mem.memInfo->rpos) 
                  << static_cast<rLenType>(mem.memInfo->memlen);
      }
    }
  }   
  return 0;

}

template <typename ReadT, typename IndexT>
inline uint32_t writeUnmappedAlignmentsToStreamSingle(
    ReadT& r, PairedAlignmentFormatter<IndexT>& formatter,
    std::vector<util::QuasiAlignment>& jointHits, fmt::MemoryWriter& sstream, bool writeOrphans, bool justMappings) {
  (void) writeOrphans;
  (void) justMappings;

	//auto& read1Temp = formatter.read1Temp;
	//auto& read2Temp = formatter.read2Temp;
	// auto& qual1Temp = formatter.qual1Temp;
	// auto& qual2Temp = formatter.qual2Temp;
	auto& cigarStr1 = formatter.cigarStr1;
	//auto& cigarStr2 = formatter.cigarStr2;

	cigarStr1.clear();
	//cigarStr2.clear();
	cigarStr1.write("*");//"{}M", r.first.seq.length());
	//cigarStr2.write("*");//"{}M", r.second.seq.length());
	//std::cerr << cigarStr1.c_str() << "\n";

	//uint16_t flags1;
	auto& readName = r.name;

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

    std::string* readSeq1 = &(r.seq);

    std::string numHitFlag = fmt::format("NH:i:0", jointHits.size());

    sstream << readName.c_str() << '\t' // QNAME
            << 4 << '\t'               // FLAGS
            << "*\t"                    // RNAME
            << 0 << '\t'                // POS (1-based)
            << 255 << '\t'              // MAPQ
            << cigarStr1.c_str()
            << '\t' // CIGAR
            << '=' << '\t'       // RNEXT
            << 0 << '\t'         // PNEXT
            << 0 << '\t'         // TLEN
            << *readSeq1 << '\t' // SEQ
            << "*\t"             // QUAL
            << numHitFlag << '\n';
    return 0;
}

template <typename ReadPairT, typename IndexT>
inline uint32_t writeUnmappedAlignmentsToStream(
    ReadPairT& r, PairedAlignmentFormatter<IndexT>& formatter,
    std::vector<util::QuasiAlignment>& jointHits, fmt::MemoryWriter& sstream, bool writeOrphans, bool justMappings) {
  (void) writeOrphans;
  (void) justMappings;
	//auto& read1Temp = formatter.read1Temp;
	//auto& read2Temp = formatter.read2Temp;
	// auto& qual1Temp = formatter.qual1Temp;
	// auto& qual2Temp = formatter.qual2Temp;
	auto& cigarStr1 = formatter.cigarStr1;
	auto& cigarStr2 = formatter.cigarStr2;

	cigarStr1.clear();
	cigarStr2.clear();
	cigarStr1.write("*");//"{}M", r.first.seq.length());
	cigarStr2.write("*");//"{}M", r.second.seq.length());
	//std::cerr << cigarStr1.c_str() << "\n";
	//uint16_t flags1, flags2;

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
      std::string* readSeq1 = &(r.first.seq);
      std::string* readSeq2 = &(r.second.seq);

     std::string numHitFlag = fmt::format("NH:i:0", jointHits.size());

      sstream << readName.c_str() << '\t'                    // QNAME
              << 77 << '\t'                              // FLAGS
              << "*\t"                             // RNAME
              << 0 << '\t'                          // POS (1-based)
              << 255 << '\t'                                   // MAPQ
              << cigarStr1.c_str() << '\t'                   // CIGAR
              //<< qa.cigar << '\t'                   // CIGAR
              << '=' << '\t'                                 // RNEXT
              << 0 << '\t'                      // PNEXT
              << 0 << '\t' // TLEN
              << *readSeq1 << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\n';

      sstream << mateName.c_str() << '\t'                    // QNAME
              << 141 << '\t'                              // FLAGS
              << "*\t"                             // RNAME
              << 0 << '\t'                      // POS (1-based)
              << 255 << '\t'                                   // MAPQ
              << cigarStr2.c_str() << '\t'                   // CIGAR
	        //<< qa.mateCigar << '\t'                   // CIGAR
              << '=' << '\t'                                 // RNEXT
              << 0 << '\t'                          // PNEXT
              << 0 << '\t' // TLEN
              << *readSeq2 << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\n';

   return 0;
}


template <typename ReadT, typename IndexT>
inline uint32_t writeAlignmentsToStreamSingle(
    ReadT& r, PairedAlignmentFormatter<IndexT>& formatter,
    std::vector<util::QuasiAlignment>& jointHits, fmt::MemoryWriter& sstream, bool writeOrphans, bool justMappings) {
  (void) writeOrphans;
  auto& read1Temp = formatter.read1Temp;
  //auto& read2Temp = formatter.read2Temp;
  // auto& qual1Temp = formatter.qual1Temp;
  // auto& qual2Temp = formatter.qual2Temp;
  auto& cigarStr1 = formatter.cigarStr1;
  //auto& cigarStr2 = formatter.cigarStr2;

  cigarStr1.clear();
  //cigarStr2.clear();
  cigarStr1.write("{}M", r.seq.length());
  //cigarStr2.write("{}M", r.second.seq.length());
  //std::cerr << cigarStr1.c_str() << "\n";
  uint16_t flags1;

  auto& readName = r.name;
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

  std::string numHitFlag = fmt::format("NH:i:{}", jointHits.size());
  uint32_t alnCtr{0};
  // uint32_t trueHitCtr{0};
  // util::QuasiAlignment* firstTrueHit{nullptr};
  bool haveRev1{false};

  size_t i{0};
  for (auto& qa : jointHits) {
    ++i;
    auto& refName = formatter.index->refName(qa.tid);
    uint32_t txpLen = formatter.index->refLength(qa.tid);
    // === SAM
      getSamFlags(qa, flags1);
      if (alnCtr != 0) {
        flags1 |= 0x100;
      }

      adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr1);

      // Reverse complement the read and reverse
      // the quality string if we need to
      std::string* readSeq1 = &(r.seq);
      // std::string* qstr1 = &(r.first.qual);
      if (!qa.fwd) {
        if (!haveRev1) {
          util::reverseRead(*readSeq1, read1Temp);
          haveRev1 = true;
        }
        readSeq1 = &(read1Temp);
        // qstr1 = &(qual1Temp);
      }

      // If the fragment overhangs the right end of the reference
      // adjust fragLen (overhanging the left end is already handled).
      int32_t read1Pos = qa.pos;
      const bool read1First{true};

      const int32_t minPos = read1Pos;
      if (minPos + qa.fragLen > txpLen) { qa.fragLen = txpLen - minPos; }

      // get the fragment length as a signed int
      const int32_t fragLen = static_cast<int32_t>(qa.fragLen);

      sstream << readName.c_str() << '\t'                    // QNAME
              << flags1 << '\t'                              // FLAGS
              << refName << '\t'                             // RNAME
              << qa.pos + 1 << '\t'                          // POS (1-based)
              << 1 << '\t'                                   // MAPQ
              << (justMappings ? cigarStr1.c_str() : qa.cigar) << '\t'                   // CIGAR
              << '=' << '\t'                                 // RNEXT
              << 0 << '\t'                      // PNEXT
              << ((read1First) ? fragLen : -fragLen) << '\t' // TLEN
              << *readSeq1 << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\n';

    ++alnCtr;
  }
  return 0;
}

template <typename ReadPairT, typename IndexT>
inline uint32_t writeAlignmentsToStream(
    ReadPairT& r, PairedAlignmentFormatter<IndexT>& formatter,
    std::vector<util::QuasiAlignment>& jointHits, fmt::MemoryWriter& sstream, 
    bool writeOrphans, 
    bool justMappings) {

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
  bool* haveRev = nullptr;

  size_t i{0};
  for (auto& qa : jointHits) {
    ++i;
    auto& refName = formatter.index->refName(qa.tid);
    uint32_t txpLen = formatter.index->refLength(qa.tid);
    // === SAM
    if (qa.isPaired) {
      getSamFlags(qa, true, flags1, flags2);
      if (alnCtr != 0) {
        flags1 |= 0x100;
        flags2 |= 0x100;
      }

      /** NOTE : WHY IS txpLen 100 here --- we should store and read this from the index. **/
      adjustOverhang(qa, txpLen, cigarStr1, cigarStr2);
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
      const int32_t minPos = read1First ? read1Pos : read2Pos;
      if (minPos + qa.fragLen > txpLen) { qa.fragLen = txpLen - minPos; }

      // get the fragment length as a signed int
      const int32_t fragLen = static_cast<int32_t>(qa.fragLen);
      sstream << readName.c_str() << '\t'                    // QNAME
              //<< qa.numHits << '\t'
              << flags1 << '\t'                              // FLAGS
              << refName << '\t'                             // RNAME
              << qa.pos + 1 << '\t'                          // POS (1-based)
              << 1 << '\t'                                   // MAPQ
              << (justMappings ? cigarStr1.c_str() : qa.cigar) << '\t'                   // CIGAR
              //<< qa.cigar << '\t'                   // CIGAR
              << '=' << '\t'                                 // RNEXT
              << qa.matePos + 1 << '\t'                      // PNEXT
              << ((read1First) ? fragLen : -fragLen) << '\t' // TLEN
              << *readSeq1 << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\n';

      sstream << mateName.c_str() << '\t'                    // QNAME
              //<< qa.numHits << '\t'
              << flags2 << '\t'                              // FLAGS
              << refName << '\t'                             // RNAME
              << qa.matePos + 1 << '\t'                      // POS (1-based)
              << 1 << '\t'                                   // MAPQ
              << (justMappings ? cigarStr2.c_str() : qa.mateCigar) << '\t'                   // CIGAR
        //<< qa.mateCigar << '\t'                   // CIGAR
              << '=' << '\t'                                 // RNEXT
              << qa.pos + 1 << '\t'                          // PNEXT
              << ((read1First) ? -fragLen : fragLen) << '\t' // TLEN
              << *readSeq2 << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\n';

    }

    else if(writeOrphans) {
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

      std::string* alignedName{nullptr} ;
      std::string* unalignedName{nullptr} ;
      std::string* readTemp{nullptr} ;

      auto& cigarStr = formatter.cigarStr1;
      cigarStr.clear();
      cigarStr.write("{}M", r.first.seq.length());

      //logic for assigning orphans
      if(qa.mateStatus == util::MateStatus::PAIRED_END_LEFT){ //left read
        alignedName = &readName ;
        unalignedName = &mateName ;

        readSeq = &(r.first.seq) ;
        unalignedSeq = &(r.second.seq) ;

        flags = flags1 ;
        unalignedFlags = flags2 ;

        haveRev = &haveRev1 ;
        readTemp = &read1Temp ;
      } else {
        alignedName = &mateName ;
        unalignedName = &readName ;

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
          util::reverseRead(*readSeq, *readTemp);
          *haveRev = true;
        }
        readSeq = readTemp;
      }

      // If the fragment overhangs the right end of the reference
      // adjust fragLen (overhanging the left end is already handled).
      int32_t read1Pos = qa.pos;
      int32_t read2Pos = qa.matePos;
      const bool read1First{read1Pos < read2Pos};

      // TODO : We don't have access to the txp len yet
      const int32_t minPos = read1First ? read1Pos : read2Pos;
      if (minPos + qa.fragLen > txpLen) { qa.fragLen = txpLen - minPos; }

      // get the fragment length as a signed int
      const int32_t fragLen = static_cast<int32_t>(qa.fragLen);

      sstream << alignedName->c_str() << '\t'                    // QNAME
              << flags << '\t'                              // FLAGS
              << refName << '\t'                             // RNAME
              << qa.pos + 1 << '\t'                          // POS (1-based)
              << 1 << '\t'                                   // MAPQ
              << (justMappings ? cigarStr1.c_str() : qa.cigar) << '\t'                   // CIGAR
              << '=' << '\t'                                 // RNEXT
              << /* qa.matePos */ qa.pos + 1 << '\t'                      // PNEXT
              << ((read1First) ? fragLen : -fragLen) << '\t' // TLEN
              << *readSeq << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\n';

      sstream << unalignedName->c_str() << '\t'                    // QNAME
              << unalignedFlags << '\t'                              // FLAGS
              << refName << '\t'                             // RNAME
              << qa.pos + 1 << '\t'                      // POS (1-based)
              << 0 << '\t'                                   // MAPQ
              << "*\t"                   // CIGAR
              << '=' << '\t'                                 // RNEXT
              << qa.pos + 1 << '\t'                          // PNEXT
              << ((read1First) ? -fragLen : fragLen) << '\t' // TLEN
              << *unalignedSeq << '\t'                           // SEQ
              << "*\t"                                       // QUAL
              << numHitFlag << '\n';

    }
    ++alnCtr;
  }

  return 0;
}

#endif
