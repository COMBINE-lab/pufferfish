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

int32_t addCigar(pufferfish::util::cigarGenerator &cigarGen, ksw_extz_t ez, bool beginGap) {
    if (!beginGap) {
        int32_t insertionDeletion = 0;
        for (int i = 0; i < ez.n_cigar; ++i) {
            std::string cigar_type = "";
            if ((ez.cigar[i] & 0xf) == 0) cigar_type = "M";
            if ((ez.cigar[i] & 0xf) == 1) {
                cigar_type = "I";
                insertionDeletion += ez.cigar[i] >> 4;
            }
            if ((ez.cigar[i] & 0xf) == 2) {
                cigar_type = "D";
                insertionDeletion -= ez.cigar[i] >> 4;
            }
            cigarGen.add_item(ez.cigar[i] >> 4, cigar_type);
        }
        return insertionDeletion;
    } else {
        uint32_t gapSize = 0;
        for (int i = ez.n_cigar - 1; i >= 0; --i) {
            std::string cigar_type = "";
            if ((ez.cigar[i] & 0xf) == 0) {
                cigar_type = "M";
                gapSize += ez.cigar[i] >> 4;
            }
            if ((ez.cigar[i] & 0xf) == 1) { cigar_type = "I"; }
            if ((ez.cigar[i] & 0xf) == 2) {
                cigar_type = "D";
                gapSize += ez.cigar[i] >> 4;
            }
            cigarGen.add_item(ez.cigar[i] >> 4, cigar_type);
        }
        return gapSize;
    }
}

std::string getRefSeq(compact::vector<uint64_t, 2> &refseq, uint64_t refAccPos, size_t tpos, uint32_t memlen) {
    if (memlen == 0) return "";
    std::string tseq = "";
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

AlignmentResult PuffAligner::alignRead(std::string read, std::vector<pufferfish::util::MemInfo> &mems, bool perfectChain,
                             bool isFw, size_t tid, AlnCacheMap &alnCache, HitCounters &hctr, bool verbose) {
//    verbose = true;
    uint32_t refExtLength = mopts->refExtendLength;
    bool firstMem = true;
    int32_t lastHitEnd_read = -1;
    int32_t currHitStart_read = 0;
    int64_t lastHitEnd_ref = -1;
    int64_t currHitStart_ref = 0;
    int32_t alignmentScore{std::numeric_limits<int32_t>::lowest()};
    int32_t alignment{0};
    uint32_t openGapLen{0};
    auto& cigarGen = cigarGen_;
    cigarGen.clear();

    std::string cigar = "";
    memset(&ez, 0, sizeof(ksw_extz_t));

    if (verbose) {
        for (size_t x = 0; x < 50; ++x)
            std::cerr << "-";
        std::cerr << " start mems\n";
    }

    uint64_t refAccPos = tid > 0 ? refAccumLengths[tid - 1] : 0;
    uint64_t refTotalLength = refAccumLengths[tid] - refAccPos;

    std::string original_read = read;
    if (!isFw)
        read = pufferfish::util::reverseComplement(read);

    int32_t keyLen = 0;
    std::unique_ptr<char[]> refSeq;
    std::string tseq = "";

    auto rpos = mems[0].rpos;
    auto memlen = mems[0].extendedlen;
    auto readLen = read.length();
    auto tpos = mems[0].tpos;

    if (mems.size() == 0)
        return AlignmentResult();
    currHitStart_read = isFw ? rpos : readLen - (rpos + memlen);
    if (currHitStart_read < 0 or currHitStart_read >= (int32_t) readLen)
        std::cerr << "Should not happen!\n";

    currHitStart_ref = tpos;
    uint32_t refStart;
    if (currHitStart_ref >= currHitStart_read) {
        refStart = currHitStart_ref - currHitStart_read;
    } else { // -- this is always true : else if (currHitStart_ref < currHitStart_read) {
        refStart = 0;
    }
    if (verbose) std::cerr<< "tpos: " << refStart << "\n";
    keyLen = (refStart + readLen < refTotalLength) ? readLen : refTotalLength - refStart;
    tseq = getRefSeq(allRefSeq, refAccPos, refStart, keyLen);

    refSeq.reset(new char[tseq.length() + 1]);
    strcpy(refSeq.get(), tseq.c_str());
    uint64_t hashKey{0};
    bool didHash{false};
    if (!alnCache.empty()) {
        // hash the reference string
        MetroHash64::Hash(reinterpret_cast<uint8_t *>(refSeq.get()), keyLen, reinterpret_cast<uint8_t *>(&hashKey), 0);
        didHash = true;
        // see if we have this hash
        auto hit = alnCache.find(hashKey);
        // if so, we know the alignment score
        if (hit != alnCache.end() and refStart + readLen + refExtLength < refTotalLength) {
            alignmentScore = hit->second.score;
            cigar = hit->second.cigar;
            openGapLen = hit->second.openGapLen;
            alignment = alignmentScore;
        }
    }

    if (alignmentScore == std::numeric_limits<int32_t>::lowest()) {
        if (perfectChain) {
            alignmentScore = readLen * mopts->matchScore;
            alignment = alignmentScore;
            cigarGen.add_item(readLen, "M");
            hctr.skippedAlignments_byCov += 1;
        } else if (mopts->fullAlignment) {
            uint32_t refStart, readStart;
            if (currHitStart_ref > currHitStart_read) {
                refStart = currHitStart_ref - currHitStart_read;
                readStart = 0;
            } else if (currHitStart_ref < currHitStart_read) {
                readStart = currHitStart_read - currHitStart_ref;
                refStart = 0;
            } else {
                readStart = currHitStart_read;
                refStart = currHitStart_ref;
            }
            auto refLength = (refStart + readLen * 2 < refTotalLength) ? readLen * 2 : refTotalLength - refStart;

            tseq = getRefSeq(allRefSeq, refAccPos, refStart, refLength);
            auto readSeq = readStart > 0 ? extractReadSeq(read, readStart, readLen, 1) : read;

            aligner(readSeq.c_str(), readSeq.length(), tseq.c_str(), tseq.length(), &ez,
                    ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
            alignmentScore = ez.mqe; //std::max(ez.mqe, ez.mte);
            openGapLen = addCigar(cigarGen, ez, false);
            if (verbose) {
                std::cerr << "Original read seq:\t" << original_read << "\n";
                std::cerr << "Total alignment with the score\t" << alignmentScore << "\t from position\t" << readStart
                          << "\t on the read:\n" << readSeq << "\n";
                std::cerr << "At\t" << refStart << "\t for length of\t" << refLength << "\t on reference:\n" << tseq
                          << "\n";
                for (int i = 0; i < ez.n_cigar; ++i) // print CIGAR
                    std::cerr<< (ez.cigar[i]>>4) << "MID"[ez.cigar[i] & 0xf];
                std::cerr << "\n";
            }
        } else {

            alignmentScore = 0;
            for (auto mem : mems) {
                rpos = mem.rpos;
                memlen = mem.extendedlen;
                tpos = mem.tpos;

                // To work around seg fault from mems
                if (!firstMem and tpos + memlen - refStart > readLen) {
                    break;
                }
                int score = 0;
                currHitStart_read = isFw ? rpos : readLen - (rpos + memlen);
                currHitStart_ref = tpos;
                if (verbose) {
                    std::cerr << "currHitStart_read:" << currHitStart_read << "\tcurrHitStart_ref:" << currHitStart_ref
                              << "\n";
                    std::cerr << "lastHitStart_read:" << lastHitEnd_read << "\tlastHitStart_ref:" << lastHitEnd_ref
                              << "\n";
                }
                if (firstMem) {
                    lastHitEnd_ref =
                            currHitStart_ref - currHitStart_read >= 0 ? currHitStart_ref - currHitStart_read - 1 : -1;
                } else {
                    if (currHitStart_read <= lastHitEnd_read - (int32_t) k + 1) {
                        if (verbose)
                            std::cerr << "The mems on read and reference are not consistent(1)\t" << original_read
                                      << "\n";
                        break;
                    }
                    if (currHitStart_ref <= lastHitEnd_ref - (int32_t) k + 1) {
                        if (verbose)
                            std::cerr << "The mems on read and reference are not consistent(2)\t" << original_read
                                      << "\n";
                        break;
                    }
                    if (currHitStart_ref > lastHitEnd_ref + (int32_t) readLen) {
                        if (verbose)
                            std::cerr << "The mems on read and reference are not consistent(3)\t" << original_read
                                      << "\n";
                        break;
                    }
                    if (currHitStart_read <= lastHitEnd_read and
                        lastHitEnd_read - currHitStart_read != lastHitEnd_ref - currHitStart_ref) {
                        if (verbose)
                            std::cerr << "to align reads with indels\t" << original_read << "\n";
                        break;
                    }
                    if (currHitStart_ref <= lastHitEnd_ref and
                        lastHitEnd_read - currHitStart_read != lastHitEnd_ref - currHitStart_ref) {
                        if (verbose)
                            std::cerr << "to align reads with indels\t" << original_read << "\n";
                        break;
                    }
                }

                // Performing full alignment to validate the socres
                if (verbose and firstMem) {
                    uint32_t refStart, readStart;
                    if (currHitStart_ref > currHitStart_read) {
                        refStart = currHitStart_ref - currHitStart_read;
                        readStart = 0;
                    } else if (currHitStart_ref < currHitStart_read) {
                        readStart = currHitStart_read - currHitStart_ref;
                        refStart = 0;
                    } else {
                        readStart = currHitStart_read;
                        refStart = currHitStart_ref;
                    }
                    auto refLength = (refStart + readLen * 2 < refTotalLength) ? readLen * 2 : refTotalLength -
                                                                                               refStart;

                    tseq = getRefSeq(allRefSeq, refAccPos, refStart, refLength);
                    auto readSeq = readStart > 0 ? extractReadSeq(read, readStart, readLen, 1) : read;

                    aligner(readSeq.c_str(), readSeq.length(), tseq.c_str(), tseq.length(), &ez,
                            ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
                    alignment = ez.mqe; //std::max(ez.mqe, ez.mte);

                    {
                        std::cerr << "Original read seq:\t" << original_read << "\n";
                        std::cerr << "Total alignment with the score\t" << alignment << "\t from position\t"
                                  << readStart << "\t on the read:\n" << readSeq << "\n";
                        std::cerr << "At\t" << refStart << "\t for length of\t" << refLength << "\t on reference:\n"
                                  << tseq << "\n";
                        for (int i = 0; i < ez.n_cigar; ++i) // print CIGAR
                            std::cerr<< (ez.cigar[i]>>4) << "MID"[ez.cigar[i] & 0xf];
                        std::cerr << "\n";
                    }
                }
                // If there is a gap in the read or reference or both
                bool gapAligned = false;
                if (currHitStart_read > lastHitEnd_read + 1 or currHitStart_ref > lastHitEnd_ref + 1) {
                    gapAligned = true;
                    auto readGapLength = currHitStart_read - lastHitEnd_read - 1;
                    auto refGapLength = currHitStart_ref - lastHitEnd_ref - 1;
                    if (verbose)
                        std::cerr << "readGapLength " << readGapLength << "\trefGapLength " << refGapLength << "\n";
                    auto readSeq = extractReadSeq(read, lastHitEnd_read + 1, lastHitEnd_read + 1 + readGapLength, 1);
                    if (verbose) {
                        for (size_t x = 0; x < 50; ++x)
                            std::cerr << "+";
                        std::cerr << " start gap alignment\n";
                        std::cerr << "read from " << lastHitEnd_read + 1 << " for length " << readGapLength
                                  << " of the read:\t " << readSeq << "\n";
                        std::cerr << "currHitStart_read: " << currHitStart_read << " lastHitEnd_read: "
                                  << lastHitEnd_read << " currHitStart_ref: " << currHitStart_ref << " lastHitEnd_ref: "
                                  << lastHitEnd_ref << "\n";
                        std::cerr << mems[0].rpos << "\t" << readLen << "\t" << rpos << "\t" << memlen << "\n";
                    }
                    if (firstMem) {
                        /*if (readGapLength > refGapLength) {
                          uint32_t startHang = readGapLength - refGapLength;
                          alignmentScore += (-1)*mopts->gapOpenPenalty + (-1)*mopts->gapExtendPenalty*(startHang);
                          cigarGen.add_item(startHang, "I");
                        }
                        if (refGapLength > 0) {
                          if (readGapLength > refGapLength) {
                            readGapLength = refGapLength;
                          }
                                  readSeq = extractReadSeq(read, lastHitEnd_read + 1, lastHitEnd_read + 1 + readGapLength, 1);
                                      std::reverse(readSeq.begin(), readSeq.end());
                                    auto refStartPos = lastHitEnd_ref > refExtLength + refGapLength ?  lastHitEnd_ref + 1 - refExtLength - refGapLength : 0;
                          refGapLength = lastHitEnd_ref - refStartPos;
                                      tseq = getRefSeq(allRefSeq, refAccPos, refStartPos, refGapLength);
                                      std::reverse(tseq.begin(), tseq.end());
                                      aligner(readSeq.c_str(), readSeq.length(), tseq.c_str(), tseq.length(), &ez, ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
                                      score = std::max(ez.mqe, ez.mte);
                          addCigar(cigarGen, ez, true);
                        }*/
                        ///// no need because of line 512: lastHitEnd_ref = currHitStart_ref > currHitStart_read ? currHitStart_ref - currHitStart_read -1: -1;
                        // Not extending the reference beyond its beginning
                        auto refStartPos = lastHitEnd_ref > refExtLength ? lastHitEnd_ref + 1 - refExtLength : 0;
                        refGapLength = currHitStart_ref - refStartPos;
                        /*if (refGapLength < readGapLength) {
                            AlignmentResult aln;
                            aln.score = -1;
                            return aln;
                        }*/
                        // SOFT CLIPPING for reads mapping with a left hanger at the beginning of a transcript
                        if (refGapLength > 0) {
                            // We reverse the strings because of ksw force alignment from the beginning of the sequences
                            std::reverse(readSeq.begin(), readSeq.end());
                            tseq = getRefSeq(allRefSeq, refAccPos, refStartPos, refGapLength);
                            std::reverse(tseq.begin(), tseq.end());
                            ksw_reset_extz(&ez);
                            aligner(readSeq.c_str(), readSeq.length(), tseq.c_str(), tseq.length(), &ez,
                                    ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
                            score = ez.mqe; //std::max(ez.mqe, ez.mte);
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
                            /*if (readGapLength > refGapLength + insertionDeletion) {
                                int32_t startHang = readGapLength - refGapLength - insertionDeletion;
                                alignmentScore +=
                                        (-1) * mopts->gapOpenPenalty + (-1) * mopts->gapExtendPenalty * (startHang);
                                cigarGen.add_item(startHang, "I");
                            }*/
                        } else {
                            score = 0;
                            int32_t startHang = readGapLength - refGapLength;
                            alignmentScore +=
                                    (-1) * mopts->gapOpenPenalty + (-1) * mopts->gapExtendPenalty * (startHang);
                            cigarGen.add_item(startHang, "I");
                        }
                    } else {
                        //Considering the case where the gap in either the read or reference is empty
                        if (readGapLength <= 0) {
                            if (verbose)
                                std::cerr << "read gap lengths was of length 0 for the ref gap length of \t"
                                          << refGapLength << "\n";
                            score = (-1) * mopts->gapOpenPenalty +
                                    (-1) * mopts->gapExtendPenalty * (refGapLength + std::abs(readGapLength));
                            cigarGen.add_item(refGapLength, "D");
                        } else if (refGapLength <= 0) {
                            if (verbose)
                                std::cerr << "ref gap lengths was of length 0 for the read gap length of \t"
                                          << readGapLength << "\n";
                            score = (-1) * mopts->gapOpenPenalty +
                                    (-1) * mopts->gapExtendPenalty * (readGapLength + std::abs(refGapLength));
                            cigarGen.add_item(readGapLength, "I");
                        } else {
                            //tseq = getRefSeq(allRefSeq, refAccPos, lastHitEnd_ref + 1, refGapLength);
                            //auto score1 = aligner(readSeq.c_str(), readSeq.length(), tseq.c_str(), tseq.length(), &ez, ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::GLOBAL>());
                            if (lastHitEnd_ref + 1 - refStart < 0 or
                                lastHitEnd_ref + 1 - refStart >= (int32_t) std::strlen(refSeq.get()))
                                std::cerr << "Should not happen: lastHitEnd_ref is " << lastHitEnd_ref
                                          << "and refStart is " << refStart << ", but refSeq length is "
                                          << std::strlen(refSeq.get()) << "\n";
                            char *refSeq1 = refSeq.get() + lastHitEnd_ref + 1 - refStart;
                            score = aligner(readSeq.c_str(), readSeq.length(), refSeq1, refGapLength, &ez,
                                            ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::GLOBAL>());
                            addCigar(cigarGen, ez, false);
                            if (verbose) {
                                tseq = getRefSeq(allRefSeq, refAccPos, lastHitEnd_ref + 1, refGapLength);
                                std::cerr << "read from " << lastHitEnd_read + 1 << " with the gap length of "
                                          << readGapLength << " is \t" << readSeq << "\n";
                                std::cerr << "reference from " << lastHitEnd_ref + 1 << " with the length "
                                          << refGapLength << " is \t" << tseq << "\n";
                                for (int i = 0; i < ez.n_cigar; ++i) // print CIGAR
                                    std::cerr<< (ez.cigar[i]>>4) << "MID"[ez.cigar[i] & 0xf];
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

                auto longest_match_extension =
                        currHitStart_ref + memlen - lastHitEnd_ref > currHitStart_read + memlen - lastHitEnd_read ?
                        currHitStart_ref + memlen - lastHitEnd_ref : currHitStart_read + memlen - lastHitEnd_read;
                score = (currHitStart_ref <= lastHitEnd_ref) ? mopts->matchScore * (longest_match_extension - 1) :
                        mopts->matchScore * memlen;
                auto num_match = (currHitStart_ref <= lastHitEnd_ref) ? (longest_match_extension - 1) : memlen;
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
                        if (currHitStart_ref - lastHitEnd_ref > currHitStart_read - lastHitEnd_read)
                            cigarGen.add_item(gapLength, "D");
                        else
                            cigarGen.add_item(gapLength, "I");
                    }
                    if (verbose) {
                        if (!gapAligned)
                            std::cerr << "penalty\t" << penalty << "\n";
                        auto readSeq = extractReadSeq(read, currHitStart_read, currHitStart_read + memlen, 1);
                        auto tseq1 = getRefSeq(allRefSeq, refAccPos, currHitStart_ref, memlen);
                        std::cerr << "read from " << currHitStart_read << "\t with the sequence \n" << readSeq << "\n";
                        std::cerr << "orientation of\t" << isFw << " at " << currHitStart_ref << " on reference " << tid
                                  << " for the length of " << memlen << "\t with sequence:\n" << tseq1 << "\n";
                        std::cerr << "with the score1:\t" << score << "\n";
                    }
                }
                alignmentScore += score;
                cigarGen.add_item(num_match, "M");

                lastHitEnd_read = currHitStart_read + memlen - 1;
                lastHitEnd_ref = tpos + memlen - 1;
                if (lastHitEnd_ref - refStart + 1 > (int32_t) std::strlen(refSeq.get()) + 1)
                    std::cerr << "Should not happen: lastHitEnd_ref is " << lastHitEnd_ref << " and refStart is "
                              << refStart << ", but refSeq length is " << std::strlen(refSeq.get()) << "\n";
                firstMem = false;
            }

            // Try and align any remaining sequence on the read
            if (lastHitEnd_read < (int32_t) readLen - 1) {
                auto readGapLength = readLen - 1 - lastHitEnd_read;
                auto refGapLength = lastHitEnd_ref + 1 + readGapLength + refExtLength < refTotalLength ? readGapLength +
                                                                                                         refExtLength :
                                    refTotalLength - lastHitEnd_ref - 1;
                if (lastHitEnd_ref + readGapLength > refTotalLength) {
                    AlignmentResult aln;
                    aln.score = -1;
                    aln.cigar = cigar;
                    aln.openGapLen = openGapLen;
                    return aln;
                }

                // SOFT CLIPPING for reads mapping with a right hanger at the beginning of a transcript
                if (refGapLength != 0) {
                    auto readSeq = extractReadSeq(read, lastHitEnd_read + 1, lastHitEnd_read + 1 + readGapLength, 1);
                    auto refStartPos = lastHitEnd_ref + 1;
                    tseq = getRefSeq(allRefSeq, refAccPos, refStartPos, refGapLength);
                    aligner(readSeq.c_str(), readSeq.length(), tseq.c_str(), refGapLength, &ez,
                            ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
                    auto score = ez.mqe; //std::max(ez.mqe, ez.mte);
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
                        cigarGen.add_item(endHang, "I");
                    }
                } else {
                    int32_t endHang = readGapLength;
                    alignmentScore += (-1) * mopts->gapOpenPenalty + (-1) * mopts->gapExtendPenalty * (endHang);
                    cigarGen.add_item(endHang, "I");
                }
            }
        }
        if (verbose)
            std::cerr << "alignmentScore: " << alignmentScore << "\n";
//        alignmentScore+=readLen*2; //addedByFatemeh
        bool cigar_fixed{false};
        cigar = cigarGen.get_cigar(readLen, cigar_fixed);
        /*if (cigar_fixed) {
          std::cerr<<read<<"\n";
        } */
        //std::cerr << cigar << " " << readLen << "\n";
        if (cigar_fixed) hctr.cigar_fixed_count++;
        if (multiMapping) { // don't bother to fill up a cache unless this is a multi-mapping read
            if (!didHash)
              MetroHash64::Hash(reinterpret_cast<uint8_t *>(refSeq.get()), keyLen, reinterpret_cast<uint8_t *>(&hashKey),
                                  0);
            AlignmentResult aln;
            aln.score = alignmentScore;
            aln.cigar = cigar;
            aln.openGapLen = openGapLen;
            alnCache[hashKey] = aln;
        }
    } else {
        hctr.skippedAlignments_byCache += 1;
    }
    //delete refSeq;
    if (verbose) {
        std::cerr << "alignmentScore\t" << alignmentScore << "\talignmment\t" << alignment << "\n";
        std::cerr << "cigar\t" << cigar << "\n";
    }
    //if (alignmentScore > 150 and alignment > alignmentScore)// != and std::abs(alignmentScore - alignment) > 5)
    //  std::cerr<< alignmentScore << " " << alignment << " " << original_read << "\n";

    return AlignmentResult{alignmentScore, cigar, openGapLen};
}


int32_t PuffAligner::calculateAlignments(std::string& read_left, std::string& read_right, pufferfish::util::JointMems& jointHit, HitCounters& hctr, bool verbose) {

    auto tid = jointHit.tid;
    double optFrac{mopts->minScoreFraction};

    auto threshold = [&, optFrac] (uint64_t len) -> double {
        return (mopts->mimicBt2Default or !mopts->matchScore)?(-0.6+-0.6*len):optFrac*mopts->matchScore*len;
    };

    if (jointHit.isOrphan()) {
      if (verbose) {std::cerr << "orphan\n";}
        hctr.totalAlignmentAttempts += 1;

        std::string read_orphan = jointHit.isLeftAvailable() ? read_left : read_right;
//        int32_t maxScore = mopts->matchScore * read_orphan.length();
        AlignmentResult ar = alignRead(read_orphan, jointHit.orphanClust()->mems, jointHit.orphanClust()->perfectChain,
                                       jointHit.orphanClust()->isFw, tid, alnCacheLeft, hctr, verbose);
        jointHit.orphanClust()->coverage =
                ar.score > threshold(read_orphan.length())  ? ar.score : std::numeric_limits<decltype(ar.score)>::min();
        jointHit.orphanClust()->cigar = ar.cigar;
        jointHit.orphanClust()->openGapLen = ar.openGapLen;
        jointHit.alignmentScore = static_cast<int32_t >(jointHit.orphanClust()->coverage);
        if (jointHit.alignmentScore < 0 and verbose) {
          std::cerr << read_orphan.length() << " " << threshold(read_orphan.length()) << " " << ar.score << "\n";
        }
        return jointHit.alignmentScore;
    } else {
        hctr.totalAlignmentAttempts += 2;
//        int32_t maxLeftScore = mopts->matchScore * read_left.length();
//        int32_t maxRightScore = mopts->matchScore * read_right.length();
        if (verbose) { std::cerr << "left\n"; }
        AlignmentResult ar_left = alignRead(read_left, jointHit.leftClust->mems, jointHit.leftClust->perfectChain,
                                            jointHit.leftClust->isFw, tid, alnCacheLeft, hctr, verbose);
        if (verbose) { std::cerr << "right\n"; }
        AlignmentResult ar_right = alignRead(read_right, jointHit.rightClust->mems, jointHit.rightClust->perfectChain,
                                             jointHit.rightClust->isFw, tid, alnCacheRight, hctr, verbose);
        auto score_left = ar_left.score > threshold(read_left.length()) ? ar_left.score
                                                                        : std::numeric_limits<decltype(ar_left.score)>::min();
        auto score_right = ar_right.score > threshold(read_right.length()) ? ar_right.score
                                                                           : std::numeric_limits<decltype(ar_right.score)>::min();
        jointHit.alignmentScore = (score_left == std::numeric_limits<decltype(ar_left.score)>::min() or
                score_right == std::numeric_limits<decltype(ar_right.score)>::min())?
                                  std::numeric_limits<decltype(ar_right.score)>::min() : score_left + score_right;
        jointHit.leftClust->cigar = ar_left.cigar;
        jointHit.leftClust->coverage = score_left;
        jointHit.leftClust->openGapLen = ar_left.openGapLen;
        jointHit.rightClust->coverage = score_right;
        jointHit.rightClust->cigar = ar_right.cigar;
        jointHit.rightClust->openGapLen = ar_right.openGapLen;
        return jointHit.alignmentScore;
    }
}

int32_t PuffAligner::calculateAlignments(pufferfish::util::JointMems &jointHit, HitCounters &hctr, bool verbose) {
  return calculateAlignments(read_left_, read_right_, jointHit, hctr, verbose);
}

bool PuffAligner::recoverSingleOrphan(std::string& read_left, std::string& read_right, pufferfish::util::MemCluster& clust, std::vector<pufferfish::util::MemCluster> &recoveredMemClusters, uint32_t tid, bool anchorIsLeft, bool verbose) {
  int32_t anchorLen = anchorIsLeft ? read_left_.length() : read_right_.length();
  auto tpos = clust.mems[0].tpos;
  auto anchorStart = clust.mems[0].isFw ? clust.mems[0].rpos : anchorLen - (clust.mems[0].rpos + clust.mems[0].extendedlen);
  uint32_t anchorPos = tpos >= anchorStart ? tpos - anchorStart : 0;

  bool recovered_fwd;
  uint32_t recovered_pos=-1;

  auto* r1 = read_left_.data();
  auto* r2 = read_right_.data();
  auto l1 = static_cast<int32_t>(read_left_.length());
  auto l2 = static_cast<int32_t>(read_right_.length());
  const char* rptr{nullptr};
  bool anchorFwd{clust.isFw};
  int32_t startPos = -1, maxDist = -1, otherLen = -1, rlen = -1;
  std::string* otherReadPtr{nullptr};
  const char* otherRead{nullptr};
  char* otherReadRC{nullptr};
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
    otherReadPtr = &read_right_;
    otherRead = r2;
    otherReadRC = r2rc;
  } else {
    anchorLen = l2;
    otherLen = l1;
    maxDist = maxDistLeft;
    otherReadPtr = &read_left_;
    otherRead = r1;
    otherReadRC = r1rc;
  }

  uint64_t refAccPos = tid > 0 ? refAccumLengths[tid - 1] : 0;
  uint64_t refLength = refAccumLengths[tid] - refAccPos;

  if (anchorFwd) {
    if (!otherReadRC){
      auto read = pufferfish::util::reverseComplement(*otherReadPtr);
      otherReadRC = const_cast<char*>(read.data());
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
  auto tseq = getRefSeq(allRefSeq, refAccPos, startPos, windowLength);
  windowSeq.reset(new char[tseq.length() + 1]);
  strcpy(windowSeq.get(), tseq.c_str());

  EdlibAlignResult result = edlibAlign(rptr, rlen, windowSeq.get(), windowLength,
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


bool PuffAligner::recoverSingleOrphan(pufferfish::util::MemCluster& clust, std::vector<pufferfish::util::MemCluster> &recoveredMemClusters, uint32_t tid, bool anchorIsLeft, bool verbose) {
  recoverSingleOrphan(read_left_, read_right_, clust, recoveredMemClusters, tid, anchorIsLeft, verbose);
}
