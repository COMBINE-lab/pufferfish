//I am not sure what are the things I would need soon
//let's go with what we have
#include <iostream>
#include <mutex>
#include <vector>
#include <random>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <thread>
#include <tuple>
#include <sstream>
#include <fstream>
#include <iostream>
#include <tuple>
#include <memory>
#include <cstring>
#include <queue>

//we already have timers

//cereal include
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <sparsepp/spp.h>


#include "spdlog/spdlog.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/sinks/stdout_sinks.h"
#include "spdlog/sinks/ansicolor_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"

#include "FastxParser.hpp"

//index header
#include "PufferfishAligner.hpp"
#include "ProgOpts.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "PufferfishLossyIndex.hpp"
#include "Kmer.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "SpinLock.hpp"
#include "MemCollector.hpp"
#include "SAMWriter.hpp"
#include "RefSeqConstructor.hpp"
#include "KSW2Aligner.hpp"
#include "zstr/zstr.hpp"

#include "metro/metrohash64.h"
#include "tsl/hopscotch_map.h"


#define START_CONTIG_ID ((uint32_t)-1)
#define END_CONTIG_ID ((uint32_t)-2)

#define MATCH_SCORE 1
#define MISMATCH_SCORE -1
#define GAP_SCORE -1

#define EPS 5

using paired_parser = fastx_parser::FastxParser<fastx_parser::ReadPair>;
using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

using HitCounters = util::HitCounters;
using QuasiAlignment = util::QuasiAlignment;
using MateStatus = util::MateStatus;

using MutexT = std::mutex;

template<typename PufferfishIndexT>
void printMemInfo(size_t tid,
                  const std::vector<util::MemCluster>::iterator lclust,
                  const std::vector<util::MemCluster>::iterator rclust,
                  PufferfishIndexT &pfi) {
    std::cerr << "\ntid:" << tid << "\t" << pfi.refName(tid) << "\n";
    std::cerr << "left:" << lclust->isFw << " size:" << lclust->mems.size() << " cov:" << lclust->coverage << "\n";
    for (size_t i = 0; i < lclust->mems.size(); i++) {
        std::cerr << "--- t" << lclust->mems[i].tpos << " r"
                  << lclust->mems[i].memInfo->rpos << " cid:"
                  << lclust->mems[i].memInfo->cid << " cpos: "
                  << lclust->mems[i].memInfo->cpos << " len:"
                  << lclust->mems[i].memInfo->memlen << " fw:"
                  << lclust->mems[i].memInfo->cIsFw << "\n";
    }
    std::cerr << "\nright:" << rclust->isFw << " size:" << rclust->mems.size() << " cov:" << rclust->coverage << "\n";
    for (size_t i = 0; i < rclust->mems.size(); i++) {
        std::cerr << "--- t" << rclust->mems[i].tpos << " r"
                  << rclust->mems[i].memInfo->rpos << " cid:"
                  << rclust->mems[i].memInfo->cid << " cpos: "
                  << rclust->mems[i].memInfo->cpos << " len:"
                  << rclust->mems[i].memInfo->memlen << " fw:"
                  << rclust->mems[i].memInfo->cIsFw << "\n";
    }
}


template<typename PufferfishIndexT>
void joinReadsAndFilter(spp::sparse_hash_map<size_t, std::vector<util::MemCluster>> &leftMemClusters,
                        spp::sparse_hash_map<size_t, std::vector<util::MemCluster>> &rightMemClusters,
                        std::vector<util::JointMems> &jointMemsList,
                        uint32_t maxFragmentLength,
                        uint32_t perfectCoverage,
                        double coverageRatio,
                        bool noDiscordant,
                        bool noOrphans,
                        PufferfishIndexT &pfi,
                        bool verbose = false) {
    if (verbose)
        std::cerr << "\n[JOINREADSANDFILTER]\n";

    uint64_t maxLeft{0}, maxRight{0}, maxLeftCnt{0}, maxRightCnt{0};
    bool isMaxLeftAndRight = false;
    if (!noOrphans) {
        for (auto &kv : leftMemClusters) {
            if (verbose)
                std::cerr << "\ntid:" << kv.first << "\t" << pfi.refName(kv.first) << "\n";
            auto &lClusts = kv.second;
            for (auto clust = lClusts.begin(); clust != lClusts.end(); clust++) {
                if (verbose) {
                    std::cerr << "left: " << clust->isFw << " cov " << clust->coverage << "\n";
                    for (size_t i = 0; i < clust->mems.size(); i++) {
                        std::cerr << "--- t" << clust->mems[i].tpos << " r"
                                  << clust->mems[i].memInfo->rpos << " cid:"
                                  << clust->mems[i].memInfo->cid << " cpos: "
                                  << clust->mems[i].memInfo->cpos << " len:"
                                  << clust->mems[i].memInfo->memlen << " fw:"
                                  << clust->mems[i].memInfo->cIsFw << "\n";
                    }
                }
                if (maxLeft == clust->coverage) {
                    maxLeftCnt += 1;
                } else if (maxLeft < clust->coverage) {
                    maxLeft = clust->coverage;
                    maxLeftCnt = 1;
                }
            }
        }
        for (auto &kv : rightMemClusters) {
            if (verbose)
                std::cerr << "\ntid:" << kv.first << "\t" << pfi.refName(kv.first) << "\n";

            auto &rClusts = kv.second;
            for (auto clust = rClusts.begin(); clust != rClusts.end(); clust++) {
                if (verbose) {
                    std::cerr << "right: " << clust->isFw << " cov " << clust->coverage << "\n";
                    for (size_t i = 0; i < clust->mems.size(); i++) {
                        std::cerr << "--- t" << clust->mems[i].tpos << " r"
                                  << clust->mems[i].memInfo->rpos << " cid:"
                                  << clust->mems[i].memInfo->cid << " cpos: "
                                  << clust->mems[i].memInfo->cpos << " len:"
                                  << clust->mems[i].memInfo->memlen << " fw:"
                                  << clust->mems[i].memInfo->cIsFw << "\n";
                    }
                }
                if (maxRight == clust->coverage) {
                    maxRightCnt += 1;
                } else if (maxRight < clust->coverage) {
                    maxRight = clust->coverage;
                    maxRightCnt = 1;
                }
            }
        }
        //std::cerr << "\nnow after joining " << round << "\n";

    }
    //orphan reads should be taken care of maybe with a flag!
    //uint32_t perfectCoverage{2*readLen};
    uint32_t maxCoverage{0};
    short round{0};

    while (round == 0 || (round == 1 && !jointMemsList.size() && !noDiscordant)) {
        for (auto &leftClustItr : leftMemClusters) {
            // reference id
            size_t tid = leftClustItr.first;
            // left mem clusters
            auto &lClusts = leftClustItr.second;
            // right mem clusters for the same reference id
            auto &rClusts = rightMemClusters[tid];
            // Compare the left clusters to the right clusters to filter by positional constraints

            for (auto lclust = lClusts.begin(); lclust != lClusts.end(); lclust++) {

                for (auto rclust = rClusts.begin(); rclust != rClusts.end(); rclust++) {

                    // if both the left and right clusters are oriented in the same direction, skip this pair
                    // NOTE: This should be optional as some libraries could allow this.
                    if (round == 0 && lclust->isFw == rclust->isFw) { // if priority 0, ends should be concordant
                        continue;
                    }
                    // FILTER 1
                    // filter read pairs based on the fragment length which is approximated by the distance between the left most start and right most hit end
                    size_t fragmentLen = rclust->lastRefPos() + rclust->lastMemLen() - lclust->firstRefPos();
                    if (lclust->firstRefPos() > rclust->firstRefPos()) {
                        fragmentLen = lclust->lastRefPos() + lclust->lastMemLen() - rclust->firstRefPos();
                    }
                    // FILTERING fragments with size smaller than maxFragmentLength
                    // FILTER just in case of priority 0 (round 0)
                    /*if (verbose) {
                        std::cerr << "fraglen: " << fragmentLen << " " << maxFragmentLength << "\n";
                    }*/
                    if (fragmentLen < maxFragmentLength || round > 0) {
                        // This will add a new potential mapping. Coverage of a mapping for read pairs is left->coverage + right->coverage
                        // If we found a perfect coverage, we would only add those mappings that have the same perfect coverage
                        auto totalCoverage = lclust->coverage + rclust->coverage;
                        if (totalCoverage >= coverageRatio * maxCoverage or
                            totalCoverage == perfectCoverage) {//}|| (lclust->coverage + rclust->coverage) == ) {
                            jointMemsList.emplace_back(tid, lclust, rclust, fragmentLen);
                            /* if (verbose) {
                              printMemInfo(tid, lclust, rclust, pfi);
                            } */
                            uint32_t currCoverage = jointMemsList.back().coverage();
                            if (maxCoverage < currCoverage) {
                                maxCoverage = currCoverage;
                                if (lclust->coverage < maxLeft || rclust->coverage < maxRight) {
                                    isMaxLeftAndRight = false;
                                } else {
                                    isMaxLeftAndRight = true;
                                }
                            }
                        }
                    }
                }
            }
        }
        round++;
    }
    // If we couldn't find any pair and we are allowed to add orphans
    if (verbose) {
        std::cerr << "isMaxLeftAndRight:" << isMaxLeftAndRight << "\n";
    }
    if (!noOrphans && (!jointMemsList.size() || !isMaxLeftAndRight || maxLeftCnt > 1 || maxRightCnt > 1)) {
        //std::cerr << "Orphans\n";
        auto orphanFiller = [&jointMemsList, &maxCoverage]
                (spp::sparse_hash_map<size_t, std::vector<util::MemCluster>> &memClusters,
                 bool isLeft) {
            // fragmentLen is set to 0
            std::vector<util::MemCluster> fake;
            for (auto &clustItr : memClusters) {
                // reference id
                size_t tid = clustItr.first;
                // left mem clusters
                auto &Clusts = clustItr.second;
                for (auto clust = Clusts.begin(); clust != Clusts.end(); clust++) {
                    if (isLeft) {
                        jointMemsList.emplace_back(tid, clust, fake.end(), 0, MateStatus::PAIRED_END_LEFT);
                    } else {
                        jointMemsList.emplace_back(tid, fake.end(), clust, 0, MateStatus::PAIRED_END_RIGHT);
                    }
                    uint32_t currCoverage = jointMemsList.back().coverage();
                    if (maxCoverage < currCoverage) {
                        maxCoverage = currCoverage;
                    }
                }
            }
        };
        //if (leftMemClusters.size() == 0 or rightMemClusters.size() == 0) {
        orphanFiller(leftMemClusters, true);
        orphanFiller(rightMemClusters, false);
        //}
        /* else {
          std::cerr << "discarded " << leftMemClusters.size() << "," << rightMemClusters.size() << "\n";
        } */

    }
    if (verbose) {
        std::cerr << "\nBefore filter " << jointMemsList.size() << " maxCov:" << maxCoverage << "\n";
        std::cerr << "\n" << jointMemsList.size() << " maxCov:" << maxCoverage << "\n";
    }
    // FILTER 2
    // filter read pairs that don't have enough base coverage
    // (i.e. their coverage is less than coverage ratio times the max coverage for this read)
    // if we've found a perfect match, we will erase any match that is not perfect
    /* if (maxCoverage == perfectCoverage) {
      jointMemsList.erase(std::remove_if(jointMemsList.begin(), jointMemsList.end(),
                                          [&maxCoverage](util::JointMems& pairedReadMems) -> bool {
                                            return pairedReadMems.coverage() != maxCoverage;
                                          }),
                           jointMemsList.end());
    }
    else  */{// ow, we will go with the heuristic that just keep those mappings that have at least
        // the required fraction of the best coverage.
        jointMemsList.erase(std::remove_if(jointMemsList.begin(), jointMemsList.end(),
                                           [&maxCoverage, coverageRatio](util::JointMems &pairedReadMems) -> bool {
                                               return pairedReadMems.coverage() < coverageRatio * maxCoverage;
                                           }),
                            jointMemsList.end());


        /* uint64_t minCoverage = 33;//0.8 * (perfectCoverage / 2);
        jointMemsList.erase(std::remove_if(jointMemsList.begin(), jointMemsList.end(),
                                         [&minCoverage](util::JointMems& pairedReadMems) -> bool {
                                           return pairedReadMems.coverage() < minCoverage ;
                                         }),
                          jointMemsList.end()); */

    }
    //TODO you should take care of this case later. Don't forget
    /*if (coverageRatio == 1) {
        std::set<std::vector<util::UniMemInfo>::iterator> lvalidUnimemSet, rvalidUnimemSet;
        if (jointMemsList.size() > 0) {
            auto &first = jointMemsList[0];
            if (first.isLeftAvailable()) {
                for (auto &mem : first.leftClust->mems) {
                    lvalidUnimemSet.insert(mem.memInfo);
                }
            }
            if (first.isRightAvailable()) {
                for (auto &mem : first.rightClust->mems) {
                    rvalidUnimemSet.insert(mem.memInfo);
                }
            }
            for (auto &jointMem : jointMemsList) {
                if (jointMem.isLeftAvailable()) {
                    for (auto &mem : jointMem.leftClust->mems) {
                        if (lvalidUnimemSet.find(mem.memInfo) == lvalidUnimemSet.end()) {
                            jointMemsList.clear();
                            return;
                        }
                    }
                }
                if (jointMem.isRightAvailable()) {
                    for (auto &mem : jointMem.rightClust->mems) {
                        if (rvalidUnimemSet.find(mem.memInfo) == rvalidUnimemSet.end()) {
                            jointMemsList.clear();
                            return;
                        }
                    }
                }
            }
        }
    }*/
    if (verbose) {
        std::cerr << "\nFinal stat: " << jointMemsList.size() << " had coverage >= " << coverageRatio * maxCoverage
                  << "\n";
        std::cerr << jointMemsList.size() << "\n";
        for (auto &jointMems : jointMemsList) {
            std::cerr << " -->" << pfi.refName(jointMems.tid) << "\n\t\t";
            for (auto &cl : jointMems.leftClust->mems) {
                std::cerr << cl.tpos << " ";
            }
            std::cerr << "\n\t\t";
            for (auto &cl : jointMems.rightClust->mems) {
                std::cerr << cl.tpos << " ";
            }
            std::cerr << "\n";
        }

    }
    if (verbose)
        std::cerr << "\n[END OF JOINREADSANDFILTER]\n";

    /* std::cerr << "\nFinal stat: " << jointMemsList.size() << " had coverage >= " << coverageRatio*maxCoverage << "\n";
    for (auto jointMems : jointMemsList) {
      std::cerr << pfi.refName(jointMems.tid) << "\t";
    }
    std::cerr << "\n"; */
}

std::string extractReadSeq(const std::string readSeq, uint32_t rstart, uint32_t rend, bool isFw) {
    std::string subseq = readSeq.substr(rstart, rend - rstart);
    if (isFw)
        return subseq;
    return util::reverseComplement(subseq); //reverse-complement the substring
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

char* getRefSeq(sdsl::int_vector<2>& refseq, uint64_t refAccPos, size_t tpos, uint32_t memlen, std::string& tseq) {
	//std::string tseq = "";
	uint64_t bucket_offset = (refAccPos+tpos)*2;
	auto len_on_vector = memlen * 2;
	for (uint32_t w = 0; w <= len_on_vector / 64; w++) {
		uint32_t len = std::min((uint32_t)64, len_on_vector - w*64);
		uint64_t word = refseq.get_int(bucket_offset, len);
		for (uint32_t i=0; i < len; i += 2) {
			uint8_t next_bits = ( (word >> i) & 0x03);
			char next = 'A';
			if (next_bits == 1)
				next = 'C';
			else if(next_bits == 2)
				next = 'G';
			else if(next_bits == 3)
				next = 'T';
			tseq += next;
			//result += next;
		}
		bucket_offset += len;
	}
	char *cstr = new char[tseq.length() + 1];
	strcpy(cstr, tseq.c_str());

	return cstr;
}

using AlnCacheMap = tsl::hopscotch_map<uint64_t, int32_t, PassthroughHash>;

int32_t PufferfishAligner::alignRead(std::string read, std::vector<util::MemInfo> mems, bool perfectChain, bool isFw, size_t tid, AlnCacheMap& alnCache, bool verbose) {
	uint32_t refExtLength = mopts->refExtendLength;
	bool firstMem = true;
	int32_t lastHitEnd_read = -1;
	int32_t currHitStart_read = 0;
	int64_t lastHitEnd_ref = -1;
	int64_t currHitStart_ref = 0;
	int32_t alignmentScore {std::numeric_limits<int32_t>::lowest()};
	int alignment = 0;

	if (verbose){
		for(size_t x=0; x<50; ++x)
			std::cerr<<"-";
		std::cerr<<" start mems\n";
	}

	auto refAccPos = refAccumLengths[tid-1];
	auto refTotalLength = refAccumLengths[tid] - refAccumLengths[tid-1];

	std::string original_read = read;
	if ( !isFw )
		read = util::reverseComplement(read);

	int32_t keyLen;
	char* refSeq;
	std::string tseq = "";

	auto mem = mems[0];
	currHitStart_read = mem.isFw ? mem.memInfo->rpos : read.length() - (mem.memInfo->rpos + mem.memInfo->memlen);
	currHitStart_ref = mem.tpos;
	uint32_t refStart, readStart;
	if ( currHitStart_ref > currHitStart_read) {
		refStart = currHitStart_ref - currHitStart_read;
	} else if (currHitStart_ref < currHitStart_read) {
		refStart = 0;
	} else{
		refStart = currHitStart_ref;
	}
	keyLen = (refStart + read.length() < refTotalLength ) ? read.length() : refTotalLength - refStart;
	refSeq = getRefSeq(allRefSeq, refAccPos, refStart, keyLen, tseq);

    uint64_t hashKey{0};
    bool didHash{false};
    if (!alnCache.empty()) {
    	// hash the reference string
      	MetroHash64::Hash(reinterpret_cast<uint8_t*>(refSeq), keyLen, reinterpret_cast<uint8_t*>(&hashKey), 0);
      	didHash = true;
      	// see if we have this hash
      	auto hit = alnCache.find(hashKey);
      	// if so, we know the alignment score
      	if (hit != alnCache.end()) {
      		alignmentScore = hit->second;
      	}
    }

	if (alignmentScore == std::numeric_limits<int32_t>::lowest()) {

		if (false) { //perfectChain) {
			alignmentScore = read.length() * mopts->matchScore; 
		} else {

		alignmentScore = 0;
		for (auto mem : mems) {
			int score = 0;
			currHitStart_read = mem.isFw ? mem.memInfo->rpos : read.length() - (mem.memInfo->rpos + mem.memInfo->memlen);
			currHitStart_ref = mem.tpos;
			if (verbose)
				std::cerr<<"currHitStart_read:" << currHitStart_read << "\tcurrHitStart_ref:" << currHitStart_ref << "\n";
			if (firstMem) {
				lastHitEnd_ref = currHitStart_ref - currHitStart_read >= 0 ? currHitStart_ref - currHitStart_read - 1 : -1;
			}
			// To work around a possible bug in the chaining algorithm, next kmer match occuring earlier
			// Example -> CGGGCATGGTGGCTCACACCTGTAATCCCAGCACTTTGGGAGGCCAAGGTGGGTGGATCATGAGGTCAGGAATTCGAGAATAGCCTGGCCAACATGGTGA
			if (currHitStart_read < lastHitEnd_read - k and !firstMem) {
				break;
			}
			// To work around a possible bug in the chaining algorithm, next kmer match far from the current match
			// Example -> GATGCAGTGGCTCATGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGCAGGCAGATCACTTGAGATCAGGAGTTCGAGACAAGCCTGGCTAAAATGGTGA
			if (currHitStart_ref > lastHitEnd_ref + read.length()) {
				break;
			}

			// Performing full alignment to validate the socres
			if (firstMem and verbose) {
				uint32_t refStart, readStart;
				if ( currHitStart_ref > currHitStart_read) {
					refStart = currHitStart_ref - currHitStart_read;
					readStart = 0;
				} else if (currHitStart_ref < currHitStart_read) {
					readStart = currHitStart_read - currHitStart_ref;
					refStart = 0;
				} else{
					readStart = currHitStart_read;
					refStart = currHitStart_ref;
				}
				auto refLength = (refStart + read.length()*2 < refTotalLength ) ? read.length()*2 : refTotalLength - refStart;
	
				tseq = "";
				getRefSeq(allRefSeq, refAccPos, refStart, refLength, tseq);
				auto readSeq = readStart > 0 ? extractReadSeq(read, readStart, read.length(), 1) : read;
				
				aligner(readSeq.c_str(), readSeq.length(), tseq.c_str(), tseq.length(), &ez, ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
				alignment = std::max(ez.mqe, ez.mte);
				if (verbose) {
					std::cerr << "Original read seq:\t" << original_read << "\n";
					std::cerr << "Total alignment with the score\t"<< alignment <<"\t from position\t" << readStart << "\t on the read:\n" << readSeq <<"\n";
					std::cerr << "At\t" << refStart << "\t for length of\t" << refLength << "\t on reference:\n" << tseq << "\n";
					for (int i = 0; i < ez.n_cigar; ++i) // print CIGAR
						printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
					std::cerr<<"\n";
				}
			}
			// If there is a gap in the read or reference or both
			if (currHitStart_read > lastHitEnd_read + 1 or currHitStart_ref > lastHitEnd_ref + 1) {
				auto readGapLength = currHitStart_read - lastHitEnd_read - 1;
				auto refGapLength = currHitStart_ref - lastHitEnd_ref - 1;
				if(verbose)
					std::cerr<<readGapLength << "\t" << refGapLength << "\n";
				auto readSeq = extractReadSeq(read, lastHitEnd_read + 1, lastHitEnd_read + 1 + readGapLength, 1);
				if (verbose) {
					for (size_t x=0; x<50; ++x)
						std::cerr<<"+";
					std::cerr<<" start gap alignment\n";
					std::cerr<<"read from " << lastHitEnd_read+1 << " for length " << readGapLength << " of the read:\t " << readSeq << "\n";
					std::cerr << "currHitStart_read: " << currHitStart_read << " lastHitEnd_read: " << lastHitEnd_read << " currHitStart_ref: " << currHitStart_ref << " lastHitEnd_ref: " << lastHitEnd_ref<< "\n";
				}

				if (firstMem) {
					lastHitEnd_ref = currHitStart_ref > currHitStart_read ? currHitStart_ref - currHitStart_read -1: -1;
					// Not extending the reference beyond its beginning
					auto refStartSeq = lastHitEnd_ref > refExtLength ?  lastHitEnd_ref + 1 - refExtLength : 0;
					refGapLength = currHitStart_ref - refStartSeq;
					// SOFT CLIPPING for reads mapping with a left hanger at the beginning of a transcript
					if (refGapLength > 0) {
						// We reverse the strings because of ksw force alignment from the beginning of the sequences
						std::reverse(readSeq.begin(), readSeq.end());
						tseq="";
						getRefSeq(allRefSeq, refAccPos, refStartSeq, refGapLength, tseq);
						std::reverse(tseq.begin(), tseq.end());
	
						aligner(readSeq.c_str(), readSeq.length(), tseq.c_str(), tseq.length(), &ez, ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
						score = std::max(ez.mqe, ez.mte);
						if (verbose){
							std::cerr<<"It's a gap at the beginning, so the reverse read from "<<lastHitEnd_read + 1<< " is\t"<<readSeq<<"\n";
							std::cerr<<"The reverse reference from " << refStartSeq << " with the length " << refGapLength << " is\t" << tseq << "\n";
							for (int i = 0; i < ez.n_cigar; ++i) // print CIGAR
								printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
							std::cout<<"\n";
						}
					} else {
						score = 0;
					}
				} else {
					//Considering the case where the gap in either the read or reference is empty
					//Example -> AAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCA
					if (readGapLength <= 0) {
						if (verbose)
							std::cerr<< "read gap lengths was of length 0 for the ref gap length of \t" <<refGapLength << "\n";
						score = (-1)*mopts->gapOpenPenalty + (-1)*mopts->gapExtendPenalty*(refGapLength-1);
					} else if (refGapLength <= 0) {
						if (verbose)
							std::cerr<< "ref gap lengths was of length 0 for the read gap length of \t" <<readGapLength << "\n";
						score = (-1)*mopts->gapOpenPenalty + (-1)*mopts->gapExtendPenalty*(readGapLength-1);
					} else {
						//tseq = "";
						//getRefSeq(allRefSeq, refAccPos, lastHitEnd_ref + 1, refGapLength, tseq);
						//auto score1 = aligner(readSeq.c_str(), readSeq.length(), tseq.c_str(), tseq.length(), &ez, ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::GLOBAL>());
						char* refSeq1 = refSeq + lastHitEnd_ref + 1 - refStart;
						score = aligner(readSeq.c_str(), readSeq.length(), refSeq1, refGapLength, &ez, ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::GLOBAL>());
						if (verbose) {
							std::cerr<<"read from "<<lastHitEnd_read + 1<< " with the gap length of " << readGapLength << " is \t"<< readSeq <<"\n";
							std::cerr<<"reference from " << lastHitEnd_ref + 1 << " with the length " << refGapLength << " is \t" << tseq << "\n";
						}
					}
				}
				alignmentScore += score;
				if (verbose) {
					std::cerr<<mem.isFw<<" at "<<currHitStart_ref<<" on "<<tid<<" for "<<refGapLength<<"\n";
					std::cerr<<"curr score2:\t" << score << "\n";
					for (size_t x=0; x<50; ++x)
						std::cerr<<"+";
					std::cerr<<" end gap alignments\n";
				}
			}
	
			score = (currHitStart_ref <= lastHitEnd_ref) ? mopts->matchScore * (currHitStart_ref + mem.memInfo->memlen - lastHitEnd_ref - 1) : mopts->matchScore * mem.memInfo->memlen;
			if (lastHitEnd_ref - currHitStart_ref == lastHitEnd_read - currHitStart_read or firstMem) {
				alignmentScore += score;
				if ( verbose ) {
					auto readSeq = extractReadSeq(read, currHitStart_read, currHitStart_read + mem.memInfo->memlen, 1);
					tseq = "";
					getRefSeq(allRefSeq, refAccPos, currHitStart_ref, mem.memInfo->memlen, tseq);
					std::cerr<<"read from "<<currHitStart_read<< "\t with the sequence \n"<<readSeq<<"\n";
					std::cerr<<"orientation of\t" <<mem.isFw<<" at "<<currHitStart_ref<<" on reference "<<tid<<" for the length of "<<mem.memInfo->memlen<<"\t with sequence:\n"<<tseq<<"\n";
					std::cerr<<"with the score1:\t" << score << "\n";
				}
			} else {
				int32_t gapLength = std::abs(lastHitEnd_ref - currHitStart_ref - lastHitEnd_read + currHitStart_read);
				int32_t penalty = -mopts->gapOpenPenalty - mopts->gapExtendPenalty*(gapLength-1);
				score += penalty - gapLength*mopts->matchScore;
				if (verbose){
					std::cerr<<"penalty\t" << penalty << "\n";
					auto readSeq = extractReadSeq(read, currHitStart_read, currHitStart_read + mem.memInfo->memlen, 1);
					tseq = "";
					getRefSeq(allRefSeq, refAccPos, currHitStart_ref, mem.memInfo->memlen, tseq);
					std::cerr<<"read from "<<currHitStart_read<< "\t with the sequence \n"<<readSeq<<"\n";
					std::cerr<<"orientation of\t" <<mem.isFw<<" at "<<currHitStart_ref<<" on reference "<<tid<<" for the length of "<<mem.memInfo->memlen<<"\t with sequence:\n"<<tseq<<"\n";
					std::cerr<<"with the score1:\t" << score << "\n";
				}
				alignmentScore += score;
			}

			lastHitEnd_read = currHitStart_read + mem.memInfo->memlen - 1;
			lastHitEnd_ref = mem.tpos + mem.memInfo->memlen - 1;
			firstMem = false;
		}

		// Try and align any remaining sequence on the read
		if (lastHitEnd_read  < read.length() - 1) {
			auto readGapLength = read.length() - 1 - lastHitEnd_read;
			auto refGapLength = lastHitEnd_ref + 1 + readGapLength + refExtLength < refTotalLength ? readGapLength + refExtLength : refTotalLength - lastHitEnd_ref - 1;

			// SOFT CLIPPING for reads mapping with a right hanger at the beginning of a transcript
			if (refGapLength != 0) {
				auto readSeq = extractReadSeq(read, lastHitEnd_read + 1, lastHitEnd_read + 1 + readGapLength, 1);
				//std::string tseq = "";
				//getRefSeq(allRefSeq, refAccPos, lastHitEnd_ref + 1, refGapLength, tseq);
				//aligner(readSeq.c_str(), readSeq.length(), tseq.c_str(), tseq.length(), &ez, ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
				//auto score1 = std::max(ez.mqe, ez.mte);
				char* refSeq1 = refSeq + lastHitEnd_ref + 1 - refStart;
				aligner(readSeq.c_str(), readSeq.length(), refSeq1, refGapLength, &ez, ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
				auto score = std::max(ez.mqe, ez.mte);
				alignmentScore += score;
				if (verbose) {
					std::cerr<<"read from "<<lastHitEnd_read + 1<< "\t\t\t"<<readSeq<<"\n";
					std::cerr<<"at "<<lastHitEnd_ref<<" for "<<refGapLength<<"\t"<<tseq<<"\n";
					std::cerr<<"curr score3:\t"<<score<<"\n";
				}
			}
		}

		}

    	if (multiMapping) { // don't bother to fill up a cache unless this is a multi-mapping read
	    	if (!didHash) 
				MetroHash64::Hash(reinterpret_cast<uint8_t*>(refSeq), keyLen, reinterpret_cast<uint8_t*>(&hashKey), 0);
        	alnCache[hashKey] = alignmentScore;
    	} // was a multi-mapper
	}

	if (verbose)
		std::cerr<<"alignmentScore\t"<<alignmentScore<< "\talignmment\t" << alignment <<"\n";
	return alignmentScore;
}

int32_t PufferfishAligner::calculateAlignments(std::string& read_left, std::string& read_right, util::JointMems& jointHit, bool verbose) {
	auto tid = jointHit.tid;
	double optFrac{mopts->minScoreFraction};

	if (jointHit.isOrphan()) {
		if (verbose)
			std::cerr<<"orphan\n";
		int32_t maxScore = mopts->matchScore * read_left.length();
		std::string read_orphan = jointHit.isLeftAvailable() ? read_left : read_right;
		int32_t score = alignRead(read_orphan, jointHit.orphanClust()->mems, jointHit.orphanClust()->perfectChain, jointHit.orphanClust()->isFw, tid, alnCacheLeft, verbose);
		jointHit.orphanClust()->coverage = score > (optFrac * maxScore) ? score : std::numeric_limits<decltype(score)>::min();
		return jointHit.orphanClust()->coverage;
	}
	else {
		int32_t maxLeftScore = mopts->matchScore * read_left.length();
		int32_t maxRightScore = mopts->matchScore * read_right.length();
		if (verbose)
			std::cerr<<"left\n";
		auto score_left = alignRead(read_left, jointHit.leftClust->mems, jointHit.leftClust->perfectChain, jointHit.leftClust->isFw, tid, alnCacheLeft, verbose);
		if (verbose)
			std::cerr<<"right\n";
		auto score_right = alignRead(read_right, jointHit.rightClust->mems, jointHit.rightClust->perfectChain, jointHit.rightClust->isFw, tid, alnCacheRight, verbose);
		jointHit.leftClust->coverage = score_left;// score_left > (optFrac * maxLeftScore) ? score_left : std::numeric_limits<decltype(score_left)>::min();
		jointHit.rightClust->coverage = score_right;//score_right > (optFrac * maxRightScore) ? score_right : std::numeric_limits<decltype(score_right)>::min();
		auto total_score = (score_left + score_right) > (optFrac * (maxLeftScore + maxRightScore)) ? score_left + score_right : std::numeric_limits<decltype(score_left)>::min();
		if (total_score < 0) {
			//std::cerr<< jointHit.coverage() << "\t" << jointHit.leftClust->coverage << "\t" << jointHit.rightClust->coverage << "\n";
		}
		return total_score;
	}
}

std::string calculateCigar(std::string &read,
                           std::string &ref,
                           ksw2pp::KSW2Aligner &aligner,
                           std::vector<util::MemCluster>::iterator clust,
                           bool verbose = true) {
    if (read.empty()) {
        clust->cigar += (std::to_string(ref.length()) + "D");
        clust->score += ref.length() * GAP_SCORE;
    } else if (ref.empty()) {
        clust->cigar += (std::to_string(read.length()) + "I");
        clust->score += read.length() * GAP_SCORE;
    } else {
        auto score = aligner(read.c_str(),
                             read.length(),
                             ref.c_str(),
                             ref.length(),
                             ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::GLOBAL>());
        clust->score += score;
        if (verbose) {
            std::cerr << "read str " << read << "\nref str " << ref << "\n";
            std::cerr << "cigar before : " << clust->cigar << "\n";
        }
        clust->cigar += cigar2str(&aligner.result());
        if (verbose) std::cerr << "cigar after : " << clust->cigar << "\n";
    }
    if (verbose) std::cerr << clust->cigar << "\n";
    return clust->cigar;
}

template<typename PufferfishIndexT>
util::ContigBlock getContigBlock(std::vector<util::UniMemInfo>::iterator memInfo,
                                 PufferfishIndexT *pfi,
                                 spp::sparse_hash_map<uint32_t, util::ContigBlock> *contigSeqCache) {
    if (contigSeqCache->find(memInfo->cid) != contigSeqCache->end()) {
        return (*contigSeqCache)[memInfo->cid];
    } else
        (*contigSeqCache)[memInfo->cid] = {memInfo->cid,
                                           memInfo->cGlobalPos,
                                           memInfo->clen,
                                           pfi->getSeqStr(memInfo->cGlobalPos, memInfo->clen)};

    return (*contigSeqCache)[memInfo->cid];
}


template<typename PufferfishIndexT>
void createSeqPairs(PufferfishIndexT *pfi,
                    std::vector<util::MemCluster>::iterator clust,
                    fastx_parser::ReadSeq &read,
                    RefSeqConstructor<PufferfishIndexT> &refSeqConstructor,
                    spp::sparse_hash_map<uint32_t, util::ContigBlock> &contigSeqCache,
                    uint32_t tid,
                    ksw2pp::KSW2Aligner &aligner,
                    bool verbose,
                    bool naive = true) {

    //(void)verbose;
    (void) naive;
    //std::string& readName = read.name ;
    std::string &readSeq = read.seq;
    auto clustSize = clust->mems.size();
    clust->score = 0;
    auto readLen = readSeq.length();

    //@debug
    if (verbose) std::cerr << "Clust size " << clustSize << "\n";

    clust->cigar = "";
    auto prevTPos = clust->mems[0].tpos;
    size_t it = 0;
    //FIXME why we have failures in graph!! if valid, discard the whole hit, otherwise, find the bug
    for (it = 0; it < clust->mems.size() - 1; ++it) {
        auto mmTstart = clust->mems[it].tpos + clust->mems[it].memInfo->memlen;
        auto mmTend = clust->mems[it + 1].tpos;

        // NOTE: Assume that consistency condition is met, merging the mems
        if (mmTend <
            mmTstart) { // this means the not matched seq is empty. i.e. we have exact match in transition from txp_i to txp_(i+1)
            continue;
        }

        uint32_t rstart = (clust->mems[it].memInfo->rpos + clust->mems[it].memInfo->memlen);
        uint32_t rend = clust->mems[it + 1].memInfo->rpos;

        if (!clust->isFw) {
            rstart = (clust->mems[it + 1].memInfo->rpos + clust->mems[it + 1].memInfo->memlen);
            rend = clust->mems[it].memInfo->rpos;
        }
        if (mmTend == mmTstart && rstart == rend) {
            continue;
        }

        clust->cigar += (std::to_string(mmTstart - prevTPos) + "M");
        clust->score += ((mmTstart - prevTPos) * MATCH_SCORE);
        prevTPos = clust->mems[it + 1].tpos;

        if (mmTend == mmTstart) { //Insertion in read
            if (rend - rstart > 0) { //validity check TODO if passed should be removed for final version
                //std::string tmp = extractReadSeq(readSeq, rstart, rend, clust->isFw) ;
                std::string readSubstr = extractReadSeq(readSeq, rstart, rend, clust->isFw);
                std::string refSeq = "";
                calculateCigar(readSubstr, refSeq, aligner, clust, verbose);
            } else {
                std::cerr << "ERROR: in pufferfishAligner tstart = tend while rend < rstart\n" << read.name << "\n";
            }
        } else if (rstart == rend) { //deletion in read
            if (mmTend - mmTstart > 0) { //validity check TODO if passed should be removed for final version
                //std::string tmp = extractReadSeq(readSeq, rstart, rend, clust->isFw) ;
                clust->cigar += (std::to_string(mmTend - mmTstart) + "D");
                clust->score += (mmTend - mmTstart) * GAP_SCORE;
                //calculateCigar(extractReadSeq(readSeq, rstart, rend, clust->isFw), "", aligner, clust);
            } else {
                std::cerr << "ERROR: in pufferfishAligner rstart = rend while tend < tstart\n" << read.name << "\n";
            }
        } else { // complex case
            // Need to fetch reference string --> Doing graph traversal
            std::string refSeq = "";
            bool firstContigDirWRTref = clust->mems[it].memInfo->cIsFw == clust->isFw; // true ~ (ref == contig)
            bool secondContigDirWRTref = clust->mems[it + 1].memInfo->cIsFw == clust->isFw;
            auto distOnTxp = clust->mems[it + 1].tpos - (clust->mems[it].tpos + clust->mems[it].memInfo->memlen);
            uint32_t cstart = firstContigDirWRTref ? (clust->mems[it].memInfo->cpos + clust->mems[it].memInfo->memlen -
                                                      1) : clust->mems[it].memInfo->cpos;
            uint32_t cend = secondContigDirWRTref ? clust->mems[it + 1].memInfo->cpos : (
                    clust->mems[it + 1].memInfo->cpos + clust->mems[it + 1].memInfo->memlen - 1);

            // heuristic : length of mismatched parts of txp and read are both one??? It's a definite mismatch!!
            if (distOnTxp == 1 && (rend - rstart) == 1) {
                clust->cigar += "1M";
                clust->score += MISMATCH_SCORE;
            } else {

                util::ContigBlock scb = getContigBlock(clust->mems[it].memInfo, pfi, &contigSeqCache);
                util::ContigBlock ecb = getContigBlock(clust->mems[it + 1].memInfo, pfi, &contigSeqCache);

                //NOTE In simplest form, assuming both start and end contigs are forward wrt the reference, then
                // cstart and cend point to the last matched base in the start contig and first matched base in last contig
                if (verbose) {
                    std::cerr << " start contig: \ncid" << scb.contigIdx_ << " relpos" << clust->mems[it].memInfo->cpos
                              << " len" << scb.contigLen_ << " ori" << firstContigDirWRTref
                              << " hitlen" << clust->mems[it].memInfo->memlen << " start pos:" << cstart << "\n"
                              << scb.seq << "\n";
                    std::cerr << " last contig: \ncid" << ecb.contigIdx_ << " relpos"
                              << clust->mems[it + 1].memInfo->cpos << " len" << ecb.contigLen_ << " ori"
                              << secondContigDirWRTref
                              << " hitlen" << clust->mems[it + 1].memInfo->memlen << " end pos:" << cend << "\n"
                              << ecb.seq << "\n";
                }
                Task res = refSeqConstructor.fillSeq(tid,
                                                     clust->mems[it].tpos + clust->mems[it].memInfo->memlen - 1,
                                                     firstContigDirWRTref,
                                                     scb, cstart, ecb, cend,
                                                     secondContigDirWRTref,
                                                     distOnTxp,
                                                     refSeq);
                if (res == Task::SUCCESS) {
                    //std::cerr << "SUCCESS\n";
                    //std::cerr << " part of read "<<extractReadSeq(readSeq, rstart, rend, clust->isFw)<<"\n"
                    //         << " part of ref  " << refSeq << "\n";
                    std::string readSubstr = extractReadSeq(readSeq, rstart, rend, clust->isFw);
                    calculateCigar(readSubstr, refSeq, aligner, clust, verbose);
                } else {
                    //FIXME discard whole hit!!!
                    clust->score = std::numeric_limits<int>::min();
                    //std::cerr << "Graph searched FAILED \n" ;
                }
            }
        }
    }
    // update cigar and add matches for last unimem(s)
    // after the loop, "it" is pointing to the last unimem in the change
    clust->cigar += (std::to_string(clust->mems[it].tpos + clust->mems[it].memInfo->memlen - prevTPos) + "M");
    clust->score += ((clust->mems[it].tpos + clust->mems[it].memInfo->memlen - prevTPos) * MATCH_SCORE);


    // Take care of left and right gaps/mismatches
    util::ContigBlock dummy = {std::numeric_limits<uint64_t>::max(), 0, 0, "", true};

    bool lastContigDirWRTref = clust->mems[it].memInfo->cIsFw == clust->isFw;
    bool firstContigDirWRTref = clust->mems[0].memInfo->cIsFw == clust->isFw;


    auto endRem = readLen - (clust->mems[it].memInfo->rpos + clust->mems[it].memInfo->memlen);
    auto startRem = clust->mems[0].memInfo->rpos;

    std::string startReadSeq = "";
    std::string endReadSeq = "";
    if (clust->isFw) {
        startReadSeq = extractReadSeq(readSeq, 0, clust->mems[0].memInfo->rpos, clust->isFw);
        endReadSeq = extractReadSeq(readSeq, (clust->mems[it].memInfo->rpos + clust->mems[it].memInfo->memlen), readLen,
                                    clust->isFw);
    }

    if (!clust->isFw) {
        endRem = clust->mems[it].memInfo->rpos;
        startRem = readLen - (clust->mems[0].memInfo->rpos + clust->mems[0].memInfo->memlen);
        startReadSeq = extractReadSeq(readSeq, (clust->mems[0].memInfo->rpos + clust->mems[0].memInfo->memlen), readLen,
                                      clust->isFw);
        endReadSeq = extractReadSeq(readSeq, 0, clust->mems[it].memInfo->rpos, clust->isFw);
    }

    /*if (verbose) std::cerr << firstContigDirWRTref << " " << startRem << " " << startReadSeq
                           << " last " << lastContigDirWRTref << " " << endRem << " " << endReadSeq
                         << "\n";
    */
    if (endRem > 0) {
        //std::cerr << "RIIIGHT hangover\n";
        std::string refSeq = "";
        util::ContigBlock scb = getContigBlock(clust->mems[it].memInfo, pfi, &contigSeqCache);
        uint32_t cstart = lastContigDirWRTref ? (clust->mems[it].memInfo->cpos + clust->mems[it].memInfo->memlen - 1)
                                              : clust->mems[it].memInfo->cpos;

        Task res = refSeqConstructor.fillSeq(tid,
                                             clust->mems[it].tpos + clust->mems[it].memInfo->memlen - 1,
                                             lastContigDirWRTref,
                                             scb, cstart, dummy, 0,
                                             lastContigDirWRTref,
                                             endRem,
                                             refSeq,
                                             verbose);
        if (res == Task::SUCCESS) {
            if (verbose)
                std::cerr << "\n\n\nSUCCESS   end extension::: " << read.name << "\n" << "read: " << endReadSeq
                          << " ref: " << refSeq << "\n"
                          << "tid " << tid << " starting from "
                          << clust->mems[it].tpos + clust->mems[it].memInfo->memlen << " endRem " << endRem << " cid "
                          << scb.contigIdx_ << "\n";

            //std::cerr << " part of read "<<endReadSeq<<"\n"
            //          << " part of ref  " << refSeq << "\n";
            calculateCigar(endReadSeq, refSeq, aligner, clust, verbose);
        } else {
            if (verbose)
                std::cerr << "\n\n\nFAILURE   end extension::: " << read.name << "\n" << "read: " << endReadSeq
                          << " ref: " << refSeq << "\n"
                          << "tid " << tid << " starting from "
                          << clust->mems[it].tpos + clust->mems[it].memInfo->memlen << " endRem " << endRem << " cid "
                          << scb.contigIdx_ << "\n";
            // discard whole hit!!!
            clust->score = std::numeric_limits<int>::min();
            //std::cerr << "Graph searched FAILED \n" ;
        }
    }
    if (startRem > 0) {
        //std::cerr << read.name << "\n";
        //std::cerr << "LEEEFT hangover\n";
        std::string refSeq = "";
        util::ContigBlock ecb = getContigBlock(clust->mems[0].memInfo, pfi, &contigSeqCache);
        if (verbose) {
            std::cerr << firstContigDirWRTref << " " << clust->mems[0].memInfo->cpos
                      << " " << (clust->mems[0].memInfo->cpos + clust->mems[0].memInfo->memlen - 1) << "\n";
            std::cerr << ecb.seq.length() << "   " << ecb.seq << "\n";
        }
        uint32_t cend = firstContigDirWRTref ? clust->mems[0].memInfo->cpos : (clust->mems[0].memInfo->cpos +
                                                                               clust->mems[0].memInfo->memlen - 1);

        Task res = refSeqConstructor.fillSeq(tid,
                                             clust->mems[0].tpos + clust->mems[0].memInfo->memlen - 1,
                                             firstContigDirWRTref,
                                             dummy, 0, ecb, cend,
                                             firstContigDirWRTref,
                                             startRem,
                                             refSeq);
        //if (verbose) std::cerr << "got here\n";
        if (res == Task::SUCCESS) {
            //std::cerr << " part of read "<<startReadSeq<<"\n"
            //          << " part of ref  " << refSeq << "\n";
            calculateCigar(endReadSeq, refSeq, aligner, clust, verbose);
        } else {
            std::cerr << "beginning extension::: " << read.name << "\n" << "read: " << endReadSeq << " ref: " << refSeq
                      << "\n";
            // discard whole hit!!!
            clust->score = std::numeric_limits<int>::min();
            //std::cerr << "Graph searched FAILED \n" ;
        }
    }
    clust->isVisited = true;
    if (verbose)
        std::cerr << read.name << " " << clust->isFw << " " << clust->mems.size() << "\n" << read.seq << "\nSCORE: "
                  << clust->score << "\nCIGAR: " << clust->cigar << "\n";
}

template<typename ReadPairT, typename PufferfishIndexT>
void traverseGraph(ReadPairT &rpair,
                   util::JointMems &hit,
                   PufferfishIndexT &pfi,
                   RefSeqConstructor<PufferfishIndexT> &refSeqConstructor,
                   spp::sparse_hash_map<uint32_t, util::ContigBlock> &contigSeqCache,
                   ksw2pp::KSW2Aligner &aligner,
                   bool verbose,
                   bool naive = false) {

    size_t tid = hit.tid;
    auto readLen = rpair.first.seq.length();
    if (verbose) std::cerr << rpair.first.name << "\n";

    if (hit.isOrphan()) {
        if (!hit.orphanClust()->isVisited && hit.orphanClust()->coverage < readLen)
            createSeqPairs(&pfi, hit.orphanClust(), rpair.first, refSeqConstructor, contigSeqCache, tid, aligner,
                           verbose, naive);
        else if (hit.orphanClust()->coverage == readLen) {
            hit.orphanClust()->score += hit.orphanClust()->coverage * MATCH_SCORE;
            hit.orphanClust()->cigar = std::to_string(hit.orphanClust()->coverage) + "M";
        }
    } else {
        if (!hit.leftClust->isVisited && hit.leftClust->coverage < readLen)
            createSeqPairs(&pfi, hit.leftClust, rpair.first, refSeqConstructor, contigSeqCache, tid, aligner, verbose,
                           naive);
        else if (hit.leftClust->coverage == readLen) {
            hit.leftClust->score += hit.leftClust->coverage * MATCH_SCORE;
            hit.leftClust->cigar = std::to_string(hit.leftClust->coverage) + "M";
        }
        //goOverClust(pfi, hit.leftClust, rpair.first, contigSeqCache, tid, verbose) ;
        if (!hit.rightClust->isVisited && hit.rightClust->coverage < readLen)
            createSeqPairs(&pfi, hit.rightClust, rpair.second, refSeqConstructor, contigSeqCache, tid, aligner, verbose,
                           naive);
        else if (hit.rightClust->coverage == readLen) {
            hit.rightClust->score += hit.rightClust->coverage * MATCH_SCORE;
            hit.rightClust->cigar = std::to_string(hit.rightClust->coverage) + "M";
        }
    }
    //goOverClust(pfi, hit.rightClust, rpair.second, contigSeqCache, tid, verbose) ;
}


bool foundTruth(std::string &refName, std::string &readName) {
    int cntr = 0;
    //gi|456351576|dbj|AP012603.1|
    size_t idx = 3;
    char refchr = refName[idx];
    char readchr = readName[idx];
    while (cntr != 3 || idx == refName.size() || idx == readName.size()) {
        if (refchr != readchr)
            return false;
        idx++;
        refchr = refName[idx];
        readchr = readName[idx];
        if (refchr == '|')
            cntr++;
    }
    if (readchr == '|')
        return true;
    return false;
}

/**
 * Main function for processing paired-end reads.
 **/
template<typename PufferfishIndexT>
void processReadsPair(paired_parser *parser,
                      PufferfishIndexT &pfi,
                      MutexT *iomutex,
                      std::shared_ptr<spdlog::logger> outQueue,
                      HitCounters &hctr,
                      AlignmentOpts *mopts) {
    MemCollector<PufferfishIndexT> memCollector(&pfi);

    //create aligner
    spp::sparse_hash_map<uint32_t, util::ContigBlock> contigSeqCache;
    RefSeqConstructor<PufferfishIndexT> refSeqConstructor(&pfi, &contigSeqCache);

    std::vector<std::string> refBlocks;

    auto logger = spdlog::get("stderrLog");
    fmt::MemoryWriter sstream;
    BinWriter bstream;
    //size_t batchSize{2500} ;
    size_t readLen{0};
    //size_t totLen{0};

    spp::sparse_hash_map<size_t, std::vector<util::MemCluster>> leftHits;
    spp::sparse_hash_map<size_t, std::vector<util::MemCluster>> rightHits;
    std::vector<util::JointMems> jointHits;
    PairedAlignmentFormatter<PufferfishIndexT *> formatter(&pfi);
    util::QueryCache qc;
    std::vector<util::MemCluster> all;

    //@fatemeh Initialize aligner ksw
    ksw2pp::KSW2Config config;
    ksw2pp::KSW2Aligner aligner(MATCH_SCORE, MISMATCH_SCORE);

    config.gapo = -1 * GAP_SCORE;
    config.gape = -1 * GAP_SCORE;
    config.bandwidth = -1;
    config.flag = KSW_EZ_RIGHT;
    aligner.config() = config;

    auto rg = parser->getReadGroup();

	//For filtering reads
	auto& txpNames = pfi.getRefNames();
    while (parser->refill(rg)) {
        for (auto &rpair : rg) {
            readLen = rpair.first.seq.length();
            //totLen = readLen + rpair.second.seq.length();
            bool verbose = mopts->verbose;
			//bool verbose = (rpair.first.seq == "AAAATGTGAAGGCCAAGATCCAAGATAAAGAAGGCATCCCCCCCGACCAGCAGAGGCTCATCTTTGCAGGCAAGCAGCTGGAAGATGGCCGCACTCTTTC" or
			//				rpair.second.seq == "AAAATGTGAAGGCCAAGATCCAAGATAAAGAAGGCATCCCCCCCGACCAGCAGAGGCTCATCTTTGCAGGCAAGCAGCTGGAAGATGGCCGCACTCTTTC");
			if (verbose) std::cerr << rpair.first.name << "\n";
            //std::cerr << "read: " << rpair.first.name << "\n\n";

            ++hctr.numReads;

            jointHits.clear();
            leftHits.clear();
            rightHits.clear();
            memCollector.clear();

            //help me to debug, will deprecate later
            //std::cerr << "\n first seq in pair " << rpair.first.seq << "\n" ;
            //std::cerr << "\n second seq in pair " << rpair.second.seq << "\n" ;

            //std::cerr << "\n going inside hit collector \n" ;
            //readLen = rpair.first.seq.length() ;
            bool lh = memCollector(rpair.first.seq,
                                   leftHits,
                                   mopts->maxSpliceGap,
                                   MateStatus::PAIRED_END_LEFT,
                                   qc,
                                   verbose

                    /*
                    mopts->consistentHits,
                    refBlocks*/);
            bool rh = memCollector(rpair.second.seq,
                                   rightHits,
                                   mopts->maxSpliceGap,
                                   MateStatus::PAIRED_END_RIGHT,
                                   qc,
                                   verbose
                    /*,
                    mopts->consistentHits,
                    refBlocks*/);

            all.clear();

			// Second one, described in slack for Fatemeh verbose = "read10059738/ENST00000247866.8:ENSG00000090266.12;mate1:102-201;mate2:257-356" == rpair.first.name;
			//First one 302 verbose = rpair.first.name == "read9998077/ENST00000223215.8:ENSG00000106484.14;mate1:757-856;mate2:892-991" or rpair.second.name=="read9998077/ENST00000223215.8:ENSG00000106484.14;mate1:757-856;mate2:892-991";
            memCollector.findOptChainAllowingOneJumpBetweenTheReadEnds/*findBestChain*/(jointHits,
                                                                                        all,
                                                                                        mopts->maxSpliceGap,
                                                                                        mopts->maxFragmentLength,
																						rpair.first.seq.length(),
                                                                                        verbose);

            hctr.numMappedAtLeastAKmer += jointHits.size() ? 1 : 0; //(leftHits.size() || rightHits.size()) ? 1 : 0;
            //do intersection on the basis of
            //performance, or going towards selective alignment
            //otherwise orphan

            // We also handle orphans inside this function
            //if(lh && rh){
            //TODO uncomment this part if you want to get back to the original joinReadsAndFilter
            /* joinReadsAndFilter(leftHits, rightHits, jointHits,
                                mopts->maxFragmentLength,
                                totLen,
                                mopts->scoreRatio,
                                mopts->noDiscordant,
                                mopts->noOrphan,
                                pfi,
                     *//* rpair.first.name == "NC_009641.fna:1:1:2066:284#0/1" *//*
                               verbose);*/
            /* } else{
              //ignore orphans for now
            } */


            //jointHits is a vector
            //this can be used for BFS
            //NOTE sanity check
            //void traverseGraph(std::string& leftReadSeq, std::string& rightReadSeq, util::JointMems& hit, PufferfishIndexT& pfi,   std::map<uint32_t, std::string>& contigSeqCache){
            int maxScore = std::numeric_limits<int>::min();
            
            std::vector<QuasiAlignment> jointAlignments;

			if (mopts->validateMappings) {		
				ksw2pp::KSW2Aligner aligner(mopts->matchScore, mopts->missMatchScore);
				PufferfishAligner puffaligner(pfi.refseq_, pfi.refAccumLengths_, pfi.k(), mopts, aligner, jointHits.size()>1);

				int32_t bestScore = std::numeric_limits<int32_t>::min() ;
				std::vector<decltype(bestScore)> scores(jointHits.size(), bestScore);
				size_t idx{0};

				std::map<int32_t, std::vector<int32_t>> transcript_set;
    	        for (auto &jointHit : jointHits) {
					// chaining bug bool verbose = rpair.first.seq == "CCAGCAGAGAGTAGTGACACAGGAGTTCTGGAGGGCTGTGCCGGGCTGCAGCTTGGAGGGCAGGGCGGGGCTGCAGCTTGGAGGGCAGGGCGGGGCTGCA" and jointHit.tid==17815;
					//bool verbose = rpair.first.seq == "CCATTTTATTTTATTTTATTTTATTTTATTTTNTTTTATTTTATTTTTGAGAAAGGGTCTCACTCTGTCACCCAGGCTGAAGTGCAGTGGTGCCATCATA" and jointHit.tid == 25088;
					//bool verbose = (rpair.first.seq == "ACTGGGAGGCAGGAGGAGCTGGGCCTGGAGAGGCTGACTCGAGGAAGTTTTGCACCTGGAGAGGCCGTCGAGAGGACGGAGCTGGGCCCAGGGAGGCCGA" or
					//			   rpair.second.seq == "ACTGGGAGGCAGGAGGAGCTGGGCCTGGAGAGGCTGACTCGAGGAAGTTTTGCACCTGGAGAGGCCGTCGAGAGGACGGAGCTGGGCCCAGGGAGGCCGA") and
					//				jointHit.tid == 44366;
					//bool verbose = rpair.first.name == "read9998077/ENST00000223215.8:ENSG00000106484.14;mate1:757-856;mate2:892-991" or rpair.second.name=="read9998077/ENST00000223215.8:ENSG00000106484.14;mate1:757-856;mate2:892-991";
					auto hitScore = puffaligner.calculateAlignments(rpair.first.seq, rpair.second.seq, jointHit, verbose);
					if (hitScore < 0)
						hitScore = std::numeric_limits<int32_t>::min();
					scores[idx] = hitScore;
					++idx;
					bestScore = (hitScore > bestScore) ? hitScore : bestScore;
				
					if (transcript_set.find(jointHit.tid) == transcript_set.end()) {
						std::vector<int32_t> coverage_hit;
						transcript_set[jointHit.tid] = coverage_hit;
						transcript_set[jointHit.tid].push_back(jointHit.coverage());
						transcript_set[jointHit.tid].push_back(idx-1);
					} else if (jointHit.coverage() > transcript_set[jointHit.tid][0] ) {
						transcript_set[jointHit.tid][0] = jointHit.coverage();
						scores[transcript_set[jointHit.tid][1]] = std::numeric_limits<int32_t>::min();
						transcript_set[jointHit.tid][1] = idx-1;
					} else {
						scores[idx-1] = std::numeric_limits<int32_t>::min();
					}
				}

				if (mopts->strictFilter) {
					int myctr = 0;
					for(auto& jointHit : jointHits){
						if(jointHit.orphanClust()->coverage < bestScore)
							scores[myctr] = std::numeric_limits<int32_t>::min();
						myctr++;
					}
				}
			
				// Filter out these alignments with low scores
				uint32_t ctr{0};
				if (bestScore > std::numeric_limits<int32_t>::min()) {
					jointHits.erase(
								std::remove_if(jointHits.begin(), jointHits.end(),
									[&ctr, &scores, bestScore] (util::JointMems& ja) -> bool{
										bool rem = (scores[ctr] == std::numeric_limits<int32_t>::min());
										++ctr;
										return rem;
									}),
								jointHits.end()
								);
            	} else {
					// There is no alignment with high quality for this read, so we skip this reads' alignments
					jointHits.clear();
            	}

			}
	        //if (mopts->salmonUnmapped.find(rpair.first.name) == mopts->salmonUnmapped.end()) {
            //	jointHits.clear();
	        //} else {
         	//} 

			hctr.totHits += jointHits.size();
            hctr.peHits += jointHits.size();


            hctr.numMapped += !jointHits.empty() ? 1 : 0;

            if (mopts->noOrphan)
                hctr.numOfOrphans += jointHits.empty() && (lh || rh);
            else
                hctr.numOfOrphans += !jointHits.empty() && (jointHits.back().isOrphan()) ? 1 : 0;

            if (jointHits.size() > hctr.maxMultimapping) {
                hctr.maxMultimapping = jointHits.size();
            }
			bool found = false;
			for (auto &jointHit : jointHits) {
				//if (rpair.first.name.find(txpNames[jointHit.tid]) != std::string::npos){
				//:	hctr.validHits++;
				//	std::cout<<"#" << rpair.first.name << "\t" << txpNames[jointHit.tid] << "\t" << jointHit.coverage() <<"\n";
				//}
					
                // FIXME : This part needs to be taken care of
                // If graph returned failure for one of the ends --> should be investigated more.
                if (!mopts->justMap and
                    ((!jointHit.isOrphan() and (jointHit.leftClust->score == std::numeric_limits<int>::min()
                                                or jointHit.rightClust->score == std::numeric_limits<int>::min()))
                     or
                     (jointHit.isLeftAvailable() and jointHit.leftClust->score == std::numeric_limits<int>::min())
                     or
                     (jointHit.isRightAvailable() and jointHit.rightClust->score == std::numeric_limits<int>::min()))
                        ) {
                    //std::cerr << "Failed: " << rpair.first.name << "\n";
                    continue;
                }
                if (mopts->justMap or (static_cast<int>(jointHit.coverage()) == maxScore)) {
                    if (jointHit.isOrphan()) {
                        if (verbose) {
                            std::cerr << "orphan\t";
                            std::cerr << jointHit.tid << "\t";
                            std::cerr << jointHit.orphanClust()->getTrFirstHitPos() << "\t";     // reference pos
                            std::cerr << (uint32_t) jointHit.orphanClust()->isFw << "\t";     // fwd direction
                            std::cerr << readLen << "\t"; // read length
                            std::cerr << jointHit.orphanClust()->cigar << "\t"; // cigar string
                            std::cerr << jointHit.fragmentLen << "\n";
                        }
                        jointAlignments.emplace_back(jointHit.tid,           // reference id
                                                     jointHit.orphanClust()->getTrFirstHitPos(),     // reference pos
                                                     jointHit.orphanClust()->isFw,     // fwd direction
                                                     readLen, // read length
                                                     jointHit.orphanClust()->cigar, // cigar string
                                                     jointHit.fragmentLen,       // fragment length
                                                     false);
                        auto &qaln = jointAlignments.back();
						qaln.numHits =jointHit.orphanClust()->coverage;
                        qaln.mateStatus = jointHit.mateStatus;
                    } else {
                        if (verbose) {
                            std::cerr << "paired\t";
                            std::cerr << jointHit.tid << "\t";
                            std::cerr << jointHit.leftClust->getTrFirstHitPos() << "\t";     // reference pos
                            std::cerr << (uint32_t) jointHit.leftClust->isFw << "\t";     // fwd direction
                            std::cerr << readLen << "\t"; // read length
                            std::cerr << jointHit.leftClust->cigar << "\t"; // cigar string
                            std::cerr << jointHit.fragmentLen << "\n";
                        }

                        jointAlignments.emplace_back(jointHit.tid,           // reference id
                                                     jointHit.leftClust->getTrFirstHitPos(),     // reference pos
                                                     jointHit.leftClust->isFw,     // fwd direction
                                                     readLen, // read length
                                                     jointHit.leftClust->cigar, // cigar string
                                                     jointHit.fragmentLen,       // fragment length
                                                     true);         // properly paired
                        // Fill in the mate info
                        auto &qaln = jointAlignments.back();
						qaln.mateLen = readLen;
						qaln.mateCigar = jointHit.rightClust->cigar;
						qaln.matePos = jointHit.rightClust->getTrFirstHitPos();
						qaln.mateIsFwd = jointHit.rightClust->isFw;
						qaln.mateStatus = MateStatus::PAIRED_END_PAIRED;
						qaln.numHits = jointHit.coverage();
					}
                }
            }

            hctr.totAlignment += jointAlignments.size();
            if (verbose) {
                std::cerr << " passed filling jointAlignments:" << hctr.totAlignment << "\n";
            }
            if (!mopts->noOutput) {
                if (mopts->krakOut) {
                    writeAlignmentsToKrakenDump(rpair, /* formatter,  */jointHits, bstream);
                } else if (mopts->salmonOut) {
                    writeAlignmentsToKrakenDump(rpair, /* formatter,  */jointHits, bstream, false);
                } else if (jointAlignments.size() > 0) {
                    writeAlignmentsToStream(rpair, formatter, jointAlignments, sstream, !mopts->noOrphan,
                                            mopts->justMap);
                } else if (jointAlignments.size() == 0) {
                    writeUnmappedAlignmentsToStream(rpair, formatter, jointAlignments, sstream, !mopts->noOrphan,
                                                    mopts->justMap);
                }
            }

            // write them on cmd
            if (hctr.numReads > hctr.lastPrint + 1000000) {
                hctr.lastPrint.store(hctr.numReads.load());
                if (!mopts->quiet and iomutex->try_lock()) {
                    if (hctr.numReads > 0) {
                        std::cerr << "\r\r";
                    }
                    std::cerr << "saw " << hctr.numReads << " reads : "
                              << "pe / read = " << hctr.peHits / static_cast<float>(hctr.numReads)
                              << " : se / read = " << hctr.seHits / static_cast<float>(hctr.numReads) << ' ';
                    iomutex->unlock();
                }
            }
        } // for all reads in this job



        // dump output
        if (!mopts->noOutput) {

            // Get rid of last newline
            if (mopts->salmonOut) {
                //std::cerr << "bytes: " << bstream.getBytes() << "\n";
                BinWriter sbw(sizeof(uint64_t));
                sbw << bstream.getBytes();
                outQueue->info("{}{}", sbw, bstream);
            } else if (mopts->krakOut) {
                outQueue->info("{}", bstream);
            } else {
                std::string outStr(sstream.str());
                if (!outStr.empty()) {
                    outStr.pop_back();
                    outQueue->info(std::move(outStr));
                }
            }
            sstream.clear();
            bstream.clear();
        }

    } // processed all reads
}

//===========
// SINGLE END
//============
template<typename PufferfishIndexT>
void processReadsSingle(single_parser *parser,
                        PufferfishIndexT &pfi,
                        MutexT *iomutex,
                        std::shared_ptr<spdlog::logger> outQueue,
                        HitCounters &hctr,
                        AlignmentOpts *mopts) {
    MemCollector<PufferfishIndexT> memCollector(&pfi);

    //create aligner
    spp::sparse_hash_map<uint32_t, util::ContigBlock> contigSeqCache;
    RefSeqConstructor<PufferfishIndexT> refSeqConstructor(&pfi, &contigSeqCache);
    std::vector<std::string> refBlocks;

    auto logger = spdlog::get("stderrLog");
    fmt::MemoryWriter sstream;
    BinWriter bstream;
    //size_t batchSize{2500} ;
    size_t readLen{0};
    //size_t totLen{0};

    spp::sparse_hash_map<size_t, std::vector<util::MemCluster>> leftHits;
    spp::sparse_hash_map<size_t, std::vector<util::MemCluster>> rightHits;
    std::vector<util::JointMems> jointHits;
    PairedAlignmentFormatter<PufferfishIndexT *> formatter(&pfi);
    util::QueryCache qc;
    std::vector<util::MemCluster> all;

    //@fatemeh Initialize aligner ksw
    ksw2pp::KSW2Config config;
    ksw2pp::KSW2Aligner aligner(MATCH_SCORE, MISMATCH_SCORE);

    config.gapo = -1 * GAP_SCORE;
    config.gape = -1 * GAP_SCORE;
    config.bandwidth = -1;
    config.flag = KSW_EZ_RIGHT;
    aligner.config() = config;

	auto& txpNames = pfi.getRefNames();

    auto rg = parser->getReadGroup();
    while (parser->refill(rg)) {
        for (auto &read : rg) {
            readLen = read.seq.length();
            //totLen = readLen;
            bool verbose = false;
            if (verbose) std::cerr << read.name << "\n";
            ++hctr.numReads;

            jointHits.clear();
            leftHits.clear();
            memCollector.clear();


            bool lh = memCollector(read.seq,
                                   leftHits,
                                   mopts->maxSpliceGap,
                                   MateStatus::SINGLE_END,
                                   qc,
                                   false //verbose
                    /*
                    mopts->consistentHits,
                    refBlocks*/);
            (void) lh;
            all.clear();
            memCollector.findOptChainAllowingOneJumpBetweenTheReadEnds/*findBestChain*/(jointHits,
                                                                                        all,
                                                                                        mopts->maxSpliceGap,
                                                                                        mopts->maxFragmentLength,
																						read.seq.length(),
                                                                                        verbose);

            // Filter left hits
/*            uint32_t maxCoverage{0};
            uint32_t perfectCoverage{static_cast<uint32_t>(totLen)};

            std::vector<std::pair<uint32_t, decltype(leftHits)::mapped_type::iterator>> validHits;
            validHits.reserve(2 * leftHits.size());
            for (auto &l : leftHits) {
                auto &lclust = l.second;
                for (auto clustIt = lclust.begin(); clustIt != lclust.end(); ++clustIt) {
                    if (clustIt->coverage > maxCoverage) { maxCoverage = clustIt->coverage; }
                    if (clustIt->coverage >= mopts->scoreRatio * maxCoverage or clustIt->coverage == perfectCoverage) {
                        validHits.emplace_back(static_cast<uint32_t>(l.first), clustIt);
                        */
/*
                        if (read.name == "spj_1622951_1623126_68819/2") {
                          std::cerr << "\ntid" << l.first << " " << pfi.refName(l.first) << " coverage:" << clustIt->coverage << "\n";
                           for (auto& clus : lclust) {
                             std::cerr << "mem size: " << clus.mems.size() << "\n";
                             for (auto& mem : clus.mems) {
                             std::cerr << "t" << mem.tpos << " r" << mem.memInfo->rpos << " cid" << mem.memInfo->cid << " -- ";
                             }
                             std::cerr << "\n";
                             }
                        }
                        *//*

                    }
                }
                double thresh = mopts->scoreRatio * maxCoverage;
                // sort the hits by coverage
                std::sort(
                        validHits.begin(), validHits.end(),
                        [thresh](
                                std::pair<uint32_t, decltype(leftHits)::mapped_type::iterator> &
                                e1, std::pair<uint32_t, decltype(leftHits)::mapped_type::iterator> &e2) -> bool {
                            return e1.second->coverage > e2.second->coverage;
                        });
                // remove those that don't achieve the threshold
                validHits.erase(std::remove_if(validHits.begin(), validHits.end(),
                                               [thresh](
                                                       std::pair<uint32_t, decltype(leftHits)::mapped_type::iterator> &e) -> bool {
                                                   return static_cast<double>(e.second->coverage) < thresh;
                                               }), validHits.end());
            }*/

            std::vector<QuasiAlignment> jointAlignments;
            std::vector<std::pair<uint32_t, std::vector<util::MemCluster>::iterator>> validHits;

			// Removing reads which contain N
			/*if (read.seq.find('N') != std::string::npos){
				jointHits.clear();
				continue;
			}*/

			if (mopts->validateMappings) {
				ksw2pp::KSW2Aligner aligner(mopts->matchScore, mopts->missMatchScore);
				PufferfishAligner puffaligner(pfi.refseq_, pfi.refAccumLengths_, pfi.k(), mopts, aligner, jointHits.size()>1);

				int32_t bestScore = std::numeric_limits<int32_t>::min();
				std::vector<decltype(bestScore)> scores(jointHits.size(), bestScore);
				size_t idx{0};
		
				std::map<int32_t, std::vector<int32_t>> transcript_set;
    	        for (auto &jointHit : jointHits) {
					//bool verbose = read.name == "read100154/ENST00000358779.9:ENSG00000008128.22;mate1:1983-2082;mate2:2149-2248";
					//bool verbose = read.name == "read16599250/ENST00000540351.1:ENSG00000150991.14;mate1:210-309;mate2:371-470" and txpNames[jointHit.tid] == "ENST00000541645.1:ENSG00000150991.14";

					int32_t hitScore = puffaligner.calculateAlignments(read.seq, read.seq, jointHit, verbose);
					scores[idx] = hitScore;
					++idx;
					bestScore = (hitScore > bestScore) ? hitScore : bestScore;

					// removing dupplicate hits from a read to the same transcript
					if (transcript_set.find(jointHit.tid) == transcript_set.end()) {
						std::vector<int32_t> coverage_hit;
						transcript_set[jointHit.tid] = coverage_hit;
						transcript_set[jointHit.tid].push_back(jointHit.orphanClust()->coverage);
						transcript_set[jointHit.tid].push_back(idx-1);
					} else if (jointHit.orphanClust()->coverage > transcript_set[jointHit.tid][0] ) {
						transcript_set[jointHit.tid][0] = jointHit.orphanClust()->coverage;
						scores[transcript_set[jointHit.tid][1]] = std::numeric_limits<int32_t>::min();
						transcript_set[jointHit.tid][1] = idx-1;
					} else {
						scores[idx-1] = std::numeric_limits<int32_t>::min();
					}
				}

				if (mopts->strictFilter) {
					int myctr = 0;
					for(auto& jointHit : jointHits){
						if(jointHit.orphanClust()->coverage < bestScore)// and validScore != std::numeric_limits<int32_t>::min())// and read.name.find(txpNames[jointHit.tid]) == std::string::npos)
							scores[myctr] = std::numeric_limits<int32_t>::min();
						myctr++;
					}
				}

				// Filter out mappings with low alignment scores
				uint32_t ctr{0};
				if (bestScore > std::numeric_limits<int32_t>::min()) {	
					jointHits.erase(
								std::remove_if(jointHits.begin(), jointHits.end(),
									[&ctr, &scores, bestScore] (util::JointMems& ja) -> bool{
										bool rem = (scores[ctr] == std::numeric_limits<int32_t>::min());
										++ctr;
										return rem;
									}),
								jointHits.end()
								);
	            } else {
					// There is no alignment with high quality for this read, so we skip this reads' alignments
					jointHits.clear();
        	    }

			}

            hctr.numMappedAtLeastAKmer += jointHits.size() > 0 ? 1:0;
            hctr.totHits += jointHits.size();
            hctr.seHits += jointHits.size();
            hctr.numMapped += !jointHits.empty() ? 1 : 0;
            if (jointHits.size() > hctr.maxMultimapping) {
                hctr.maxMultimapping = jointHits.size();
            }

            for (auto &jointHit : jointHits) {
				//if (read.name.find(txpNames[jointHit.tid]) != std::string::npos) {
				//} else {
				//}

                jointAlignments.emplace_back(jointHit.tid,           // reference id
                                                     jointHit.orphanClust()->getTrFirstHitPos(),     // reference pos
                                                     jointHit.orphanClust()->isFw,     // fwd direction
                                                     readLen, // read length
                                                     jointHit.orphanClust()->cigar, // cigar string
                                                     readLen,       // fragment length
                                                     false);

                auto &qaln = jointAlignments.back();
                qaln.mateLen = readLen;
                qaln.mateCigar = "";
                qaln.matePos = 0;       // jointHit.rightClust->getTrFirstHitPos();
                qaln.mateIsFwd = false; // jointHit.rightClust->isFw;
                qaln.mateStatus = MateStatus::SINGLE_END;
				qaln.numHits = jointHit.orphanClust()->coverage;
                validHits.emplace_back(jointHit.tid, jointHit.orphanClust());
            }

            hctr.totAlignment += jointHits.size();

            // write puffkrak format output
            if (mopts->krakOut) {
                writeAlignmentsToKrakenDump(read, /* formatter,  */validHits, bstream);
            } else if (mopts->salmonOut) {
                writeAlignmentsToKrakenDump(read, /* formatter,  */validHits, bstream, false);
            } else if (jointHits.size() > 0 and !mopts->noOutput) {
                // write sam output for mapped reads
                writeAlignmentsToStreamSingle(read, formatter, jointAlignments, sstream,
                                              !mopts->noOrphan, mopts->justMap);
            } else if (jointHits.size() == 0 and !mopts->noOutput) {
                // write sam output for un-mapped reads
                writeUnmappedAlignmentsToStreamSingle(read, formatter, jointAlignments,
                                                      sstream, !mopts->noOrphan,
                                                      mopts->justMap);
            }

            // write them on cmd
            if (hctr.numReads > hctr.lastPrint + 1000000) {
                hctr.lastPrint.store(hctr.numReads.load());
                if (!mopts->quiet and iomutex->try_lock()) {
                    if (hctr.numReads > 0) {
                        std::cerr << "\r\r";
                    }
                    std::cerr << "saw " << hctr.numReads << " reads : "
                              << "pe / read = "
                              << hctr.peHits / static_cast<float>(hctr.numReads)
                              << " : se / read = "
                              << hctr.seHits / static_cast<float>(hctr.numReads) << ' ';
                    iomutex->unlock();
                }
            }
        } // for all reads in this job


        // dump output
        if (!mopts->noOutput) {
            // Get rid of last newline
            if (mopts->krakOut || mopts->salmonOut) {
                if (mopts->salmonOut && bstream.getBytes() > 0) {
                    BinWriter sbw(64);
                    sbw << bstream.getBytes();
                    outQueue->info("{}{}", sbw, bstream);
                }else if (mopts->krakOut) {
                    outQueue->info("{}", bstream);
                }
                bstream.clear();
            } else {
                std::string outStr(sstream.str());
                // Get rid of last newline
                if (!outStr.empty()) {
                    outStr.pop_back();
                    outQueue->info(std::move(outStr));
                }
                sstream.clear();
            }
        }

    } // processed all reads
}

template<typename PufferfishIndexT>
bool spawnProcessReadsthreads(
        uint32_t nthread,
        paired_parser *parser,
        PufferfishIndexT &pfi,
        MutexT &iomutex,
        std::shared_ptr<spdlog::logger> outQueue,
        HitCounters &hctr,
        AlignmentOpts *mopts) {

    std::vector<std::thread> threads;

    for (size_t i = 0; i < nthread; ++i) {

        threads.emplace_back(processReadsPair<PufferfishIndexT>,
                             parser,
                             std::ref(pfi),
                             &iomutex,
                             outQueue,
                             std::ref(hctr),
                             mopts);
    }
    for (auto &t : threads) { t.join(); }

    return true;
}

template<typename PufferfishIndexT>
bool spawnProcessReadsthreads(
        uint32_t nthread,
        single_parser *parser,
        PufferfishIndexT &pfi,
        MutexT &iomutex,
        std::shared_ptr<spdlog::logger> outQueue,
        HitCounters &hctr,
        AlignmentOpts *mopts) {

    std::vector<std::thread> threads;

    for (size_t i = 0; i < nthread; ++i) {

        threads.emplace_back(processReadsSingle<PufferfishIndexT>,
                             parser,
                             std::ref(pfi),
                             &iomutex,
                             outQueue,
                             std::ref(hctr),
                             mopts);
    }
    for (auto &t : threads) { t.join(); }

    return true;
}

void printAlignmentSummary(HitCounters &hctrs, std::shared_ptr<spdlog::logger> consoleLog) {
    consoleLog->info("Done mapping reads.");
    consoleLog->info("\n\n");
    consoleLog->info("=====");
    consoleLog->info("Observed {} reads", hctrs.numReads);
    consoleLog->info("Rate of Fragments with at least one found k-mer: {:03.2f}%",
                     (100.0 * static_cast<float>(hctrs.numMappedAtLeastAKmer)) / hctrs.numReads);
    consoleLog->info("Discordant Rate: {:03.2f}%",
                     (100.0 * static_cast<float>(hctrs.numOfOrphans)) / hctrs.numReads);
    consoleLog->info("Total reads Mapped: {}", (hctrs.numMapped));
    consoleLog->info("Mapping rate : {:03.2f}%", (100.0 * static_cast<float>(hctrs.numMapped)) / hctrs.numReads);
    consoleLog->info("Average # hits per read : {}", hctrs.totHits / static_cast<float>(hctrs.numReads));
    consoleLog->info("Total # of alignments : {}", hctrs.totAlignment);
    consoleLog->info("Total # of orphans : {}", hctrs.numOfOrphans);
    consoleLog->info("Total # of pe hits : {}", hctrs.peHits);
    consoleLog->info("Total # of total Hits : {}", hctrs.totHits);
    //consoleLog->info("Total # of valid hits : {}", hctrs.validHits);
    consoleLog->info("Max multimapping group : {}", hctrs.maxMultimapping);
    consoleLog->info("=====");
}

template<typename PufferfishIndexT>
bool alignReads(
        PufferfishIndexT &pfi,
        std::shared_ptr<spdlog::logger> consoleLog,
        AlignmentOpts *mopts) {

	/*std::string output;
	std::ifstream unmapped;
	//unmapped.open("/mnt/scratch4/mohsen/puffalign_test/quants/polyester_0.001/salmon_puf_single_noChange/aux_info/unmapped_names.txt");
	//unmapped.open("/mnt/scratch4/mohsen/puffalign_test/allMappings_strict_Equal_mohsen");
	unmapped.open("/mnt/scratch4/mohsen/puffalign_test/allMappings_equal_pe");
 	while (!unmapped.eof()) {
    	unmapped >> output;
		mopts->salmonUnmapped.insert(output);
	}
	unmapped.close();*/

    std::streambuf *outBuf;
    std::ofstream outFile;
    std::unique_ptr<std::ostream> outStream{nullptr};
    //bool haveOutputFile{false} ;
    std::shared_ptr<spdlog::logger> outLog{nullptr};
    if (!mopts->noOutput) {
        if (mopts->outname == "-") {
            outBuf = std::cout.rdbuf();
        } else {
            outFile.open(mopts->outname);
            outBuf = outFile.rdbuf();
            //haveOutputFile = true ;
        }

        // out stream to the buffer
        // it can be std::cerr or a file
        if (mopts->compressedOutput) {
            outStream.reset(new zstr::ostream(outBuf));
        } else {
            outStream.reset(new std::ostream(outBuf));
        }
        // the async queue size must be a power of 2
        size_t queueSize{268435456};
        spdlog::set_async_mode(queueSize);

        if (mopts->krakOut || mopts->salmonOut) {
            auto outputSink = std::make_shared<ostream_bin_sink_mt>(*outStream);
            outLog = std::make_shared<spdlog::logger>("puffer::outLog", outputSink);
            outLog->set_pattern("");
        } else {
            auto outputSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(*outStream);
            outLog = std::make_shared<spdlog::logger>("puffer::outLog", outputSink);
            outLog->set_pattern("%v");
        }
        // write the SAM Header
        // If nothing gets printed by this time we are in trouble
        if (mopts->krakOut || mopts->salmonOut) {
            writeKrakOutHeader(pfi, outLog, mopts);
        } else {
            writeSAMHeader(pfi, outLog);
        }
    }

    uint32_t nthread = mopts->numThreads;
    std::unique_ptr<paired_parser> pairParserPtr{nullptr};
    std::unique_ptr<single_parser> singleParserPtr{nullptr};

    size_t chunkSize{10000};
    MutexT iomutex;

    if (!mopts->singleEnd) {
        ScopedTimer timer(!mopts->quiet);
        HitCounters hctrs;
        consoleLog->info("mapping reads ... \n\n\n");
        std::vector<std::string> read1Vec = util::tokenize(mopts->read1, ',');
        std::vector<std::string> read2Vec = util::tokenize(mopts->read2, ',');

        if (read1Vec.size() != read2Vec.size()) {
            consoleLog->error("The number of provided files for"
                              "-1 and -2 are not same!");
            std::exit(1);
        }

        uint32_t nprod = (read1Vec.size() > 1) ? 2 : 1;
        pairParserPtr.reset(new paired_parser(read1Vec, read2Vec, nthread, nprod, chunkSize));
        pairParserPtr->start();

        spawnProcessReadsthreads(nthread, pairParserPtr.get(), pfi, iomutex,
                                 outLog, hctrs, mopts);

        pairParserPtr->stop();
        consoleLog->info("flushing output queue.");
        printAlignmentSummary(hctrs, consoleLog);
        if (outLog) { outLog->flush(); }
    } else {
        ScopedTimer timer(!mopts->quiet);
        HitCounters hctrs;
        consoleLog->info("mapping reads ... \n\n\n");
        std::vector<std::string> readVec = util::tokenize(mopts->unmatedReads, ',');

        uint32_t nprod = (readVec.size() > 1) ? 2 : 1;
        singleParserPtr.reset(new single_parser(readVec, nthread, nprod, chunkSize));
        singleParserPtr->start();

        spawnProcessReadsthreads(nthread, singleParserPtr.get(), pfi, iomutex,
                                 outLog, hctrs, mopts);

        singleParserPtr->stop();
        consoleLog->info("flushing output queue.");
        printAlignmentSummary(hctrs, consoleLog);
        if (outLog) { outLog->flush(); }
    }
    return true;
}

int pufferfishAligner(AlignmentOpts &alnargs) {

    auto consoleLog = spdlog::stderr_color_mt("console");
    bool success{false};
    auto indexDir = alnargs.indexDir;

    std::string indexType;
    {
        std::ifstream infoStream(indexDir + "/info.json");
        cereal::JSONInputArchive infoArchive(infoStream);
        infoArchive(cereal::make_nvp("sampling_type", indexType));
        std::cerr << "Index type = " << indexType << "\n";
        infoStream.close();
    }

    if (indexType == "dense") {
        PufferfishIndex pfi(indexDir);
        success = alignReads(pfi, consoleLog, &alnargs);

    } else if (indexType == "sparse") {
        PufferfishSparseIndex pfi(indexDir);
        success = alignReads(pfi, consoleLog, &alnargs);
    } else if (indexType == "lossy") {
        PufferfishLossyIndex pfi(indexDir);
        success = alignReads(pfi, consoleLog, &alnargs);
    }

    if (!success) {
        consoleLog->warn("Problem mapping.");
    }
    return 0;
}
