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
#include "parallel_hashmap/phmap.h"

#include "spdlog/spdlog.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/sinks/stdout_sinks.h"
#include "spdlog/sinks/ansicolor_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"

#include "FastxParser.hpp"

//index header
#include "PuffAligner.hpp"
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


#define MATCH_SCORE 1
#define MISMATCH_SCORE -1
#define GAP_SCORE -1

#define EPS 5

#define ALLOW_VERBOSE 0

using paired_parser = fastx_parser::FastxParser<fastx_parser::ReadPair>;
using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

using HitCounters = pufferfish::util::HitCounters;
using QuasiAlignment = pufferfish::util::QuasiAlignment;
using MateStatus = pufferfish::util::MateStatus;

using MutexT = std::mutex;

void joinReadsAndFilterSingle( pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& leftMemClusters,
                               //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>> &leftMemClusters,
                              std::vector<pufferfish::util::JointMems> &jointMemsList,
                              uint32_t perfectCoverage,
                              double coverageRatio,
                              bool verbose) {
    uint32_t maxCoverage{0};
    for (auto &leftClustItr : leftMemClusters) {
        // reference id
        size_t tid = leftClustItr.first;
        // left mem clusters
        auto &lClusts = *(leftClustItr.second);
#if ALLOW_VERBOSE
        if (verbose)
            std::cerr << "\ntid:" << leftClustItr.first << "\n";
#endif // ALLOW_VERBOSE
        // Compare the left clusters to the right clusters to filter by positional constraints
        std::vector<pufferfish::util::MemCluster> dummyCluster;
        for (auto lclust = lClusts.begin(); lclust != lClusts.end(); lclust++) {
            auto totalCoverage = lclust->coverage;
            if (totalCoverage >= coverageRatio * maxCoverage or totalCoverage == perfectCoverage) {
                jointMemsList.emplace_back(tid, lclust, dummyCluster.end(), 0, MateStatus::PAIRED_END_LEFT);
                uint32_t currCoverage = jointMemsList.back().coverage();
                if (maxCoverage < currCoverage) {
                    maxCoverage = currCoverage;
                }
            }
        }
    }
}


pufferfish::util::MergeResult joinReadsAndFilter(
                                                 pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& leftMemClusters,
                                                 pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& rightMemClusters,
                                                 //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>> &leftMemClusters,
                                                 //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>> &rightMemClusters,
                        std::vector<pufferfish::util::JointMems> &jointMemsList,
                        uint32_t maxFragmentLength,
                        uint32_t perfectCoverage,
                        double coverageRatio,
                        bool noDiscordant,
                        bool noOrphans,
                        bool verbose) {
    using pufferfish::util::MergeResult;
    MergeResult mergeRes{MergeResult::HAD_NONE};

#if ALLOW_VERBOSE
    if (verbose) { std::cerr << "\n[JOINREADSANDFILTER]\n";}
#endif // ALLOW_VERBOSE

    uint64_t maxLeft{0}, maxRight{0}, maxLeftCnt{0}, maxRightCnt{0};
    bool isMaxLeftAndRight = false;
    if (!noOrphans) {
        for (auto &kv : leftMemClusters) {
#if ALLOW_VERBOSE
          if (verbose) {std::cerr << "\ntid:" << kv.first << "\n"; }
#endif // ALLOW_VERBOSE
            auto &lClusts = *(kv.second);
            for (auto clust = lClusts.begin(); clust != lClusts.end(); clust++) {
                if (maxLeft == clust->coverage) {
                    maxLeftCnt += 1;
                } else if (maxLeft < clust->coverage) {
                    maxLeft = clust->coverage;
                    maxLeftCnt = 1;
                }
            }
        } // leftMemClusters
        for (auto &kv : rightMemClusters) {
#if ALLOW_VERBOSE
          if (verbose) { std::cerr << "\ntid:" << kv.first << "\n"; }
#endif // ALLOW_VERBOSE
          auto &rClusts = *(kv.second);
          for (auto clust = rClusts.begin(); clust != rClusts.end(); clust++) {
            if (maxRight == clust->coverage) {
              maxRightCnt += 1;
            } else if (maxRight < clust->coverage) {
              maxRight = clust->coverage;
              maxRightCnt = 1;
            }
          }
        } // rightMemClusters
    } // !noOrphans


    auto maxLeftOrRight = maxLeft > maxRight ? maxLeft : maxRight;

    //orphan reads should be taken care of maybe with a flag!
    uint32_t maxCoverage{0};
    uint8_t round{0};
    uint32_t sameTxpCount{0};

    while (round == 0 or (round == 1 and !jointMemsList.size() and !noDiscordant)) {
      for (auto &leftClustItr : leftMemClusters) {
//      for (auto leftClustItr = leftMemClusters.key_begin(); leftClustItr != leftMemClusters.key_end(); ++ leftClustItr) {
            // reference id
            size_t tid = leftClustItr.first;
//            size_t tid = leftClustItr->first;
            // left mem clusters
            auto &lClusts = *(leftClustItr.second);
//            auto &lClusts = leftMemClusters.cache_index(leftClustItr->second);
            // right mem clusters for the same reference id
            auto &rClusts = rightMemClusters[tid];

            // Compare the left clusters to the right clusters to filter by positional constraints
            for (auto lclust = lClusts.begin(); lclust != lClusts.end(); lclust++) {

                for (auto rclust = rClusts.begin(); rclust != rClusts.end(); rclust++) {
                    // if both the left and right clusters are oriented in the same direction, skip this pair
                    // NOTE: This should be optional as some libraries could allow this.
                    if (round == 0 and lclust->isFw == rclust->isFw) { // if priority 0, ends should be concordant
                        continue;
                    }

                    // FILTER 1
                    // filter read pairs based on the fragment length which is approximated by the distance between the left most start and right most hit end
                    int32_t fragmentLen = rclust->lastRefPos() + rclust->lastMemLen() - lclust->firstRefPos();
                    if (lclust->firstRefPos() > rclust->firstRefPos()) {
                        fragmentLen = lclust->lastRefPos() + lclust->lastMemLen() - rclust->firstRefPos();
                    }
                    if (fragmentLen < 0) { // @fatemeh : should we even be checking for this?
                        std::cerr << "Fragment length cannot be smaller than zero!\n";
                        exit(1);
                    }

                    // FILTERING fragments with size smaller than maxFragmentLength
                    // FILTER just in case of priority 0 (round 0)
                    if ((fragmentLen < maxFragmentLength) or (round > 0)) {
                        // This will add a new potential mapping. Coverage of a mapping for read pairs is left->coverage + right->coverage
                        // If we found a perfect coverage, we would only add those mappings that have the same perfect coverage
                        auto totalCoverage = lclust->coverage + rclust->coverage;
                        if ( (totalCoverage >= coverageRatio * maxCoverage) or
                              (totalCoverage == perfectCoverage) ) {
                            ++sameTxpCount;
                            jointMemsList.emplace_back(tid, lclust, rclust, fragmentLen);
                            uint32_t currCoverage = jointMemsList.back().coverage();
                            if (maxCoverage < currCoverage) {
                                maxCoverage = currCoverage;
                                if ( (lclust->coverage < maxLeft) or (rclust->coverage < maxRight)) {
                                    isMaxLeftAndRight = false;
                                } else {
                                    isMaxLeftAndRight = true;
                                }
                            }
                        }
                    }
                }
            }
        } // @fatemeh : this nesting just seems too many levels deep.  Can we re-work the logic here to make things simpler?
        round++;
    }

#if ALLOW_VERBOSE
    // If we couldn't find any pair and we are allowed to add orphans
    if (verbose) {
        std::cerr << "isMaxLeftAndRight:" << isMaxLeftAndRight << "\n";
    }
#endif // ALLOW_VERBOSE
    
    bool leftOrphan = false; bool rightOrphan = false;
    if (!noOrphans and (!jointMemsList.size() or !isMaxLeftAndRight or maxLeftCnt > 1 or maxRightCnt > 1)) {
        auto orphanFiller = [&jointMemsList, &maxCoverage, &coverageRatio, &maxLeftOrRight, &leftOrphan, &rightOrphan]
        (pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>> &memClusters,
                 bool isLeft) {
            // fragmentLen is set to 0
            std::vector<pufferfish::util::MemCluster> dummyCluster;
            for (auto &clustItr : memClusters) {
                // reference id
                size_t tid = clustItr.first;
                // left mem clusters
                auto &Clusts = *(clustItr.second);
                for (auto clust = Clusts.begin(); clust != Clusts.end(); clust++) {
                    if (clust->coverage >= coverageRatio * maxLeftOrRight) {
                        if (isLeft) {
                            leftOrphan = true;
                            jointMemsList.emplace_back(tid, clust, dummyCluster.end(), 0, MateStatus::PAIRED_END_LEFT);
                        } else {
                            rightOrphan = true;
                            jointMemsList.emplace_back(tid, dummyCluster.end(), clust, 0, MateStatus::PAIRED_END_RIGHT);
                        }
                        uint32_t currCoverage = jointMemsList.back().coverage();
                        if (maxCoverage < currCoverage) {
                            maxCoverage = currCoverage;
                        }
                    }
                }
            }
        };
        orphanFiller(leftMemClusters, true);
        orphanFiller(rightMemClusters, false);
    }
    if (sameTxpCount == 0) {
      if (leftOrphan and !rightOrphan) {
        mergeRes = MergeResult::HAD_ONLY_LEFT;
      } else if (!leftOrphan and rightOrphan) {
        mergeRes = MergeResult::HAD_ONLY_RIGHT;
      } else if (leftOrphan and rightOrphan) {
        mergeRes = MergeResult::HAD_EMPTY_INTERSECTION;
      } else {
        mergeRes = MergeResult::HAD_NONE;
      }
    } else {
      mergeRes = (round == 0) ? MergeResult::HAD_CONCORDANT : MergeResult::HAD_DISCORDANT;
    }

#if ALLOW_VERBOSE
    if (verbose) {
        std::cerr << "\nBefore filter " << jointMemsList.size() << " maxCov:" << maxCoverage << "\n";
    }
#endif // ALLOW_VERBOSE

    jointMemsList.erase(std::remove_if(jointMemsList.begin(), jointMemsList.end(),
                                       [&maxCoverage, coverageRatio](pufferfish::util::JointMems &pairedReadMems) -> bool {
                                           return pairedReadMems.coverage() < coverageRatio * maxCoverage;
                                       }),
                        jointMemsList.end());

#if ALLOW_VERBOSE
    if (verbose) {
        std::cerr << "\nAfter:" << jointMemsList.size() << " maxCov:" << maxCoverage << "\n";
        std::cerr << "\n[END OF JOINREADSANDFILTER]\n";
    }
#endif // ALLOW_VERBOSE

    return mergeRes;
}

void recoverOrphans(std::vector<pufferfish::util::MemCluster> &recoveredMemClusters, std::vector<pufferfish::util::JointMems> &jointMemsList, PuffAligner& puffaligner, bool verbose) {
  puffaligner.orphanRecoveryMemCollection.reserve(2 * jointMemsList.size() + 1);
  recoveredMemClusters.reserve(2 * jointMemsList.size() + 1);
  for (auto& jointMem : jointMemsList) {
    auto& leftClust = jointMem.leftClust;
    auto& rightClust = jointMem.rightClust;
    auto tid = jointMem.tid;
    if (jointMem.mateStatus == MateStatus::PAIRED_END_LEFT) {
      bool recovered = puffaligner.recoverSingleOrphan(*leftClust, recoveredMemClusters, tid, true, verbose);
     if (recovered) {
       jointMem.rightClust = recoveredMemClusters.begin() + recoveredMemClusters.size() - 1;
       jointMem.recovered = true;
       jointMem.mateStatus = MateStatus::PAIRED_END_PAIRED;
    }
  } else if (jointMem.mateStatus == MateStatus::PAIRED_END_RIGHT) {
      bool recovered = puffaligner.recoverSingleOrphan(*rightClust, recoveredMemClusters, tid, false, verbose);
      if (recovered) {
        jointMem.leftClust = recoveredMemClusters.begin() + recoveredMemClusters.size() - 1;
        jointMem.recovered = true;
        jointMem.mateStatus = MateStatus::PAIRED_END_PAIRED;
      }
    }
  }
}


//===========
// PAIRED END
//============
template<typename PufferfishIndexT>
void processReadsPair(paired_parser *parser,
                      PufferfishIndexT &pfi,
                      MutexT *iomutex,
                      std::shared_ptr<spdlog::logger> outQueue,
                      HitCounters &hctr,
                      phmap::flat_hash_set<std::string>& gene_names,
                      AlignmentOpts *mopts) {
    MemCollector<PufferfishIndexT> memCollector(&pfi);
    memCollector.configureMemClusterer(mopts->maxAllowedRefsPerHit);

    auto logger = spdlog::get("stderrLog");
    fmt::MemoryWriter sstream;
    BinWriter bstream;

    //size_t batchSize{2500} ;
    size_t readLen{0}, mateLen{0};
    size_t totLen{0};


    pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>> leftHits;
    pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>> rightHits;
    //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>> leftHits;
    //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>> rightHits;

    std::vector<pufferfish::util::MemCluster> recoveredHits;
    std::vector<pufferfish::util::JointMems> jointHits;
    PairedAlignmentFormatter<PufferfishIndexT *> formatter(&pfi);
    pufferfish::util::QueryCache qc;
    std::vector<pufferfish::util::MemCluster> all;

    //@fatemeh Initialize aligner ksw
    ksw2pp::KSW2Aligner aligner(mopts->matchScore, mopts->missMatchScore);
    ksw2pp::KSW2Config config;

    config.dropoff = -1;
    config.gapo = mopts->gapOpenPenalty;
    config.gape = mopts->gapExtendPenalty;
    config.bandwidth = 15;
    config.flag = 0;
    config.flag |= KSW_EZ_RIGHT;
    config.flag |= KSW_EZ_SCORE_ONLY;
    aligner.config() = config;

    constexpr const int32_t invalidScore = std::numeric_limits<int32_t>::min();

    // don't think we should have these, they should
    // be baked into the index I think.
    bool filterGenomics = mopts->filterGenomics;
    bool filterMicrobiom = mopts->filterMicrobiom;

    auto rg = parser->getReadGroup();


    PuffAligner puffaligner(pfi.refseq_, pfi.refAccumLengths_, pfi.k(),
                            mopts, aligner, true);


    //For filtering reads
    bool verbose = mopts->verbose;
    auto &txpNames = pfi.getRefNames();
    while (parser->refill(rg)) {
        for (auto &rpair : rg) {
            readLen = rpair.first.seq.length();
            mateLen = rpair.second.seq.length();
            totLen = readLen + mateLen;

            ++hctr.numReads;

            jointHits.clear();
            leftHits.clear();
            rightHits.clear();
            recoveredHits.clear();
            memCollector.clear();

            // There is no way to revocer the following case other than aligning indels
            //verbose = rpair.first.seq == "CAGTGAGCCAAGATGGCGCCACTGCACTCCAGCCTGGGCAAAAAGAAACTCCATCTAAAAAAAAAAAAAAAAAAAAAAAAAAGAGAAAACCCTGGTCCCT" or
            //          rpair.second.seq == "CAGTGAGCCAAGATGGCGCCACTGCACTCCAGCCTGGGCAAAAAGAAACTCCATCTAAAAAAAAAAAAAAAAAAAAAAAAAAGAGAAAACCCTGGTCCCT";

            // The only way to get it right is with non-heuristic chaining
            // verbose = rpair.first.seq == "AGCAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGTGGTGGGGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTAGAGAGGCACCAGCA" or
            //           rpair.second.seq == "AGCAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGTGGTGGGGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTAGAGAGGCACCAGCA";

            //verbose = rpair.first.name == "mason_sample5_primary_1M_random.fasta.000050010/1";
            bool lh = memCollector(rpair.first.seq,
                                   qc,
                                   true, // isLeft
                                   verbose);
            bool rh = memCollector(rpair.second.seq,
                                   qc,
                                   false, // isLeft
                                   verbose);
            memCollector.findChains(rpair.first.seq,
                                   leftHits,
                                   mopts->maxSpliceGap,
                                   MateStatus::PAIRED_END_LEFT,
                                   mopts->heuristicChaining,
                                   true, // isLeft
                                   verbose);
            memCollector.findChains(rpair.second.seq,
                                   rightHits,
                                   mopts->maxSpliceGap,
                                   MateStatus::PAIRED_END_RIGHT,
                                   mopts->heuristicChaining,
                                   false, // isLeft
                                   verbose);

            all.clear();
            hctr.numMappedAtLeastAKmer += (leftHits.size() > 0 || rightHits.size() > 0) ? 1 : 0;
            //do intersection on the basis of
            //performance, or going towards selective alignment
            //otherwise orphan

            // We also handle orphans inside this function
            //TODO uncomment this part if you want to get back to the original joinReadsAndFilter

            auto mergeRes = joinReadsAndFilter(leftHits, rightHits, jointHits,
                               mopts->maxFragmentLength,
                               totLen,
                               mopts->scoreRatio,
                               mopts->noDiscordant,
                               mopts->noOrphan,
                               verbose);

            PuffAligner puffaligner(pfi.refseq_, pfi.refAccumLengths_, pfi.k(),
                                    rpair.first.seq, rpair.second.seq, mopts, aligner, true);


            bool mergeStatusOR = (mergeRes == pufferfish::util::MergeResult::HAD_EMPTY_INTERSECTION or
                                  mergeRes == pufferfish::util::MergeResult::HAD_ONLY_LEFT or
                                  mergeRes == pufferfish::util::MergeResult::HAD_ONLY_RIGHT);

            if ( mergeStatusOR and mopts->recoverOrphans) {
              recoverOrphans(recoveredHits, jointHits, puffaligner, verbose);
            }

            std::vector<QuasiAlignment> jointAlignments;

            hctr.peHits += jointHits.size();

#if ALLOW_VERBOSE
            if (verbose)
                std::cerr<<"Number of hits: "<<jointHits.size()<<"\n";
#endif // ALLOW_VERBOSE

            if (!mopts->justMap) {
              puffaligner.clear();
                int32_t bestScore = invalidScore;
                std::vector<decltype(bestScore)> scores(jointHits.size(), bestScore);
                size_t idx{0};

                bool bestScoreGenomic{false};
                bool bestScoreTxpomic{false};
                std::map<int32_t, std::vector<int32_t>> transcript_set;
                for (auto &jointHit : jointHits) {
                  auto hitScore = puffaligner.calculateAlignments(rpair.first.seq, rpair.second.seq, jointHit, hctr, verbose);
                    /*if (hitScore < 0)
                        hitScore = invalidScore;*/
                    scores[idx] = hitScore;
                    const std::string& ref_name = txpNames[jointHit.tid];
                    if (filterGenomics or filterMicrobiom) {
                        if (hitScore > bestScore) {
                            if (gene_names.find(ref_name) != gene_names.end()) {
                                bestScoreGenomic = true;
                                bestScoreTxpomic = false;
                            } else {
                                bestScoreGenomic = false;
                                bestScoreTxpomic = true;
                            }
                        } else if (hitScore == bestScore) {
                            if (gene_names.find(ref_name) != gene_names.end())
                                bestScoreGenomic = true;
                            else
                                bestScoreTxpomic = true;
                        }
                        // Genomic alignments should not be reported in this case
                        if (gene_names.find(ref_name) != gene_names.end())
                            scores[idx] = invalidScore;
                    }

                    bestScore = (hitScore > bestScore) ? hitScore : bestScore;

#if ALLOW_VERBOSE
                    if (verbose)
                        std::cerr << txpNames[jointHit.tid] << " " << jointHit.coverage() << " "
                                  << jointHit.leftClust->mems[0].tpos << "\n";
#endif // ALLOW_VERBOSE

                    if (!mopts->genomicReads) {
                        // removing dupplicate hits from a read to the same transcript
                        if (transcript_set.find(jointHit.tid) == transcript_set.end()) {
                            std::vector<int32_t> coverage_hit;
                            transcript_set[jointHit.tid] = coverage_hit;
                            transcript_set[jointHit.tid].push_back(jointHit.coverage());
                            transcript_set[jointHit.tid].push_back(idx);
                        } else if (jointHit.coverage() > transcript_set[jointHit.tid][0]) {
                            transcript_set[jointHit.tid][0] = jointHit.coverage();
                            scores[transcript_set[jointHit.tid][1]] = invalidScore;
                            transcript_set[jointHit.tid][1] = idx;
                        } else {
                            scores[idx] = invalidScore;
                        }
                    }
                    ++idx;
                }
                if (filterGenomics and bestScoreGenomic and !bestScoreTxpomic) {
                    // This read is likely come from decoy sequence and should be discarded from reference alignments
                    continue;
                }
                if (filterMicrobiom and bestScoreGenomic) {
                    continue;
                }

                if (mopts->bestStrata or mopts->primaryAlignment) {
                    bool primaryDetected = false;
                    int myctr = 0;
                    for (auto &jointHit : jointHits) {
                        if (primaryDetected and mopts->primaryAlignment) {
                            scores[myctr] = invalidScore;
                        }
                        if (jointHit.coverage() < bestScore)
                            scores[myctr] = invalidScore;
                        else {
                            primaryDetected = true;
                        }
                        myctr++;
                    }
                }

                // Filter out these alignments with low scores
                uint32_t ctr{0};
                if (bestScore > invalidScore) {
                    jointHits.erase(
                            std::remove_if(jointHits.begin(), jointHits.end(),
                                           [&ctr, &scores, bestScore](pufferfish::util::JointMems &) -> bool {
                                               bool rem = (scores[ctr] == invalidScore);
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

            if (jointHits.size() > 200) {
                std::sort(jointHits.begin(), jointHits.end(),
                          [](const auto &lhs, const auto &rhs) {
                              return lhs.alignmentScore < rhs.alignmentScore;
                          });
                jointHits.erase(jointHits.begin() + 200, jointHits.end());
            }

            hctr.totHits += !jointHits.empty() && !jointHits.back().isOrphan() ? 1 : 0;;
            hctr.numMapped += !jointHits.empty() ? 1 : 0;
            if (mopts->noOrphan)
                hctr.numOfOrphans += jointHits.empty() && (lh || rh);
            else
                hctr.numOfOrphans += !jointHits.empty() && (jointHits.back().isOrphan()) ? 1 : 0;

            if (jointHits.size() > hctr.maxMultimapping)
                hctr.maxMultimapping = jointHits.size();

            for (auto &jointHit : jointHits) {
                //FIXME : This part needs to be taken care of
                //If graph returned failure for one of the ends --> should be investigated more.
                /*if (!mopts->justMap and ((!jointHit.isOrphan() and (jointHit.leftClust->score == invalidScore
                             or jointHit.rightClust->score == invalidScore ))
                             or (jointHit.isLeftAvailable() and jointHit.leftClust->score == invalidScore )
                             or (jointHit.isRightAvailable() and jointHit.rightClust->score == invalidScore ))
                                ) {
                  if (verbose)
                    std::cerr << "Failed: " << rpair.first.name << "\n";
                  continue;
                }*/
                if (jointHit.isOrphan()) {
#if ALLOW_VERBOSE
                    if (verbose) {
                        std::cerr << "orphan\t";
                        std::cerr << jointHit.tid << "\t";
                        std::cerr << jointHit.orphanClust()->getTrFirstHitPos() << "\t";     // reference pos
                        std::cerr << (uint32_t) jointHit.orphanClust()->isFw << "\t";     // fwd direction
                        std::cerr << readLen << "\t"; // read length
                        std::cerr << jointHit.orphanClust()->cigar << "\t"; // cigar string
                        std::cerr << jointHit.fragmentLen << "\n";
                    }
#endif // ALLOW_VERBOSE
                    jointAlignments.emplace_back(jointHit.tid,           // reference id
                                                 0, //jointHit.orphanClust()->getTrFirstHitPos(),     // reference pos
                                                 jointHit.orphanClust()->isFw,     // fwd direction
                                                 readLen, // read length
                                                 jointHit.orphanClust()->cigar, // cigar string
                                                 jointHit.fragmentLen,       // fragment length
                                                 false);
                    auto &qaln = jointAlignments.back();
                    qaln.numHits = jointHit.orphanClust()->coverage;
                    qaln.mateStatus = jointHit.mateStatus;
                } else {
#if ALLOW_VERBOSE
                    if (verbose) {
                        std::cerr << "paired\t";
                        std::cerr << jointHit.tid << "\t";
                        std::cerr << jointHit.leftClust->getTrFirstHitPos() << "\t";     // reference pos
                        std::cerr << (uint32_t) jointHit.leftClust->isFw << "\t";     // fwd direction
                        std::cerr << readLen << "\t"; // read length
                        std::cerr << jointHit.leftClust->cigar << "\t"; // cigar string
                        std::cerr << jointHit.fragmentLen << "\n";
                    }
#endif // ALLOW_VERBOSE
                    jointAlignments.emplace_back(jointHit.tid,           // reference id
                                                 jointHit.leftClust->getTrFirstHitPos(),     // reference pos
                                                 jointHit.leftClust->isFw,     // fwd direction
                                                 readLen, // read length
                                                 jointHit.leftClust->cigar, // cigar string
                                                 jointHit.fragmentLen,       // fragment length
                                                 true);         // properly paired
                    // Fill in the mate info
                    auto &qaln = jointAlignments.back();
                    qaln.mateLen = mateLen;
                    qaln.mateCigar = jointHit.rightClust->cigar;
                    qaln.matePos = jointHit.rightClust->getTrFirstHitPos();
                    qaln.mateIsFwd = jointHit.rightClust->isFw;
                    qaln.mateStatus = MateStatus::PAIRED_END_PAIRED;
                    qaln.numHits = jointHit.coverage();
                    qaln.score = jointHit.leftClust->coverage;
                    qaln.mateScore = jointHit.rightClust->coverage;
                }
            }

            hctr.totAlignment += jointAlignments.size();
#if ALLOW_VERBOSE
            if (verbose) {
                std::cerr << " passed filling jointAlignments:" << hctr.totAlignment << "\n";
            }
#endif // ALLOW_VERBOSE
            if (!mopts->noOutput) {
                if (mopts->krakOut) {
                  // writeAlignmentsToKrakenDump(rpair, /* formatter,  */jointHits, bstream);
                } else if (mopts->salmonOut) {
                  // writeAlignmentsToKrakenDump(rpair, /* formatter,  */jointHits, bstream, false);
                } else if (jointAlignments.size() > 0) {
                    writeAlignmentsToStream(rpair, formatter, jointAlignments, sstream, !mopts->noOrphan,
                                            mopts->justMap);
                } else if (jointAlignments.size() == 0) {
                    writeUnmappedAlignmentsToStream(rpair, formatter, jointAlignments, sstream, !mopts->noOrphan,
                                                    mopts->justMap);
                }
            }

            // write them on cmd
            if (hctr.numReads > hctr.lastPrint + 100000) {
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
            //puffaligner.clear();
        } // for all reads in this job
        // dump output
        if (!mopts->noOutput) {
            // Get rid of last newline
            if (mopts->salmonOut) {
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
                        phmap::flat_hash_set<std::string>& gene_names,
                        AlignmentOpts *mopts) {
    MemCollector<PufferfishIndexT> memCollector(&pfi);
    memCollector.configureMemClusterer(mopts->maxAllowedRefsPerHit);

    auto logger = spdlog::get("stderrLog");
    fmt::MemoryWriter sstream;
    BinWriter bstream;
    //size_t batchSize{2500} ;
    size_t readLen{0};
    //size_t totLen{0};

    pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>> leftHits;
    std::vector<pufferfish::util::JointMems> jointHits;
    PairedAlignmentFormatter<PufferfishIndexT *> formatter(&pfi);
    pufferfish::util::QueryCache qc;
    std::vector<pufferfish::util::MemCluster> all;

    //@fatemeh Initialize aligner ksw
    ksw2pp::KSW2Aligner aligner(mopts->matchScore, mopts->missMatchScore);
    ksw2pp::KSW2Config config;

    config.dropoff = -1;
    config.gapo = mopts->gapOpenPenalty;
    config.gape = mopts->gapExtendPenalty;
    config.bandwidth = 10;
    config.flag = 0;
    config.flag |= KSW_EZ_RIGHT;
    aligner.config() = config;

    constexpr const int32_t invalidScore = std::numeric_limits<int32_t>::min();

    auto &txpNames = pfi.getRefNames();
    auto rg = parser->getReadGroup();
    while (parser->refill(rg)) {
        for (auto &read : rg) {
            readLen = read.seq.length();
            auto totLen = readLen;
            bool verbose = false;
            //if (verbose) std::cerr << read.name << "\n";
            ++hctr.numReads;

            jointHits.clear();
            leftHits.clear();
            memCollector.clear();


            bool filterGenomics = mopts->filterGenomics;
            bool filterMicrobiom = mopts->filterMicrobiom;

            bool lh = memCollector(read.seq,
                                   qc,
                                   true, // isLeft
                                   verbose);
            memCollector.findChains(read.seq,
                                   leftHits,
                                   mopts->maxSpliceGap,
                                   MateStatus::SINGLE_END,
                                   mopts->heuristicChaining,
                                   true, // isLeft
                                   verbose);
            (void) lh;
            all.clear();
            joinReadsAndFilterSingle(leftHits, jointHits,
                                     totLen,
                                     mopts->scoreRatio,
                                     verbose);

            std::vector<QuasiAlignment> jointAlignments;
            std::vector<std::pair<uint32_t, std::vector<pufferfish::util::MemCluster>::iterator>> validHits;

            if (mopts->justMap) {
//                ksw2pp::KSW2Aligner aligner(mopts->matchScore, mopts->missMatchScore);
                PuffAligner puffaligner(pfi.refseq_, pfi.refAccumLengths_, pfi.k(),
                                        read.seq, NULL, mopts, aligner, jointHits.size() > 1);

                int32_t bestScore = invalidScore;
                std::vector<decltype(bestScore)> scores(jointHits.size(), bestScore);
                size_t idx{0};

                bool bestScoreGenomic{false};
                bool bestScoreTxpomic{false};
                std::map<int32_t, std::vector<int32_t>> transcript_set;
                for (auto &jointHit : jointHits) {
                    int32_t hitScore = puffaligner.calculateAlignments(jointHit, hctr, verbose);
                    scores[idx] = hitScore;

                    std::string ref_name = txpNames[jointHit.tid];
                    if (filterGenomics or filterMicrobiom) {
                        if (hitScore > bestScore) {
                            if (gene_names.find(ref_name) != gene_names.end()) {
                                bestScoreGenomic = true;
                                bestScoreTxpomic = false;
                            } else {
                                bestScoreGenomic = false;
                                bestScoreTxpomic = true;
                            }
                        } else if (hitScore == bestScore) {
                            if (gene_names.find(ref_name) != gene_names.end())
                                bestScoreGenomic = true;
                            else
                                bestScoreTxpomic = true;
                        }
                        // Genomic alignments should not be reported in this case
                        if (gene_names.find(ref_name) != gene_names.end())
                            scores[idx] = invalidScore;
                    }

                    bestScore = (hitScore > bestScore) ? hitScore : bestScore;

                    if (verbose)
                        std::cerr << txpNames[jointHit.tid] << " " << jointHit.coverage() << " "
                                  << jointHit.leftClust->mems[0].tpos << "\n";
                    if (!mopts->genomicReads) {
                        // removing dupplicate hits from a read to the same transcript
                        if (transcript_set.find(jointHit.tid) == transcript_set.end()) {
                            std::vector<int32_t> coverage_hit;
                            transcript_set[jointHit.tid] = coverage_hit;
                            transcript_set[jointHit.tid].push_back(jointHit.orphanClust()->coverage);
                            transcript_set[jointHit.tid].push_back(idx);
                        } else if (jointHit.orphanClust()->coverage > transcript_set[jointHit.tid][0]) {
                            transcript_set[jointHit.tid][0] = jointHit.orphanClust()->coverage;
                            scores[transcript_set[jointHit.tid][1]] = invalidScore;
                            transcript_set[jointHit.tid][1] = idx;
                        } else {
                            scores[idx] = invalidScore;
                        }
                    }
                    ++idx;
                }

                if (filterGenomics and bestScoreGenomic and !bestScoreTxpomic) {
                    // This read is likely come from the genome and should be discarded from txptomic alignments
                    continue;
                }
                if (filterMicrobiom and bestScoreGenomic) {
                    continue;
                }

                if (mopts->bestStrata or mopts->primaryAlignment) {
                    bool primaryDetected = false;
                    int myctr = 0;
                    for (auto &jointHit : jointHits) {
                        if (primaryDetected and mopts->primaryAlignment) {
                            scores[myctr] = invalidScore;
                        }
                        if (jointHit.coverage() < bestScore)
                            scores[myctr] = invalidScore;
                        else {
                            primaryDetected = true;
                        }
                        myctr++;
                    }
                }
                // Filter out mappings with low alignment scores
                uint32_t ctr{0};
                if (bestScore > invalidScore) {
                    jointHits.erase(
                            std::remove_if(jointHits.begin(), jointHits.end(),
                                           [&ctr, &scores, bestScore](pufferfish::util::JointMems &) -> bool {
                                               bool rem = (scores[ctr] == invalidScore);
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

            if (jointHits.size() > 200) {
                std::sort(jointHits.begin(), jointHits.end(),
                          [](const auto &lhs, const auto &rhs) {
                              return lhs.alignmentScore < rhs.alignmentScore;
                          });
                jointHits.erase(jointHits.begin() + 200, jointHits.end());
            }

            hctr.numMappedAtLeastAKmer += jointHits.size() > 0 ? 1 : 0;
            hctr.totHits += jointHits.size();
            hctr.seHits += jointHits.size();
            hctr.numMapped += !jointHits.empty() ? 1 : 0;
            if (jointHits.size() > hctr.maxMultimapping) {
                hctr.maxMultimapping = jointHits.size();
            }

            for (auto &jointHit : jointHits) {
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
              //writeAlignmentsToKrakenDump(read, // formatter,  
              //                            validHits, bstream);
            } else if (mopts->salmonOut) {
              //writeAlignmentsToKrakenDump(read, / formatter,  
              //                              validHits, bstream, false);
            } else if (jointHits.size() > 0 and !mopts->noOutput) {
                // write sam output for mapped reads
                writeAlignmentsToStreamSingle(read, formatter, jointAlignments, sstream,
                                              !mopts->noOrphan, mopts->justMap);
            } else if (jointHits.size() == 0 and !mopts->noOutput) {
                // write sam output for un-mapped reads
                writeUnmappedAlignmentsToStreamSingle(read, formatter, jointAlignments,
                                                      sstream, !mopts->noOrphan, mopts->justMap);
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
                } else if (mopts->krakOut) {
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

//===========
// PAIRED END
//============
template<typename PufferfishIndexT>
bool spawnProcessReadsThreads(
        uint32_t nthread,
        paired_parser *parser,
        PufferfishIndexT &pfi,
        MutexT &iomutex,
        std::shared_ptr<spdlog::logger> outQueue,
        HitCounters &hctr,
        phmap::flat_hash_set<std::string>& gene_names,
        AlignmentOpts *mopts) {

    std::vector<std::thread> threads;

    for (size_t i = 0; i < nthread; ++i) {

        threads.emplace_back(processReadsPair<PufferfishIndexT>,
                             parser,
                             std::ref(pfi),
                             &iomutex,
                             outQueue,
                             std::ref(hctr),
                             std::ref(gene_names),
                             mopts);
    }
    for (auto &t : threads) { t.join(); }

    return true;
}

//===========
// SINGLE END
//============
template<typename PufferfishIndexT>
bool spawnProcessReadsThreads(
        uint32_t nthread,
        single_parser *parser,
        PufferfishIndexT &pfi,
        MutexT &iomutex,
        std::shared_ptr<spdlog::logger> outQueue,
        HitCounters &hctr,
        phmap::flat_hash_set<std::string>& gene_names,
        AlignmentOpts *mopts) {

    std::vector<std::thread> threads;

    for (size_t i = 0; i < nthread; ++i) {

        threads.emplace_back(processReadsSingle<PufferfishIndexT>,
                             parser,
                             std::ref(pfi),
                             &iomutex,
                             outQueue,
                             std::ref(hctr),
                             std::ref(gene_names),
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
    consoleLog->info("Average # hits per read : {}", hctrs.totAlignment / static_cast<float>(hctrs.numReads));
    consoleLog->info("Total # of alignments : {}", hctrs.totAlignment);
    consoleLog->info("Total # of orphans : {}", hctrs.numOfOrphans);
    consoleLog->info("Total # of pe hits : {}", hctrs.peHits);
    consoleLog->info("Total # of total Hits : {}", hctrs.totHits);
    //consoleLog->info("Total # of valid hits : {}", hctrs.validHits);
    consoleLog->info("Max multimapping group : {}", hctrs.maxMultimapping);
    consoleLog->info("Total number of alignment attempts : {}", hctrs.totalAlignmentAttempts);
    consoleLog->info("Number of skipped alignments because of cache hits : {}", hctrs.skippedAlignments_byCache);
    consoleLog->info("Number of skipped alignments because of perfect chains : {}", hctrs.skippedAlignments_byCov);

    consoleLog->info("Number of cigar strings which are fixed: {}", hctrs.cigar_fixed_count);
    consoleLog->info("=====");
}

template<typename PufferfishIndexT>
bool alignReads(
        PufferfishIndexT &pfi,
        std::shared_ptr<spdlog::logger> consoleLog,
        AlignmentOpts *mopts) {

    std::streambuf *outBuf;
    std::ofstream outFile;
    std::unique_ptr<std::ostream> outStream{nullptr};
    //bool haveOutputFile{false} ;
    std::shared_ptr<spdlog::logger> outLog{nullptr};

    phmap::flat_hash_set<std::string> gene_names;
    if (mopts->filterGenomics or mopts->filterMicrobiom) {
        std::ifstream gene_names_file;
        std::string gene_name;
        gene_names_file.open(mopts->genesNamesFile);
        if (!gene_names_file) {
            std::cerr << "Genomic file does not exist\n";
            exit(1);
        }
        std::string sepStr = " \t";
        while (gene_names_file >> gene_name) {
            auto processedName =
                    gene_name.substr(0, gene_name.find_first_of(sepStr));
            gene_names.insert(processedName);
        }
        gene_names_file.close();
    }


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
            writeSAMHeader(pfi, outLog, mopts->filterGenomics or mopts->filterMicrobiom, gene_names);
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
        std::vector<std::string> read1Vec = pufferfish::util::tokenize(mopts->read1, ',');
        std::vector<std::string> read2Vec = pufferfish::util::tokenize(mopts->read2, ',');

        if (read1Vec.size() != read2Vec.size()) {
            consoleLog->error("The number of provided files for"
                              "-1 and -2 are not same!");
            std::exit(1);
        }

        uint32_t nprod = (read1Vec.size() > 1) ? 2 : 1;
        pairParserPtr.reset(new paired_parser(read1Vec, read2Vec, nthread, nprod, chunkSize));
        pairParserPtr->start();
        spawnProcessReadsThreads(nthread, pairParserPtr.get(), pfi, iomutex,
                                 outLog, hctrs, gene_names, mopts);
        pairParserPtr->stop();
        consoleLog->info("flushing output queue.");
        printAlignmentSummary(hctrs, consoleLog);
        if (outLog) { outLog->flush(); }
    } else {
        ScopedTimer timer(!mopts->quiet);
        HitCounters hctrs;
        consoleLog->info("mapping reads ... \n\n\n");
        std::vector<std::string> readVec = pufferfish::util::tokenize(mopts->unmatedReads, ',');

        uint32_t nprod = (readVec.size() > 1) ? 2 : 1;
        singleParserPtr.reset(new single_parser(readVec, nthread, nprod, chunkSize));
        singleParserPtr->start();

        spawnProcessReadsThreads(nthread, singleParserPtr.get(), pfi, iomutex,
                                 outLog, hctrs, gene_names, mopts);

        singleParserPtr->stop();
        consoleLog->info("flushing output queue.");
        printAlignmentSummary(hctrs, consoleLog);
        if (outLog) { outLog->flush(); }
    }
    return true;
}

template<typename PufferfishIndexT>
bool alignReadsWrapper(
        PufferfishIndexT &pfi,
        std::shared_ptr<spdlog::logger> consoleLog,
        AlignmentOpts *mopts) {
    bool res = true;
    if (mopts->listOfReads) {
        if (mopts->singleEnd) {
            std::string unmatedReadsFile = mopts->unmatedReads;
            std::ifstream unmatedReadsF(unmatedReadsFile);
            while (unmatedReadsF.good()) {
                unmatedReadsF >> mopts->unmatedReads;
                uint64_t start = mopts->unmatedReads.find_last_of('/');
                uint64_t end = mopts->unmatedReads.find_last_of('.');
                mopts->outname += mopts->unmatedReads.substr(start + 1, end - start - 1);
                res &= alignReads(pfi, consoleLog, mopts);
            }
        } else {
            std::string readFile1 = mopts->read1;
            std::string readFile2 = mopts->read2;
            std::ifstream readF1(readFile1);
            std::ifstream readF2(readFile2);
            std::string outname = mopts->outname;
            while (readF1.good() and readF2.good()) {
                readF1 >> mopts->read1;
                readF2 >> mopts->read2;
                uint64_t start = mopts->read1.find_last_of('/');
                uint64_t end = mopts->read1.find_last_of('_');
                mopts->outname = outname + mopts->read1.substr(start + 1, end - start - 1);
//                std::cerr << mopts->read1 << "\n";
//                std::cerr << mopts->read2 << "\n";
//                std::cerr << mopts->outname << "\n";
                res &= alignReads(pfi, consoleLog, mopts);
            }
        }
    } else res &= alignReads(pfi, consoleLog, mopts);
    return res;
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
        success = alignReadsWrapper(pfi, consoleLog, &alnargs);
    } else if (indexType == "sparse") {
        PufferfishSparseIndex pfi(indexDir);
        success = alignReadsWrapper(pfi, consoleLog, &alnargs);
    } else if (indexType == "lossy") {
        PufferfishLossyIndex pfi(indexDir);
        success = alignReadsWrapper(pfi, consoleLog, &alnargs);
    }

    if (!success) {
        consoleLog->warn("Problem mapping.");
    }
    return 0;
}
