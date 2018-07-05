#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib> // std::abs
#include <cmath> // std::abs
#include <memory>
#include <unordered_set>

#include "cedar.hpp"

#include "clipp.h"
#include "EquivalenceClassBuilder.hpp"
#include "CLI/Timer.hpp"
#include "PufferFS.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"

// trim from both ends (in place)
static inline std::string& trim(std::string &s) {
    std::string chars = "\t\n\v\f\r ";
    s.erase(0, s.find_first_not_of(chars));
    s.erase(s.find_last_not_of(chars) + 1);
    return s;
}

struct CedarOpts {
    std::string taxonomyTree_filename;
    std::string refId2TaxId_filename;
    std::string mapperOutput_filename;
    std::string output_filename;
    std::string level = "species";
    size_t maxIter = 1000;
    double eps = 0.001;
    double filterThreshold = 0;
    bool flatAbund{false};
    bool requireConcordance{false};
    bool isPuffOut{false};
    bool isSAM{false};
    bool onlyUniq{false};
    bool onlyPerfect{false};
};

template<class ReaderType>
Cedar<ReaderType>::Cedar(std::string& taxonomyTree_filename, 
                 std::string& refId2TaxId_filename, 
                 std::string pruneLevelIn,
                 double filteringThresholdIn,
                 bool flatAbundIn,
                 std::shared_ptr<spdlog::logger> loggerIn) {
    flatAbund = flatAbundIn;
    logger = loggerIn;
    logger->info("Cedar: Construct ..");

    
    // map rank string values to enum values
    filteringThreshold = filteringThresholdIn;
    if (!flatAbund) {
        pruningLevel = TaxaNode::str2rank(pruneLevelIn);
        std::ifstream tfile;
        uint32_t id, pid;
        std::string rank, name;
        // load the reference id (name) to its corresponding taxonomy id
        tfile.open(refId2TaxId_filename);
        while(!tfile.eof()) {
            tfile >> name >> id;
            refId2taxId[name] = id;
        }
        tfile.close();

        // load the taxonomy child-parent tree and the rank of each node
        tfile.open(taxonomyTree_filename);
        std::string tmp, line;
        uint64_t lcntr{0};
        //std::cerr << "\n\nreading taxonomy tree\n\n";
        while (!tfile.eof()) {
            try {
                lcntr++;
                std::getline(tfile, line);
                uint64_t first = line.find_first_of('|');
                std::string fstr = line.substr(0, first);
                id = (uint32_t) std::stoul(trim(fstr));
                uint64_t second = line.find_first_of('|', first + 1);
                fstr = line.substr(first + 1, second - first - 1);
                pid = (uint32_t) std::stoul(trim(fstr));
                uint64_t third = line.find_first_of('|', second + 1);
                fstr = line.substr(second + 1, third - second - 1);
                rank = trim(fstr);
                if (rank.find("HS") != std::string::npos)
                    std::cerr << "here:" << id << " " << pid << " " << tmp << "\n";
                //if (lcntr == 1369188) {
                //std::cerr << lcntr << " " /*<< line << "\n"*/ << id << " " << pid << " " << rank << "\n";
                //}

                taxaNodeMap[id] = TaxaNode(id, pid, TaxaNode::str2rank(rank));
                if (taxaNodeMap[id].isRoot()) {
                    rootId = id;
                    logger->info("Root Id : {}", id);
                }
            } catch (const std::invalid_argument& ia) {
                std::cerr << "Invalid argument: " << ia.what() << '\n';
                continue;
            }
            
        }

        tfile.close();
        //std::cerr << "\n\ndone reading taxonomy tree\n\n";
    } 
}

template<class ReaderType>
void Cedar<ReaderType>::loadMappingInfo(std::string mapperOutput_filename,
                            bool requireConcordance,
                            bool onlyUniq,
                            bool onlyPerfect) {
    int32_t rangeFactorization{4};
    uint64_t totalReadCnt{0}, seqNotFound{0}, totalUnmappedReads{0}, tid, totalMultiMappedReads{0};
    logger->info("Cedar: Load Mapping File ..");
    logger->info("Mapping Output File: {}", mapperOutput_filename);
    mappings.load(mapperOutput_filename, logger);
    logger->info("is dataset paired end? {}", mappings.isMappingPaired());
    ReadInfo readInfo;
    TaxaNode* prevTaxa{nullptr};
    size_t conflicting{0};
    size_t discordantMappings{0};

    constexpr const bool getReadName = true;
    bool checkTruthExists = true;
    //size_t numCorrectHits{0};
    size_t numMapped{0};
    bool wasMapped = false;
    //size_t gid{0};
    //bool foundTruth = false;
    // std::ofstream covDist("uniqCovDist.txt");
    while(mappings.nextRead(readInfo, getReadName)) {
        if (checkTruthExists) {
            //foundTruth = false;
            wasMapped = false;
            /*
            size_t start = readInfo.rid.find_first_of('|');
            size_t finish = readInfo.rid.find_first_of('|', start+1);
            std::cout <<readInfo.rid<< "\t";
                        std::cout << start+1 << " " << finish-(start+1) << "\n";
                        //std::cout << std::stoul(readInfo.rid.substr(start+1, finish-(start+1))) << "\n";
            gid = std::stoul(readInfo.rid.substr(start+1, finish-(start+1)));
            */
        }
        totalReadCnt++;
        if (totalReadCnt % 1000000 == 0) {
          logger->info("Processed {} reads",totalReadCnt);
        }
        activeTaxa.clear();
        double readMappingsScoreSum = 0;
        std::vector<std::pair<uint64_t, double>> readPerStrainProbInst;
        readPerStrainProbInst.reserve(readInfo.cnt);
        bool isConflicting = true;
        if (readInfo.cnt != 0) {
            //std::cout << "read len " << readInfo.len << "\n";
            if (!wasMapped) { wasMapped = true; ++numMapped; }
            std::set<uint64_t> seen;
            prevTaxa = nullptr;
            for (auto& mapping : readInfo.mappings) {
                // first condition: Ignore those references that we don't have a
                // taxaId for 
                // second condition: Ignore repeated exactly identical
                // mappings (FIXME thing)
                // note: Fatemeh!! don't change the or and and again!!
                // if we are on flatAbund, we want to count for multiple occurence of a reference
                if( (flatAbund or 
                     (refId2taxId.find(mappings.refName(mapping.getId())) != refId2taxId.end())) 
                  and
                  activeTaxa.find(mapping.getId()) == activeTaxa.end()) {

                    if (prevTaxa != nullptr and
                        activeTaxa.find(mapping.getId()) == activeTaxa.end() and 
                        !prevTaxa->compareIntervals(mapping)) {
                        isConflicting = false;
                    }
                    prevTaxa = &mapping;
                    /*
                        if (activeTaxa.find(mapping.getId()) == activeTaxa.end() and
                            (flatAbund or 
                            refId2taxId.find(mappings.refName(mapping.getId())) != refId2taxId.end())) {
                    */
                    if (requireConcordance && mappings.isMappingPaired() && (!mapping.isConcordant() || mapping.isFw(ReadEnd::LEFT) == mapping.isFw(ReadEnd::RIGHT) ) ) {
                        discordantMappings++;
                        continue;
                    }

                    tid = flatAbund ? mapping.getId() : refId2taxId[mappings.refName(mapping.getId())];
                    seqToTaxMap[mapping.getId()] = tid;
                    activeTaxa.insert(mapping.getId());
                    if (cov.find(refId2taxId[mappings.refName(mapping.getId())]) != cov.end()) {
                        cov[refId2taxId[mappings.refName(mapping.getId())]] += mapping.getScore();
                    }
                    else {
                        cov[refId2taxId[mappings.refName(mapping.getId())]] = mapping.getScore();
                    }
                    readPerStrainProbInst.emplace_back(mapping.getId(), static_cast<double>( mapping.getScore()) / static_cast<double>(mappings.refLength(mapping.getId())) /* / static_cast<double>(tlen) */);
                    readMappingsScoreSum += readPerStrainProbInst.back().second;
                    
                    /* notMappedReadsFile  << readInfo.rid << "\t" 
                                        << mappings.refName(mapping.getId()) 
                                        << "\t" << mapping.getScore() << "\n"; */
                    /* if (!foundTruth and checkTruthExists) {
                      const auto& refName = mappings.refName(tid);
                      auto start = refName.find_first_of('|');
                      auto finish = refName.find_first_of('|', start+1);
                      auto refId = std::stoul(refName.substr(start+1, finish-(start+1)));
                      if (refId == gid) { ++numCorrectHits; foundTruth = true; }
                    } */

                }
                else {
                   // std::cerr << mappings.refName(mapping.getId()) << "\n";
                }
            } 
            if (activeTaxa.size() == 0) {
                seqNotFound++;
            } else if ( (!onlyUniq and !onlyPerfect) 
                        or (onlyUniq and activeTaxa.size() == 1)
                        or (onlyPerfect and activeTaxa.size() == 1 
                            and 
                            readInfo.mappings[0].getScore() >= readInfo.len) ) {
                if (!isConflicting) {conflicting++;}
                //if (!foundTruth) { notMappedReadsFile << readInfo.rid << "\n";}
                // bool isUnique = (readPerStrainProbInst.size() == 1);
                // it->first : strain id
                // it->second : prob of current read comming from this strain id
                for (auto it = readPerStrainProbInst.begin(); it != readPerStrainProbInst.end(); it++) {
                    it->second = it->second/readMappingsScoreSum; // normalize the probabilities for each read
                    // strain[it->first].first : read count for strainCnt
                    // strain[it->first].second : strain length
                    strain[it->first] += 1.0/static_cast<double>(readPerStrainProbInst.size());
                }
                // SAVE MEMORY, don't push this
                //readPerStrainProb.push_back(readPerStrainProbInst);

                // construct the range factorized eq class here 
                std::sort(readPerStrainProbInst.begin(), readPerStrainProbInst.end(), 
                [](std::pair<uint64_t, double>& a, std::pair<uint64_t, double>& b) {
                    return a.first < b.first;
                });
                std::vector<uint32_t> genomeIDs; genomeIDs.reserve(2*readPerStrainProbInst.size());
                std::vector<double> probs; probs.reserve(readPerStrainProbInst.size());
                for (auto it = readPerStrainProbInst.begin(); it != readPerStrainProbInst.end(); it++) {
                    genomeIDs.push_back(static_cast<uint32_t>(it->first));
                    probs.push_back(it->second);
                } 
                if (rangeFactorization > 0) {
                    int genomeSize = genomeIDs.size();
                    int rangeCount = std::sqrt(genomeSize) + rangeFactorization;
                    for (int i = 0; i < genomeSize; i++) {
                        int rangeNumber = probs[i] * rangeCount;
                        genomeIDs.push_back(static_cast<uint32_t>(rangeNumber));
                    }
                }
                readCnt++;
                TargetGroup tg(genomeIDs);
                eqb.addGroup(std::move(tg), probs);
            }
        } else {
                totalUnmappedReads++;
        }
    }
    //notMappedReadsFile.close();
    logger->info("Total # of conflicting reads reported: {}", conflicting); 
    logger->info("Total # of multimapped reads: {}", totalMultiMappedReads);
    if (requireConcordance)
        logger->info("Discarded {} discordant mappings.", discordantMappings);
}

template<class ReaderType>
bool Cedar<ReaderType>::basicEM(size_t maxIter, double eps) {
    eqb.finish();
    auto& eqvec = eqb.eqVec();
    int64_t maxSeqID{-1};
    for (auto& kv : strain) { 
        maxSeqID = (static_cast<int64_t>(kv.first) > maxSeqID) ? static_cast<int64_t>(kv.first) : maxSeqID;
    }

    std::vector<double> newStrainCnt(maxSeqID+1,0.0); 
    std::vector<double> strainCnt(maxSeqID+1);
    for (auto& kv : strain) { 
      strainCnt[kv.first] = kv.second;
    }

    logger->info("maxSeqID : {}", maxSeqID);
    logger->info("found : {} equivalence classes",eqvec.size()); 
    size_t totCount{0};
    for (auto& eqc : eqvec) {
        totCount += eqc.second.count;
    }
    logger->info("total staring count {}", totCount);
    logger->info("mapped readCnt {}", readCnt);

    size_t cntr = 0;
    bool converged = false;
    while (cntr++ < maxIter && !converged) {
        // M step
        // Find the best (most likely) count assignment
        for (auto& eqc : eqvec) {
            auto& tg = eqc.first;
            auto& v = eqc.second;
            auto csize = v.weights.size();
            std::vector<double> tmpReadProb(csize, 0.0);
            double denom{0.0};
            for (size_t readMappingCntr = 0; readMappingCntr < csize; ++readMappingCntr) {
                auto& tgt = tg.tgts[readMappingCntr];
               tmpReadProb[readMappingCntr] = v.weights[readMappingCntr] * strainCnt[tgt];// * (1.0/refLengths[tgt]);
               denom += tmpReadProb[readMappingCntr];
            }
            for (size_t readMappingCntr = 0; readMappingCntr < csize; ++readMappingCntr) {
                auto& tgt = tg.tgts[readMappingCntr];
                newStrainCnt[tgt] += v.count * (tmpReadProb[readMappingCntr] / denom);
            }
        }
    
        // E step
        // normalize strain probabilities using the denum : p(s) = (count(s)/total_read_cnt) 
        double readCntValidator = 0;
        converged = true;   
        double maxDiff={0.0};
        for (size_t i = 0; i < strainCnt.size(); ++i) {
            readCntValidator += newStrainCnt[i];
            auto adiff = std::abs(newStrainCnt[i] - strainCnt[i]);
            if ( adiff > eps) {
                converged = false;
            }
            maxDiff = (adiff > maxDiff) ? adiff : maxDiff;
            strainCnt[i] = newStrainCnt[i];
            newStrainCnt[i] = 0.0;
        }

        if (std::abs(readCntValidator - readCnt) > 10) {
            logger->error("Total read count changed during the EM process");
            logger->error("original: {}, current : {}, diff : {}", readCnt, 
                           readCntValidator, std::abs(readCntValidator - readCnt));
            std::exit(1);
        }
        
        if (cntr > 0 and cntr % 100 == 0) {
            logger->info("max diff : {}", maxDiff);
        }
    }
    logger->info( "iterator cnt: {}", cntr); 

    // We have done the EM in the space of sequence / reference IDs
    // but we need to output results in terms of taxa IDs.  Here, we 
    // will map our reference IDs back to taxa IDs, and put the resulting
    // computed abundances in the "strain" member variable.
    decltype(strain) outputMap;
    outputMap.reserve(strain.size());
    for (auto& kv : strain) { 
        outputMap[seqToTaxMap[kv.first]] += strainCnt[kv.first];
    }
    std::swap(strain, outputMap);
    
    return cntr < maxIter;
}

template<class ReaderType>
void Cedar<ReaderType>::serialize(std::string& output_filename) {
    logger->info("Write results into the file: {}", output_filename);
    logger->info("# of strains: {}", strain.size());
    std::ofstream ofile(output_filename);
    ofile << "taxaId\ttaxaRank\tcount\n";
    spp::sparse_hash_map<uint64_t, double> validTaxa;
    std::cerr << "strain size: " << strain.size() << "\n";
    for (auto& kv : strain) {
        if (taxaNodeMap.find(kv.first) != taxaNodeMap.end()) {
            TaxaNode *walker = &taxaNodeMap[kv.first];
            //std::cerr << "s" << walker->getId() << " ";
            while (!walker->isRoot() && walker->getRank() != pruningLevel) {
                //std::cerr << "p" << walker->getParentId() << " ";
                walker = &taxaNodeMap[walker->getParentId()];
                //std::cerr << walker->getId() << " ";
                if (walker->getId() > 18000000000000) std::exit(1);
            }
            //std::cerr << "\n";
            if (!walker->isRoot()) {
                if (validTaxa.find(walker->getId()) == validTaxa.end()) {
                    validTaxa[walker->getId()] = kv.second;
                } else {
                    validTaxa[walker->getId()] += kv.second;
                }
            }
        } else {
            std::cerr << "taxa not found: " << kv.first << "\n";
        }
    }
    for (auto& kv : validTaxa) { 
        ofile << kv.first << "\t" 
              << TaxaNode::rank2str(taxaNodeMap[kv.first].getRank()) 
              << "\t" << kv.second << "\n";
    }
    ofile.close();

    std::ofstream covOfile(output_filename + ".coverage");
    for (auto& kv: cov) {
        covOfile << kv.first << "\t" << kv.second << "\n";
    }
    covOfile.close();
    
}

template<class ReaderType>
void Cedar<ReaderType>::serializeFlat(std::string& output_filename) {
    logger->info("[FlatAbund]");
    // validate final count:
    uint64_t finalMappedReadCnt = 0;
    for (auto &kv : strain) {
        finalMappedReadCnt += kv.second;
    }
    logger->info("Before writing results in the file, total # of mapped reads is {}", finalMappedReadCnt);
    logger->info("Write results in the file: {}", output_filename);
    std::ofstream ofile(output_filename);
    ofile << "taxaId\ttaxaRank\tcount\n";
    std::cerr << "NUMREFS: " << mappings.numRefs() << "\n";
    for (uint32_t i = 0; i < mappings.numRefs(); ++i) { 
        //for (auto& kv : strain) {
        auto it = strain.find(i);
        double abund = 0.0;
        if (it != strain.end()){
        abund = it->second;
        }
        ofile << mappings.refName(i) << "\t" 
            << "flat" 
            << "\t" << abund << "\n";
    }
    ofile.close();
}

template class Cedar<PuffMappingReader>;
template class Cedar<SAMReader>;

/**
 * "How to run" example:
 * make Pufferfish!
 * In the Pufferfish build directory run the following command:
 * /usr/bin/time src/cedar 
 * -t /mnt/scratch2/avi/meta-map/kraken/KrakenDB/taxonomy/nodes.dmp  
 * -s /mnt/scratch2/avi/meta-map/kraken/KrakenDB/seqid2taxid.map 
 * -m /mnt/scratch2/avi/meta-map/kraken/puff/dmps/HC1.dmp 
 * -o HC1.out 
 * -l genus (optional)
 * -f 0.8 (optional)
 **/
int main(int argc, char* argv[]) {
  (void)argc;
  using namespace clipp;
  CedarOpts kopts;
  bool showHelp{false};

  auto checkLevel = [](const char* lin) -> void {
    std::string l(lin);
    std::unordered_set<std::string> valid{"species", "genus", "family", "order", "class", "phylum"};
    if (valid.find(l) == valid.end()) {
      std::string s = "The level " + l + " is not valid.";
      throw std::range_error(s);
    }
  };

  auto cli = (
            (required("--flat").set(kopts.flatAbund, true) % "estimate flat abundance (i.e. there is no taxonomy given)"
            | (
              required("--taxtree", "-t") & value("taxtree", kopts.taxonomyTree_filename) % "path to the taxonomy tree file",
              required("--seq2taxa", "-s") & value("seq2taxa", kopts.refId2TaxId_filename) % "path to the refId 2 taxId file "
            )),
              ( (required("--puffMapperOut", "-p").set(kopts.isPuffOut, true) & value("mapout", kopts.mapperOutput_filename) % "path to the pufferfish mapper output file")
              |
              (required("--sam").set(kopts.isSAM, true) & value("mapout", kopts.mapperOutput_filename) % "path to the SAM file")
              ),
              required("--output", "-o") & value("output", kopts.output_filename) % "path to the output file to write results",
              option("--maxIter", "-i") & value("iter", kopts.maxIter) % "maximum number of EM iteratons (default : 100)",
              option("--eps", "-e") & value("eps", kopts.eps) % "epsilon for EM convergence (default : 0.001)",
              option("--level", "-l") & value("level", kopts.level).call(checkLevel) % "choose between (species, genus, family, order, class, phylum). (default : species)",
              option("--filter", "-f") & value("filter", kopts.filterThreshold) % "choose the threshold [0,1] below which to filter out mappings (default : no filter)",
              option("--noDiscordant").set(kopts.requireConcordance, true) % "ignore orphans for paired end reads",
              option("--unique").set(kopts.onlyUniq, true) % "report abundance based on unique reads",
              option("--perfect").set(kopts.onlyPerfect, true) % "report abundance based on perfect reads (unique and with complete coverage)",
              option("--help", "-h").set(showHelp, true) % "show help",
              option("-v", "--version").call([]{std::cout << "version 0.1.0\n\n";}).doc("show version")
              );
   
  //Multithreaded console logger(with color support)
  auto console = spdlog::stderr_color_mt("console");

  decltype(parse(argc, argv, cli)) res;
  try {
    res = parse(argc, argv, cli);
    if (showHelp){
      std::cout << make_man_page(cli, "cedar"); 
      return 0;
    }
  } catch (std::exception& e) {
    std::cout << "\n\nparsing command line failed with exception: " << e.what() << "\n";
    std::cout << "\n\n";
    std::cout << make_man_page(cli, "cedar");
    return 1;
  }

  if(res) {
    if (kopts.isSAM) {
        Cedar<SAMReader> cedar(kopts.taxonomyTree_filename, kopts.refId2TaxId_filename, kopts.level, kopts.filterThreshold, kopts.flatAbund, console);
        cedar.run(kopts.mapperOutput_filename, 
                  kopts.requireConcordance,
                  kopts.maxIter, 
                  kopts.eps,
                  kopts.output_filename,
                  kopts.onlyUniq,
                  kopts.onlyPerfect);
        //std::cout << "inside the scope\n";
    }
    else {
        Cedar<PuffMappingReader> cedar(kopts.taxonomyTree_filename, kopts.refId2TaxId_filename, kopts.level, kopts.filterThreshold, kopts.flatAbund, console);
        cedar.run(kopts.mapperOutput_filename, 
            kopts.requireConcordance,
            kopts.maxIter, 
            kopts.eps,
            kopts.output_filename,
            kopts.onlyUniq,
            kopts.onlyPerfect);
    }
    //std::cout << "We even got here!!\n";
    return 0;
  } else {
    std::cout << usage_lines(cli, "cedar") << '\n';
    return 1;
  }
}

