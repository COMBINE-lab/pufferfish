#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib> // std::abs
#include <cmath> // std::abs
#include "CLI/CLI.hpp"
#include "cedar.hpp"

#define LEFT true
#define RIGHT true
#define SCALE_FACTOR 1000000

struct CedarOpts {
    std::string taxonomyTree_filename;
    std::string refId2TaxId_filename;
    std::string mapperOutput_filename;
    std::string output_filename;
    std::string level = "species";
    double filterThreshold = 0;
    double eps = 0.001;
    size_t maxIter = 100;
};

Cedar::Cedar(std::string& taxonomyTree_filename, 
                 std::string& refId2TaxId_filename, 
                 std::string pruneLevelIn,
                 double filteringThresholdIn) {

    std::cerr << "KrakMap: Construct ..\n";
    // map rank string values to enum values
    filteringThreshold = filteringThresholdIn;
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
    std::string tmp;
    while (!tfile.eof()) {
        tfile >> id >> tmp >> pid >> tmp >> rank >> tmp;
        size_t i = 0;
        while (i < tmp.size()-1 && isspace(tmp[i]))
            i++;
        if (tmp[i] != '|') {
            rank += " " + tmp;
        }
        taxaNodeMap[id] = TaxaNode(id, pid, TaxaNode::str2rank(rank));
        if (taxaNodeMap[id].isRoot()) {
            rootId = id;
            std::cerr << "Root Id : " << id << "\n";
        }
        std::getline(tfile, tmp);
        
    }

    tfile.close();  
}

bool Cedar::readHeader(std::ifstream& mfile) {
    std::string tmp, readType;
    mfile >> tmp >> readType;
    if (tmp != "#")
        return false;
    if (readType == "LT:S") 
        isPaired = false;
    else if (readType == "LT:P")
        isPaired = true;
    else
        return false;
    std::getline(mfile, tmp);
    return true;
}

void Cedar::loadMappingInfo(std::string mapperOutput_filename) {
    std::string rid, tname, tmp;// read id, taxa name, temp
    uint64_t lcnt, rcnt, tid, puff_tid, tlen, ibeg, ilen;
    std::cerr << "Cedar: Load Mapping File ..\n";
    std::cerr << "\tMapping Output File: " << mapperOutput_filename << "\n";
    std::ifstream mfile(mapperOutput_filename);
    uint64_t rlen, mcnt; // taxa id, read mapping count, # of interals, interval start, interval length
    uint64_t totalReadCnt = 0, totalUnmappedReads = 0, seqNotFound = 0;
    if (!readHeader(mfile)) {
        std::cerr << "ERROR: Invalid header for mapping output file.\n";
        std::exit(1);
    }
    std::cout<< "is dataset paired end? " << isPaired << "\n";
    while (!mfile.eof()) {
        mfile >> rid >> mcnt;
        totalReadCnt++;
        activeTaxa.clear();
        float readMappingsScoreSum = 0;
        std::vector<std::pair<uint64_t, float>> readPerStrainProbInst;
        readPerStrainProbInst.reserve(10);
        //std::cout << "r" << rid << " " << mcnt << "\n";
        if (mcnt != 0) {
            if (isPaired) {
                uint64_t rllen, rrlen;
                mfile >> rllen >> rrlen;
                rlen = rllen + rrlen;
            }
            else {
                mfile >> rlen;
            }
            std::set<uint64_t> seen;
            for (size_t mappingCntr = 0; mappingCntr < mcnt; mappingCntr++) {
                mfile >> puff_tid >> tname >> tlen; //txp_id, txp_name, txp_len
                // first condition: Ignore those references that we don't have a taxaId for
                // secon condition: Ignore repeated exactly identical mappings (FIXME thing)
                if (refId2taxId.find(tname) != refId2taxId.end() &&
                    activeTaxa.find(refId2taxId[tname]) == activeTaxa.end()) { 
                    tid = refId2taxId[tname];
                    activeTaxa.insert(tid);
                    
                    // fetch the taxon from the map
                    TaxaNode taxaPtr;
                    mfile >> lcnt;
                    if (isPaired)
                        mfile >> rcnt;
                    for (size_t i = 0; i < lcnt; ++i) {
                        mfile >> ibeg >> ilen;
                        taxaPtr.addInterval(ibeg, ilen, LEFT);
                    }
                    if (isPaired)
                        for (size_t i = 0; i < rcnt; ++i) {
                            mfile >> ibeg >> ilen;
                            taxaPtr.addInterval(ibeg, ilen, RIGHT);
                        }
                    taxaPtr.cleanIntervals(LEFT);
                    taxaPtr.cleanIntervals(RIGHT);
                    taxaPtr.updateScore();
                    readPerStrainProbInst.emplace_back(tid, static_cast<float>(taxaPtr.getScore())/static_cast<float>(tlen));
                    readMappingsScoreSum += readPerStrainProbInst.back().second;
                }
                else { // otherwise we have to read till the end of the line and throw it away
                    std::getline(mfile, tmp);
                }
            } 
            if (activeTaxa.size() == 0) {
                seqNotFound++;
            }
            else {
                readCnt++;
                // it->first : strain id
                // it->second : prob of current read comming from this strain id
                for (auto it = readPerStrainProbInst.begin(); it != readPerStrainProbInst.end(); it++) {
                    it->second = it->second/readMappingsScoreSum; // normalize the probabilities for each read
                    // strain[it->first].first : read count for strainCnt
                    // strain[it->first].second : strain length
                    if (strain.find(it->first) == strain.end()) {
                        strain[it->first] = 1.0/static_cast<float>(readPerStrainProbInst.size());
                    }
                    else {
                        strain[it->first] += 1.0/static_cast<float>(readPerStrainProbInst.size());
                    }
                }
                readPerStrainProb.push_back(readPerStrainProbInst);
            }
        } else {
            totalUnmappedReads++;
            std::getline(mfile, tmp);
        }
    }  
}

bool Cedar::basicEM(size_t maxIter, double eps) {
    size_t cntr = 0;
    bool converged = false;
    while (cntr++ < maxIter && !converged) {
        spp::sparse_hash_map<uint64_t, float> newStrainCnt;
        // M step
        // Find the best (most likely) count assignment
        for (auto readIt = readPerStrainProb.begin(); readIt != readPerStrainProb.end(); readIt++) {
            float denum = 0;
            std::vector<float> tmpReadProb(readIt->size());
            size_t readMappingCntr = 0;
            for (auto strainIt = readIt->begin(); strainIt != readIt->end(); strainIt++, readMappingCntr++) {
                tmpReadProb[readMappingCntr] = strainIt->second*strain[strainIt->first]/readCnt;
                denum += tmpReadProb[readMappingCntr];
            }
            readMappingCntr = 0;
            for (auto strainIt = readIt->begin(); strainIt != readIt->end(); strainIt++, readMappingCntr++) {
                if (newStrainCnt.find(strainIt->first) == newStrainCnt.end())
                    newStrainCnt[strainIt->first] = tmpReadProb[readMappingCntr]/denum;
                else
                    newStrainCnt[strainIt->first] += tmpReadProb[readMappingCntr]/denum;
                //std::cout << strain[strainIt->first] << "-" 
                 //         << strainIt->second << "-" << readCnt << "-"
                   //       << strain[strainIt->first]*strainIt->second/readCnt << " ";
            }
        }

        // E step
        // normalize strain probabilities using the denum : p(s) = (count(s)/total_read_cnt) 
        float readCntValidator = 0;
        converged = true;   
        for (auto it = strain.begin(); it != strain.end(); it++) {
            readCntValidator += it->second;
            if (std::abs(newStrainCnt[it->first] - it->second) > eps) {
                converged = false;
            }
            it->second = newStrainCnt[it->first];
        }
        if (std::abs(readCntValidator - readCnt) > 10) {
            std::cerr << "ERROR: Total read count changed during the EM process\n";
            std::cerr << "original: " << readCnt << " current: " << readCntValidator << " diff: " 
                      << std::abs(readCntValidator - readCnt) << "\n";
            std::exit(1);
        }
    }
    std::cout << "iterator cnt: " << cntr << "\n";
    return cntr < maxIter;
}

void Cedar::serialize(std::string& output_filename) {
    std::cerr << "Write results in the file:\n" << output_filename << "\n";
    std::ofstream ofile(output_filename);
    ofile << "taxaId\ttaxaRank\tcount\n";
    for (auto& kv : strain) {
        ofile << kv.first << "\t" << TaxaNode::rank2str(taxaNodeMap[kv.first].getRank()) << "\t" << kv.second << "\n";
    }
    ofile.close();
}

/**
 * "How to run" example:
 * make Pufferfish!
 * In the Pufferfish build directory run the following command:
 * /usr/bin/time src/krakmap 
 * -t /mnt/scratch2/avi/meta-map/kraken/KrakenDB/taxonomy/nodes.dmp  
 * -s /mnt/scratch2/avi/meta-map/kraken/KrakenDB/seqid2taxid.map 
 * -m /mnt/scratch2/avi/meta-map/kraken/puff/dmps/HC1.dmp 
 * -o HC1.out 
 * -l genus (optional)
 * -f 0.8 (optional)
 **/
int main(int argc, char* argv[]) {
  (void)argc;

  CedarOpts kopts;
  CLI::App app{"krakMap : Taxonomy identification based on the output of Pufferfish mapper through the same process as Kraken."};
  app.add_option("-t,--taxtree", kopts.taxonomyTree_filename,
                 "path to the taxonomy tree file")
      ->required();
  app.add_option("-s,--seq2taxa", kopts.refId2TaxId_filename, "path to the refId 2 taxaId file")
      ->required();
  app.add_option("-m,--mapperout", kopts.mapperOutput_filename, "path to the pufferfish mapper output file")
      ->required();
  app.add_option("-o,--output", kopts.output_filename, "path to the output file to write results")
      ->required();
  app.add_option("-i,--maxIter", kopts.maxIter, "maximum allowed number of iterations for EM")
      ->required(false);
  app.add_option("-e,--eps", kopts.eps, "epsilon for EM convergence condition")
      ->required(false);
  app.add_option("-l,--level", kopts.level, "choose between (species, genus, family, order, class, phylum). Default:species")
      ->required(false);
  app.add_option("-f,--filter", kopts.filterThreshold, "choose the threshold (0-1) to filter out mappings with a score below that. Default: no filter")
      ->required(false);

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }
  Cedar cedar(kopts.taxonomyTree_filename, kopts.refId2TaxId_filename, kopts.level, kopts.filterThreshold);
  cedar.loadMappingInfo(kopts.mapperOutput_filename);
  cedar.basicEM(kopts.maxIter, kopts.eps);
  cedar.serialize(kopts.output_filename);
  return 0;
}
