#include <fstream>
#include <iostream>
#include <string>
#include <algorithm> // std::sort
#include "CLI/CLI.hpp"
#include "krakMap.h"

struct krakMapOpts {
    std::string taxonomyTree_filename;
    std::string refId2TaxId_filename;
    std::string mapperOutput_filename;
};

void TaxaNode::addInterval(uint64_t begin, uint64_t len) {
    intervals.emplace_back(begin, begin+len);
}
/**
 * Sorts intervals
 * Merge intervals if possible
 * Calculates score
**/
void TaxaNode::cleanIntervalsAndCalcScore() {
    // if we were writing intervals wrt read position we wouldn't need this part
    std::sort(intervals.begin(), intervals.end(),
    [](Interval& i1, Interval& i2){
        return i1.begin != i2.begin?i1.begin < i2.begin:i1.end < i2.end;
    });

    // merge intervals if necessary
    // calculate score / coverage !! this whole process is repetition of mapping coverage calc!!
    for (auto it=intervals.begin(); it != intervals.end();) {
        // start from next item and merge (and erase) as much as possible
        bool merged = true;
        auto next = it+1;
        while (next != intervals.end() && merged) {
            // if they overlap, merge them
            if (it->end >= next->begin) {
                if (it->end < next->end)
                    it->end = next->end;
                intervals.erase(next); // erase next after merging with it
                next = it+1; // it pointer and its next are always valid
            }
            else {// no overlap anymore and no merging. Update score and leave the loop
                score += (it->end-it->begin);
                merged = false;
            }
        }
        // there is nothing to merge with the current it, so increase it
        it++;
    }
}

bool TaxaNode::addChild(TaxaNode* child) { 
    children.push_back(child);
    return true;
}


KrakMap::KrakMap(std::string& taxonomyTree_filename, std::string& refId2TaxId_filename) {
    // map rank string values to enum values
    initializeRanks();

    std::ifstream tfile;
    uint32_t id, pid;
    std::string rank, name;
    // load the reference id (name) to its corresponding taxonomy id
    tfile.open(refId2TaxId_filename);
    while(!tfile.eof()) {
        tfile >> name >> id;
        /* if (refId2taxId.find(name) != refId2taxId.end()) {
            std::cerr << "MULTIPLE REF_NAMES: " << name << " " << id << " --> ";
            std::cerr << refId2taxId[name] << " " << refId2taxId.size() << "\n";
        } */
        refId2taxId[name] = id;
    }
    tfile.close();

    // load the taxonomy child-parent tree and the rank of each node
    tfile.open(taxonomyTree_filename);
    std::string tmp;
    while (!tfile.eof()) {
        tfile >> id >> tmp >> pid >> tmp >> rank;
        // std::cout << id << tmp << pid << tmp << rank << "\n";
        /* if (taxonomyC2P.find(id) != taxonomyC2P.end()) {
            std::cerr << "MULTIPLE NODE IDs: " << id << " " << rank << " --> ";
            std::cerr << taxonomyC2P[id].first << " " << taxonomyC2P.size() << "\n";
        } */
        taxaNodeMap[id] = TaxaNode(id, str2Rank[rank], pid);
        std::getline(tfile, tmp);
    }
    tfile.close();

    // add pointers to children for each taxa node
    for (auto& kv : taxaNodeMap) {
        TaxaNode& taxaNode = kv.second;
        taxaNodeMap[taxaNode.getParentId()].addChild(&taxaNode);
    }     
}

bool KrakMap::classify(std::string& mapperOutput_filename) {
    std::ifstream mfile(mapperOutput_filename);
    std::string rid, tname; // read id, taxa name temp
    uint64_t rlen, tid, mcnt, icnt, ibeg, ilen; // taxa id, read mapping count, # of interals, interval start, interval length
    while (!mfile.eof()) {
        mfile >> rid >> mcnt >> rlen;
        // reset everything we've done for previous read
        // leaves.clear();
        leaves = std::priority_queue <TaxaNode*, std::vector<TaxaNode*>, TaxaNodePtrComparator>();
        resetIntervals();

        // construct intervals for leaves
        for (size_t mappingCntr = 0; mappingCntr < mcnt; mappingCntr++) {
            mfile >> tname >> icnt;
            tid = refId2taxId[tname];
            activeTaxa.insert(tid);
            // fetch the taxum from the map
            TaxaNode* taxaPtr = &taxaNodeMap[tid];
            for (size_t i = 0; i < icnt; ++i) {
                mfile >> ibeg >> ilen;
                taxaPtr->addInterval(ibeg, ilen);
            }
            taxaPtr->cleanIntervalsAndCalcScore();
            leaves.push(taxaPtr);
        }

        // calculate score and intervals for all internal nodes

        // find best path for this read

        // update counter of the corresponding level in the taxonomy tree
    }
    return true;
}

void KrakMap::resetIntervals() {
    for (auto & activeTaxum : activeTaxa) {
        taxaNodeMap[activeTaxum].clearIntervals();
    }
    activeTaxa.clear();
}


int main(int argc, char* argv[]) {
  (void)argc;

  krakMapOpts kopts;
  CLI::App app{"krakMap : Taxonomy identification based on the output of Pufferfish mapper through the same process as Kraken."};
  app.add_option("-t,--taxtree", kopts.taxonomyTree_filename,
                 "path to the taxonomy tree file")
      ->required();
  app.add_option("-s,--seq2taxa", kopts.refId2TaxId_filename, "path to the refId 2 taxaId file")
      ->required();
  app.add_option("-m,--mapperout", kopts.mapperOutput_filename, "path to the pufferfish mapper output file")
      ->required();

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }
  KrakMap krakMap(kopts.taxonomyTree_filename, kopts.refId2TaxId_filename);
  return 0;
}
