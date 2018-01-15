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
    std::string output_filename;
    std::string level;
};

void TaxaNode::addInterval(uint64_t begin, uint64_t len) {
    intervals.emplace_back(begin, begin+len);
}

void TaxaNode::updateIntervals(TaxaNode* child) {

    // reduce the unready child counter
    notIncorporatedChildrenCounter--;

    // merge two sorted interval lists into parent
    // update parent score
    score = 0;
    auto& childIntervals = child->getIntervals();
    std::vector<Interval> parentIntervals(intervals.size());
    std::copy(intervals.begin(), intervals.end(), parentIntervals.begin());
    
    intervals.clear();
    intervals.reserve(parentIntervals.size()+childIntervals.size());
    
    std::vector<Interval>::iterator pit = parentIntervals.begin();
    std::vector<Interval>::iterator cit = childIntervals.begin();
    std::vector<Interval>::iterator fit = intervals.begin();
    std::cerr << parentIntervals.size() << " " << childIntervals.size() << " ";

    // add the smallest interval as the first interval
    if ( (pit == parentIntervals.end() && cit != parentIntervals.end()) || cit->begin < pit->begin) {
        intervals.emplace_back(cit->begin, cit->end);
        cit++;
    }
    else if (pit != parentIntervals.end()) {
        intervals.emplace_back(pit->begin, pit->end);
        pit++;
    }
    else {
        std::cerr << "ERROR!! Both parent an child intervals were empty.\n";
    }
    std::vector<Interval>::iterator cur;
    while (pit != parentIntervals.end() || cit != childIntervals.end()) {
        // find the smallest interval between the heads of the two lists
        if (pit == parentIntervals.end() || cit->begin < pit->begin) {
            cur = cit;
            cit++;
        }
        else {
            cur = pit;
            pit++;  
        }
        // merge the new interval
        // Note: since both lists are sorted
        // the new interval's begin is always >= the last inserted interval's
        if (fit->end >= cur->begin) { // if the new interval has an overlap with the last inserted one
            fit->end = cur->end; // merge them
        } else { // insert the interval as a separate one and move fit forward
            intervals.emplace_back(cur->begin, cur->end);
            fit++;
        }
    }
    for (auto& it : intervals) {
        score += it.end - it.begin;
    }
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
    // if this child was not already added
    if (activeChildren.find(child->getId()) == activeChildren.end()) {
        activeChildren.insert(child->getId()); // add it
        notIncorporatedChildrenCounter++; // increase not incorprated children counter
        return true;
    }
    return false;
}

void TaxaNode::reset() {
     intervals.clear();
     activeChildren.clear();
     notIncorporatedChildrenCounter = 0;
}


KrakMap::KrakMap(std::string& taxonomyTree_filename, 
                 std::string& refId2TaxId_filename, 
                 std::string pruneLevelIn="species") {

    std::cerr << "KrakMap: Construct ..\n";
    // map rank string values to enum values
    initializeRanks();
    pruningLevel = str2Rank[pruneLevelIn];
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
        if (taxaNodeMap[id].isRoot())
            root = &taxaNodeMap[id];
        std::getline(tfile, tmp);
    }
    tfile.close();  
}

bool KrakMap::classify(std::string& mapperOutput_filename) {
    std::cerr << "KrakMap: Classify ..\n";
    std::cerr << "\tMapping Output File: " << mapperOutput_filename << "\n";
    std::ifstream mfile(mapperOutput_filename);
    std::string rid, tname; // read id, taxa name temp
    uint64_t rlen, tid, mcnt, icnt, ibeg, ilen; // taxa id, read mapping count, # of interals, interval start, interval length
    while (!mfile.eof()) {
        mfile >> rid >> mcnt >> rlen;
        // reset everything we've done for previous read
        clearReadSubTree();
        // construct intervals for leaves
        // std::cerr << "# of hits: " << mcnt << "...\n";
        for (size_t mappingCntr = 0; mappingCntr < mcnt; mappingCntr++) {
            mfile >> tname >> icnt;
            tid = refId2taxId[tname];
            activeTaxa.insert(tid);
            // fetch the taxum from the map
            TaxaNode* taxaPtr = &taxaNodeMap[tid];
            walk2theRoot(taxaPtr);
            for (size_t i = 0; i < icnt; ++i) {
                mfile >> ibeg >> ilen;
                taxaPtr->addInterval(ibeg, ilen);
            }
            taxaPtr->cleanIntervalsAndCalcScore();
            hits.push_front(taxaPtr);
        }

        // propagate score and intervals to all internal nodes
        std::cerr << "Update intervals and scores of internal nodes ..\n";
        propagateInfo();
        // find best path for this read
        std::cerr << "Assign Read ..\n";
        assignRead();
    }
    return true;
}

void KrakMap::walk2theRoot(TaxaNode* child) {
    while (!child->isRoot()) {
        TaxaNode* parent = &taxaNodeMap[child->getParentId()]; // fetch parent
        activeTaxa.insert(parent->getId());
        parent->addChild(child); // add current node as its child
        child = parent; // walk to parent --> parent is the new child
    }
}

//TODO don't like it 
// It's not as beautiful as when the root is ripe, we are done, but it does the same thing :/
// will we pass being ripe and become overripe? (negative # of incorporated children .. )
// we shouldn't .. but it just doesn't happen because it doesn't .. :/
void KrakMap::propagateInfo() {
    while (!hits.empty()) {
        TaxaNode* taxaPtr = hits.back();
        hits.pop_back();
        // if the hit itself is not ripe, no need to push it back to the queue
        // when it's the turn for one of the other hits, it'll get ripe and updated
        while (!taxaPtr->isRoot() && taxaPtr->isRipe()) {
            TaxaNode* parentPtr = &taxaNodeMap[taxaPtr->getParentId()];
            parentPtr->updateIntervals(taxaPtr);
            taxaPtr = parentPtr;
        }
    }
}

void KrakMap::assignRead() {
    TaxaNode* walker = root;
    
    while (walker->getRank() != pruningLevel) {
        uint64_t maxScore=0, maxId = -1, maxCntr;
        for (auto childId : walker->getActiveChildren()) {
            TaxaNode& child = taxaNodeMap[childId];
            if (child.getScore() == maxScore) {
                maxCntr++;
            }
            else if (child.getScore() > maxScore) {
                maxId = childId;
                maxScore = child.getScore();
                maxCntr = 1;
            }
        }
        if (maxCntr != 1) { // zero --> no children (it's a leaf) || > 1 --> more than one child with max score
            break;
        }
        walker = &taxaNodeMap[maxId];
    }
    if (mappedReadCntr.find(walker->getId()) == mappedReadCntr.end())
        mappedReadCntr[walker->getId()] = std::make_pair(1, walker->getRank());
    else
        mappedReadCntr[walker->getId()].first += 1;
}

void KrakMap::serialize(std::string& output_filename) {
    std::ofstream ofile(output_filename);
    ofile << "taxaId\ttaxaRank\tcount\n";
    for (auto& kv : mappedReadCntr) {
        ofile << kv.first << "\t" << rankToStr(kv.second.second) << "\t" << kv.second.first << "\n";
    }
    ofile.close();
}
    
void KrakMap::clearReadSubTree() {
    for (auto & activeTaxum : activeTaxa) {
        taxaNodeMap[activeTaxum].reset();
    }
    activeTaxa.clear();
    hits.clear();
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
  app.add_option("-o,--output", kopts.output_filename, "path to the output file to write results")
      ->required();
  app.add_option("-l,--level", kopts.level, "choose between (species, genus, family, order, class, phylum)")
      ->required();

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }
  KrakMap krakMap(kopts.taxonomyTree_filename, kopts.refId2TaxId_filename, kopts.level);
  krakMap.classify(kopts.mapperOutput_filename);
  krakMap.serialize(kopts.output_filename);
  return 0;
}
