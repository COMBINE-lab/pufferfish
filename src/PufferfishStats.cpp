#include "FastxParser.hpp"
#include <cmath>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>
#include <bitset>
#include <cereal/archives/binary.hpp>
#include <set>

#include "ProgOpts.hpp"
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferFS.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "jellyfish/mer_dna.hpp"
#include "spdlog/spdlog.h"

#include "PufferfishIndex.hpp"
#include "Util.hpp"
#include "PufferfishBinaryGFAReader.hpp"

struct PrefixTree {
    std::string txps;
    uint32_t exactCnt;
    uint32_t prefixCnt;
    uint32_t prefixTotalCnt;
    uint32_t numTxps;

    PrefixTree(const std::string &txps, uint32_t exactCnt, uint32_t prefixCnt, uint32_t prefixTotalCnt)
            : txps(txps), exactCnt(exactCnt), prefixCnt(prefixCnt), prefixTotalCnt(prefixTotalCnt) {
        numTxps = std::count(txps.begin(), txps.end(), '-') + 1;
    }
};

struct NodeInfo {
    uint32_t cid;
    std::set<uint32_t> refid;

    NodeInfo(uint32_t cid) : cid(cid) {}
};

void doFindMotifs(std::string &indexDir) {
    PufferfishIndex pfi(indexDir);
    compact::vector<uint64_t, 1> visited(pfi.numContigs()-1);
    visited.clear_mem();
    std::cerr << "\nTotal contigs: " << visited.size() << "\n";
    uint64_t prevV{0};
    std::vector<char> chars{'A', 'C', 'G', 'T'};
    for (uint64_t v = 0; v < visited.size(); v++) {
        if (not visited[v]) {
            if (v - prevV > 1000) {
                std::cerr << "\r" << v;
                prevV = v;
            }
            visited[v] = 1;
            std::vector<CanonicalKmer> cks;
            cks.push_back(pfi.getStartKmer(v));
            auto end = pfi.getEndKmer(v);
            if (cks[0].getCanonicalWord() != end.getCanonicalWord()) {
                cks.push_back(end);
            }
            auto chits = pfi.getRefPos(cks[0]);
            NodeInfo self(chits.contigIdx_); // should be equal to v
            for (auto &r: chits.refRange) {
                self.refid.insert(r.transcript_id());
            }
            uint32_t intersect{0}, x{0}, y{0}, z{0};
            bool found = false;
            for (auto &ck : cks) {
                std::vector<NodeInfo> oneEndNeighbors;
                uint32_t cntr = 0;
                while (not found and cntr < 2) {
                    for (auto c: chars) {
                        auto cck(ck);
                        if (cntr) {
                            cck.shiftFw(c);
                        } else {
                            cck.shiftBw(c);
                        }
                        chits = pfi.getRefPos(cck);
                        if (not chits.empty() and not visited[chits.contigIdx_]) {
                            oneEndNeighbors.emplace_back(chits.contigIdx_);
                            for (auto &r: chits.refRange) {
                                oneEndNeighbors.back().refid.insert(r.transcript_id());
                            }
                        }

                        if (oneEndNeighbors.size() >= 2) {
                            for (auto r : self.refid) {
                                if (oneEndNeighbors[0].refid.find(r) != oneEndNeighbors[0].refid.end() and
                                    oneEndNeighbors[1].refid.find(r) != oneEndNeighbors[1].refid.end()) {
                                    intersect++;
                                }
                            }
                            x = self.refid.size() - intersect;
                            y = oneEndNeighbors[0].refid.size() - intersect;
                            z = oneEndNeighbors[1].refid.size() - intersect;
                            found = true;
                            visited[oneEndNeighbors[0].cid] = 1;
                            visited[oneEndNeighbors[1].cid] = 1;
                            break;
                        }
                    }
                    cntr++;
                }
            }
            if (found) {
                std::cout << intersect << "\t" << x << "\t" << y << "\t" << z << "\n";
            } else {
                std::cout << self.refid.size() << "\t-1\t-1\t-1\n";
            }
        }
    }
    std::cerr << "\n";
}

void doCtabStats(std::string &indexDir) {
    std::vector<std::string> refNames_;
    std::vector<uint32_t> refExt_;
    std::ifstream contigTableStream(indexDir + "/" + pufferfish::util::CTABLE);
    cereal::BinaryInputArchive contigTableArchive(contigTableStream);
    contigTableArchive(refNames_);
    contigTableArchive(refExt_);
    refNames_.clear();
    refNames_.shrink_to_fit();
    refExt_.clear();
    refExt_.shrink_to_fit();

    std::vector<pufferfish::util::Position> rawContigTable;
    contigTableArchive(rawContigTable);
    contigTableStream.close();
    std::string pfile = indexDir + "/" + pufferfish::util::CONTIG_OFFSETS;
    compact::vector<uint64_t> rawOffsets{16};
    rawOffsets.deserialize(pfile, false);
    std::vector<uint64_t> offsets(rawOffsets.size());
    for (size_t i = 0; i < rawOffsets.size(); ++i) {
      offsets[i] = rawOffsets[i];
    }
    pufferfish::ContigTable contigTable;
    contigTable.build(std::move(rawContigTable), offsets);
    std::cerr << "contigTable size: " << contigTable.numPositions()
              << " numContigs: " << contigTable.numContigs() << "\n";
    std::vector<std::string> txps;
    for (uint64_t i = 0; i < contigTable.numContigs(); i++) {
        auto range = contigTable.contigRange(i);
        auto it = range.begin();
        auto txp = it->transcript_id();
        std::string txpstr = std::to_string(txp);
        ++it;
        while (it != range.end()) {
            if (it->transcript_id() != txp) {
                txp = it->transcript_id();
                txpstr += ("-" + std::to_string(txp));
            }
            ++it;
        }
        txps.push_back(txpstr);
    }
    std::sort(txps.begin(), txps.end());
    auto prevt = txps[0];
    std::vector<PrefixTree> tree;
    tree.emplace_back(prevt, 0, 0, 0);
    for (auto &t : txps) {
        if (t == prevt) {
            tree.back().exactCnt++;
        } else {
            prevt = t;
            tree.emplace_back(t, 1, 0, 0);
        }
    }
    std::cerr << "original cnt= " << txps.size() << " distinct cnt=" << tree.size() << "\n";
    for (uint64_t i = 0; i < tree.size(); i++) {
        std::string prefix = tree[i].txps;
        uint64_t j = i;
        while (j < tree.size()) {
            bool isPrefix = prefix.size() <= tree[j].txps.size() &&
                            std::mismatch(prefix.begin(), prefix.end(),
                                          tree[j].txps.begin(), tree[j].txps.end()).first == prefix.end();
            if (isPrefix) {
                tree[i].prefixCnt++;
                tree[i].prefixTotalCnt += tree[j].exactCnt;
            } else break;
            j++;
        }
    }
    std::cout << "numTxps exactCnt prefixCnt multiPrefixCnt\n";
    for (auto &t:tree) {
        std::cout << t.numTxps << " "
                  << t.exactCnt << " " << t.prefixCnt << " " << t.prefixTotalCnt << "\n";
    }
};

int pufferfishStats(pufferfish::StatsOptions &statsOpts) {
    auto indexDir = statsOpts.indexDir;
    switch (statsOpts.statType) {
        case pufferfish::StatType::ctab:
            doCtabStats(indexDir);
            break;
        case pufferfish::StatType::motif:
            doFindMotifs(indexDir);
            break;
    }
    return 0;
}


