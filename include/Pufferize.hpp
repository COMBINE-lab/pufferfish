#ifndef OUR_GFA_READER_H
#define OUR_GFA_READER_H

#include "FatPufferGraph.hpp"
#include "Util.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/vector.hpp"
#include "sparsepp/spp.h"
#include "sdsl/int_vector.hpp"

#include "string_view.hpp"
#include "Kmer.hpp"
#include "zstr/zstr.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

constexpr const uint64_t HighBitMask = 1uLL << (sizeof(uint64_t) * 8 - 1);
constexpr const uint64_t LowBitsMask = HighBitMask - 1;

struct SubseqView {
    uint64_t contigId;
    uint32_t startp;
    uint32_t len;

    SubseqView() {};

    SubseqView(uint64_t contigIdIn,
               uint32_t startpIn,
               uint32_t lenIn) : contigId(contigIdIn),
                                 startp(startpIn),
                                 len(lenIn) {}
};

struct ContigOri {
    uint64_t wrd;

    inline bool ori() {
        return (wrd & HighBitMask);
    }

    inline uint64_t id() {
        return (wrd & LowBitsMask);
    }

    ContigOri() {}

    ContigOri(uint64_t contig, bool orientation) {
        if (contig & HighBitMask) {
            std::cerr << "ERROR!! highest bit was already set. shouldn't happen";
            std::exit(1);
        }
        wrd = contig | (((uint64_t) orientation) << 63);
    }

    ContigOri(uint64_t wrdIn) : wrd(wrdIn) {}

    bool const operator<(const ContigOri &o) {
        return wrd < o.wrd;
    }
};

enum class Side : uint8_t {
    LEFT = 0,
    RIGHT = 1,
    LEFT_RIGHT = 2
};

class Pufferize {
private:
    std::string gfaFileName_;
    std::unique_ptr<zstr::ifstream> file;
    size_t k;
    spp::sparse_hash_map<uint64_t, SubseqView> id2seq;
    std::vector<uint64_t> ksizeContig;
    // path maps each transcript_id to a pair of <contig_id, orientation>
    // orientation : +/true main, -/false reverse

    spp::sparse_hash_map<uint64_t, std::vector<std::string>> start2path;
    spp::sparse_hash_map<std::string, ContigOri> path2start;
    spp::sparse_hash_map<uint64_t, Side> donTouch;
    spp::sparse_hash_map<uint64_t, std::vector<ContigOri> > old2newids;
    pufg::Graph semiCG;

    std::vector<sdsl::int_vector<2>> contigs;

    void encodeSeq(sdsl::int_vector<2> &contig, std::string &str);

    void
    processContigSeq(uint64_t &contigId, std::string &contigSeq,
                     spp::sparse_hash_map<uint64_t, uint64_t> &seq2id,
                     uint64_t &idCntr);

    std::vector<std::pair<uint64_t, bool> > &convertOld2newPath(std::vector<std::pair<uint64_t, bool> > &oContigVec);

    void setDonTouchSide(std::vector<std::pair<uint64_t, bool> > &contigVec);

    void updateDonTouchSide(uint64_t &id, Side side);

    void mergeDonTouchSides(uint64_t &mergedId, bool mergedOri, uint64_t &mergingId, bool mergingOri);

    bool dont_touch(uint64_t &nodeId, Side side);

    void setPathStart(std::string &pathId, std::pair<uint64_t, bool> &startNode);

    void updatePathStart(uint64_t &oStartNodeId, bool oStartNodeOri, uint64_t &nStartNodeId, bool nStartNodeOri);

    bool isCornerCase(pufg::Node &n, bool mergeLeft);

    void mergeLeft(pufg::Node &n);

    void mergeRight(pufg::Node &n);

    void updateGraph(std::vector<std::pair<uint64_t, bool> > &contigVec);

    std::string convertSubview2str(SubseqView ss);

public:
    Pufferize(const char *gfaFileName, size_t input_k);

    void parseFile();

//  void buildGraph();
    void randomWalk();

    void reconstructPathAndWrite(const char *gfaOutFileName, const char *refFileName);
};

#endif
