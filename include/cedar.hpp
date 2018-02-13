#include <vector>
#include <set> // std::set
#include <deque> // std::deque
//#include <queue> // std::priority_queue
#include "sparsepp/spp.h"
#include "taxa.h"


class Cedar {
    public:
        Cedar(std::string& taxonomyTree_filename, std::string& refId2TaxId_filename, std::string pruneLevelIn, double filteringThresholdIn);
        void loadMappingInfo(std::string mapperOutput_filename);
        bool basicEM();
        void serialize(std::string& output_filename);
    private:
        bool readHeader(std::ifstream& mfile);
        spp::sparse_hash_map<uint32_t, TaxaNode> taxaNodeMap;
        spp::sparse_hash_map<std::string, uint32_t> refId2taxId;
        Rank pruningLevel = Rank::SPECIES;
        std::set<uint64_t> activeTaxa;
        uint64_t rootId = 1;
        double filteringThreshold = 0;
        spp::sparse_hash_map<uint64_t, std::pair<float, uint64_t>> strain;
        uint64_t readCntr = 0;
        bool isPaired = true;
        std::vector<std::vector<std::pair<uint64_t, float>>> readPerStrainProb;
};