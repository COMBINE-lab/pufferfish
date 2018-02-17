#include <vector>
#include <set> // std::set
#include <deque> // std::deque
//#include <queue> // std::priority_queue
#include <memory>
#include "spdlog/spdlog.h"
#include "sparsepp/spp.h"
#include "taxa.h"
#include "EquivalenceClassBuilder.hpp"


class Cedar {
    public:
        Cedar(std::string& taxonomyTree_filename, std::string& 
              refId2TaxId_filename, std::string pruneLevelIn, double filteringThresholdIn,
              bool flatAbund,
              std::shared_ptr<spdlog::logger> loggerIn);
        void loadMappingInfo(std::string mapperOutput_filename);
        bool basicEM(size_t maxIter, double eps);
        void serialize(std::string& output_filename);
        void serializeFlat(std::string& output_filename);
    private:
        bool readHeader(std::ifstream& mfile);
        spp::sparse_hash_map<uint32_t, TaxaNode> taxaNodeMap;
        spp::sparse_hash_map<std::string, uint32_t> refId2taxId;
        spp::sparse_hash_map<uint32_t, uint32_t> seqToTaxMap;
        Rank pruningLevel = Rank::SPECIES;
        std::set<uint64_t> activeTaxa;
        uint64_t rootId = 1;
        double filteringThreshold = 0;
        spp::sparse_hash_map<uint64_t, double> strain;
        uint64_t readCnt = 0;
        bool isPaired = true;
        bool flatAbund = false;
        std::vector<std::vector<std::pair<uint64_t, double>>> readPerStrainProb;
        EquivalenceClassBuilder eqb;
        std::vector<uint32_t> refLengths;
        std::vector<std::string> refNames;
        std::shared_ptr<spdlog::logger> logger;
};