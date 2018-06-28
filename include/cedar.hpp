#include <vector>
#include <set> // std::set
#include <deque> // std::deque
//#include <queue> // std::priority_queue
#include <memory>
#include "spdlog/spdlog.h"
#include "sparsepp/spp.h"
#include "EquivalenceClassBuilder.hpp"

#include "taxa.h"

struct ReadInfo {
    std::string rid;
    uint32_t cnt = 0;
    uint32_t len = 0;
    std::vector<TaxaNode> mappings;
};

#include "PuffMappingReader.hpp"
#include "SAMReader.hpp"

template<class ReaderType>
class Cedar {
    public:
        Cedar(std::string& taxonomyTree_filename, std::string& 
              refId2TaxId_filename, std::string pruneLevelIn, double filteringThresholdIn,
              bool flatAbund,
              std::shared_ptr<spdlog::logger> loggerIn);
        void run(std::string mapperOutput_filename, 
                 bool requireConcordance,
                 size_t maxIter, 
                 double eps,
                 std::string& output_filename,
                 bool onlyUniq,
                 bool onlyPerf) {
                    loadMappingInfo(mapperOutput_filename, requireConcordance, onlyUniq, onlyPerf);
                    basicEM(maxIter, eps);
                    if (!flatAbund) {
                        serialize(output_filename);
                    } else {
                        serializeFlat(output_filename);
                    }
                    //std::cout << "I guess that's it\n";
                 };
    private:
        bool readHeader(std::ifstream& mfile);
        void loadMappingInfo(std::string mapperOutput_filename, bool requireConcordance, bool onlyUniq, bool onlyPerfect);
        bool basicEM(size_t maxIter, double eps);
        void serialize(std::string& output_filename);
        void serializeFlat(std::string& output_filename);
        
        spp::sparse_hash_map<uint32_t, TaxaNode> taxaNodeMap;
        spp::sparse_hash_map<std::string, uint32_t> refId2taxId;
        spp::sparse_hash_map<uint32_t, uint32_t> seqToTaxMap;        
        Rank pruningLevel = Rank::SPECIES;
        std::set<uint64_t> activeTaxa;
        uint64_t rootId = 1;
        double filteringThreshold = 0;
        spp::sparse_hash_map<uint64_t, double> strain;
        uint64_t readCnt = 0;
        bool flatAbund = false;
        std::vector<std::vector<std::pair<uint64_t, double>>> readPerStrainProb;
        EquivalenceClassBuilder eqb;
        ReaderType mappings;
        std::shared_ptr<spdlog::logger> logger;
        
        spp::sparse_hash_map<uint32_t, uint64_t> cov;


};

