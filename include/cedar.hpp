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
    Cedar(std::string &taxonomyTree_filename, std::string &
    refId2TaxId_filename, std::string pruneLevelIn, double filteringThresholdIn,
          bool flatAbund,
          std::shared_ptr<spdlog::logger> loggerIn);

    void run(std::string mapperOutput_filename,
             bool requireConcordance,
             size_t maxIter,
             double eps,
             double minCnt,
             std::string &output_filename,
             bool onlyUniq,
             bool onlyPerf,
             uint32_t segmentSize,
             uint32_t rangeFactorizationBins);

private:
    bool readHeader(std::ifstream &mfile);

    void loadMappingInfo(std::string mapperOutput_filename, bool requireConcordance, bool onlyUniq, bool onlyPerfect,
                         uint32_t segmentSize, uint32_t rangeFactorizationBins);

    bool applySetCover(std::vector<double> &strainCnt, std::vector<bool> &strainValid,
                       std::vector<bool> &strainPotentiallyRemovable, double minCnt, bool canHelp, bool verbose = false);

    bool basicEM(size_t maxIter, double eps, double minCnt, bool verbose = false);

    void serialize(std::string &output_filename);

    void serializeFlat(std::string &output_filename);

    void calculateCoverage();

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

    spp::sparse_hash_map<uint64_t, double> strain_coverage;
    spp::sparse_hash_map<uint64_t, std::vector<uint32_t> > strain_coverage_bins;
};

