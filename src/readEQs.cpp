#include <string>
#include <string.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"
// #include "cereal/types/string.hpp"


int main (int argc, char* argv[]) {
    bool isIndex = false;
    if (argc == 3) {
        isIndex = strcmp(argv[2],"index")==0;
    }
    std::ifstream eqf(argv[1]);
    cereal::BinaryInputArchive eq(eqf);
    std::vector<uint32_t> ids;
    std::vector<std::vector<uint32_t>> tgs;
    std::vector<uint32_t> cnts;
    eq(ids);
    eq(tgs);
    eq(cnts);
    std::string outf = argv[1];
    outf += ".mat";
    std::ofstream eqMat(outf);
    uint64_t strain_count = 0;
    std::unordered_map<std::string, std::pair<int, std::vector<uint32_t>>> eqm; 
    for (uint64_t i = 0; i < tgs.size(); i++) {
        std::string colName;
        std::sort(tgs[i].begin(), tgs[i].end());
        // std::cerr << tgs[i].size() << " to ";
        tgs[i].erase(std::unique(tgs[i].begin(), tgs[i].end()), tgs[i].end());
        // std::cerr << tgs[i].size() << "\n";
        if (tgs[i].empty()) {
            std::cerr << " eq " << i << " is empty. Discarded\n";
            continue; 
        }
        for (auto c : tgs[i]) {
            if (strain_count < c) strain_count = c;
            colName += std::to_string(c) + "_";
        }
        colName = colName.substr(0, colName.find_last_of("_"));
        if (eqm.find(colName) == eqm.end()) {
            eqm[colName] = std::make_pair(0, tgs[i]);
        }
        eqm[colName].first += cnts[i];
    }
    strain_count++;
    std::cerr << "total number of strains: " << strain_count << " total number of eqs: " << eqm.size() << "\n";
    std::vector<std::vector<float>> strains(strain_count, std::vector<float>(tgs.size()));
    for (auto &kv : eqm) {
        float cnt = kv.second.first;
        if (not isIndex) {
            cnt = (((kv.second.first) * 100) / kv.second.second.size())/100.0;
        }
        for (auto c : kv.second.second) {
            eqMat << c << "," << kv.first << "," << cnt << "\n";
        }
    }
    eqMat.close(); 
}