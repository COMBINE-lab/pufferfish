#include "FastxParser.hpp"
#include <cmath>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>
#include <bitset>

#include "ProgOpts.hpp"
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferFS.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "jellyfish/mer_dna.hpp"
#include "spdlog/spdlog.h"

#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "Util.hpp"
#include "PufferfishBinaryGFAReader.hpp"

void doCtabStats(std::string& indexDir) {
    std::cerr << "We're here!!\n\n" << indexDir << "\n";
};

int pufferfishStats(pufferfish::StatsOptions& statsOpts) {
  auto indexDir = statsOpts.indexDir;
  switch (statsOpts.statType) {
      case pufferfish::StatType::ctab:
        doCtabStats(indexDir);
        break;
  }
    return 0;
}


