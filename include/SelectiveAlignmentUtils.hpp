#ifndef __SELECTIVE_ALIGNMENT_UTILS__
#define __SELECTIVE_ALIGNMENT_UTILS__

#include "PuffAligner.hpp"

namespace selective_alignment {
  namespace utils {


void recoverOrphans(std::string& leftRead, std::string& rightRead, std::vector<pufferfish::util::MemCluster> &recoveredMemClusters, std::vector<pufferfish::util::JointMems> &jointMemsList, PuffAligner& puffaligner, bool verbose) {
  using MateStatus = pufferfish::util::MateStatus;

  puffaligner.orphanRecoveryMemCollection.reserve(2 * jointMemsList.size() + 1);
  recoveredMemClusters.reserve(2 * jointMemsList.size() + 1);
  for (auto& jointMem : jointMemsList) {
    auto& leftClust = jointMem.leftClust;
    auto& rightClust = jointMem.rightClust;
    auto tid = jointMem.tid;
    if (jointMem.mateStatus == MateStatus::PAIRED_END_LEFT) {
      bool recovered = puffaligner.recoverSingleOrphan(leftRead, rightRead, *leftClust, recoveredMemClusters, tid, true, verbose);
     if (recovered) {
       jointMem.rightClust = recoveredMemClusters.begin() + recoveredMemClusters.size() - 1;
       jointMem.recovered = true;
       jointMem.mateStatus = MateStatus::PAIRED_END_PAIRED;
    }
  } else if (jointMem.mateStatus == MateStatus::PAIRED_END_RIGHT) {
      bool recovered = puffaligner.recoverSingleOrphan(leftRead, rightRead, *rightClust, recoveredMemClusters, tid, false, verbose);
      if (recovered) {
        jointMem.leftClust = recoveredMemClusters.begin() + recoveredMemClusters.size() - 1;
        jointMem.recovered = true;
        jointMem.mateStatus = MateStatus::PAIRED_END_PAIRED;
      }
    }
  }
}

  }
}

#endif // __SELECTIVE_ALIGNMENT_UTILS__
