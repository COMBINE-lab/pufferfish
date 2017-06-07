#ifndef _PUFFERFISH_INDEX_HPP_
#define _PUFFERFISH_INDEX_HPP_

#include <vector>

#include "core/range.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"
#include "cereal/archives/json.hpp"

#include "CanonicalKmer.hpp"
#include "BooPHF.h"
#include "Util.hpp"

class PufferfishIndex {
  typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
  typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;
private:
  uint32_t k_{0};
  uint64_t numKmers_{0};
  std::vector<std::string> refNames_;
  std::vector<std::vector<util::Position>> contigTable_;
  sdsl::bit_vector contigBoundary_;
  sdsl::bit_vector::rank_1_type contigRank_;
  sdsl::bit_vector::select_1_type contigSelect_;
  sdsl::int_vector<> seq_;
  sdsl::int_vector<> pos_;
  std::unique_ptr<boophf_t> hash_{nullptr};
public:
  PufferfishIndex();
  PufferfishIndex(const std::string& indexPath);

  const uint32_t k();
  const std::vector<util::Position>& refList(uint64_t contigRank);
  const std::string& refName(uint64_t refRank);
  //uint64_t getRawPos(uint64_t km);
  bool contains(CanonicalKmer& mer);
  uint64_t getRawPos(CanonicalKmer& mer);
  bool isValidPos(uint64_t pos);
  auto getRefPos(CanonicalKmer& mer) -> util::ProjectedHits;
};

#endif // _PUFFERFISH_INDEX_HPP_
