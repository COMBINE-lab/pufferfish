#ifndef _PUFFERFISH_INDEX_HPP_
#define _PUFFERFISH_INDEX_HPP_

#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"
#include "cereal/archives/json.hpp"

#include "BooPHF.h"

class PufferfishIndex {
public:
  typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
  typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

  PufferfishIndex();
  PufferfishIndex(const std::string& indexPath);

private:
  uint32_t k_{1};
  sdsl::bit_vector contigBoundary_;
  sdsl::bit_vector::rank_1_type contigRank_;
  sdsl::bit_vector::select_1_type contigSelect_;
  sdsl::int_vector<> seq_;
  sdsl::int_vector<> pos_;
  std::unique_ptr<boophf_t> hash_{nullptr};
};

#endif // _PUFFERFISH_INDEX_HPP_
