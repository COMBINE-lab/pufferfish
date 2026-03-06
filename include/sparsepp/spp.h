#ifndef SPARSEPP_SPP_H
#define SPARSEPP_SPP_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <type_traits>
#include <utility>

#include "parallel_hashmap/phmap.h"

namespace spp {

template <typename T>
using spp_hash = std::hash<T>;

template <typename T>
using libc_allocator = std::allocator<T>;

template <typename T>
using spp_allocator = std::allocator<T>;

template <typename T>
inline void hash_combine(std::size_t& seed, const T& v) {
  seed ^= spp_hash<T>{}(v) + 0x9e3779b9U + (seed << 6U) + (seed >> 2U);
}

template <
    class Key,
    class T,
    class HashFcn = spp_hash<Key>,
    class EqualKey = std::equal_to<Key>,
    class Alloc = std::allocator<std::pair<const Key, T>>>
using sparse_hash_map = phmap::flat_hash_map<Key, T, HashFcn, EqualKey, Alloc>;

template <
    class Value,
    class HashFcn = spp_hash<Value>,
    class EqualKey = std::equal_to<Value>,
    class Alloc = std::allocator<Value>>
using sparse_hash_set = phmap::flat_hash_set<Value, HashFcn, EqualKey, Alloc>;

}  // namespace spp

#endif  // SPARSEPP_SPP_H
