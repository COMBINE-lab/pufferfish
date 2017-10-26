#ifndef __PUFFERFISH_COMMON_TYPES_HPP__
#define __PUFFERFISH_COMMON_TYPES_HPP__

#include "core/range.hpp"

namespace pufferfish {
  namespace common_types {
    using ReferenceID = size_t;
    template <typename T>
    using IterRange = core::range<typename T::iterator>;
  }
}

#endif // __PUFFERFISH_COMMON_TYPES_HPP__
