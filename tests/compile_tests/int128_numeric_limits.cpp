/*  This file is part of Jellyfish.

    This work is dual-licensed under 3-Clause BSD License or GPL 3.0.
    You can choose between one of them if you use this work.

`SPDX-License-Identifier: BSD-3-Clause OR  GPL-3.0`
*/

#include <unistd.h>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <limits>
#include <cstring>

using namespace std;

namespace std {
template<>
class numeric_limits<__int128> {
public:
  static const bool is_specialized = true;
  static __int128 max() { return (unsigned __int128)-1 >> 1; }
  static __int128 min() { return  max() + 1; }
  static const int  digits     = 127;
  static const int  digits10   = 38;
#define NLS64 numeric_limits<int64_t>
  static const bool is_signed  = NLS64::is_signed;
  static const bool is_integer = NLS64::is_integer;
  static const bool is_exact   = NLS64::is_exact;
  static const int  radix      = NLS64::radix;
  static __int128 epsilon() { return NLS64::epsilon(); }
  static __int128 round_error() { return NLS64::round_error(); }
  static const int                min_exponent      = NLS64::min_exponent;
  static const int                min_exponent10    = NLS64::min_exponent10;
  static const int                max_exponent      = NLS64::max_exponent;
  static const int                max_exponent10    = NLS64::max_exponent10;
  static const bool               has_infinity      = NLS64::has_infinity;
  static const bool               has_quiet_NaN     = NLS64::has_quiet_NaN;
  static const bool               has_signaling_NaN = NLS64::has_signaling_NaN;
  static const float_denorm_style has_denorm        = NLS64::has_denorm;
  static const bool               has_denorm_loss   = NLS64::has_denorm_loss;
  static __int128 infinity() { return NLS64::infinity(); }
  static __int128 quiet_NaN() { return NLS64::quiet_NaN(); }
  static __int128 signaling_NaN() { return NLS64::signaling_NaN(); }
  static __int128 denorm_min() { return NLS64::denorm_min(); }
  static const bool              is_iec559       = NLS64::is_iec559;
  static const bool              is_bounded      = NLS64::is_bounded;
  static const bool              is_modulo       = NLS64::is_modulo;
  static const bool              traps           = NLS64::traps;
  static const bool              tinyness_before = NLS64::tinyness_before;
  static const float_round_style round_style     = NLS64::round_style;
};

template<>
class numeric_limits<unsigned __int128> {
public:
  static const bool is_specialized = true;
  static __int128 max() { return (unsigned __int128)-1; }
  static __int128 min() { return  0; }
  static const int  digits     = 128;
  static const int  digits10   = 39;
#define NLU64 numeric_limits<uint64_t>
  static const bool is_signed  = NLU64::is_signed;
  static const bool is_integer = NLU64::is_integer;
  static const bool is_exact   = NLU64::is_exact;
  static const int  radix      = NLU64::radix;
  static __int128 epsilon() { return NLU64::epsilon(); }
  static __int128 round_error() { return NLU64::round_error(); }
  static const int                min_exponent      = NLU64::min_exponent;
  static const int                min_exponent10    = NLU64::min_exponent10;
  static const int                max_exponent      = NLU64::max_exponent;
  static const int                max_exponent10    = NLU64::max_exponent10;
  static const bool               has_infinity      = NLU64::has_infinity;
  static const bool               has_quiet_NaN     = NLU64::has_quiet_NaN;
  static const bool               has_signaling_NaN = NLU64::has_signaling_NaN;
  static const float_denorm_style has_denorm        = NLU64::has_denorm;
  static const bool               has_denorm_loss   = NLU64::has_denorm_loss;
  static __int128 infinity() { return NLU64::infinity(); }
  static __int128 quiet_NaN() { return NLU64::quiet_NaN(); }
  static __int128 signaling_NaN() { return NLU64::signaling_NaN(); }
  static __int128 denorm_min() { return NLU64::denorm_min(); }
  static const bool              is_iec559       = NLU64::is_iec559;
  static const bool              is_bounded      = NLU64::is_bounded;
  static const bool              is_modulo       = NLU64::is_modulo;
  static const bool              traps           = NLU64::traps;
  static const bool              tinyness_before = NLU64::tinyness_before;
  static const float_round_style round_style     = NLU64::round_style;
};
} // namespace std

int main(int argc, char* argv[]) {
    __int128 a = std::numeric_limits<__int128>::max();
    std::cout << ((a > static_cast<__int128>(0)) ? 'T' : 'F') << '\n';
    return 0;
}
