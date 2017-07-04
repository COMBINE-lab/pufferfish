#ifndef __CANONICAL_KMER_HPP__
#define __CANONICAL_KMER_HPP__

#include "jellyfish/mer_dna.hpp"

// NO_MATCH => two k-mers k1, and k2 are distinct such that k1 != k2 and rc(k1)
// != k2
// IDENTITY_MATCH => k1 = k2
// TWIN_MATCH => rc(k1) = k2
enum class KmerMatchType : uint8_t { NO_MATCH = 0, IDENTITY_MATCH, TWIN_MATCH };

using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

/**
 * This class wraps a pair of jellifish k-mers
 * (i.e., mer_dna objects).  It maintains both a
 * k-mer and its reverse complement at all times
 * to make the operation of retreiving the canonical
 * k-mer efficent.
 */
class CanonicalKmer {
private:
  my_mer fw_;
  my_mer rc_;

public:
  CanonicalKmer() = default;
  CanonicalKmer(CanonicalKmer&& other) = default;
  CanonicalKmer(CanonicalKmer& other) = default;
  CanonicalKmer(const CanonicalKmer& other) = default;
  CanonicalKmer& operator=(CanonicalKmer& other) = default;

  static inline void k(int kIn) { my_mer::k(kIn); }
  static inline int k() { return my_mer::k(); }

  inline bool fromStr(const std::string& s) {
    auto k = my_mer::k();
    if (s.length() < k) {
      return false;
    }
    for (size_t i = 0; i < k; ++i) {
      fw_.shift_right(s[i]);
      rc_.shift_left(my_mer::complement(s[i]));
    }
    return true;
  }

  inline bool fromStr(const char* s) {
    auto k = my_mer::k();
    // if (s.length() < k) { return false; }
    for (size_t i = 0; i < k; ++i) {
      fw_.shift_right(s[i]);
      rc_.shift_left(my_mer::complement(s[i]));
    }
    return true;
  }

  inline void fromNum(uint64_t w) {
    fw_.word__(0) = w;
    rc_ = fw_.get_reverse_complement();
  }

  inline void swap(){
    std::swap(fw_, rc_);
    //my_mer tmp = fw_ ;
    //fw_ = rc_ ;
    //rc_ = tmp ;
  }

  inline bool isFwCanonical() const { return fw_ < rc_; }

  inline auto shiftFw(int c) -> decltype(this->fw_.shift_right(c)) {
    rc_.shift_left(my_mer::complement(c));
    return fw_.shift_right(c);
  }

  inline auto shiftBw(int c) -> decltype(this->fw_.shift_left(c)) {
    rc_.shift_right(my_mer::complement(c));
    return fw_.shift_left(c);
  }

  inline auto shiftFw(char c) -> decltype(this->fw_.shift_right(c)) {
    int x = my_mer::code(c);
    if (x == -1)
      return 'N';
    rc_.shift_left(my_mer::complement(x));
    return my_mer::rev_code(fw_.shift_right(x));
  }

  inline auto shiftBw(char c) -> decltype(this->fw_.shift_left(c)) {
    int x = my_mer::code(c);
    if (x == -1)
      return 'N';
    rc_.shift_right(my_mer::complement(x));
    return my_mer::rev_code(fw_.shift_left(x));
  }

  inline uint64_t getCanonicalWord() const {
    return (fw_.word(0) < rc_.word(0)) ? fw_.word(0) : rc_.word(0);
  }

  inline const my_mer& getCanonical() const {
    return (fw_.word(0) < rc_.word(0)) ? fw_ : rc_;
  }

  inline const my_mer& fwMer() const { return fw_; }

  inline const my_mer& rcMer() const { return rc_; }

  inline uint64_t fwWord() const { return fw_.word(0); }

  inline uint64_t rcWord() const { return rc_.word(0); }

  inline KmerMatchType isEquivalent(const my_mer& m) const {
    return m.word(0) == fwWord()
               ? KmerMatchType::IDENTITY_MATCH
               : (m.word(0) == rcWord() ? KmerMatchType::TWIN_MATCH
                                        : KmerMatchType::NO_MATCH);
  }
  inline KmerMatchType isEquivalent(uint64_t m) const {
    return m == fwWord() ? KmerMatchType::IDENTITY_MATCH
                         : (m == rcWord() ? KmerMatchType::TWIN_MATCH
                                          : KmerMatchType::NO_MATCH);
  }

  inline std::string to_str() const {
    std::string s = fw_.to_str();
    std::reverse(s.begin(), s.end());
    return s;
  }

  bool operator==(const CanonicalKmer& rhs) const {
    return this->fw_ == rhs.fw_;
  }
  bool operator!=(const CanonicalKmer& rhs) const {
    return !this->operator==(rhs);
  }
  bool operator<(const CanonicalKmer& rhs) const { return this->fw_ < rhs.fw_; }
  bool operator<=(const CanonicalKmer& rhs) const {
    return *this < rhs || *this == rhs;
  }
  bool operator>(const CanonicalKmer& rhs) const { return !(*this <= rhs); }
  bool operator>=(const CanonicalKmer& rhs) const { return !(*this < rhs); }
  bool is_homopolymer() const { return fw_.is_homopolymer(); }
};

#endif // __CANONICAL_KMER_HPP__
