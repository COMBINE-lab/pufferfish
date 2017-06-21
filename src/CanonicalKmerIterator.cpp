// This code is derived and modified from Pall Melsted's
// KmerIterator class (part of bfgraph) : https://github.com/pmelsted/bfgraph/blob/master/src/KmerIterator.cpp

#include <iterator>
#include <utility>
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"

namespace pufferfish {
/* Note: That an iter is exhausted means that (iter._invalid == true) */

// use:  ++iter;
// pre:
// post: *iter is now exhausted
//       OR *iter is the next valid pair of kmer and location
CanonicalKmerIterator& CanonicalKmerIterator::operator++() {
  int pos_ = p_.second;
  if (!invalid_) {
    if (pos_+k_ >= s_.length()) {
      invalid_ = true;
      return *this;
    } else {
      find_next(pos_,pos_+k_-1,true);
      return *this;
    }
  }
  return *this;
}


// use:  iter++;
// pre:
// post: iter has been incremented by one
CanonicalKmerIterator CanonicalKmerIterator::operator++(int) {
  CanonicalKmerIterator tmp(*this);
  operator++();
  return tmp;
}


// use:  val = (a == b);
// pre:
// post: (val == true) if a and b are both exhausted
//       OR a and b are in the same location of the same string.
//       (val == false) otherwise.
bool CanonicalKmerIterator::operator==(const CanonicalKmerIterator& o) {
  if (invalid_  || o.invalid_) {
    return invalid_ && o.invalid_;
  } else {
    return (s_ == o.s_) && (p_.second == o.p_.second);
  }
}


// use:  p = *iter;
// pre:
// post: p is NULL or a pair of Kmer and int
std::pair<CanonicalKmer, int>& CanonicalKmerIterator::operator*() {
  return p_;
}


// use:  example 1: km = iter->first;
//       example 2:  i = iter->second;
// pre:  *iter is not NULL
// post: km will be (*iter).first, i will be (*iter).second
std::pair<CanonicalKmer, int> *CanonicalKmerIterator::operator->() {
  return &(operator*());
}


// use:  iter.raise(km, rep);
// post: iter has been incremented by one
//       if iter is not invalid, km is iter->first and rep is km.rep()
/*
void CanonicalKmerIterator::raise(Kmer& km, Kmer& rep) {
  operator++();
  if (!invalid_) {
    km = p_.first;
    rep = km.rep();
  }
}
*/

// use:  find_next(i,j, last_valid);
// pre:
// post: *iter is either invalid or is a pair of:
//       1) the next valid kmer in the string that does not have any 'N'
//       2) the location of that kmer in the string
void CanonicalKmerIterator::find_next(int i, int j, bool last_valid) {
  ++i;
  ++j;
  bool valid{false};
  // j is the last nucleotide in the k-mer we're building
  while (j < s_.length()) {
    // get the code for the last nucleotide, save it as c
    int c = my_mer::code(s_[j]);
    // c is a valid code if != -1
    if (c != -1) {
      p_.first.shiftFw(c);
      valid = (j - lastinvalid_ >= k_);
    } else {
      // if c is not a valid code, then j is the last invalid position
      lastinvalid_ = j;
      // the start position is the next (potentially) valid position
      i = j+1;
      // this k-mer is clearly not valid
      valid = false;
    }
    if (valid) {
      p_.second = i;
      return;
    }
    ++j;
  }
  invalid_ = true;
}

// use:  find_next(i,j, last_valid);
// pre:
// post: *iter is either invalid or is a pair of:
//       1) the next valid kmer in the string that does not have any 'N'
//       2) the location of that kmer in the string
  /*
void CanonicalKmerIterator::find_next(size_t i, size_t j, bool last_valid) {

  ++i;
  ++j;

  while (s_[j] != 0) {
    //auto c = 
    //char c = s_[j] & 0xDF; // mask lowercase bit
    int c = my_mer::code(s_[j]);
    if (c != -1) {//} == 'A' || c == 'C' || c == 'G' || c == 'T') {
      if (last_valid) {
        p_.first.shiftFw(c);
        break; // default case,
      } else {
        if (i + k_ - 1 == j) {
          p_.first.fromStr(s_+i);
          last_valid = true;
          break; // create k-mer from scratch
        } else {
          ++j;
        }
      }
    } else {
      ++j;
      i = j;
      last_valid = false;
    }
  }
  if (i+k_-1 == j && s_[j] != 0) {
    p_.second = i;
  } else {
    invalid_ = true;
  }
}
  */
}
