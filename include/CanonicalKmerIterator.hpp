// This code is derived and modified from Pall Melsted's
// KmerIterator class (part of bfgraph) : https://github.com/pmelsted/bfgraph/blob/master/src/KmerIterator.cpp

#ifndef MER_ITERATOR_HPP
#define MER_ITERATOR_HPP

#include <iterator>
#include "string_view.hpp"
#include "CanonicalKmer.hpp"

namespace pufferfish {
  class CanonicalKmerIterator : public std::iterator<std::input_iterator_tag, std::pair<CanonicalKmer, int>, int> {
    stx::string_view s_;
    std::pair<CanonicalKmer, int> p_;
    bool invalid_;
    int lastinvalid_;
    int k_;
  public:
		typedef CanonicalKmer value_type;
		typedef value_type& reference;
		typedef value_type* pointer;
		typedef std::input_iterator_tag iterator_category;
		typedef int64_t difference_type;
    CanonicalKmerIterator() : s_(), p_(), invalid_(true), lastinvalid_(-1), k_(CanonicalKmer::k()) {}
    CanonicalKmerIterator(const std::string& s) : s_(s), p_(), invalid_(false), lastinvalid_(-1), k_(CanonicalKmer::k()) { find_next(-1,-1,false);}
    CanonicalKmerIterator(const CanonicalKmerIterator& o) : s_(o.s_), p_(o.p_), invalid_(o.invalid_), lastinvalid_(o.lastinvalid_), k_(o.k_) {}

    CanonicalKmerIterator& operator++();
    CanonicalKmerIterator operator++(int);
    //void raise(CanonicalKmer& km, CanonicalKmer& rep);

    bool operator==(const CanonicalKmerIterator& o);
    bool operator!=(const CanonicalKmerIterator& o) { return !this->operator==(o);}

    std::pair<CanonicalKmer, int>& operator*();
    std::pair<CanonicalKmer, int>* operator->();

  private:
    void find_next(int i, int j, bool last_valid);
  };
}

#endif // MER_ITERATOR_HPP

