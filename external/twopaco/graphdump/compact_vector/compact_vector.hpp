#ifndef __COMPACT_VECTOR_H__
#define __COMPACT_VECTOR_H__

#include <new>
#include <stdexcept>
#include <fstream>
#include <cstring>
#include <system_error> // for std::error_code
#include "mio.hpp"
#include "compact_iterator.hpp"
#include <bitset>
#include <iostream>

namespace compact {

inline uint64_t get_bits_per_element(const std::string& fname) {
	// load the vector by reading from file
	std::ifstream ifile(fname, std::ios::binary);
	uint64_t static_flag{0};
	ifile.read(reinterpret_cast<char*>(&static_flag), sizeof(static_flag));
	uint64_t bits_per_element;
	ifile.read(reinterpret_cast<char*>(&bits_per_element), sizeof(bits_per_element));
	ifile.close();
	return bits_per_element;
}

namespace vector_imp {
inline int clz(unsigned int x) { return __builtin_clz(x); }
inline int clz(unsigned long x) { return __builtin_clzl(x); }
inline int clz(unsigned long long x) { return __builtin_clzll(x); }

// XXX TODO Missing copy and move constructors
template<class Derived,
         typename IDX, unsigned BITS, typename W, typename Allocator, unsigned UB, bool TS>
class vector {
  Allocator m_allocator;
  size_t    m_size;             // Size in number of elements
  size_t    m_capacity;         // Capacity in number of elements
  W*        m_mem;
  mio::mmap_source ro_mmap;

public:
  // Number of bits required for indices/values in the range [0, s).
  static unsigned required_bits(size_t s) {
    unsigned res = bitsof<size_t>::val - 1 - clz(s);
    res += (s > ((size_t)1 << res)) + (std::is_signed<IDX>::value ? 1 : 0);
    return res;
  }

  static size_t elements_to_words(size_t size, unsigned bits) {
    size_t total_bits = size * bits;
    return total_bits / UB + (total_bits % UB != 0);
  }

  typedef compact::iterator<IDX, BITS, W, TS, UB>   iterator;
  typedef compact::const_iterator<IDX, BITS, W, UB> const_iterator;
  typedef compact::iterator<IDX, BITS, W, true, UB> mt_iterator; // Multi thread safe version
  typedef std::reverse_iterator<iterator>        reverse_iterator;
  typedef std::reverse_iterator<const_iterator>  const_reverse_iterator;

  vector(size_t s, size_t mem, Allocator allocator = Allocator())
    : m_allocator(allocator)
    , m_size(s)
    , m_capacity(s)
    , m_mem(m_allocator.allocate(mem))
  {
    static_assert(UB <= bitsof<W>::val, "used_bits must be less or equal to the number of bits in the word_type");
    static_assert(BITS <= UB, "number of bits larger than usable bits");
  }
  explicit vector(Allocator allocator = Allocator())
    : vector(0, 0, allocator)
  { }
  ~vector() {
    if (!ro_mmap.is_mapped()) {
      m_allocator.deallocate(m_mem, elements_to_words(m_capacity, bits()));
    }
  }

  const_iterator begin() const { return const_iterator(m_mem, bits(), 0); }
  iterator begin() { return iterator(m_mem, bits(), 0); }
  const_iterator end() const { return begin() + m_size; }
  iterator end() { return begin() + m_size; }
  const_iterator cbegin() const { return begin(); }
  const_iterator cend() const { return end(); }
  const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
  reverse_iterator rbegin() { return reverse_iterator(end()); }
  const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }
  reverse_iterator rend() { return reverse_iterator(begin()); }
  const_reverse_iterator crbegin() const { return const_reverse_iterator(end()); }
  const_reverse_iterator crend() const { return const_reverse_iterator(begin()); }

  // Multi thread safe iterator
  mt_iterator mt_begin() { return mt_iterator(m_mem, bits(), 0); }
  mt_iterator mt_end() { return begin() + m_size; }

  IDX operator[](size_t i) const {
    return BITS
      ? *const_iterator(m_mem + (i * BITS) / UB, BITS, (i * BITS) % UB)
      : *const_iterator(m_mem + (i * bits()) / UB, bits(), (i * bits()) % UB);
    // return cbegin()[i];
  }
  typename iterator::lhs_setter_type operator[](size_t i) {
    return BITS
      ? typename iterator::lhs_setter_type(m_mem + (i * BITS) / UB, BITS, (i * BITS) % UB)
      : typename iterator::lhs_setter_type(m_mem + (i * bits()) / UB, bits(), (i * bits()) % UB);
    //  return begin()[i];
  }
  IDX front() const { return *cbegin(); }
  typename iterator::lhs_setter_type front() { return *begin(); }
  IDX back() const { return *(cbegin() + (m_size - 1)); }
  typename iterator::lhs_setter_type back() { return *(begin() + (m_size - 1)); }

  size_t size() const { return m_size; }
  bool empty() const { return m_size == 0; }
  size_t capacity() const { return m_capacity; }
  uint64_t * get_words() const { return m_mem; }

  void push_back(IDX x) {
    if(m_size == m_capacity)
      enlarge();
    *end() = x;
    ++m_size;
  }

  void pop_back() { --m_size; }
  void clear() { m_size = 0; }
  void emplace_back(IDX x) { push_back(x); }

  W* get() { return m_mem; }
  const W* get() const { return m_mem; }
  size_t bytes() const { return sizeof(W) * elements_to_words(m_size, bits()); }
  inline unsigned bits() const { return static_cast<const Derived*>(this)->bits(); }
  static constexpr unsigned static_bits() { return BITS; }
  static constexpr unsigned used_bits() { return UB; }
  static constexpr bool thread_safe() { return TS; }

  void serialize(std::ofstream& of) {
    uint64_t static_flag = (static_bits() == bits()) ? 1 : 0;
    of.write(reinterpret_cast<char*>(&static_flag), sizeof(static_flag));
    if (static_flag != 0) {
      uint64_t bits_per_element = static_bits();
      of.write(reinterpret_cast<char*>(&bits_per_element), sizeof(bits_per_element));
    } else {
      uint64_t bits_per_element = bits();
      of.write(reinterpret_cast<char*>(&bits_per_element), sizeof(bits_per_element));
    }
    uint64_t w_size = m_size;
    std::cerr << "size: " << w_size << "\n";
    of.write(reinterpret_cast<char*>(&w_size), sizeof(w_size));
    uint64_t w_capacity = m_capacity;
    of.write(reinterpret_cast<char*>(&w_capacity), sizeof(w_capacity));
    of.write(reinterpret_cast<char*>(m_mem), bytes());
    //std::cerr << "wrote " << bytes() << " bytes of data at the end\n";
  }

  void deserialize(const std::string& fname, bool mmap) {
    std::error_code error;
    if (mmap) {
      // load the vector *read only* by mmap
      ro_mmap.map(fname, error);
      if (error) { std::cerr << "error = " << error << "\n"; }
      const char* data = ro_mmap.data();
      data += sizeof(uint64_t);
      uint64_t bits_per_element;
      std::memcpy(reinterpret_cast<void*>(&bits_per_element), reinterpret_cast<void*>(const_cast<char*>(data)), sizeof(bits_per_element));
      //std::cerr<< "bits / element = " << bits_per_element << "\n";
      data += sizeof(W);
      uint64_t w_size{0};
      std::memcpy(reinterpret_cast<void*>(&w_size), reinterpret_cast<void*>(const_cast<char*>(data)), sizeof(w_size));
      m_size = w_size;
      std::cerr<< "size = " << m_size << "\n";
      data += sizeof(w_size);
      uint64_t w_capacity{0};
      std::memcpy(reinterpret_cast<void*>(&w_capacity), reinterpret_cast<void*>(const_cast<char*>(data)), sizeof(w_capacity));
      m_capacity = w_capacity;
      //std::cerr<< "capacity = " << m_capacity << "\n";
      data += sizeof(w_capacity);
      m_allocator.deallocate(m_mem, elements_to_words(m_capacity, bits()));
      m_mem = reinterpret_cast<W*>(const_cast<char*>(data));
    } else {
      // load the vector by reading from file
      std::ifstream ifile(fname, std::ios::binary);
      uint64_t static_flag{0};
      ifile.read(reinterpret_cast<char*>(&static_flag), sizeof(static_flag));

      uint64_t bits_per_element;
      ifile.read(reinterpret_cast<char*>(&bits_per_element), sizeof(bits_per_element));

      //std::cerr<< "bits / element = " << bits_per_element << "\n";

      uint64_t w_size{0};
      ifile.read(reinterpret_cast<char*>(&w_size), sizeof(w_size));
      m_size = w_size;
      std::cerr<< "size = " << m_size << "\n";

      uint64_t w_capacity{0};
      ifile.read(reinterpret_cast<char*>(&w_capacity), sizeof(w_capacity));
      m_capacity = w_capacity;
      //std::cerr<< "capacity = " << m_capacity << "\n";

      m_allocator.deallocate(m_mem, elements_to_words(m_capacity, bits()));
      m_mem = m_allocator.allocate(elements_to_words(m_capacity, bits()));
      if(m_mem == nullptr) throw std::bad_alloc();
      ifile.read(reinterpret_cast<char*>(m_mem), sizeof(W) * elements_to_words(m_size, bits()));
    }

  }
  void touch_all_pages(uint64_t bits_per_element) {                                                                                                        
    uint64_t sum = 0;
    std::cerr<<"number of elements:"<<this->size() << "\n";
    std::cerr<<"page size:"<<mio::page_size() << "\n";
    auto elements_per_page = 8 * mio::page_size() / bits_per_element;
    std::cerr<<"elements per page:"<<elements_per_page<<"\n";
    std::cerr<<"bits per element:"<<bits_per_element<<"\n";
    for ( size_t i = 0; i<this->size(); i+=elements_per_page ) {
      sum += (*this)[i];
    }
    std::cerr<<sum << "\n";
  }
  uint64_t get_int(uint64_t from, uint64_t len) {
    uint64_t result = 0;
    if (len == 0) {
      std::cerr<<"get_int with length of zero is meaningless!";
      exit(1);
    }
    for(uint64_t i=0; i<=(len-1); i++) {
      uint64_t current = (*this)[from+i];
      auto shift_bits = (i)*BITS;
      result += (current<<shift_bits);
    }
    return result;
  }
  void reserve(size_t m) { m_capacity = m; m_mem = m_allocator.allocate(elements_to_words(m,BITS)); }

  void resize(size_t m) { m_size = m; m_capacity = m; m_mem = m_allocator.allocate(elements_to_words(m,BITS)); }

  vector& operator=(vector &vec){
    m_allocator = vec.m_allocator;
    m_size = vec.m_size;
    m_capacity = vec.m_capacity;
    m_mem = vec.m_mem;
    return *this;
  }

protected:
  void enlarge() {
    const size_t new_capacity = std::max(m_capacity * 2, (size_t)1);
    W* new_mem = m_allocator.allocate(new_capacity);
    if(new_mem == nullptr) throw std::bad_alloc();
    std::copy(m_mem, m_mem + elements_to_words(m_capacity, bits()), new_mem);
    m_allocator.deallocate(m_mem, m_capacity);
    m_mem      = new_mem;
    m_capacity = new_capacity;
  }
};

template<typename IDX, typename W, typename Allocator, unsigned UB, bool TS>
class vector_dyn
  : public vector_imp::vector<vector_dyn<IDX, W, Allocator, UB, TS>, IDX, 0, W, Allocator, UB, TS>
{
  typedef vector_imp::vector<vector_dyn<IDX, W, Allocator, UB, TS>, IDX, 0, W, Allocator, UB, TS> super;
  unsigned m_bits;    // Number of bits in an element

public:
  typedef typename super::iterator              iterator;
  typedef typename super::const_iterator        const_iterator;
  typedef IDX                                   value_type;
  typedef Allocator                             allocator_type;
  typedef typename iterator::lhs_setter_type    reference;
  typedef const reference                       const_reference;
  typedef iterator                              pointer;
  typedef const_iterator                        const_pointer;
  typedef std::reverse_iterator<iterator>       reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef ptrdiff_t                             difference_type;
  typedef size_t                                size_type;
  typedef W                                     word_type;

  vector_dyn(unsigned b, size_t s, Allocator allocator = Allocator())
    : super(s, super::elements_to_words(s, b), allocator)
    , m_bits(b)
  { }
  vector_dyn(unsigned b, Allocator allocator = Allocator())
    : super(allocator)
    , m_bits(b)
  { }

  inline unsigned bits() const { return m_bits; }
  void set_m_bits(size_t m) { m_bits = m; }
};

} // namespace vector_imp

template<typename IDX, unsigned BITS = 0, typename W = uint64_t, typename Allocator = std::allocator<W>>
class vector
  : public vector_imp::vector<vector<IDX, BITS, W, Allocator>, IDX, BITS, W, Allocator, bitsof<W>::val, false>
{
  typedef vector_imp::vector<vector<IDX, BITS, W, Allocator>, IDX, BITS, W, Allocator, bitsof<W>::val, false> super;

public:
  typedef typename super::iterator              iterator;
  typedef typename super::const_iterator        const_iterator;
  typedef IDX                                   value_type;
  typedef Allocator                             allocator_type;
  typedef typename iterator::lhs_setter_type    reference;
  typedef const reference                       const_reference;
  typedef iterator                              pointer;
  typedef const_iterator                        const_pointer;
  typedef std::reverse_iterator<iterator>       reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef ptrdiff_t                             difference_type;
  typedef size_t                                size_type;
  typedef W                                     word_type;

  vector(size_t s, Allocator allocator = Allocator())
    : super(s, super::elements_to_words(s, BITS), allocator)
  { }
  vector(Allocator allocator = Allocator())
    : super(allocator)
  { }

  static constexpr unsigned bits() { return BITS; }
};

template<typename IDX, typename W, typename Allocator>
class vector<IDX, 0, W, Allocator>
  : public vector_imp::vector_dyn<IDX, W, Allocator, bitsof<W>::val, false>
{
  typedef vector_imp::vector_dyn<IDX, W, Allocator, bitsof<W>::val, false> super;

public:
  typedef typename super::iterator              iterator;
  typedef typename super::const_iterator        const_iterator;
  typedef IDX                                   value_type;
  typedef Allocator                             allocator_type;
  typedef typename iterator::lhs_setter_type    reference;
  typedef const reference                       const_reference;
  typedef iterator                              pointer;
  typedef const_iterator                        const_pointer;
  typedef std::reverse_iterator<iterator>       reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef ptrdiff_t                             difference_type;
  typedef size_t                                size_type;
  typedef W                                     word_type;

  vector(unsigned b, size_t s, Allocator allocator = Allocator())
    : super(b, s, allocator)
  {
    if(b > bitsof<W>::val)
      throw std::out_of_range("Number of bits larger than usable bits");
  }
  vector(unsigned b, Allocator allocator = Allocator())
    : super(b, allocator)
  { }
};

template<typename IDX, unsigned BITS = 0, typename W = uint64_t, typename Allocator = std::allocator<W>>
class ts_vector
  : public vector_imp::vector<ts_vector<IDX, BITS, W, Allocator>, IDX, BITS, W, Allocator, bitsof<W>::val, true>
{
  typedef vector_imp::vector<ts_vector<IDX, BITS, W, Allocator>, IDX, BITS, W, Allocator, bitsof<W>::val, true> super;

public:
  typedef typename super::iterator              iterator;
  typedef typename super::const_iterator        const_iterator;
  typedef IDX                                   value_type;
  typedef Allocator                             allocator_type;
  typedef typename iterator::lhs_setter_type    reference;
  typedef const reference                       const_reference;
  typedef iterator                              pointer;
  typedef const_iterator                        const_pointer;
  typedef std::reverse_iterator<iterator>       reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef ptrdiff_t                             difference_type;
  typedef size_t                                size_type;
  typedef W                                     word_type;

  ts_vector(size_t s, Allocator allocator = Allocator())
    : super(s, super::elements_to_words(s, BITS), allocator)
  { }
  ts_vector(Allocator allocator = Allocator())
    : super(allocator)
  { }

  static constexpr unsigned bits() { return BITS; }
};


template<typename IDX, typename W, typename Allocator>
class ts_vector<IDX, 0, W, Allocator>
  : public vector_imp::vector_dyn<IDX, W, Allocator, bitsof<W>::val, true>
{
  typedef vector_imp::vector_dyn<IDX, W, Allocator, bitsof<W>::val, true> super;
public:
  typedef typename super::iterator              iterator;
  typedef typename super::const_iterator        const_iterator;
  typedef IDX                                   value_type;
  typedef Allocator                             allocator_type;
  typedef typename iterator::lhs_setter_type    reference;
  typedef const reference                       const_reference;
  typedef iterator                              pointer;
  typedef const_iterator                        const_pointer;
  typedef std::reverse_iterator<iterator>       reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef ptrdiff_t                             difference_type;
  typedef size_t                                size_type;
  typedef W                                     word_type;

  ts_vector(unsigned b, size_t s, Allocator allocator = Allocator())
    : super(b, s, allocator)
  {
    if(b > bitsof<W>::val)
      throw std::out_of_range("Number of bits larger than usable bits");
  }
  ts_vector(unsigned b, Allocator allocator = Allocator())
    : super(b, allocator)
  { }
};

template<typename IDX, unsigned BITS = 0, typename W = uint64_t, typename Allocator = std::allocator<W>>
class cas_vector
  : public vector_imp::vector<cas_vector<IDX, BITS, W, Allocator>, IDX, BITS, W, Allocator, bitsof<W>::val-1, true>
{
  typedef vector_imp::vector<cas_vector<IDX, BITS, W, Allocator>, IDX, BITS, W, Allocator, bitsof<W>::val-1, true> super;

public:
  typedef typename super::iterator              iterator;
  typedef typename super::const_iterator        const_iterator;
  typedef IDX                                   value_type;
  typedef Allocator                             allocator_type;
  typedef typename iterator::lhs_setter_type    reference;
  typedef const reference                       const_reference;
  typedef iterator                              pointer;
  typedef const_iterator                        const_pointer;
  typedef std::reverse_iterator<iterator>       reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef ptrdiff_t                             difference_type;
  typedef size_t                                size_type;
  typedef W                                     word_type;

  cas_vector(size_t s, Allocator allocator = Allocator())
    : super(s, super::elements_to_words(s, BITS), allocator)
  { }
  cas_vector(Allocator allocator = Allocator())
    : super(allocator)
  { }

  static constexpr unsigned bits() { return BITS; }
};

template<typename IDX, typename W, typename Allocator>
class cas_vector<IDX, 0, W, Allocator>
  : public vector_imp::vector_dyn<IDX, W, Allocator, bitsof<W>::val - 1, true>
{
  typedef vector_imp::vector_dyn<IDX, W, Allocator, bitsof<W>::val - 1, true> super;
public:
  typedef typename super::iterator              iterator;
  typedef typename super::const_iterator        const_iterator;
  typedef IDX                                   value_type;
  typedef Allocator                             allocator_type;
  typedef typename iterator::lhs_setter_type    reference;
  typedef const reference                       const_reference;
  typedef iterator                              pointer;
  typedef const_iterator                        const_pointer;
  typedef std::reverse_iterator<iterator>       reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef ptrdiff_t                             difference_type;
  typedef size_t                                size_type;
  typedef W                                     word_type;

  cas_vector(unsigned b, size_t s, Allocator allocator = Allocator())
    : super(b, s, allocator)
  {
    if(b > bitsof<W>::val - 1)
      throw std::out_of_range("Number of bits larger than usable bits");
  }
  cas_vector(unsigned b, Allocator allocator = Allocator())
    : super(b, allocator)
  { }
};

} // namespace compact

#endif /* __COMPACT_VECTOR_H__ */
