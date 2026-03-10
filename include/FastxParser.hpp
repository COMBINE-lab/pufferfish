#ifndef __FASTX_PARSER__
#define __FASTX_PARSER__

#include "concurrentqueue.h"
#include "fcntl.h"
#include "kseq++.hpp"
#include "unistd.h"
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <thread>
#include <utility>
#include <vector>

using std::make_unique;

namespace fastx_parser {

// configuration for a FastxParser
struct ParserConfig {
  uint32_t numConsumers{1};
  uint32_t numParsers{1};
  uint32_t chunkSize{1000};
  bool parallelParsing{false};

  static ParserConfig with_consumers_single(uint32_t numConsumers) {
    return {numConsumers, 1, 1000, false};
  }

  static ParserConfig with_consumers_multi(uint32_t numConsumers) {
    return {numConsumers, 1, 1000, true};
  }
};

class ParserConfigBuilder {
public:
  ParserConfigBuilder() = default;

  ParserConfigBuilder& with_consumers(uint32_t numConsumers) {
    c_.numConsumers = numConsumers;
    return *this;
  }

  ParserConfigBuilder& with_parsers(uint32_t numParsers) {
    c_.numParsers = numParsers;
    return *this;
  }

  ParserConfigBuilder& within_set_parallelism(bool parallelParsing) {
    c_.parallelParsing = parallelParsing;
    return *this;
  }

  ParserConfigBuilder& with_chunk_size(uint32_t chunkSize) {
    c_.chunkSize = chunkSize;
    return *this;
  }

  ParserConfig build() { return c_; }

private:
  ParserConfig c_;
};

// holds a "set" of files that correspond to components (in different files)
// of the same fragment. For single-end reads, this is just a file, for
// paired-end reads, it is a pair of files, etc.
struct FileGroup {
  template <typename... Strings>
  explicit FileGroup(Strings&&... strs)
      : file_names{std::forward<Strings>(strs)...}, arity(sizeof...(strs)) {}

  template <typename Iterator>
  FileGroup(Iterator first, Iterator last)
      : file_names(first, last), arity(file_names.size()) {}

  std::vector<std::string> file_names;
  size_t arity{0};
};

// forward declaration of hte read trait
template <typename T> struct ReadTrait;

// Generic ReadSet that works for any arity
template <size_t N> struct ReadSet {
  std::array<klibpp::KSeq, N> reads;

  // Array-like access
  klibpp::KSeq& operator[](size_t i) { return reads[i]; }
  const klibpp::KSeq& operator[](size_t i) const { return reads[i]; }

  // Named accessors for convenience (only enabled when N is large enough)
  template <size_t M = N, typename = std::enable_if_t<(M >= 1)>>
  klibpp::KSeq& first() {
    return reads[0];
  }

  template <size_t M = N, typename = std::enable_if_t<(M >= 2)>>
  klibpp::KSeq& second() {
    return reads[1];
  }

  template <size_t M = N, typename = std::enable_if_t<(M >= 3)>>
  klibpp::KSeq& third() {
    return reads[2];
  }

  template <size_t M = N, typename = std::enable_if_t<(M >= 1)>>
  const klibpp::KSeq& first() const {
    return reads[0];
  }

  template <size_t M = N, typename = std::enable_if_t<(M >= 2)>>
  const klibpp::KSeq& second() const {
    return reads[1];
  }

  template <size_t M = N, typename = std::enable_if_t<(M >= 3)>>
  const klibpp::KSeq& third() const {
    return reads[2];
  }
};

// Full specialization for N = 1
template <> struct ReadSet<1> {
  klibpp::KSeq read;
  
  klibpp::KSeq& operator[](size_t i) { (void)i; return read; }
  const klibpp::KSeq& operator[](size_t i) const { (void)i; return read; }
  
  // Simplified accessors - no need for enable_if
  klibpp::KSeq& first() { return read; }
  const klibpp::KSeq& first() const { return read; }
  
  // You could also add a get() method for more intuitive single-element access
  klibpp::KSeq& get() { return read; }
  const klibpp::KSeq& get() const { return read; }
};

// Specialization of ReadTrait for ReadSet<N>
template <size_t N> struct ReadTrait<ReadSet<N>> {
  static constexpr size_t arity = N;
  static klibpp::KSeq& get(ReadSet<N>& t, size_t i) { return t[i]; }
};

// If you want to distinguish qual vs non-qual types:
template <size_t N> struct ReadQualSet {
  std::array<klibpp::KSeq, N> reads;

  klibpp::KSeq& operator[](size_t i) { return reads[i]; }
  const klibpp::KSeq& operator[](size_t i) const { return reads[i]; }

  template <size_t M = N, typename = std::enable_if_t<(M >= 1)>>
  klibpp::KSeq& first() {
    return reads[0];
  }

  template <size_t M = N, typename = std::enable_if_t<(M >= 2)>>
  klibpp::KSeq& second() {
    return reads[1];
  }

  template <size_t M = N, typename = std::enable_if_t<(M >= 3)>>
  klibpp::KSeq& third() {
    return reads[2];
  }
};

// Specialization of ReadTrait for ReadQualSet<N>
template <size_t N> struct ReadTrait<ReadQualSet<N>> {
  static constexpr size_t arity = N;
  static klibpp::KSeq& get(ReadQualSet<N>& t, size_t i) { return t[i]; }
};

// Specialization for KSeq (single read) - keep this for backward compatibility
template <> struct ReadTrait<klibpp::KSeq> {
  static constexpr size_t arity = 1;
  static klibpp::KSeq& get(klibpp::KSeq& t, size_t) { return t; }
};

// Type aliases for convenience and backward compatibility
using ReadSeq = ReadSet<1>;//klibpp::KSeq;
using ReadPair = ReadSet<2>;
using ReadTriple = ReadSet<3>;
/*
using ReadQuad = ReadSet<4>;
using ReadQuint = ReadSet<5>;
using ReadSextuple = ReadSet<6>;
using ReadSeptuple = ReadSet<7>;
using ReadOctuple = ReadSet<8>;
using ReadQualQuad = ReadQualSet<4>;
*/

using ReadQualPair = ReadQualSet<2>;
using ReadQualTriple = ReadQualSet<3>;

struct ChunkFragOffset {
  uint32_t file_idx{0};
  uint32_t frag_idx{0};
};

template <typename T> class ReadChunk {
public:
  ReadChunk(size_t want) : group_(want), want_(want), have_(want) {}
  inline void have(size_t num) { have_ = num; }
  inline size_t size() { return have_; }
  inline size_t want() const { return want_; }
  T& operator[](size_t i) { return group_[i]; }
  typename std::vector<T>::iterator begin() { return group_.begin(); }
  typename std::vector<T>::iterator end() { return group_.begin() + have_; }

  void set_chunk_frag_offset(uint32_t file_num, uint64_t frag_num) {
    frag_offset_.file_idx = file_num;
    frag_offset_.frag_idx = frag_num;
  }

  ChunkFragOffset chunk_frag_offset() const { return frag_offset_; }

private:
  std::vector<T> group_;
  size_t want_;
  size_t have_;
  ChunkFragOffset frag_offset_;
};

template <typename T> class ReadGroup {
public:
  ReadGroup(moodycamel::ProducerToken&& pt, moodycamel::ConsumerToken&& ct)
      : pt_(std::move(pt)), ct_(std::move(ct)) {}
  moodycamel::ConsumerToken& consumerToken() { return ct_; }
  moodycamel::ProducerToken& producerToken() { return pt_; }
  // get a reference to the chunk this ReadGroup owns
  std::unique_ptr<ReadChunk<T>>& chunkPtr() { return chunk_; }
  // get a *moveable* reference to the chunk this ReadGroup owns
  std::unique_ptr<ReadChunk<T>>&& takeChunkPtr() { return std::move(chunk_); }
  inline void have(size_t num) { chunk_->have(num); }
  inline size_t size() { return chunk_->size(); }
  inline size_t want() const { return chunk_->want(); }
  T& operator[](size_t i) { return (*chunk_)[i]; }
  typename std::vector<T>::iterator begin() { return chunk_->begin(); }
  typename std::vector<T>::iterator end() {
    return chunk_->begin() + chunk_->size();
  }
  void setChunkEmpty() { chunk_.release(); }
  bool empty() const { return chunk_.get() == nullptr; }
  ChunkFragOffset chunk_frag_offset() const {
    return chunk_->chunk_frag_offset();
  }

private:
  std::unique_ptr<ReadChunk<T>> chunk_{nullptr};
  moodycamel::ProducerToken pt_;
  moodycamel::ConsumerToken ct_;
};

template <typename T> class FastxParser {
public:
  template <typename... FileVectors>
  FastxParser(fastx_parser::ParserConfig& c, FileVectors&&... fileVectors)
      : inputStreamSets_{std::forward<FileVectors>(fileVectors)...},
        numParsing_(0), parallelParsing_(c.parallelParsing),
        blockSize_(c.chunkSize) {

    constexpr size_t arity = sizeof...(fileVectors);

    // Static assert to ensure arity matches T
    static_assert(arity == ReadTrait<T>::arity,
                  "Number of file vectors must match read type arity");

    // Validate that all vectors have the same size
    if (inputStreamSets_.empty()) {
      throw std::invalid_argument("Must provide at least one file vector");
    }

    size_t numFiles = inputStreamSets_[0].size();
    for (size_t i = 1; i < arity; ++i) {
      if (inputStreamSets_[i].size() != numFiles) {
        throw std::invalid_argument(
            "All file vectors must have the same number of files");
      }
    }

    // Validate no duplicate files across the same file set
    for (size_t fileIdx = 0; fileIdx < numFiles; ++fileIdx) {
      for (size_t i = 0; i < arity; ++i) {
        for (size_t j = i + 1; j < arity; ++j) {
          if (inputStreamSets_[i][fileIdx] == inputStreamSets_[j][fileIdx]) {
            std::cerr << "[WARNING]: Same file provided for multiple reads: "
                      << inputStreamSets_[i][fileIdx] << "\n";
          }
        }
      }
    }

    // Adjust numParsers if needed
    if (!parallelParsing_ && c.numParsers > numFiles) {
      auto limit = numFiles;
      std::cerr
          << "[INFO]: In serial-within-set mode, can't make use of more parsing threads than file sets (" << limit << "); "
             "setting # of parsing threads to "
          << limit << '\n';
      c.numParsers = limit;
    } else if (parallelParsing_ && c.numParsers > (numFiles + 1) * arity) {
      auto limit = (numFiles + 1) * arity;
      std::cerr
          << "[INFO]: In parallel-within-set mode, can't make use of more parsing threads than arity (" << arity << " + 1) * file sets (" << numFiles << "); "
             "setting # of parsing threads to "
          << limit << '\n';
      c.numParsers = limit;
    }

    numParsers_ = c.numParsers;
    numParsing_ = 0;

    // Initialize concurrent queues
    readQueue_ = moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>(
        4 * c.numConsumers, c.numParsers, 0);

    seqContainerQueue_ =
        moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>(
            4 * c.numConsumers, 1 + c.numConsumers, 0);

    workQueue_ = moodycamel::ConcurrentQueue<uint32_t>(numParsers_);

    // Push all file indices onto the work queue
    for (size_t i = 0; i < numFiles; ++i) {
      workQueue_.enqueue(i);
    }

    // Create tokens for each parsing thread
    for (size_t i = 0; i < numParsers_; ++i) {
      consumeContainers_.emplace_back(
          new moodycamel::ConsumerToken(seqContainerQueue_));
      produceReads_.emplace_back(new moodycamel::ProducerToken(readQueue_));
    }

    // Pre-allocate read chunks
    moodycamel::ProducerToken produceContainer(seqContainerQueue_);
    for (size_t i = 0; i < 4 * c.numConsumers; ++i) {
      auto chunk = make_unique<ReadChunk<T>>(blockSize_);
      seqContainerQueue_.enqueue(produceContainer, std::move(chunk));
    }
  }
  ~FastxParser();
  bool start();
  bool stop();
  ReadGroup<T> getReadGroup();
  bool refill(ReadGroup<T>& rg);
  void finishedWithGroup(ReadGroup<T>& s);

private:
  moodycamel::ProducerToken getProducerToken_();
  moodycamel::ConsumerToken getConsumerToken_();

  std::vector<std::string> inputStreams_;
  std::vector<std::string> inputStreams2_;
  std::vector<std::string> inputStreams3_; // For triplet files
  std::vector<std::vector<std::string>> inputStreamSets_;

  uint32_t numParsers_;
  alignas(64) std::atomic<uint32_t> numParsing_;
  bool parallelParsing_{true}; // Enable parallel parsing for multi-file modes

  // NOTE: Would like to use std::future<int> here instead, but that
  // solution doesn't seem to work.  It's unclear exactly why
  // see (https://twitter.com/nomad421/status/917748383321817088)
  std::vector<std::unique_ptr<std::thread>> parsingThreads_;

  // holds the results of the parsing threads, which is simply equal to
  // the return value of kseq_read() for the last call to that function.
  // A value < -1 signifies some sort of error.
  std::vector<int> threadResults_;

  size_t blockSize_;
  moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>> readQueue_,
      seqContainerQueue_;

  // holds the indices of files (file-pairs) to be processed
  moodycamel::ConcurrentQueue<uint32_t> workQueue_;

  std::vector<std::unique_ptr<moodycamel::ProducerToken>> produceReads_;
  std::vector<std::unique_ptr<moodycamel::ConsumerToken>> consumeContainers_;
  bool isActive_{false};

  // Helper for parallel parsing of N-way read sets
  template <size_t N> bool start_parallel_parsing_impl();
};
} // namespace fastx_parser

#endif // __FASTX_PARSER__
