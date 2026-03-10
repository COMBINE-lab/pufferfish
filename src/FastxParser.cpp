#include "FastxParser.hpp"
#include "FastxParserThreadUtils.hpp"

#include "fcntl.h"
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <poll.h>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>
#include <zlib.h>

namespace fastx_parser {

// ============================================================================
// Parallel parsing functions for multi-file modes
// ============================================================================

// Parse a single file and push ReadChunks to an intermediate queue
template <typename SingleReadT>
int parse_single_file(
    const std::string& filename, uint32_t file_idx,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<SingleReadT>>>&
        outputQueue,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<SingleReadT>>>&
        recycleQueue,
    uint32_t chunkSize = 1000) {
  using namespace klibpp;

  gzFile fp = gzopen(filename.c_str(), "r");
  if (!fp) {
    // Signal end-of-file with nullptr
    thread_utils::simple_wait([&]() { 
      return outputQueue.try_enqueue(nullptr);
    });
    return -4;
  }

  auto seq = make_kstream(fp, gzread, mode::in);
  uint32_t recordsInChunk = 0;

  // helper to allocate or recycle a chunk
  auto get_chunk = [&](size_t size) {
    std::unique_ptr<ReadChunk<SingleReadT>> chunk;
    if (recycleQueue.try_dequeue(chunk)) {
      chunk->have(0); // Reset count
      return chunk;
    }
    return std::make_unique<ReadChunk<SingleReadT>>(size);
  };

  // allocate initial chunk
  auto currentChunk = get_chunk(chunkSize);

  while (seq >> (*currentChunk)[recordsInChunk]) {
    recordsInChunk++;
    if (recordsInChunk == chunkSize) {
      currentChunk->have(recordsInChunk);

      thread_utils::simple_wait([&]() { 
        return outputQueue.try_enqueue(std::move(currentChunk));
      });

      // allocate next chunk
      currentChunk = get_chunk(chunkSize);
      recordsInChunk = 0;
    }
  }

  int result = 0;
  if (seq.err()) {
    result = -3;
  } else if (seq.tqs()) {
    result = -2;
  }

  // flush remaining read in last chunk
  if (recordsInChunk > 0) {
    currentChunk->have(recordsInChunk);
    thread_utils::simple_wait([&]() { 
      return outputQueue.try_enqueue(std::move(currentChunk));
    });
  }

  // Signal end-of-file with nullptr
  thread_utils::simple_wait([&]() { 
      return outputQueue.try_enqueue(nullptr);
  });

  gzclose(fp);
  return result;
}

template <typename T> ReadGroup<T> FastxParser<T>::getReadGroup() {
  return ReadGroup<T>(getProducerToken_(), getConsumerToken_());
}

template <typename T>
moodycamel::ProducerToken FastxParser<T>::getProducerToken_() {
  return moodycamel::ProducerToken(seqContainerQueue_);
}

template <typename T>
moodycamel::ConsumerToken FastxParser<T>::getConsumerToken_() {
  return moodycamel::ConsumerToken(readQueue_);
}

template <typename T> FastxParser<T>::~FastxParser() {
  if (isActive_ or numParsing_ > 0) {
    // Think about if this is too noisy --- but the user really shouldn't do
    // this.
    std::cerr
        << "\n\nEncountered FastxParser destructor while parser was still "
           "marked active (or while parsing threads were still active). "
        << "Be sure to call stop() before letting FastxParser leave scope!\n";
    try {
      stop();
    } catch (const std::exception& e) {
      // Should exiting here be a user-definable behavior?
      // What is the right mechanism for that.
      std::cerr << "\n\nParser encountered exception : " << e.what() << "\n";
      std::exit(-1);
    }
  }
  // Otherwise, we are good to go (i.e., destruct)
}

template <typename T> bool FastxParser<T>::stop() {
  bool ret{false};
  if (isActive_) {
    for (auto& t : parsingThreads_) {
      t->join();
    }
    isActive_ = false;
    for (auto& res : threadResults_) {
      if (res == -3) {
        throw std::range_error("Error reading from the FASTA/Q stream. Make "
                               "sure the file is valid.");
      } else if (res < -1) {
        std::stringstream ss;
        ss << "Error reading from the FASTA/Q stream. Minimum return code for "
              "left and right read was ("
           << res << "). Make sure the file is valid.";
        throw std::range_error(ss.str());
      }
    }
    ret = true;
  } else {
    // Is this being too loud?  Again, if this triggers, the user has violated
    // the API.
    std::cerr << "stop() was called on a FastxParser that was not marked "
                 "active. Did you remember "
              << "to call start() on this parser?\n";
  }
  return ret;
}

template <typename T>
int parse_reads(
    std::vector<std::string>& inputStreams, std::atomic<uint32_t>& numParsing,
    moodycamel::ConsumerToken* cCont, moodycamel::ProducerToken* pRead,
    moodycamel::ConcurrentQueue<uint32_t>& workQueue,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>&
        seqContainerQueue_,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>& readQueue_) {

  using namespace klibpp;
  using fastx_parser::thread_utils::MIN_BACKOFF_ITERS;
  auto curMaxDelay = MIN_BACKOFF_ITERS;
  T* s;

  uint32_t fn{0};
  while (workQueue.try_dequeue(fn)) {
    auto& file = inputStreams[fn];
    std::unique_ptr<ReadChunk<T>> local;
    while (!seqContainerQueue_.try_dequeue(*cCont, local)) {
      fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
      // Think of a way to do this that wouldn't be loud (or would allow a
      // user-definable logging mechanism) std::cerr << "couldn't dequeue read
      // chunk\n";
    }
    size_t numObtained{local->size()};
    // open the file and init the parser
    gzFile fp = gzopen(file.c_str(), "r");

    // we start off with the 0-th fragment in this
    // file.
    uint64_t frag_id{0};
    uint64_t first_frag_of_chunk{frag_id};

    // The number of reads we have in the local vector
    size_t numWaiting{0};

    auto seq = make_kstream(fp, gzread, mode::in);

    s = &((*local)[numWaiting]);
    while (seq >> s->first()) { // ksv >= 0
      frag_id++;
      numWaiting++;
      // If we've filled the local vector, then dump to the concurrent queue
      if (numWaiting == numObtained) {
        curMaxDelay = MIN_BACKOFF_ITERS;
        local->set_chunk_frag_offset(fn, first_frag_of_chunk);
        while (!readQueue_.try_enqueue(std::move(local))) {
          fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
        }
        first_frag_of_chunk = frag_id;
        numWaiting = 0;
        numObtained = 0;
        // And get more empty reads
        curMaxDelay = MIN_BACKOFF_ITERS;
        while (!seqContainerQueue_.try_dequeue(*cCont, local)) {
          fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
        }
        numObtained = local->size();
      }
      s = &((*local)[numWaiting]);
    }

    // if we had an error in the stream
    if (seq.err()) {
      --numParsing;
      return -3;
    } else if (seq.tqs()) {
      // if we had a quality string of the wrong length
      // tqs == truncated quality string
      --numParsing;
      return -2;
    }

    // If we hit the end of the file and have any reads in our local buffer
    // then dump them here.
    if (numWaiting > 0) {
      local->have(numWaiting);
      local->set_chunk_frag_offset(fn, first_frag_of_chunk);
      curMaxDelay = MIN_BACKOFF_ITERS;
      while (!readQueue_.try_enqueue(*pRead, std::move(local))) {
        fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
      }
      numWaiting = 0;
    } else if (numObtained > 0) {
      curMaxDelay = MIN_BACKOFF_ITERS;
      local->set_chunk_frag_offset(fn, first_frag_of_chunk);
      while (!seqContainerQueue_.try_enqueue(std::move(local))) {
        fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
      }
    }
    // destroy the parser and close the file
    gzclose(fp);
  }

  --numParsing;
  return 0;
}

/// make an std::array of size N of kstreams based on the 
/// compile-time template parameter N.
template <size_t N, size_t... Is>
auto make_kstream_array(
    const std::array<gzFile, N>& fptrs,
    std::index_sequence<Is...>) {
  return std::array{klibpp::make_kstream(fptrs[Is], gzread, klibpp::mode::in)...};
}

/// generic function to parse a set (of arity N) of reads such that each parser 
/// handles all N files in one slice of the set.
template <typename T, size_t N>
int parse_read_set_serial(
    std::vector<std::vector<std::string>>& inputStreams,
    std::atomic<uint32_t>& numParsing,
    moodycamel::ConsumerToken* cCont, moodycamel::ProducerToken* pRead,
    moodycamel::ConcurrentQueue<uint32_t>& workQueue,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>&
        seqContainerQueue_,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>& readQueue_) {

  using namespace klibpp;
  using fastx_parser::thread_utils::MIN_BACKOFF_ITERS;
  T* s;

  uint32_t fn{0};
  while (workQueue.try_dequeue(fn)) {

    std::unique_ptr<ReadChunk<T>> local;
    thread_utils::simple_wait([&]() { 
      return seqContainerQueue_.try_dequeue(*cCont, local);
    });
    size_t numObtained{local->size()};

    // open the file and init the parser
    std::array<gzFile, N> fptrs{};
    for (size_t i = 0; i < N; ++i){
      fptrs[i] = gzopen(inputStreams[i][fn].c_str(), "r");
    }

    // we start off with the 0-th fragment in this
    // file.
    uint64_t frag_id{0};
    uint64_t first_frag_of_chunk{frag_id};

    // The number of reads we have in the local vector
    size_t numWaiting{0};

    auto seqs = make_kstream_array(fptrs, std::make_index_sequence<N>{});

    auto get_seqs_from_files = [&](T* s) -> bool {
      for (size_t i = 0; i < N; ++i) {
        if (!(seqs[i] >> (*s)[i])) { return false; }
      }
      return true;
    };

    s = &((*local)[numWaiting]);
    while ( get_seqs_from_files(s) ) { 
      frag_id++;
      numWaiting++;
      // If we've filled the local vector, then dump to the concurrent queue
      if (numWaiting == numObtained) {
        local->set_chunk_frag_offset(fn, first_frag_of_chunk);
        thread_utils::simple_wait([&]() { 
          return readQueue_.try_enqueue(std::move(local));
        });

        first_frag_of_chunk = frag_id;
        numWaiting = 0;
        numObtained = 0;
        // And get more empty reads
        thread_utils::simple_wait([&]() { 
          return seqContainerQueue_.try_dequeue(*cCont, local);
        });
        numObtained = local->size();
      }
      s = &((*local)[numWaiting]);
    }

    // if we had an error in the stream
    bool had_err = false;
    for (size_t i = 0; i < N; ++i) {
      had_err = had_err || seqs[i].err();
      if (had_err) { break; }
    }
    bool had_tqs = false;
    for (size_t i = 0; i < N; ++i) {
      had_tqs = had_tqs || seqs[i].tqs();
      if (had_tqs) { break; }
    }

    if (had_err) {
      --numParsing;
      return -3;
    } else if (had_tqs) {
      // if we had a quality string of the wrong length
      // tqs == truncated quality string
      --numParsing;
      return -2;
    }

    // If we hit the end of the file and have any reads in our local buffer
    // then dump them here.
    if (numWaiting > 0) {
      local->have(numWaiting);
      local->set_chunk_frag_offset(fn, first_frag_of_chunk);

      thread_utils::simple_wait([&]() { 
        return readQueue_.try_enqueue(*pRead, std::move(local));
      });

      numWaiting = 0;
    } else if (numObtained > 0) {
      local->set_chunk_frag_offset(fn, first_frag_of_chunk);
      thread_utils::simple_wait([&]() { 
        return seqContainerQueue_.try_enqueue(std::move(local));
      });
    }
    // destroy the parser and close the file
    for (size_t i = 0; i < N; ++i) {
      gzclose(fptrs[i]);
    }
  }

  --numParsing;
  return 0;
}

// Template member function implementation for parallel parsing
template <typename T>
template <size_t N>
bool FastxParser<T>::start_parallel_parsing_impl() {
  
  if (numParsing_ != 0) {
    return false;
  }

  isActive_ = true;

  size_t numFiles = inputStreamSets_[0].size();
  
  if (inputStreamSets_.size() != N) {
    throw std::logic_error("inputStreamSets_ size doesn't match template arity");
  }

  if (!parallelParsing_) {
    throw std::runtime_error(
        "Multi-file parsing with arity > 1 requires parallelParsing=true");
  }

  size_t numConcurrentFileSets = std::min(static_cast<size_t>(numParsers_), numFiles);
  // size_t totalThreads = numConcurrentFileSets * (N + 1);  // Not currently used
  
  /*
  std::cerr << "Processing " << numFiles << " file sets with " 
            << numConcurrentFileSets << " concurrent producers (" 
            << totalThreads << " total threads)\n";
  */
  // HEAP-ALLOCATE the work queue so it outlives this function
  auto fileWorkQueue = std::make_shared<moodycamel::ConcurrentQueue<uint32_t>>(numFiles);
  for (size_t i = 0; i < numFiles; ++i) {
    fileWorkQueue->enqueue(static_cast<uint32_t>(i));
  }

  threadResults_.resize(numFiles * (N + 1));
  std::fill(threadResults_.begin(), threadResults_.end(), 0);


  for (size_t producerIdx = 0; producerIdx < numConcurrentFileSets; ++producerIdx) {
    ++numParsing_;
    
    // Capture fileWorkQueue by VALUE (it's a shared_ptr, so the copy keeps the queue alive)
    auto processFileSets = [this, fileWorkQueue, producerIdx]() {

      // Allocate queues and recycle queues with simpler lifetime management
      // Use unique_ptr instead of nested shared_ptr to reduce atomic operations
      struct FileSetContext {
        std::array<moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<klibpp::KSeq>>>, N> queues;
        std::array<moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<klibpp::KSeq>>>, N> recycleQueues;
        
        FileSetContext() {
          // Queues are already default-constructed; no need for placement new
        }
      };
      
      auto ctx = std::make_unique<FileSetContext>();

      for (size_t i = 0; i < N; ++i) {
        // Pre-allocate chunks in recycle queue
        constexpr size_t NUM_PREALLOCATED = 4;
        for (size_t j = 0; j < NUM_PREALLOCATED; ++j) {
          ctx->recycleQueues[i].enqueue(
            std::make_unique<ReadChunk<klibpp::KSeq>>(this->blockSize_));
        }
      }

      // Create persistent parser worker threads (one per file in the arity)
      // These threads will process multiple file sets
      struct WorkItem {
        uint32_t fileIdx{0};
        bool done{false};
      };
      
      std::array<moodycamel::ConcurrentQueue<WorkItem>, N> parserWorkQueues;
      std::array<std::atomic<bool>, N> parserBusy;
      for (size_t i = 0; i < N; ++i) {
        parserBusy[i].store(false);
      }
      
      std::vector<std::thread> parserThreads;
      
      for (size_t i = 0; i < N; ++i) {
        parserThreads.emplace_back(
            [this, i, &ctx, &parserWorkQueues, &parserBusy]() {
              WorkItem work;
              while (true) {
                // Wait for work
                thread_utils::simple_wait([&]() { 
                  return parserWorkQueues[i].try_dequeue(work);
                });
                
                if (work.done) {
                  break; // Exit thread
                }
                
                parserBusy[i].store(true, std::memory_order_release);
                const std::string& filename = inputStreamSets_[i][work.fileIdx];
                this->threadResults_[work.fileIdx * (N + 1) + i] = 
                    parse_single_file<klibpp::KSeq>(
                        filename, work.fileIdx,  
                        ctx->queues[i], ctx->recycleQueues[i], this->blockSize_);
                parserBusy[i].store(false, std::memory_order_release);
              }
            });
      }

      // Create persistent assembler thread
      moodycamel::ConcurrentQueue<WorkItem> assemblerWorkQueue;
      std::atomic<bool> assemblerBusy{false};
      size_t tokenIdx = producerIdx;
      std::thread assemblerThread(
          [this, tokenIdx, &ctx, &assemblerWorkQueue, &assemblerBusy]() {
            WorkItem work;
            while (true) {
              // Wait for work
              thread_utils::simple_wait([&]() { 
                return assemblerWorkQueue.try_dequeue(work);
              });
              
              if (work.done) {
                break; // Exit thread
              }
              
              assemblerBusy.store(true, std::memory_order_release);
              this->threadResults_[work.fileIdx * (N + 1) + N] = 
                  thread_utils::assemble_read_set_raw<T, N>(
                      ctx->queues, ctx->recycleQueues,
                      this->consumeContainers_[tokenIdx].get(),
                      this->produceReads_[tokenIdx].get(),
                      this->seqContainerQueue_,
                      this->readQueue_,
                      work.fileIdx);
              assemblerBusy.store(false, std::memory_order_release);
            }
          });

      // Process all file sets using persistent threads
      uint32_t fn{0};
      while (fileWorkQueue->try_dequeue(fn)) {
        // Submit work to all parser threads
        for (size_t i = 0; i < N; ++i) {
          parserWorkQueues[i].enqueue(WorkItem{fn, false});
        }
        
        // Submit work to assembler thread
        assemblerWorkQueue.enqueue(WorkItem{fn, false});
        
        // Wait for all parser threads to complete
        auto all_parsers_idle = [&]() {
          for (size_t i = 0; i < N; ++i) {
            if (parserBusy[i].load(std::memory_order_acquire)) {
              return false;
            }
          }
          return true;
        };
        
        // Wait for assembler to finish processing this file set
        thread_utils::simple_wait([&]() { 
          return !assemblerBusy.load(std::memory_order_acquire) && all_parsers_idle();
        });
      }

      // Signal all threads to exit
      for (size_t i = 0; i < N; ++i) {
        parserWorkQueues[i].enqueue(WorkItem{0, true});
      }
      assemblerWorkQueue.enqueue(WorkItem{0, true});

      // Join all persistent threads
      for (auto& t : parserThreads) {
        t.join();
      }
      assemblerThread.join();

      --numParsing_;
    };

    parsingThreads_.emplace_back(new std::thread(processFileSets));
  }

  return true;
}


template <> bool FastxParser<ReadSeq>::start() {
  if (numParsing_ == 0) {
    isActive_ = true;
    threadResults_.resize(numParsers_);
    std::fill(threadResults_.begin(), threadResults_.end(), 0);
    for (size_t i = 0; i < numParsers_; ++i) {
      ++numParsing_;
      parsingThreads_.emplace_back(new std::thread([this, i]() {
        this->threadResults_[i] = parse_reads(
            this->inputStreamSets_[0], this->numParsing_,
            this->consumeContainers_[i].get(), this->produceReads_[i].get(),
            this->workQueue_, this->seqContainerQueue_, this->readQueue_);
      }));
    }
    return true;
  } else {
    return false;
  }
}

template <> bool FastxParser<ReadPair>::start() {
  if (parallelParsing_ && inputStreamSets_.size() > 0) {
    return start_parallel_parsing_impl<2>();
  } else {
    // Fall back to sequential parsing
    if (numParsing_ == 0) {
      isActive_ = true;
      if (inputStreams_.size() != inputStreams2_.size()) {
        throw std::invalid_argument("There should be the same number of files "
                                    "for the left and right reads");
      }

      threadResults_.resize(numParsers_);
      std::fill(threadResults_.begin(), threadResults_.end(), 0);

      for (size_t i = 0; i < numParsers_; ++i) {
        ++numParsing_;
        parsingThreads_.emplace_back(new std::thread([this, i]() {
          this->threadResults_[i] = parse_read_set_serial<ReadPair, 2>(
              this->inputStreamSets_,
              this->numParsing_, this->consumeContainers_[i].get(),
              this->produceReads_[i].get(), this->workQueue_,
              this->seqContainerQueue_, this->readQueue_);
        }));
      }
      return true;
    }
    return false;
  }
}

template <> bool FastxParser<ReadQualPair>::start() {
  if (parallelParsing_ && inputStreamSets_.size() > 0) {
    return start_parallel_parsing_impl<2>();
  } else {
    // Fall back to sequential
    if (numParsing_ == 0) {
      isActive_ = true;
      threadResults_.resize(numParsers_);
      std::fill(threadResults_.begin(), threadResults_.end(), 0);

      for (size_t i = 0; i < numParsers_; ++i) {
        ++numParsing_;
        parsingThreads_.emplace_back(new std::thread([this, i]() {
          this->threadResults_[i] = parse_read_set_serial<ReadQualPair, 2>(
              this->inputStreamSets_,
              this->numParsing_, this->consumeContainers_[i].get(),
              this->produceReads_[i].get(), this->workQueue_,
              this->seqContainerQueue_, this->readQueue_);
        }));
      }
      return true;
    }
    return false;
  }
}

template <> bool FastxParser<ReadTriple>::start() {
  if (parallelParsing_ && inputStreamSets_.size() > 0) {
    return start_parallel_parsing_impl<3>();
  } else {
    // Fall back to sequential
    if (numParsing_ == 0) {
      isActive_ = true;
      threadResults_.resize(numParsers_);
      std::fill(threadResults_.begin(), threadResults_.end(), 0);

      for (size_t i = 0; i < numParsers_; ++i) {
        ++numParsing_;
        parsingThreads_.emplace_back(new std::thread([this, i]() {
          this->threadResults_[i] = parse_read_set_serial<ReadTriple, 3>(
              this->inputStreamSets_,
              this->numParsing_, this->consumeContainers_[i].get(),
              this->produceReads_[i].get(), this->workQueue_,
              this->seqContainerQueue_, this->readQueue_);
        }));
      }
      return true;
    }
    return false;
  }
}

template <> bool FastxParser<ReadQualTriple>::start() {
  if (parallelParsing_ && inputStreamSets_.size() > 0) {
    return start_parallel_parsing_impl<3>();
  } else {
    // Fall back to sequential
    if (numParsing_ == 0) {
      isActive_ = true;
      threadResults_.resize(numParsers_);
      std::fill(threadResults_.begin(), threadResults_.end(), 0);

      for (size_t i = 0; i < numParsers_; ++i) {
        ++numParsing_;
        parsingThreads_.emplace_back(new std::thread([this, i]() {
          this->threadResults_[i] = parse_read_set_serial<ReadQualTriple, 3>(
            this->inputStreamSets_,
            this->numParsing_, this->consumeContainers_[i].get(),
            this->produceReads_[i].get(), this->workQueue_,
            this->seqContainerQueue_, this->readQueue_);
        }));
      }
      return true;
    }
    return false;
  }
}

template <typename T> bool FastxParser<T>::refill(ReadGroup<T>& seqs) {
  finishedWithGroup(seqs);
  auto curMaxDelay = fastx_parser::thread_utils::MIN_BACKOFF_ITERS;
  while (numParsing_ > 0) {
    if (readQueue_.try_dequeue(seqs.consumerToken(), seqs.chunkPtr())) {
      return true;
    }
    fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
  }
  return readQueue_.try_dequeue(seqs.consumerToken(), seqs.chunkPtr());
}

template <typename T> void FastxParser<T>::finishedWithGroup(ReadGroup<T>& s) {
  // If this read group is holding a valid chunk, then give it back
  if (!s.empty()) {
    seqContainerQueue_.enqueue(s.producerToken(), std::move(s.takeChunkPtr()));
    s.setChunkEmpty();
  }
}

template class FastxParser<ReadSeq>;
template class FastxParser<ReadPair>;
// template class FastxParser<ReadQual>;
template class FastxParser<ReadQualPair>;
template class FastxParser<ReadTriple>;
template class FastxParser<ReadQualTriple>;
} // namespace fastx_parser
