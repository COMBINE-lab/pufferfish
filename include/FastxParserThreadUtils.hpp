#ifndef FASTX_PARSER_THREAD_UTILS_HPP
#define FASTX_PARSER_THREAD_UTILS_HPP

#include <cassert>
#include <chrono>
#include <pthread.h>
#include <random>
#include <thread>

#if defined(__SSE2__)
 #if defined(HAVE_SIMDE)
  #include "simde/x86/sse2.h"
 #else
  #include <emmintrin.h>
 #endif
#endif

// Most of this code is taken directly from
// https://github.com/geidav/spinlocks-bench/blob/master/os.hpp. However, things
// may be renamed, modified, or randomly mangled over time.
#define ALWAYS_INLINE inline __attribute__((__always_inline__))

namespace fastx_parser {
namespace thread_utils {

static const constexpr size_t MIN_BACKOFF_ITERS = 32;
static const size_t MAX_BACKOFF_ITERS = 1024;

ALWAYS_INLINE static void cpuRelax() {
#if defined(__SSE2__)  // AMD and Intel
  #if defined(HAVE_SIMDE)
    simde_mm_pause();
  #else
    _mm_pause();
  #endif
#elif defined(__i386__) || defined(__x86_64__)
  asm volatile("pause");
#elif defined(__aarch64__)
  asm volatile("isb");
#elif defined(__armel__) || defined(__ARMEL__)
  asm volatile ("nop" ::: "memory");
#elif defined(__arm__) 
  __asm__ __volatile__ ("yield" ::: "memory");
#elif defined(__ia64__)  // IA64
  __asm__ __volatile__ ("hint @pause");
#elif defined(__powerpc__) || defined(__ppc__) || defined(__PPC__)
   __asm__ __volatile__ ("or 27,27,27" ::: "memory");
#else  // everything else.
   asm volatile ("nop" ::: "memory");
#endif
}

ALWAYS_INLINE void yieldSleep() {
  using namespace std::chrono;
  std::chrono::microseconds ytime(20);
  std::this_thread::sleep_for(ytime);
}

ALWAYS_INLINE void backoffExp(size_t& curMaxIters) {
  thread_local std::uniform_int_distribution<size_t> dist;

  // see : https://github.com/coryan/google-cloud-cpp-common/blob/a6e7b6b362d72451d6dc1fec5bc7643693dbea96/google/cloud/internal/random.cc
  #if defined(__linux) && defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
    thread_local std::random_device rd("/dev/urandom");
  #else
    thread_local std::random_device rd;
  #endif  // defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128

  thread_local std::minstd_rand gen(rd());
  const size_t spinIters =
      dist(gen, decltype(dist)::param_type{0, curMaxIters});
  curMaxIters = std::min(2 * curMaxIters, MAX_BACKOFF_ITERS);
  for (size_t i = 0; i < spinIters; i++) {
    cpuRelax();
  }
}

ALWAYS_INLINE void backoffOrYield(size_t& curMaxDelay) {
  if (curMaxDelay >= MAX_BACKOFF_ITERS) {
    yieldSleep();
    curMaxDelay = MIN_BACKOFF_ITERS;
  }
  backoffExp(curMaxDelay);
}

// Dead simple - inline everything
template<typename Func>
ALWAYS_INLINE void simple_wait(Func&& try_op) {
    auto curMaxDelay = MIN_BACKOFF_ITERS;
    while(!try_op()) {
      backoffOrYield(curMaxDelay);
    }
    return;
    /*
    // Try a few times with just a pause
    for (int i = 0; i < 32; ++i) {
        if (try_op()) return;
        cpuRelax();
    }
    
    // If that didn't work, yield to scheduler
    // (this is likely I/O bound, not contention)
    while (!try_op()) {
        std::this_thread::yield();
    }
    */
}

template <typename T, size_t N>
int assemble_read_set(
    std::array<std::shared_ptr<moodycamel::ConcurrentQueue<
        std::unique_ptr<ReadChunk<klibpp::KSeq>>>>, N>& queues,
    std::array<std::shared_ptr<moodycamel::ConcurrentQueue<
        std::unique_ptr<ReadChunk<klibpp::KSeq>>>>, N>& recycleQueues,
    moodycamel::ConsumerToken* cCont,
    moodycamel::ProducerToken* pRead,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>& seqContainerQueue,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>& readQueue,
    uint32_t file_idx) {  // Changed from numAssembling

  std::array<std::unique_ptr<ReadChunk<klibpp::KSeq>>, N> chunks;
  std::array<size_t, N> indices{};
  std::array<bool, N> fileDone{};
  
  // get initial output chunk
  std::unique_ptr<ReadChunk<T>> local;
   
  bool got_chunk = seqContainerQueue.try_dequeue(*cCont, local);
  if (!got_chunk) {
    thread_utils::simple_wait([&]() { 
      return seqContainerQueue.try_dequeue(*cCont, local);
    });
  }
  
  size_t numObtained = local->size();
  size_t numWaiting = 0;
  uint64_t gathered_count = 0;

  // Lambda to fetch chunk from a specific queue
  auto fetch_chunk = [&](size_t idx) -> bool {
    if (chunks[idx] && indices[idx] < chunks[idx]->size())
      return true;
    if (fileDone[idx])
      return false;

    std::unique_ptr<ReadChunk<klibpp::KSeq>> next_chunk;
    // Access the shared_ptr at queues[idx] and dereference it
    if (queues[idx]->try_dequeue(next_chunk)) {
      if (next_chunk == nullptr) {
        fileDone[idx] = true;
        if (chunks[idx])
          recycleQueues[idx]->enqueue(std::move(chunks[idx]));
        chunks[idx] = nullptr;
        return false;
      }
      if (chunks[idx]) {
        recycleQueues[idx]->enqueue(std::move(chunks[idx]));
      }
      chunks[idx] = std::move(next_chunk);
      indices[idx] = 0;
      return true;
    }
    return false;
  };

  // Check if all files are done
  auto all_done = [&]() {
    for (size_t i = 0; i < N; ++i) {
      if (!fileDone[i] || chunks[i])
        return false;
    }
    return true;
  };

  while (!all_done()) {
    // Try to fetch from all queues to update their done status
    std::array<bool, N> haveData;
    for (size_t i = 0; i < N; ++i) {
      haveData[i] = fetch_chunk(i);
    }

    // Check if all queues have data
    bool allHaveData = true;
    for (size_t i = 0; i < N; ++i) {
      if (!haveData[i]) {
        allHaveData = false;
        break;
      }
    }

    if (allHaveData) {
      /*
       // Ensure that the ranks in each chunk match
      size_t first_rank = std::numeric_limits<size_t>::max();
      for (size_t i = 0; i < N; ++i) {
        if (first_rank == std::numeric_limits<size_t>::max()) {
          first_rank = chunks[i]->chunk_frag_offset().frag_idx;
        }
        if (chunks[i]->chunk_frag_offset().frag_idx != first_rank) {
          std::cerr << "[ERROR]: Rank of first chunk in this set was " << first_rank << ", but part " << i << " has rank " << chunks[i]->chunk_frag_offset().frag_idx << "\n";
        }
      }
      */
      // Assemble N-tuple
      T& readSet = (*local)[numWaiting];
      for (size_t i = 0; i < N; ++i) {
        readSet[i] = std::move((*chunks[i])[indices[i]++]);
      }
      ++numWaiting;
      ++gathered_count;

      if (numWaiting == numObtained) {
        local->have(numWaiting);
        local->set_chunk_frag_offset(file_idx, gathered_count - numWaiting);

        thread_utils::simple_wait([&]() { 
          return readQueue.try_enqueue(*pRead, std::move(local));
        });

        // get next output chunk
        numWaiting = 0;

        thread_utils::simple_wait([&]() { 
          return seqContainerQueue.try_dequeue(*cCont, local);
        });
        numObtained = local->size();
      }
    } else {
      // not all queues have data, but check if we're done before backing off
      if (all_done()) {
        break;
      }
    }

  }

  // flush remaining
  if (numWaiting > 0) {
    local->have(numWaiting);
    local->set_chunk_frag_offset(file_idx, gathered_count - numWaiting);

    thread_utils::simple_wait([&]() { 
      return readQueue.try_enqueue(*pRead, std::move(local));
    });
  } else {
    thread_utils::simple_wait([&]() { 
      return seqContainerQueue.try_enqueue(std::move(local));
    });
  }

  return 0;
}

// Version that takes raw array references instead of shared_ptr arrays
template <typename T, size_t N>
int assemble_read_set_raw(
    std::array<moodycamel::ConcurrentQueue<
        std::unique_ptr<ReadChunk<klibpp::KSeq>>>, N>& queues,
    std::array<moodycamel::ConcurrentQueue<
        std::unique_ptr<ReadChunk<klibpp::KSeq>>>, N>& recycleQueues,
    moodycamel::ConsumerToken* cCont,
    moodycamel::ProducerToken* pRead,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>& seqContainerQueue,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>& readQueue,
    uint32_t file_idx) {

  std::array<std::unique_ptr<ReadChunk<klibpp::KSeq>>, N> chunks;
  std::array<size_t, N> indices{};
  std::array<bool, N> fileDone{};
  
  // get initial output chunk
  std::unique_ptr<ReadChunk<T>> local;
   
  bool got_chunk = seqContainerQueue.try_dequeue(*cCont, local);
  if (!got_chunk) {
    thread_utils::simple_wait([&]() { 
      return seqContainerQueue.try_dequeue(*cCont, local);
    });
  }
  
  size_t numObtained = local->size();
  size_t numWaiting = 0;
  uint64_t gathered_count = 0;

  // Lambda to fetch chunk from a specific queue
  auto fetch_chunk = [&](size_t idx) -> bool {
    if (chunks[idx] && indices[idx] < chunks[idx]->size())
      return true;
    if (fileDone[idx])
      return false;

    std::unique_ptr<ReadChunk<klibpp::KSeq>> next_chunk;
    // Access queue directly (not through shared_ptr)
    if (queues[idx].try_dequeue(next_chunk)) {
      if (next_chunk == nullptr) {
        fileDone[idx] = true;
        if (chunks[idx])
          recycleQueues[idx].enqueue(std::move(chunks[idx]));
        chunks[idx] = nullptr;
        return false;
      }
      if (chunks[idx]) {
        recycleQueues[idx].enqueue(std::move(chunks[idx]));
      }
      chunks[idx] = std::move(next_chunk);
      indices[idx] = 0;
      return true;
    }
    return false;
  };

  // Check if all files are done
  auto all_done = [&]() {
    for (size_t i = 0; i < N; ++i) {
      if (!fileDone[i] || chunks[i])
        return false;
    }
    return true;
  };

  while (!all_done()) {
    // Try to fetch from all queues to update their done status
    std::array<bool, N> haveData;
    for (size_t i = 0; i < N; ++i) {
      haveData[i] = fetch_chunk(i);
    }

    // Check if all queues have data
    bool allHaveData = true;
    for (size_t i = 0; i < N; ++i) {
      if (!haveData[i]) {
        allHaveData = false;
        break;
      }
    }

    if (allHaveData) {
      // Assemble N-tuple
      T& readSet = (*local)[numWaiting];
      for (size_t i = 0; i < N; ++i) {
        readSet[i] = std::move((*chunks[i])[indices[i]++]);
      }
      ++numWaiting;
      ++gathered_count;

      if (numWaiting == numObtained) {
        local->have(numWaiting);
        local->set_chunk_frag_offset(file_idx, gathered_count - numWaiting);

        thread_utils::simple_wait([&]() { 
          return readQueue.try_enqueue(*pRead, std::move(local));
        });

        // get next output chunk
        numWaiting = 0;

        thread_utils::simple_wait([&]() { 
          return seqContainerQueue.try_dequeue(*cCont, local);
        });
        numObtained = local->size();
      }
    } else {
      // not all queues have data, but check if we're done before backing off
      if (all_done()) {
        break;
      }
    }

  }

  // flush remaining
  if (numWaiting > 0) {
    local->have(numWaiting);
    local->set_chunk_frag_offset(file_idx, gathered_count - numWaiting);

    thread_utils::simple_wait([&]() { 
      return readQueue.try_enqueue(*pRead, std::move(local));
    });
  } else {
    thread_utils::simple_wait([&]() { 
      return seqContainerQueue.try_enqueue(std::move(local));
    });
  }

  return 0;
}

} // namespace thread_utils
} // namespace fastx_parser

#endif // FASTX_PARSER_THREAD_UTILS_HPP
