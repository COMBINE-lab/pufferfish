#ifndef _COMMON_H_
#define _COMMON_H_
#define __STDC_LIMIT_MACROS

#include <atomic>
#include <vector>
#include <memory>
#include <cstdint>
#include <climits>
#include <algorithm>
#include <functional>

#include <dnachar.h>
#include "oneapi/tbb/concurrent_queue.h"

namespace TwoPaCo
{

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
	TypeName(const TypeName&); \
	void operator=(const TypeName&)

	using std::max;
	using std::min;

	extern const int64_t INVALID_VERTEX;
	extern const uint32_t MAX_COUNTER;

	struct Task
	{
		bool isFinal;
		uint64_t piece;
		uint64_t start;
		uint64_t seqId;
		uint64_t offset;
		std::string str;
#ifdef _DEBUG
		static const size_t TASK_SIZE = 32;
#else
		static const size_t TASK_SIZE = 1 << 19;
#endif					
		static const size_t GAME_OVER = SIZE_MAX;
		Task() {}
		Task(uint64_t seqId, uint64_t start, uint64_t piece, uint64_t offset, bool isFinal, std::string && str) :
			seqId(seqId), start(start), piece(piece), offset(offset), isFinal(isFinal), str(std::move(str)) {}
	};

	typedef oneapi::tbb::concurrent_bounded_queue<Task> TaskQueue;
	typedef std::unique_ptr<TaskQueue> TaskQueuePtr;
}

#endif
