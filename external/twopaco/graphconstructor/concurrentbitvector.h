#ifndef _CONCURRENT_BIT_VECTOR_
#define _CONCURRENT_BIT_VECTOR_

#include <cstdlib>
#include <vector>
#include <atomic>

#include "common.h"

namespace TwoPaCo
{
	class ConcurrentBitVector
	{
	public:
		~ConcurrentBitVector();
		ConcurrentBitVector(size_t size);
		void Reset(size_t size);
		size_t Size() const;
		bool GetBit(size_t idx) const;
		void SetBitConcurrently(size_t idx);
		void MergeOr(const ConcurrentBitVector & mask, size_t size);
		void WriteToFile(const std::string & fileName) const;
		void WriteToFile(std::ofstream & maskStorage, size_t offset, size_t size) const;
		void ReadFromFile(const std::string & fileName, size_t size, bool cleanUp);
		void ReadFromFile(std::ifstream & maskStorage, size_t offset, size_t size) const;
	private:
		DISALLOW_COPY_AND_ASSIGN(ConcurrentBitVector);
		static const size_t SUCCESS = -1;
		typedef uint32_t BASIC_TYPE;
		typedef std::atomic<BASIC_TYPE> UInt;
		size_t size_;
		size_t realSize_;
		UInt * filter_;
		void GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const;
	};	

}

#endif