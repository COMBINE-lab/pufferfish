#include <cstdio>
#include <fstream>
#include <cassert>
#include <iostream>
#include <stdexcept>

#include "concurrentbitvector.h"

namespace TwoPaCo
{
	ConcurrentBitVector::ConcurrentBitVector(size_t size)
		: size_(size), realSize_(size / 32 + 1), filter_(new UInt[realSize_])
	{
		Reset(size);
	}

	void ConcurrentBitVector::Reset(size_t size)
	{
		size_t end = size / 32 + 1;
		for (size_t i = 0; i < end; i++)
		{
			filter_[i] = 0;
		}
	}

	size_t ConcurrentBitVector::Size() const
	{
		return size_;
	}

	void ConcurrentBitVector::SetBitConcurrently(size_t idx)
	{
		uint64_t bit;
		uint64_t element;
		GetCoord(idx, element, bit);
		filter_[element].fetch_or(uint32_t(1) << uint32_t(bit));
	}

	bool ConcurrentBitVector::GetBit(size_t idx) const
	{
		uint64_t bit;
		uint64_t element;
		GetCoord(idx, element, bit);
		return (filter_[element] & (uint32_t(1) << uint32_t(bit))) != 0;
	}

	void ConcurrentBitVector::GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const
	{
		bit = idx & ((uint32_t(1) << uint64_t(5)) - 1);
		element = idx >> 5;
		assert(element < size_ / 32 + 1);
	}

	ConcurrentBitVector::~ConcurrentBitVector()
	{
		delete[] filter_;
	}

	void ConcurrentBitVector::WriteToFile(std::ofstream & maskStorage, size_t offset, size_t size) const
	{
		maskStorage.seekp(offset * sizeof(BASIC_TYPE));
		size_t nowSize = size / 32 + 1;
		if (!maskStorage.write(reinterpret_cast<const char*>(filter_), nowSize * sizeof(BASIC_TYPE)))
		{
			throw std::runtime_error("Can't write to a temporary file");
		}
	}

	void ConcurrentBitVector::WriteToFile(const std::string & fileName) const
	{
		std::ofstream candidateMaskFile(fileName.c_str(), std::ios::binary);
		if (!candidateMaskFile)
		{
			throw std::runtime_error("Can't open a temporary file");
		}

		if (!candidateMaskFile.write(reinterpret_cast<const char*>(filter_), realSize_ * sizeof(BASIC_TYPE)))
		{
			throw std::runtime_error("Can't write to a temporary file");
		}
	}

	void ConcurrentBitVector::ReadFromFile(std::ifstream & maskStorage, size_t offset, size_t size) const
	{
		maskStorage.seekg(offset * sizeof(BASIC_TYPE));
		size_t nowSize = size / 32 + 1;
		if (!maskStorage.read(reinterpret_cast<char*>(filter_), nowSize * sizeof(BASIC_TYPE)))
		{
			throw std::runtime_error("Can't read from a temporary file");
		}
	}

	void ConcurrentBitVector::ReadFromFile(const std::string & fileName, size_t size, bool cleanUp)
	{
		{
			size_t nowSize = size / 32 + 1;
			std::ifstream candidateMaskFile(fileName.c_str(), std::ios::binary);
			if (!candidateMaskFile)
			{
				throw std::runtime_error("Can't open a temporary file");
			}

			if (!candidateMaskFile.read(reinterpret_cast<char*>(filter_), nowSize * sizeof(BASIC_TYPE)))
			{
				throw std::runtime_error("Can't read from a temporary file");
			}
		}

		if (cleanUp)
		{
			std::remove(fileName.c_str());
		}
	}

	void ConcurrentBitVector::MergeOr(const ConcurrentBitVector & mask, size_t size)
	{
		size_t nowSize = size / 32 + 1;
		for (size_t i = 0; i < nowSize; i++)
		{
			filter_[i] |= mask.filter_[i];
		}
	}
}
