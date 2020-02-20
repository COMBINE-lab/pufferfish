#include <set>
#include <deque>
#include <ctime>
#include <memory>
#include <bitset>
#include <numeric>
#include <cassert>
#include <fstream>
#include <iostream>

#include "vertexenumerator.h"

namespace TwoPaCo
{
	namespace
	{
		template<size_t CAPACITY>
		std::unique_ptr<VertexEnumerator> CreateEnumeratorImpl(const std::vector<std::string> & fileName,
			size_t vertexLength,
			size_t filterSize,
			size_t hashFunctions,
			size_t rounds,
			size_t threads,
			size_t abundance,
			const std::string & tmpFileName,
			const std::string & outFileName,
			std::ostream & logStream)
		{
			size_t neededCapacity = CalculateNeededCapacity(vertexLength);
			if (CAPACITY == neededCapacity)
			{
				return std::unique_ptr<VertexEnumerator>(new VertexEnumeratorImpl<CAPACITY>(fileName,
					vertexLength,
					filterSize,
					hashFunctions,
					rounds,
					threads,
					abundance,
					tmpFileName,
					outFileName,
					logStream));
			}
			
			return CreateEnumeratorImpl<CAPACITY + 1>(fileName,
				vertexLength,
				filterSize,
				hashFunctions,
				rounds,
				threads,
				abundance,
				tmpFileName,
				outFileName,
				logStream);
		}

		template<>
		std::unique_ptr<VertexEnumerator> CreateEnumeratorImpl<MAX_CAPACITY>(const std::vector<std::string> & fileName,
			size_t vertexLength,
			size_t filterSize,
			size_t hashFunctions,
			size_t rounds,
			size_t threads,
			size_t abundance,
			const std::string & tmpFileName,
			const std::string & outFileName,
			std::ostream & logStream)
		{
			throw std::runtime_error("The value of K is too big. Please refer to documentaion how to increase the max supported value of K.");
			return 0;
		}		
	}

	std::unique_ptr<VertexEnumerator> CreateEnumerator(const std::vector<std::string> & fileName,
		size_t vertexLength,
		size_t filterSize,
		size_t hashFunctions,
		size_t rounds,
		size_t threads,
		size_t abundance,
		const std::string & tmpFileName,
		const std::string & outFileName,		
		std::ostream & logStream)
	{
		return CreateEnumeratorImpl<1>(fileName,
			vertexLength,
			filterSize,
			hashFunctions,
			rounds,
			threads,
			abundance,
			tmpFileName,
			outFileName,
			logStream);
	}
}