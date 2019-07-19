#ifndef _TEST_H_
#define _TEST_H_

#include <vector>
#include <string>
#include <iostream>

namespace TwoPaCo
{
	typedef std::pair<size_t, size_t> Range;	
	bool RunTests(size_t tests, size_t filterBits, size_t length, size_t chrNumber, Range vertexSize, Range hashFunctions, Range rounds, Range threads, double changeRate, double indelRate, const std::string & temporaryDir);
}

#endif