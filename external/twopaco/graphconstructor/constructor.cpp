#include <set>
#include <ctime>
#include <string>
#include <vector>
#include <memory>
#include <cassert>
#include <cstdint>
#include <cassert>
#include <fstream>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <sstream>

#include <tclap/CmdLine.h>

#include "tbb/scalable_allocator.h"
#include "test.h"
#include "assemblyedgeconstructor.h"

size_t Atoi(const char * str)
{
	size_t ret;
	std::stringstream ss(str);
	ss >> ret;
	return ret;
}

class OddConstraint : public TCLAP::Constraint <unsigned int>
{
public:
	~OddConstraint()
	{

	}

	std::string description() const
	{
		return "value of K must be odd";
	}

	std::string shortID() const
	{
		return "oddc";
	}

	bool check(const unsigned & value) const
	{
		return (value % 2) == 1;
	}
};


int buildGraphMain(std::vector<std::string>& args)  
{
	OddConstraint constraint;
	try
	{
		TCLAP::CmdLine cmd("Program for construction of the condensed de Bruijn graph from complete genomes", ' ', "0.9.3");

		TCLAP::ValueArg<unsigned int> kvalue("k",
			"kvalue",
			"Value of k",
			false,
			25,
			&constraint,
			cmd);

		TCLAP::ValueArg<unsigned int> filterSize("f",
			"filtersize",
			"Size of the filter",
			true,
			32,
			"integer");

		TCLAP::ValueArg<double> filterMemory("",
                        "filtermemory",
                        "Memory in GBs allocated for the filter",
                        true,
                        4,
                        "float");

		TCLAP::ValueArg<unsigned int> hashFunctions("q",
			"hashfnumber",
			"Number of hash functions",
			false,
			5,
			"integer",
			cmd);

		TCLAP::ValueArg<unsigned int> rounds("r",
			"rounds",
			"Number of computation rounds",
			false,
			1,
			"integer",
			cmd);

		TCLAP::ValueArg<unsigned int> threads("t",
			"threads",
			"Number of worker threads",
			false,
			1,
			"integer",
			cmd);

		TCLAP::ValueArg<size_t> abundance("a",
			"abundance",
			"Vertex abundance threshold",
			false,
			UINT64_MAX,
			"integer",
			cmd);

		TCLAP::ValueArg<std::string> tmpDirName("",
			"tmpdir",
			"Temporary directory name",
			false,
			".",
			"directory name",
			cmd);
		
		TCLAP::SwitchArg runTests("",
			"test",
			"Run tests",				
			cmd);

		TCLAP::UnlabeledMultiArg<std::string> fileName("filenames",
			"FASTA file(s) with nucleotide sequences.",
			true,
			"fasta files with genomes",
			cmd);

		TCLAP::ValueArg<std::string> outFileName("o",
			"outfile",
			"Output file name prefix",
			false,
			"de_bruijn.bin",
			"file name",
			cmd);

		cmd.xorAdd(filterSize, filterMemory);
		cmd.parse(args);
		using TwoPaCo::Range;
		if (runTests.getValue()) {
			TwoPaCo::RunTests(10, 20, 9000, 6, Range(3, 11), Range(1, 2), Range(1, 5), Range(4, 5), 0.05, 0.1, tmpDirName.getValue());
			return 0;
		}

		int64_t filterBits = 1;
		if (filterSize.isSet()) {
			filterBits = filterSize.getValue();
		} else {
			filterBits = log2(filterMemory.getValue() * 8e+9);
		}

		std::unique_ptr<TwoPaCo::VertexEnumerator> vid = TwoPaCo::CreateEnumerator(fileName.getValue(),
			kvalue.getValue(),
			filterBits,
			hashFunctions.getValue(),
			rounds.getValue(),
			threads.getValue(),
			abundance.getValue(),
			tmpDirName.getValue(),
			outFileName.getValue(),
			std::cout);

		if (vid)
		{
			std::cout << "Distinct junctions = " << vid->GetVerticesCount() << std::endl;
			std::cout << std::endl;
		}

	}
	catch (TCLAP::ArgException & e)
	{
		std::cerr << std::endl << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	catch (std::runtime_error & e)
	{
		std::cerr << std::endl << "Error: " << e.what() << std::endl;
		return 1;
	}
    scalable_allocation_command(TBBMALLOC_CLEAN_ALL_BUFFERS, 0);
	scalable_allocation_command(TBBMALLOC_CLEAN_THREAD_BUFFERS, 0);
	return 0;
}