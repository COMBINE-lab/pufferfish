#include <set>
#include <map>
#include <deque>
#include <string> 
#include <vector>
#include <fstream>
#include <cstdint>
#include <iostream>
#include <exception>

#include <tclap/CmdLine.h>

#include <tpie/tpie.h>
#include <tpie/sort.h>
#include <tpie/file_stream.h>
#include <tpie/dummy_progress.h>

#include <streamfastaparser.h>
#include <junctionapi/junctionapi.h>

const size_t MAX_K = 128;

template <size_t k>
struct NaiveGraphConstructor
{
public:
	typedef uint32_t DnaChar;
	typedef std::vector<uint32_t> DnaString;
	DnaChar currentUnknownChar;

	NaiveGraphConstructor() : currentUnknownChar(4)
	{
	}

	struct KMerOccurence
	{
		static const size_t SIZE = k + 2;
		DnaChar body[SIZE];

		KMerOccurence()
		{

		}

		KMerOccurence(DnaString::const_iterator it)
		{
			std::copy(it, it + SIZE, body);
		}
	};

	struct KMerOccurenceComparer
	{
	public:
		bool operator()(const KMerOccurence & a, const KMerOccurence & b) const
		{
			const size_t OFF = 1;
			return std::lexicographical_compare(a.body + OFF, a.body + OFF + k, b.body + OFF, b.body + OFF + k);
		}
	};

	uint32_t MakeUpChar(char ch)
	{
		if (!TwoPaCo::DnaChar::IsDefinite(ch))
		{
			if (currentUnknownChar == UINT32_MAX)
			{
				throw std::runtime_error("Too many unknown characters");
			}

			return currentUnknownChar++;
		}

		return TwoPaCo::DnaChar::MakeUpChar(ch);
	}

	uint32_t ReverseMadeUpChar(uint32_t ch)
	{
		char before = TwoPaCo::DnaChar::UnMakeUpChar(ch);
		return MakeUpChar(TwoPaCo::DnaChar::ReverseChar(before));
	}

	template<class It>
	DnaString GetKMer(It pos)
	{
		DnaString ret;
		std::copy(pos, pos + k, ret.body);
		return ret;
	}

	void InsertJunction(std::map<DnaString, uint64_t> & junctionMap, const DnaString & newJunction)
	{
		if (junctionMap.count(newJunction) == 0)
		{
			size_t id = junctionMap.size();
			junctionMap[newJunction] = id;
		}
	}

	bool OccurenceEqual(const KMerOccurence & a, const KMerOccurence & b)
	{
		const size_t OFF = 1;
		return std::equal(a.body + OFF, a.body + OFF + k, b.body + OFF);
	}

	DnaString SubStr(const DnaString & str, size_t pos)
	{
		return DnaString(str.begin() + pos, str.begin() + pos + k);
	}

	void MakeDeBruijnGraph(const std::vector<std::string> & fileName, const std::string outFileName, const std::string & tmpDirName)
	{		
		std::vector<std::vector<DnaString> > strand(2);
		std::cout << "Parsing input..." << std::endl;
		for (auto name : fileName)
		{
			for (TwoPaCo::StreamFastaParser parser(name); parser.ReadRecord();)
			{
				strand[0].push_back(DnaString());
				strand[1].push_back(DnaString());
				for (char ch; parser.GetChar(ch);)
				{
					strand[0].back().push_back(MakeUpChar(ch));
				}

				for (auto it = strand[0].back().rbegin(); it != strand[0].back().rend(); ++it)
				{
					strand[1].back().push_back(ReverseMadeUpChar(*it));
				}
			}
		}

		uint64_t records = 0;
		std::cout << "Storing (k + 2)-mers..." << std::endl;
		std::map<DnaString, uint64_t> junctionMap;
		tpie::tpie_init(tpie::ALL);
		tpie::file_stream<KMerOccurence> tmpFile;
		std::string tmpFileName = tmpDirName + "/tmp.bin";
		tmpFile.open(tmpFileName.c_str(), tpie::access_write);
		for (auto currentStrand : strand)
		{
			for (auto str : currentStrand)
			{
				if (str.size() >= k)
				{
					InsertJunction(junctionMap, SubStr(str, 0));
					InsertJunction(junctionMap, SubStr(str, str.size() - k));
				}

				for (size_t pos = 0; pos + k + 2 <= str.size(); ++pos)
				{
					++records;
					tmpFile.write(KMerOccurence(str.begin() + pos));
				}
			}
		}

		std::cout << "Sorting..." << std::endl;
		tmpFile.close();
		tmpFile.open(tmpFileName.c_str());
		tpie::progress_indicator_null pi;
		tpie::sort(tmpFile, KMerOccurenceComparer(), pi);
		tmpFile.seek(0);

		KMerOccurence it;
		std::cout << "Checking junctions..." << std::endl;
		for (uint64_t i = 0; i < records;)
		{
			if (i == 0)
			{
				++i;
				it = tmpFile.read();
			}

			KMerOccurence jt = it;
			std::set<DnaChar> inGoing;
			std::set<DnaChar> outGoing;
			for (; OccurenceEqual(it, jt);)
			{
				inGoing.insert(jt.body[0]);
				outGoing.insert(jt.body[k + 1]);
				if (i < records)
				{
					++i;
					jt = tmpFile.read();
				}
				else
				{
					break;
				}
			}

			if (inGoing.size() > 1 || outGoing.size() > 1)
			{
				const size_t OFF = 1;
				InsertJunction(junctionMap, DnaString(it.body + OFF, it.body + OFF + k));
			}

			it = jt;
		}
		
		size_t occurences = 0;
		std::cout << "Generating edges..." << std::endl;
		TwoPaCo::JunctionPositionWriter posWriter(outFileName);
		for (size_t chr = 0; chr < strand[0].size(); chr++)
		{
			const DnaString & str = strand[0][chr];
			if (str.size() >= k)
			{
				for (size_t i = 0; i <= str.size() - k; i++)
				{
 					DnaString kmer = SubStr(str, i);
					auto it = junctionMap.find(kmer);
					if (it != junctionMap.end())
					{
						++occurences;
 						posWriter.WriteJunction(TwoPaCo::JunctionPosition(chr, i, it->second));
					}
				}
			}			
		}

		
		std::cout << "Vertices: " << junctionMap.size() << std::endl;
		std::cout << "Occurences: " << occurences << std::endl;
	}
};

template<size_t NOW_K>
void MakeDeBruijnGraph(const std::vector<std::string> & fileName, const std::string outFileName, const std::string & tmpDirName, size_t k)
{
	if (NOW_K == k)
	{
		NaiveGraphConstructor<NOW_K> ctor;
		ctor.MakeDeBruijnGraph(fileName, outFileName, tmpDirName);
	}

	MakeDeBruijnGraph<NOW_K + 1>(fileName, outFileName, tmpDirName, k);
}

template<>
void MakeDeBruijnGraph<MAX_K>(const std::vector<std::string> & fileName, const std::string outFileName, const std::string & tmpDirName, size_t k)
{

}

int main(int argc, char * argv[])
{	
	try
	{
		TCLAP::CmdLine cmd("A really naive program for condensed de Bruijn graph construction", ' ', "0");

		TCLAP::UnlabeledMultiArg<std::string> fileName("filenames",
			"FASTA file(s) with nucleotide sequences.",
			true,
			"fasta files with genomes",
			cmd);

		TCLAP::ValueArg<std::string> outFileName("o",
			"outfile",
			"Output file name prefix",
			false,
			"de_bruijn_slow",
			"file name",
			cmd);

		TCLAP::ValueArg<std::string> tmpFileName("t",
			"tmpdir",
			"Temporary directory name",
			false,
			"./",
			"file name",
			cmd);

		TCLAP::ValueArg<unsigned int> kvalue("k",
			"kvalue",
			"Value of k",
			false,
			25,
			"integer",
			cmd);

		cmd.parse(argc, argv);
		MakeDeBruijnGraph<1>(fileName.getValue(), outFileName.getValue(), tmpFileName.getValue(), kvalue.getValue());

	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	catch (std::runtime_error & e)
	{
		std::cerr << "error: " << e.what() << std::endl;
		return 1;
	}

	return 0;

}