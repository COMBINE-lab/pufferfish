#ifndef OUR_GFA_READER_H
#define OUR_GFA_READER_H

#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <memory>
#include "string_view.hpp"
#include "sparsepp/spp.h"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"
#include "Util.hpp"
#include "FatPufferGraph.hpp"
#include "zstr/zstr.hpp"

class GFAConverter {
	private:
	  	std::unique_ptr<zstr::ifstream> file;
		size_t k;
		struct Contig {
			std::string seq;
			std::string id;
		};
		spp::sparse_hash_map<std::string, std::pair<std::string, std::vector< std::pair<std::string, bool> > > > new2seqAoldids;
		spp::sparse_hash_map<std::string, std::vector<std::pair<std::string, bool> > > old2newids;
		// path maps each transcript_id to a pair of <contig_id, orientation>
		//orientation : +/true main, -/false reverse
		spp::sparse_hash_map<std::string, std::vector< std::pair<std::string, bool> > > path;
		
		std::map<std::pair<std::string, bool>, bool, util::cmpByPair> pathStart ;
		std::map<std::pair<std::string, bool>, bool, util::cmpByPair> pathEnd ;

		pufg::Graph semiCG;

	public:
		GFAConverter(const char* gfaFileName, size_t input_k);
		void parseFile();
		void processContigSeq(std::string & contigId, std::string & contigSeq, spp::sparse_hash_map<std::string, std::string> & seq2newid);
};

#endif

