#ifndef OUR_GFA_READER_H
#define OUR_GFA_READER_H

#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "string_view.hpp"
#include "sparsepp/spp.h"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"

		// For the time being, assume < 4B contigs
		// and that each contig is < 4B bases
		struct Position {
			//std::string transcript_id;
			uint32_t transcript_id;
			uint32_t pos;
			bool orien;

			Position(uint32_t tid, uint32_t tpos, bool torien) {
				transcript_id = tid;
				pos = tpos;
				orien = torien;
			}

      template <class Archive>
      void serialize(Archive& ar) {
        ar(transcript_id, pos, orien);
      }
		};


class PosFinder {
	private:
		std::ifstream file;
		size_t k;
		struct Contig {
			std::string seq;
			std::string id;
		};
    spp::sparse_hash_map<std::string, std::string> contigid2seq;//map of contig_id to # of letters in contig (contig length)
	spp::sparse_hash_map<std::string, std::string> seq2contigid;
		// path maps each transcript_id to a pair of <contig_id, orientation>
		//orientation : +/true main, -/false reverse
		spp::sparse_hash_map<std::string, std::vector< std::pair<std::string, bool> > > path;
		spp::sparse_hash_map<std::string, uint32_t> refIDs;
		
		spp::sparse_hash_map<uint32_t, std::string> refMap;
		// maps each contig to a list of positions in different transcripts
		std::vector<std::pair<std::string, bool> > explode(const stx::string_view str, const char& ch);		
		bool is_number(const std::string& s);

		// Avoiding un-necessary stream creation + replacing strings with string view
		// is a bit > than a 2x win!
		// implementation from : https://marcoarena.wordpress.com/tag/string_view/
		std::vector<stx::string_view> split(stx::string_view str, char delims);
	public:
		spp::sparse_hash_map<std::string, std::vector<Position> > contig2pos;  
		PosFinder(const char* gfaFileName, size_t input_k);
    spp::sparse_hash_map<std::string, std::string>& getContigNameMap();
    spp::sparse_hash_map<std::string, std::string>& getContigIDMap();
    spp::sparse_hash_map<uint32_t, std::string>& getRefIDs();

		void parseFile();	
		void mapContig2Pos();
    void serializeContigTable(const std::string& ofile);
		void deserializeContigTable();
};

#endif

