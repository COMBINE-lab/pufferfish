#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "string_view.hpp"
#include "sparsepp/spp.h"

class PosFinder {
	private:
		std::ifstream file;
		size_t k;
		struct Contig {
			std::string seq;
			std::string id;
		};
		spp::sparse_hash_map<std::string, uint64_t> contigid2cnt;//std::map<std::string, long> contigid2cnt;
		//orientation : + main, - reverse
		//std::map<std::pair<long, long>, short> link;
		//std::map<std::string, std::vector< std::pair<std::string, bool> > > path;
		spp::sparse_hash_map<std::string, std::vector< std::pair<std::string, bool> > > path;
		spp::sparse_hash_map<std::string, uint32_t> refIDs;

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
		};

		//std::vector<std::

		// maps each contig to a list of positions in different transcripts

		std::vector<std::pair<std::string, bool> > explode(const stx::string_view str, const char& ch) {
				std::string next;
				std::vector< std::pair<std::string, bool> > result;
				// For each character in the string
				for (auto it = str.begin(); it != str.end(); it++) {
					// If we've hit the terminal character
					if (*it == '+' or *it == '-') {
						bool orientation = true;
						// If we have some characters accumulated
						// Add them to the result vector
						if (!next.empty()) {
							if (*it == '-') {
								orientation = false;
							}
							result.emplace_back(next, orientation);
							next.clear();
						}
					}
					else if (*it != ch) {
						// Accumulate the next character into the sequence
						next += *it;
					}
				}
				if (!next.empty())
					result.emplace_back(next, true); //this case shouldn't even happen
				return result;
			}

		bool is_number(const std::string& s)
		{
				    return !s.empty() && std::find_if(s.begin(), 
									        s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
		}

		// Avoiding un-necessary stream creation + replacing strings with string view
		// is a bit > than a 2x win!
		// implementation from : https://marcoarena.wordpress.com/tag/string_view/
		std::vector<stx::string_view> split(stx::string_view str, char delims)
		{
			std::vector<stx::string_view> ret;

			stx::string_view::size_type start = 0;
			auto pos = str.find_first_of(delims, start);
			while (pos != stx::string_view::npos) {
				if (pos != start) {
					ret.push_back(str.substr(start, pos - start));
				}
				start = pos + 1;
				pos = str.find_first_of(delims, start);
			}
			if (start < str.length()) {
				ret.push_back(str.substr(start, str.length() - start));
			}
			return ret;
		}

		std::vector<std::string> tokenize(const std::string &s, char delim) {
			std::stringstream ss(s);
			std::string item;
			std::vector<std::string> elems;
			while (std::getline(ss, item, delim)) {
				elems.push_back(item);
			}
			return elems;
		}
	public:
		std::map<std::string, std::vector<Position> > contig2pos;  
		PosFinder(char* gfaFileName, size_t input_k) {
				std::cerr << "Reading GFA file " << gfaFileName << "\n";
				file.open(gfaFileName);
				k = input_k;
		}

		void parseFile() {
			std::string ln;
			std::string tag, id, value;
			size_t contig_cnt{0};
			size_t ref_cnt{0};
			while(std::getline(file, ln)) {
					char firstC = ln[0];
					if (firstC != 'S' and firstC != 'P') continue;
					stx::string_view lnview(ln);
					std::vector<stx::string_view> splited = split(lnview, '\t');
					tag = splited[0].to_string();
					id = splited[1].to_string();
					value = splited[2].to_string();
					// A segment line
					if (tag == "S") {						
						if (is_number(id)) {
							contigid2cnt[id] = value.length();
						}
						contig_cnt++;
					} 
					// A path line
					if (tag == "P") {
						auto pvalue = splited[2];
						std::vector<std::pair<std::string, bool> > contigVec = explode(pvalue, ',');
						// parse value and add all conitgs to contigVec
						path[id] = contigVec;
						refIDs[id] = ref_cnt++;
					}
			}

			std::cerr << " Total # of Contigs : " << contig_cnt << " Total # of numerical Contigs : " << contigid2cnt.size() << "\n\n";
		}
		

		void mapContig2Pos() {
			uint64_t pos = 0;
			uint64_t accumPos;
			uint64_t currContigLength = 0;
			uint64_t total_output_lines = 0;
			for (auto const &ent : path) {
				const std::string& tr = ent.first;
				const std::vector<std::pair<std::string, bool> >& contigs = ent.second;
				accumPos = 0;
				for (long i = 0; i < contigs.size(); i++) {
					if (contig2pos.find(contigs[i].first) == contig2pos.end()) {
						contig2pos[contigs[i].first] = {};
						total_output_lines += 1;
					}
					if (contigid2cnt.find(contigs[i].first) == contigid2cnt.end()) {
						std::cerr<< contigs[i].first << "\n";
					}
					pos = accumPos;				
					currContigLength = contigid2cnt[contigs[i].first];
					accumPos += currContigLength - k;
					(contig2pos[contigs[i].first]).push_back(Position(refIDs[tr], pos, contigs[i].second));
				}
			}
			std::cerr << "\nTotal # of segments we have position for : " << total_output_lines << "\n";
		}

};

int main(int argc, char* argv[]) {
	//std::vector<std::string> gfa_file = {argv[1]};
	//std::vector<std::string> read_file = {argv[2]};
	std::cerr<<" Starting ... \n";
	PosFinder pf(argv[1], 30);
	pf.parseFile();
	
	pf.mapContig2Pos();
	
	std::ofstream ct("contig_table.tsv");
	for (auto const & ent : pf.contig2pos) {
		ct << ent.first;
		for (auto const & ent2 : ent.second) {
			ct << '\t' << ent2.transcript_id << ',' << ent2.pos << ',' << ent2.orien;
		}
		ct << '\n';
	}
	ct.close();
	
	
}

