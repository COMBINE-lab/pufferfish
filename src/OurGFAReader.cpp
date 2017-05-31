#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>

class PosFinder {
	private:
		std::ifstream file;
		size_t k;
		struct Contig {
			std::string seq;
			std::string id;
		};
		std::map<std::string, long> contigid2cnt;
		//orientation : + main, - reverse
		//std::map<std::pair<long, long>, short> link;
		std::map<std::string, std::vector< std::pair<std::string, bool> > > path;

		struct Position {
			std::string transcript_id;
			uint64_t pos;
			bool orien;

			Position(std::string tid, uint64_t tpos, bool torien) {
				transcript_id = tid;
				pos = tpos;
				orien = torien;
			}
		};
		// maps each contig to a list of positions in different transcripts

		std::vector<std::pair<std::string, bool> > explode(const std::string& str, const char& ch) {
				std::string next;
				std::vector< std::pair<std::string, bool> > result;
				// For each character in the string
				for (std::string::const_iterator it = str.begin(); it != str.end(); it++) {
					// If we've hit the terminal character
					if (*it == '+' or *it == '-') {
						bool orientation = true;
						// If we have some characters accumulated
						// Add them to the result vector
						if (!next.empty()) {
							if (*it == '-')
								orientation = false;
							result.push_back(std::make_pair(next, orientation));
							next.clear();
						}
					}
					else if (*it != ch) {
						// Accumulate the next character into the sequence
						next += *it;
					}
				}
				if (!next.empty())
					result.push_back(std::make_pair(next, true)); //this case shouldn't even happen
				return result;
			}

		bool is_number(const std::string& s)
		{
				    return !s.empty() && std::find_if(s.begin(), 
									        s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
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
			size_t contig_cnt = 0;
			while(std::getline(file, ln)) {
					std::vector<std::string> splited = tokenize(ln, '\t');
					tag = splited[0];
					//id = splited[1];
					//value = splited[2];
					//std::cout<<tag<<":"<<id<<","<<value<<std::endl;
					if (tag == "S") {						
						id = splited[1];
						value = splited[2];
						contig_cnt++;
						if (is_number(id)) {
							contigid2cnt[id] = value.length();
						}
						//std::cout<<"S "<<stol(id)<<","<<value.length()<<"\n";
					} 
					if (tag == "P") {
						id = splited[1];
						value = splited[2];
						std::vector<std::pair<std::string, bool> > contigVec = explode(value, ',');
						// parse value and add all conitgs to contigVec
						path[id] = contigVec;
						//std::cout<<"P "<<id<<"\n";
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
				std::string tr = ent.first;
				std::vector<std::pair<std::string, bool> > contigs = ent.second;
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
					(contig2pos[contigs[i].first]).push_back(Position(tr, pos, contigs[i].second));
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
	for (auto const & ent : pf.contig2pos) {
		std::cout<<ent.first<< ":";
		for (auto const & ent2 : ent.second) {
			std::cout<<ent2.transcript_id<<","<<ent2.pos<< ","<<ent2.orien<<";";
		}
		std::cout<<std::endl;
	}
}

