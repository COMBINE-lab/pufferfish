
#include <cstdint>
#include <vector>
#include <iostream>
#include <fstream>

#include "sparsepp/spp.h"
#include "FastxParser.hpp"
#include "Kmer.hpp"
#include "clipp.h"

struct pconfig {
    uint32_t k=31;
    std::string ref;
    std::string gfa;
};

using KmerT = combinelib::kmers::Kmer<32,1>;

struct KmerContigMap {
    spp::sparse_hash_set<uint64_t> skmer;
    spp::sparse_hash_set<uint64_t> ekmer;
};

struct UnitigMap {
    spp::sparse_hash_map<uint64_t, std::pair<uint32_t, uint32_t>> skmer;
    spp::sparse_hash_map<uint64_t, std::pair<uint32_t, uint32_t>> ekmer;
};


KmerContigMap
getTerminalKmers(fastx_parser::FastxParser<fastx_parser::ReadSeq>& parser, uint32_t k) {
    KmerContigMap km;
    KmerT::k(k);
    KmerT start;
    KmerT end;
  // Get the read group by which this thread will
  // communicate with the parser (*once per-thread*)
  auto rg = parser.getReadGroup();

  while (parser.refill(rg)) {
    // Here, rg will contain a chunk of read pairs
    // we can process.
    for (auto& rp : rg) {
      auto& r1 = rp.seq;
      start.fromChars(rp.seq.begin());
      end.fromChars(rp.seq.end() - k);
      km.skmer.insert(start.word(0));
      km.ekmer.insert(end.word(0));
    }
  }
  return km;
}

uint32_t saveUnitig(std::string& seq, uint32_t id, std::ofstream& ofile) {
    ofile << "S\t" << id << "\t" << seq << '\n'; 
    return id + 1;
}

uint32_t createUnitig(UnitigMap& umap, uint32_t k, const std::string& header, std::string&& seq, 
                      uint32_t id, std::ofstream& ofile) {
    if (seq.length() == k) {
        KmerT ut(seq.begin());
        seq = ut.getCanonical().toStr();
    }
    auto nid = saveUnitig(seq, id, ofile);

    KmerT canonFirst(seq.begin());
    canonFirst.canonicalize();
    if (umap.skmer.contains(canonFirst.word(0)) or umap.ekmer.contains(canonFirst.word(0))) {
        std::cerr << "Error: Initial kmer " << canonFirst.toStr() << " is repeated.\n";
        std::exit(1);
    }

    KmerT canonLast(seq.end() - k);
    canonLast.canonicalize();
    if (umap.skmer.contains(canonLast.word(0)) or umap.ekmer.contains(canonLast.word(0))) {
        std::cerr << "Error: Last kmer is repeated.\n";
        std::exit(1);
    }
    umap.ekmer[canonLast.word(0)] = std::make_pair(id, seq.length());
    umap.skmer[canonFirst.word(0)] = std::make_pair(id, seq.length());
    return nid;
}

UnitigMap splitUnitigs(fastx_parser::FastxParser<fastx_parser::ReadSeq>& parser, KmerContigMap& km, uint32_t k, std::ofstream& ufile) {
  // Get the read group by which this thread will
  // communicate with the parser (*once per-thread*)
  UnitigMap umap;
  auto rg = parser.getReadGroup();
  uint32_t unitigId{0};
  KmerT mer;
  while (parser.refill(rg)) {
    for (auto& rp : rg) {
      auto& h = rp.name;
      auto& seq = rp.seq;
      uint32_t prev{0};
      bool first{true};
      mer.fromChars(seq.begin());
      for (uint32_t i = 0; i < seq.length() - k + 1; ++i) {
          if (!first) {
              mer.append(seq[i+k-1]);
          }
          first = false;
          auto rcmer = mer.getRC();
          if (km.skmer.contains(mer.word(0)) or km.ekmer.contains(rcmer.word(0))) {
              if (i + k - 1 - prev >= k) {
              uint64_t len = i + k - 1  - prev;
              unitigId = createUnitig(umap, k, h, seq.substr(prev, len), unitigId, ufile);
              prev = i;
              }
          }
          if (km.ekmer.contains(mer.word(0)) or km.skmer.contains(rcmer.word(0))) {
              uint64_t len = i + k - prev;
              unitigId = createUnitig(umap, k, h, seq.substr(prev, len), unitigId, ufile);
              prev = i + 1;
          }
      } 
      if (seq.length() - prev >= k) {
          unitigId = createUnitig(umap, k, h, seq.substr(prev), unitigId, ufile);
      }
    }
  }
  return umap;
}

bool buildPaths(fastx_parser::FastxParser<fastx_parser::ReadSeq>& parser, UnitigMap& umap, uint32_t k, std::ofstream& ufile) {
    uint32_t pctr{0};
    KmerT mer;
    char ori = '*';
    auto rg = parser.getReadGroup();
    while (parser.refill(rg)) {
        for (auto& rp : rg) { 
            auto& ref = rp.seq;
            uint32_t i{0};
            ufile << "P\t" << rp.name << '\t';
            bool first{true};
            while (i < ref.length() - k + 1) {
                mer.fromChars(ref.begin() + i); 
                auto nkmer = mer.getCanonical();
                std::decay<decltype(umap.skmer[nkmer.word(0)])>::type unitigInfo;
                if (umap.skmer.contains(nkmer.word(0)) and umap.ekmer.contains(nkmer.word(0))) {
                    unitigInfo = umap.skmer[nkmer.word(0)];
                    ori = (mer == nkmer) ? '+' : '-';
                } else if (umap.skmer.contains(nkmer.word(0))) {
                    unitigInfo = umap.skmer[nkmer.word(0)];
                    ori = '+';
                } else if (umap.ekmer.contains(nkmer.word(0))) {
                    unitigInfo = umap.ekmer[nkmer.word(0)];
                    ori = '-';
                } else {
                    std::cerr << pctr << " paths reconstructed.\n";
                    std::cerr << "ERROR: kmer is not found in the start or end of a unitig\n";
                    std::cerr << mer.toStr() << "   ,   " << nkmer.toStr() << "\n";
                    std::exit(1);
                }
                i += unitigInfo.second - k + 1;
                char comma = (i < ref.length() - k + 1) ? ',' : '\n';
                ufile << unitigInfo.first << ori << comma;
            }
            ++pctr;
        }
    }
    ufile << '\n';
}

int main(int argc, char* argv[]) {
    using namespace clipp;
    pconfig config;

    auto cli = (
        required("-r", "--ref") & value("reference file", config.ref) % "reference fasta file",
        required("-u", "--uinitigs") & value("BCALM unitig file", config.gfa) % "unitig file",
        option("-k", "--klen") & value("length of k", config.k) % "length of kmer (default : 31)"
    );

    if(!parse(argc, argv, cli)) std::cout << make_man_page(cli, argv[0]);
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser({config.ref}, 1, 1);
    parser.start();
    auto kmap = getTerminalKmers(parser, config.k);
    parser.stop();
    std::cout << "have " << kmap.skmer.size() << " contigs\n";
    
    std::cout << "splitting unitigs\n";

    std::string ufname = config.gfa + ".pufferized.gfa";
    std::ofstream ufile(ufname);
    fastx_parser::FastxParser<fastx_parser::ReadSeq> uparser({config.gfa}, 1, 1);
    uparser.start();
    auto umap = splitUnitigs(uparser, kmap, config.k, ufile);
    uparser.stop();
    std::cerr << "done\n";
    std::cerr << "umap size = " << umap.skmer.size() << "\n";
    
    std::cerr << "start reconstructing the paths\n";
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser({config.ref}, 1, 1);
    rparser.start();
    buildPaths(rparser, umap, config.k, ufile);
    rparser.stop();
    std::cerr << "done!\n";
    return 0;
}