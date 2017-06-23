//#include <iterator>
//#include <type_traits>
//#include "sdsl/int_vector.hpp"
//#include "jellyfish/mer_dna.hpp"
//#include "string_view.hpp"
//#include "BooPHF.h"
//#include "popl.hpp"
//#include "ScopedTimer.hpp"
#include "CLI/CLI.hpp"
#include "sparsepp/spp.h"

#include "Util.hpp"

int main(int argc, char* argv[]) {
  CLI::App app{"pufferizer : prepare a GFA file for pufferfish"};

  std::string gfa_file;
  app.add_option("filename", gfa_file, "the input file")->required();

  try {
    app.parse(argc, argv);
  } catch (const CLI::Error& e) {
    return app.exit(e);
  }

  std::ifstream file(gfa_file);
  short k = 31;
  std::string ln;
  std::string tag, id, value;
  size_t contig_cnt{0};
  size_t ref_cnt{0};
  spp::sparse_hash_map<std::string, bool> right_clipped;
  spp::sparse_hash_map<std::string, bool> left_clipped;
  spp::sparse_hash_map<std::string, std::string> contigid2seq;
  spp::sparse_hash_map<std::string, std::string> kmer2contigid;
  spp::sparse_hash_map<std::string, std::string> reconstructedTr;
  spp::sparse_hash_map<std::string, uint32_t> refIDs;
  spp::sparse_hash_map<uint32_t, std::string> refMap;
  spp::sparse_hash_map<std::string, std::vector<std::pair<std::string, bool>>>
      path;

  uint64_t newContigCntr = 0;
  std::string prefix = "00"; // TODO we just keep it number to pass is_number,
                             // but remember not to EVER convert contigID to
                             // number (e.g. uint64_t) b/c it'll have conflict
                             // of IDs
  while (std::getline(file, ln)) {
    char firstC = ln[0];
    if (firstC != 'S' and firstC != 'P')
      continue;
    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = util::split(lnview, '\t');
    tag = splited[0].to_string();
    id = splited[1].to_string();
    value = splited[2].to_string();
    // A segment line
    if (tag == "S") {
      if (util::is_number(id)) {
        contigid2seq[id] = value;
      }
      contig_cnt++;
    }
    // A path line
    if (tag == "P") {
      auto pvalue = splited[2];
      std::vector<std::pair<std::string, bool>> contigVec =
          util::explode(pvalue, ',');
      std::vector<std::pair<std::string, bool>> updatedContigVec = {};
      if (contigVec.size() > 1) {
        for (uint64_t i = 0; i < contigVec.size() - 1; i++) {
          updatedContigVec.push_back(contigVec[i]);
          if (!contigid2seq.contains(contigVec[i].first)) {
            std::cerr << "can't find " << contigVec[i].first << '\n';
          }
          std::string contigSeq = contigid2seq[contigVec[i].first];
          std::cerr << contigSeq.length() << '\n';
          // left + : get the right most k nucleotides
          std::string kmer =
              contigSeq.substr(contigSeq.size() - k, contigSeq.size());
          if (!contigVec[i].second) { // which means the orientation of the
                                      // contig in negative in the path
            // left - : get the left most k nucleotides
            kmer = contigSeq.substr(0, k);
          }

          const auto& curr = contigVec[i].first;
          const auto& next = contigVec[i + 1].first;
          if (contigid2seq[curr].length() <= k + 1 or
              contigid2seq[next].length() <= k + 1) {
            continue;
          }
          // First we add the left most segment id to the path with its original
          // orientation
          // We also update the contig sequence if it hasn't been updated before
          /*
           * clipping
           * left + : right
           * left - : left
           * right + : left
           * right - : right
           */
          if (contigVec[i].second) { // left +
            if (right_clipped.find(contigVec[i].first) == right_clipped.end()) {
              contigid2seq[contigVec[i].first].erase(
                  contigid2seq[contigVec[i].first].size() - 1,
                  1); // erase right most nucleotide
              right_clipped[contigVec[i].first] = true;
            }
          } else { // left -
            if (left_clipped.find(contigVec[i].first) == left_clipped.end()) {
              contigid2seq[contigVec[i].first].erase(
                  0, 1); // erase left most nucleotide
              left_clipped[contigVec[i].first] = true;
            }
          }
          if (contigVec[i + 1].second) { // right +
            if (left_clipped.find(contigVec[i + 1].first) ==
                left_clipped.end()) {
              contigid2seq[contigVec[i + 1].first].erase(
                  0, 1); // erase left most nucleotide
              left_clipped[contigVec[i + 1].first] = true;
            }
          } else { // right -
            if (right_clipped.find(contigVec[i + 1].first) ==
                right_clipped.end()) {
              contigid2seq[contigVec[i + 1].first].erase(
                  contigid2seq[contigVec[i + 1].first].size() - 1,
                  1); // erase right most nucleotide
              right_clipped[contigVec[i + 1].first] = true;
            }
          }

          // then we add the kmer segment
          std::string kmerId;
          if (kmer2contigid.find(kmer) != kmer2contigid.end()) {
            kmerId = kmer2contigid[kmer];
            updatedContigVec.emplace_back(
                kmerId, contigVec[i].second); // insert the kmer segment with
                                              // the same orientation as left
                                              // segment
          } else if (kmer2contigid.find(util::revcomp(kmer)) !=
                     kmer2contigid.end()) {
            kmerId = kmer2contigid[util::revcomp(kmer)];
            updatedContigVec.emplace_back(
                kmerId, !contigVec[i].second); // insert the kmer segment with
                                               // the REVERSE orientation as
                                               // left segment
          } else {
            kmer2contigid[kmer] =
                prefix +
                std::to_string(
                    newContigCntr++); // at the same time increase newContigCntr
            kmerId = kmer2contigid[kmer];
            updatedContigVec.emplace_back(
                kmerId, contigVec[i].second); // insert the kmer segment with
                                              // the same orientation as left
                                              // segment
          }
        }
        updatedContigVec.push_back(contigVec[contigVec.size() - 1]);
      } else
        updatedContigVec.push_back(contigVec[0]);
      // parse value and add all conitgs to contigVec
      path[id] = updatedContigVec;
      refMap[ref_cnt] = id;
      refIDs[id] = ref_cnt++;
    }
  }
  return 0;
}
