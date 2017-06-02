#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <random>
#include <type_traits>
#include <unordered_map>
#include <map>
#include <vector>
#include "FastxParser.hpp"
#include "xxhash.h"
#include "sparsepp/spp.h"
#include "jellyfish/mer_dna.hpp"
#include "spdlog/spdlog.h"
#include "popl.hpp"

using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

void fixFasta(single_parser* parser,
              //std::string& outputDir,
              bool keepDuplicates,
              uint32_t k,
              std::mutex& iomutex,
              std::shared_ptr<spdlog::logger> log,
              std::string outFile) {
              //std::shared_ptr<spdlog::logger> log) {
  // Create a random uniform distribution
  std::default_random_engine eng(271828);
  std::uniform_int_distribution<> dis(0, 3);
  // Hashers for getting txome signature
  //picosha2::hash256_one_by_one seqHasher; seqHasher.init();
  //picosha2::hash256_one_by_one nameHasher; nameHasher.init();

  uint32_t n{0};
  std::vector<std::string> transcriptNames;
  std::vector<int64_t> transcriptStarts;
  constexpr char bases[] = {'A', 'C', 'G', 'T'};
  uint32_t polyAClipLength{10};
  uint32_t numPolyAsClipped{0};
  uint32_t numNucleotidesReplaced{0};
  std::string polyA(polyAClipLength, 'A');

  using TranscriptList = std::vector<uint32_t>;
  using KmerBinT = uint64_t;

  bool clipPolyA = true; 

  struct DupInfo {
    uint64_t txId;
    uint64_t txOffset;
    uint32_t txLen;
  };

  std::string sepStr = " \t";

  // http://biology.stackexchange.com/questions/21329/whats-the-longest-transcript-known
  // longest human transcript is Titin (108861), so this gives us a *lot* of
  // leeway before
  // we issue any warning.
  size_t tooLong = 200000;
  size_t numDistinctKmers{0};
  size_t numKmers{0};
  size_t currIndex{0};
  size_t numDups{0};
  std::map<XXH64_hash_t, std::vector<DupInfo>> potentialDuplicates;
  spp::sparse_hash_map<uint64_t, std::vector<std::string>> duplicateNames;
  std::cerr << "\n[Step 1 of 4] : counting k-mers\n";

  // rsdic::RSDicBuilder rsdb;
  std::vector<uint64_t>
      onePos; // Positions in the bit array where we should write a '1'
  // remember the initial lengths (e.g., before clipping etc., of all transcripts)
  std::vector<uint32_t> completeLengths;
  // the stream of transcript sequence
  fmt::MemoryWriter txpSeqStream;
  {
    //ScopedTimer timer;
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser->getReadGroup();

    while (parser->refill(rg)) {
      for (auto& read : rg) { // for each sequence
        std::string& readStr = read.seq;
        readStr.erase(
            std::remove_if(readStr.begin(), readStr.end(),
                           [](const char a) -> bool { return !(isprint(a)); }),
            readStr.end());

        //seqHasher.process(readStr.begin(), readStr.end());

        uint32_t readLen = readStr.size();
        uint32_t completeLen = readLen;

        // get the hash to check for collisions before we change anything.
        auto txStringHash = XXH64(reinterpret_cast<void*>(const_cast<char*>(readStr.data())), readLen, 0);

        // First, replace non ATCG nucleotides
        for (size_t b = 0; b < readLen; ++b) {
          readStr[b] = ::toupper(readStr[b]);
          int c = jellyfish::mer_dna::code(readStr[b]);
          // Replace non-ACGT bases with pseudo-random bases
          if (jellyfish::mer_dna::not_dna(c)) {
            char rbase = bases[dis(eng)];
            c = jellyfish::mer_dna::code(rbase);
            readStr[b] = rbase;
            ++numNucleotidesReplaced;
          }
        }

        // Now, do Kallisto-esque clipping of polyA tails
        if (clipPolyA) {
          if (readStr.size() > polyAClipLength and
              readStr.substr(readStr.length() - polyAClipLength) == polyA) {

            auto newEndPos = readStr.find_last_not_of("Aa");
            // If it was all As
            if (newEndPos == std::string::npos) {
              log->warn("Entry with header [{}] appeared to be all A's; it "
                        "will be removed from the index!",
                        read.name);
              readStr.resize(0);
            } else {
              readStr.resize(newEndPos + 1);
            }
            ++numPolyAsClipped;
          }
        }

        readLen = readStr.size();
        // If the transcript was completely removed during clipping, don't
        // include it in the index.
        if (readStr.size() > 0) {
          // If we're suspicious the user has fed in a *genome* rather
          // than a transcriptome, say so here.
          if (readStr.size() >= tooLong) {
            log->warn("Entry with header [{}] was longer than {} nucleotides.  "
                      "Are you certain that "
                      "we are indexing a transcriptome and not a genome?",
                      read.name, tooLong);
          } else if (readStr.size() < k) {
            log->warn("Entry with header [{}], had length less than "
                      "the k-mer length of {} (perhaps after poly-A clipping)",
                      read.name, k);
          }

          uint32_t txpIndex = n++;

          // The name of the current transcript
          auto& recHeader = read.name;
          auto processedName = recHeader.substr(0, recHeader.find_first_of(sepStr));

          // Add this transcript, indexed by it's sequence's hash value
          // to the potential duplicate list.
          bool didCollide{false};
          auto dupIt = potentialDuplicates.find(txStringHash);
          if (dupIt != potentialDuplicates.end()) {
            auto& dupList = dupIt->second;
            for (auto& dupInfo : dupList) {
              // they must be of the same length
              if (readLen == dupInfo.txLen) {
                bool collision = (readStr.compare(0, readLen, txpSeqStream.data() + dupInfo.txOffset, readLen) == 0);
                if (collision) {
                  ++numDups;
                  didCollide = true;
                  duplicateNames[dupInfo.txId].push_back(processedName);
                  continue;
                } // if collision
              } // if readLen == dupInfo.txLen
            } // for dupInfo : dupList
          } // if we had a potential duplicate

          if (!keepDuplicates and didCollide) {
            // roll back the txp index & skip the rest of this loop
            n--;
            continue;
          }

          // If there was no collision, then add the transcript
          transcriptNames.emplace_back(processedName);
          //nameHasher.process(processedName.begin(), processedName.end());

          // The position at which this transcript starts
          transcriptStarts.push_back(currIndex);
          // The un-molested length of this transcript
          completeLengths.push_back(completeLen);

          // If we made it here, we were not an actual duplicate, so add this transcript
          // for future duplicate checking.
          if (!keepDuplicates or (keepDuplicates and !didCollide)) {
            potentialDuplicates[txStringHash].push_back({txpIndex, currIndex, readLen});
          }

          txpSeqStream << readStr;
          txpSeqStream << '$';
          currIndex += readLen + 1;
          onePos.push_back(currIndex - 1);
        } else {
            log->warn("Discarding entry with header [{}], since it had length 0 "
                      "(perhaps after poly-A clipping)",
                      read.name);
        }
      }
      if (n % 10000 == 0) {
        std::cerr << "\r\rcounted k-mers for " << n << " transcripts";
      }
    }
  }
  std::cerr << "\n";
  if (numDups > 0) {
    if (!keepDuplicates) {
      log->warn("Removed {} transcripts that were sequence duplicates of indexed transcripts.", numDups);
      log->warn("If you wish to retain duplicate transcripts, please use the `--keepDuplicates` flag");
    } else {
      log->warn("There were {} transcripts that would need to be removed to avoid duplicates.", numDups);
    }
  }

  /*
  std::ofstream dupClusterStream(outputDir + "duplicate_clusters.tsv");
  {
    dupClusterStream << "RetainedTxp" << '\t' << "DuplicateTxp" << '\n';
    for (auto kvIt = duplicateNames.begin(); kvIt != duplicateNames.end(); ++kvIt) {
      auto& retainedName = transcriptNames[kvIt->first];
      for (auto& droppedName : kvIt->second) {
        dupClusterStream << retainedName << '\t' << droppedName << '\n';
      }
    }
  }
  dupClusterStream.close();
  */

  std::cerr << "Replaced " << numNucleotidesReplaced
            << " non-ATCG nucleotides\n";
  std::cerr << "Clipped poly-A tails from " << numPolyAsClipped
            << " transcripts\n";

  // Put the concatenated text in a string
  std::string concatText = txpSeqStream.str();
  // And clear the stream
  txpSeqStream.clear();


  std::ofstream ffa(outFile);
  size_t prev1{0};
  for (size_t i = 0; i < transcriptNames.size(); ++i) {
    size_t next1 = onePos[i];
    size_t len = next1 - prev1;
    ffa << "> " << transcriptNames[i] << "\n";
    ffa << concatText.substr(prev1, len) << "\n";
    prev1 = next1 + 1;
  }
  ffa.close();
}


int main(int argc, char* argv[]) {
  uint32_t k;
  std::string cfile;
  std::string rfile;
  popl::Switch helpOption("h", "help", "produce help message");
  popl::Value<uint32_t> kOpt("k", "klen", "length of the k-mer with which the compacted dBG was built", 31, &k);
  popl::Value<std::string> inOpt("i", "input", "input fasta file");
  popl::Value<std::string> outOpt("o", "out", "output fasta file");

  popl::OptionParser op("Allowed options");
	op.add(helpOption)
    .add(kOpt)
    .add(inOpt)
    .add(outOpt);

  op.parse(argc, argv);
  if (helpOption.isSet()) {
    std::cout << op << '\n';
    std::exit(0);
  }

  std::string refFile = inOpt.getValue();
  std::string outFile = outOpt.getValue();

  size_t numThreads{1};
  std::unique_ptr<single_parser> transcriptParserPtr{nullptr};
  size_t numProd = 1;
  std::vector<std::string> refFiles{refFile};

  auto console = spdlog::stderr_color_mt("console");

  transcriptParserPtr.reset(new single_parser(refFiles, numThreads, numProd));
  transcriptParserPtr->start();
  std::mutex iomutex;
  bool keepDuplicates{true};
  fixFasta(transcriptParserPtr.get(), keepDuplicates, k, iomutex, console, outFile);
  return 0;
}
