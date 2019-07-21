#include "FastxParser.hpp"
#include "jellyfish/mer_dna.hpp"
#include "clipp.h"
#include "sparsepp/spp.h"
#include "spdlog/spdlog.h"
#include "xxhash.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <unordered_set>

#include "cereal/cereal.hpp"
#include "cereal/archives/json.hpp"
#include "string_view.hpp"
#include "digestpp/digestpp.hpp"
#include "ghc/filesystem.hpp"

using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

void fixFasta(single_parser* parser,
              // std::string& outputDir,
              spp::sparse_hash_set<std::string>& decoyNames,
              bool keepDuplicates, uint32_t k,
              std::string& sepStr, std::mutex& iomutex,
              std::shared_ptr<spdlog::logger> log, std::string outFile) {
  (void)iomutex;

  ghc::filesystem::path outFilePath{outFile};
  ghc::filesystem::path outDir = outFilePath.parent_path();

  // std::shared_ptr<spdlog::logger> log) {
  // Create a random uniform distribution
  std::default_random_engine eng(271828);
  std::uniform_int_distribution<> dis(0, 3);

  // Hashers for getting txome signature
  digestpp::sha256 seqHasher256;
  digestpp::sha256 nameHasher256;
  digestpp::sha512 seqHasher512;
  digestpp::sha512 nameHasher512;
  digestpp::sha256 decoySeqHasher256;
  digestpp::sha256 decoyNameHasher256;

  // Keep track of if we've seen a decoy sequence yet.
  // The index enforces that all decoy sequences are consecutive, and that
  // they come after all valid (non-decoy) sequences.  If we see a non-decoy
  // sequence after having observed a decoy, then we complain and exit.
  bool sawDecoy{false};
  uint64_t numberOfDecoys{0};
  uint64_t firstDecoyIndex{std::numeric_limits<uint64_t>::max()};

  bool firstRecord{true};
  bool hasGencodeSep = (sepStr.find('|') != std::string::npos);
  uint32_t n{0};
  std::vector<std::string> transcriptNames;
  std::unordered_set<std::string> transcriptNameSet;
  std::map<std::string, bool> shortFlag ;

  std::vector<int64_t> transcriptStarts;
  constexpr char bases[] = {'A', 'C', 'G', 'T'};
  uint32_t polyAClipLength{10};
  uint32_t numPolyAsClipped{0};
  uint32_t numNucleotidesReplaced{0};
  std::string polyA(polyAClipLength, 'A');

  //using TranscriptList = std::vector<uint32_t>;
  //using KmerBinT = uint64_t;

  bool haveDecoys = !decoyNames.empty();
  bool clipPolyA = true;

  struct DupInfo {
    uint64_t txId;
    uint64_t txOffset;
    uint32_t txLen;
  };

  // http://biology.stackexchange.com/questions/21329/whats-the-longest-transcript-known
  // longest human transcript is Titin (108861), so this gives us a *lot* of
  // leeway before
  // we issue any warning.
  size_t tooLong = 200000;
  //size_t numDistinctKmers{0};
  //size_t numKmers{0};
  size_t currIndex{0};
  size_t numDups{0};
  std::map<XXH64_hash_t, std::vector<DupInfo>> potentialDuplicates;
  spp::sparse_hash_map<uint64_t, std::vector<std::string>> duplicateNames;
  std::cerr << "\n[Step 1 of 4] : counting k-mers\n";

  // rsdic::RSDicBuilder rsdb;
  std::vector<uint64_t>
      onePos; // Positions in the bit array where we should write a '1'
  // remember the initial lengths (e.g., before clipping etc., of all
  // transcripts)
  std::vector<uint32_t> completeLengths;
  // the stream of transcript sequence
  fmt::MemoryWriter txpSeqStream;
  {
    // ScopedTimer timer;
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser->getReadGroup();
    bool tooShort{false};
    while (parser->refill(rg)) {
      for (auto& read : rg) { // for each sequence
        tooShort = false;
        std::string& readStr = read.seq;
        readStr.erase(
            std::remove_if(readStr.begin(), readStr.end(),
                           [](const char a) -> bool { return !(isprint(a)); }),
            readStr.end());

        uint32_t readLen = readStr.size();
        uint32_t completeLen = readLen;

        // get the hash to check for collisions before we change anything.
        auto txStringHash =
            XXH64(reinterpret_cast<void*>(const_cast<char*>(readStr.data())),
                  readLen, 0);
        auto& readName = read.name;

        // check if we think this is a gencode transcriptome, and the user has not passed the gencode flag
        if (firstRecord and !hasGencodeSep) {
          constexpr const size_t numGencodeSep{8};
          if ( std::count(readName.begin(), readName.end(), '|') == numGencodeSep ) {
            log->warn("It appears that this may be a GENCODE transcriptome (from analyzing the separators in the FASTA header).  However, "
                      "you have not set \'|\' as a header separator.  If this is a GENCODE transcriptome, consider passing --gencode to the "
                      "pufferfish index command.\n\n");
          }
          firstRecord = false;
        }

        bool isDecoy = (haveDecoys) ? decoyNames.contains(readName) : false;
        // If this is *not* a decoy sequence, make sure that
        // we haven't seen any decoys yet.  Otherwise we are violating
        // the condition that decoys must come last.
        if (!isDecoy and sawDecoy) {
          log->critical("Observed a non-decoy sequence [{}] after having already observed a decoy. "
                        "However, it is required that any decoy target records appear, consecutively, "
                        "at the end of the input fasta file.  Please re-format your input file so that "
                        "all decoy records appear contiguously at the end of the file, after all valid "
                        "(non-decoy) records", readName);
          log->flush();
          spdlog::drop_all();
          std::exit(1);
        }

        // If this was a decoy, add it to the decoy hash
        if (isDecoy) {
          decoySeqHasher256.absorb(readStr.begin(), readStr.end());
        } else { // otherwise the ref hash
          seqHasher256.absorb(readStr.begin(), readStr.end());
          seqHasher512.absorb(readStr.begin(), readStr.end());
        }

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
          if (readStr.size() >= tooLong and !isDecoy) {
            log->warn("Entry with header [{}] was longer than {} nucleotides.  "
                      "Are you certain that "
                      "we are indexing a transcriptome and not a genome?",
                      read.name, tooLong);
          } else if (readStr.size() < k) {
            log->warn("Entry with header [{}], had length less than "
                      "the k-mer length of {} (perhaps after poly-A clipping)",
                      read.name, k);
            tooShort = true;
          }

          uint32_t txpIndex = n++;

          // The name of the current transcript
          auto& recHeader = read.name;
          auto processedName =
              recHeader.substr(0, recHeader.find_first_of(sepStr));

          // Add this transcript, indexed by it's sequence's hash value
          // to the potential duplicate list.
          bool didCollide{false};
          auto dupIt = potentialDuplicates.find(txStringHash);
          if (dupIt != potentialDuplicates.end()) {
            auto& dupList = dupIt->second;
            for (auto& dupInfo : dupList) {
              // they must be of the same length
              if (readLen == dupInfo.txLen) {
                bool collision =
                    (readStr.compare(0, readLen,
                                     txpSeqStream.data() + dupInfo.txOffset,
                                     readLen) == 0);
                if (collision) {
                  ++numDups;
                  didCollide = true;
                  duplicateNames[dupInfo.txId].push_back(processedName);
                  continue;
                } // if collision
              }   // if readLen == dupInfo.txLen
            }     // for dupInfo : dupList
          }       // if we had a potential duplicate

          if (!keepDuplicates and didCollide) {
            // roll back the txp index & skip the rest of this loop
            n--;
            continue;
          }

          // Check for duplicate name
          if (transcriptNameSet.find(processedName) != transcriptNameSet.end()) {
            log->error("In FixFasta, two references with the same name but different sequences: {}. "
                       "We require that all input records have a unique name "
                       "up to the first whitespace character.", processedName);
            std::exit(1);
          }
          // If there was no collision, then add the transcript
          transcriptNameSet.insert(processedName);
          transcriptNames.emplace_back(processedName);
          if (isDecoy) {
            decoyNameHasher256.absorb(processedName.begin(), processedName.end());
          } else {
            nameHasher256.absorb(processedName.begin(), processedName.end());
            nameHasher512.absorb(processedName.begin(), processedName.end());
          }

          if(!tooShort) {
              shortFlag[processedName] = false ;
          } else {
              shortFlag[processedName] = true ;
          }
          // nameHasher.process(processedName.begin(), processedName.end());

          // The position at which this transcript starts
          transcriptStarts.push_back(currIndex);
          // The un-molested length of this transcript
          completeLengths.push_back(completeLen);

          if (isDecoy) {
            // if we haven't seen another decoy yet, this is the first decoy
            // index
            if (!sawDecoy) {
              firstDecoyIndex = txpIndex;
            }
            // once we see the first decoy, saw decoy is set to true
            // for the rest of the processing.
            sawDecoy = true;
            ++numberOfDecoys;
            //decoyIndices.push_back(txpIndex);
          }

          // If we made it here, we were not an actual duplicate, so add this
          // transcript
          // for future duplicate checking.
          if (!keepDuplicates or (keepDuplicates and !didCollide)) {
            potentialDuplicates[txStringHash].push_back(
                {txpIndex, currIndex, readLen});
          }

          txpSeqStream << readStr;
          currIndex += readLen;
          onePos.push_back(currIndex);
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
      log->warn("Removed {} transcripts that were sequence duplicates of "
                "indexed transcripts.",
                numDups);
      log->warn("If you wish to retain duplicate transcripts, please use the "
                "`--keepDuplicates` flag");
    } else {
      log->warn("There were {} transcripts that would need to be removed to "
                "avoid duplicates.",
                numDups);
    }
  }

  /*
  std::ofstream dupClusterStream(outputDir + "duplicate_clusters.tsv");
  {
    dupClusterStream << "RetainedTxp" << '\t' << "DuplicateTxp" << '\n';
    for (auto kvIt = duplicateNames.begin(); kvIt != duplicateNames.end();
  ++kvIt) {
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
  stx::string_view concatTextView(concatText);
  // And clear the stream
  txpSeqStream.clear();

  std::ofstream ffa(outFile);
  size_t prev1{0};
  size_t numWritten{0};
  for (size_t i = 0; i < transcriptNames.size(); ++i) {
    size_t next1 = onePos[i];
    size_t len = next1 - prev1;
    if(!shortFlag[transcriptNames[i]]){
        ffa << ">" << transcriptNames[i] << "\n";
        ffa << concatTextView.substr(prev1, len) << "\n";
        ++numWritten;
    }
    prev1 = next1;
  }
  ffa.close();
  std::cerr << "wrote " << numWritten << " contigs\n";

  // Set the hash info
  std::string seqHash256 = seqHasher256.hexdigest();
  std::string nameHash256 = nameHasher256.hexdigest();
  std::string seqHash512 = seqHasher512.hexdigest();
  std::string nameHash512 = nameHasher512.hexdigest();
  std::string decoySeqHash256 = decoySeqHasher256.hexdigest();
  std::string decoyNameHash256 = decoyNameHasher256.hexdigest();

  {
    ghc::filesystem::path sigPath = outDir / ghc::filesystem::path{"ref_sigs.json"};
    std::ofstream os(sigPath.string());
    cereal::JSONOutputArchive ar(os);
    ar( cereal::make_nvp("num_decoys", numberOfDecoys));
    ar( cereal::make_nvp("first_decoy_index", firstDecoyIndex));
    ar( cereal::make_nvp("SeqHash", seqHash256) );
    ar( cereal::make_nvp("NameHash", nameHash256) );
    ar( cereal::make_nvp("SeqHash512", seqHash512) );
    ar( cereal::make_nvp("NameHash512", nameHash512) );
    ar( cereal::make_nvp("DecoySeqHash", decoySeqHash256) );
    ar( cereal::make_nvp("DecoyNameHash", decoyNameHash256) );
  }

  std::cerr << "seqHash 256 : " << seqHash256 << "\n";
  std::cerr << "seqHash 512 : " << seqHash512 << "\n";
  std::cerr << "nameHash 256 : " << nameHash256 << "\n";
  std::cerr << "nameHash 512 : " << nameHash512 << "\n";

  /*  header.setSeqHash256(seqHash256);
  header.setNameHash256(nameHash256);
  header.setSeqHash512(seqHash512);
  header.setNameHash512(nameHash512);
  header.setDecoySeqHash256(decoySeqHash256);
  header.setDecoyNameHash256(decoyNameHash256);
  */
}

spp::sparse_hash_set<std::string> populateDecoyHash(const std::string& fname, std::shared_ptr<spdlog::logger> log) {
  spp::sparse_hash_set<std::string> dset;
  std::ifstream dfile(fname);

  std::string dname;
  while (dfile >> dname) {
    auto it = dset.insert(dname);
    if (!it.second) {
      log->warn("The decoy name {} was encountered more than once --- please be sure all decoy names and sequences are unique.", dname);
    }
  }

  dfile.close();
  return dset;
}



int fixFastaMain(std::vector<std::string>& args) {
  using namespace clipp;

  uint32_t k{31};
  std::vector<std::string> refFiles;
  std::string outFile;
  std::string decoyFile;
  bool printHelp{false};
  std::string sepStr{" \t"};

  auto cli = (
              option("--help", "-h").set(printHelp, true) % "show usage",
              required("--input", "-i") & values("input", refFiles) % "input FASTA file",
              required("--output", "-o") & value("output", outFile) % "output FASTA file",
              option("--headerSep", "-s") & value("sep_strs", sepStr) %
              "Instead of a space or tab, break the header at the first "
              "occurrence of this string, and name the transcript as the token before "
              "the first separator (default = space & tab)",
              option("--decoys", "-d") & value("decoys", decoyFile) %
              "Treat these sequences as decoys that may be sequence-similar to some known indexed reference",
              option("--klen", "-k") & value("k-mer length", k) % "length of the k-mer used to build the cDBG (default = 31)"
              );

  //  if (parse(argc, argv, cli)) {
  if (parse(args, cli)) {
    if (printHelp) {
      std::cout << make_man_page(cli, "fixFasta");
      return 0;
    }

    auto console = spdlog::stderr_color_mt("ff::console");

    spp::sparse_hash_set<std::string> decoyNames;
    if (!decoyFile.empty()) {
      bool decoyFileExists = ghc::filesystem::exists(decoyFile);
      if (!decoyFileExists) {
        console->error("The decoy file {} does not exist.", decoyFile);
        std::exit(1);
      }
      decoyNames = populateDecoyHash(decoyFile, console);
    }

    size_t numThreads{1};
    std::unique_ptr<single_parser> transcriptParserPtr{nullptr};
    size_t numProd = 1;

    transcriptParserPtr.reset(new single_parser(refFiles, numThreads, numProd));
    transcriptParserPtr->start();
    std::mutex iomutex;
    bool keepDuplicates{true};
    fixFasta(transcriptParserPtr.get(), decoyNames, keepDuplicates, k, sepStr, iomutex, console,
             outFile);
    transcriptParserPtr->stop();
    return 0;
  } else {
    std::cout << usage_lines(cli, "fixFasta") << '\n';
    return 1;
  }
}
