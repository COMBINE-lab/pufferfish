
#include "ghc/filesystem.hpp"
#include "cereal/archives/json.hpp"
#include "parallel_hashmap/phmap.h"
//index header
#include "PuffAligner.hpp"
#include "ProgOpts.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "PufferfishLossyIndex.hpp"
#include "Kmer.hpp"
#include "Util.hpp"
#include "SpinLock.hpp"

template <typename Index>
bool dump_index_fasta(Index& pi, std::string& out) {
  const auto& refNames = pi.getFullRefNames();
  const auto& refLengths = pi.getFullRefLengths();
  if (!pi.hasReferenceSequence()) {
    std::cerr << "cannot dump fasta from an index that does not contain reference sequence\n";
    return false;
  }

  std::ofstream ofile(out);

  uint32_t k = pi.k();
  size_t nr = refNames.size();
  size_t curr{0};
  for (size_t i = 0; i < nr; ++i) {
    const auto& s = refNames[i];
    const auto& l = refLengths[i];
    ofile << ">" << s << "\n";
    // if this is a reference shorter than the k-mer length, then
    // we don't store it's sequence (since nothing can ever map to it).
    std::string seq = (l <= k) ? std::string(l, 'N') : pi.getRefSeqStr(curr, static_cast<int64_t>(l));
    ofile << seq << "\n";
    curr += (l < k) ? 0 : l;
  }
  ofile.close();

  return true;
}


template <typename Index>
bool dump_kmer_freq(Index& pi, std::string& out) {
  const auto numContigs = pi.numContigs();
  std::ofstream ofile(out);
  // maps from a frequency (key) to how many k-mers have that frequency
  phmap::flat_hash_map<uint64_t, uint64_t> freqMap;

  uint32_t k = pi.k();
  size_t curr{0};
  for (size_t i = 0; i < numContigs; ++i) {
    auto refRange = pi.refList(i);
    auto contigLen = pi.getContigLen(i);
    uint32_t numKmers = contigLen - k + 1;
    freqMap[refRange.size()] += numKmers;
  }
  std::vector<std::pair<uint64_t, uint64_t>> freqList;
  freqList.reserve(freqMap.size());
  for (auto& kv : freqMap) {
    freqList.push_back(std::make_pair(kv.first, kv.second));
  }
  std::sort(freqList.begin(), freqList.end());
  for (auto&& e : freqList) {
    ofile << e.first << "\t" << e.second << "\n";
  }
  ofile.close();
  return true;
}

int pufferfishExamine(ExamineOptions& opts) {
  bool dump_fasta = !opts.fasta_out.empty();
  bool kmer_freq = !opts.kmer_freq_out.empty();

  auto indexDir = opts.index_dir;
  std::string indexType;
  {
    std::ifstream infoStream(indexDir + "/info.json");
    cereal::JSONInputArchive infoArchive(infoStream);
    infoArchive(cereal::make_nvp("sampling_type", indexType));
    std::cerr << "Index type = " << indexType << '\n';
    infoStream.close();
  }
  bool s{false};
  if (indexType == "sparse") {
    PufferfishSparseIndex pi(opts.index_dir);
    if (dump_fasta) {
      s = dump_index_fasta(pi, opts.fasta_out);
    }
    if (kmer_freq) {
      s = dump_kmer_freq(pi, opts.kmer_freq_out);
    }
  } else if (indexType == "dense") {
    PufferfishIndex pi(opts.index_dir);
    if (dump_fasta) {
      s = dump_index_fasta(pi, opts.fasta_out);
    }
    if (kmer_freq) {
      s = dump_kmer_freq(pi, opts.kmer_freq_out);
    }
  } else if (indexType == "lossy") {
    PufferfishLossyIndex pi(opts.index_dir);
    if (dump_fasta) {
      s = dump_index_fasta(pi, opts.fasta_out);
    }
    if (kmer_freq) {
      s = dump_kmer_freq(pi, opts.kmer_freq_out);
    }
  }

  if (!s) {
    std::cerr << "Didn't dump reference fasta.";
  }

  return 0;
}
