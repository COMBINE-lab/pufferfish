
#include "ghc/filesystem.hpp"
#include "cereal/archives/json.hpp"

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


int pufferfishExamine(ExamineOptions& opts) {
  bool dump_fasta = !opts.fasta_out.empty();

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
  } else if (indexType == "dense") {
    PufferfishIndex pi(opts.index_dir);
    if (dump_fasta) {
      s = dump_index_fasta(pi, opts.fasta_out);
    }
  } else if (indexType == "lossy") {
    PufferfishLossyIndex pi(opts.index_dir);
    if (dump_fasta) {
      s = dump_index_fasta(pi, opts.fasta_out);
    }
  }

  if (!s) {
    std::cerr << "Didn't dump reference fasta.";
  }

  return 0;
}
