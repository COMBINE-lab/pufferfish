
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
  const auto& refNames = pi.getRefNames();
  const auto& refLengths = pi.getRefLengths();
  if (!pi.hasReferenceSequence()) {
    std::cerr << "cannot dump fasta from an index that does not contain reference sequence\n";
    return false;
  }

  std::ofstream ofile(out);

  size_t nr = refNames.size();
  size_t curr{0};
  for (size_t i = 0; i < nr; ++i) {
    const auto& s = refNames[i];
    const auto& l = refLengths[i];
    ofile << ">" << s << "\n";
    ofile << pi.getRefSeqStr(curr, static_cast<int64_t>(l)) << "\n";
    curr += l;
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

  if (indexType == "sparse") {
    PufferfishSparseIndex pi(opts.index_dir);
    if (dump_fasta) {
      bool s = dump_index_fasta(pi, opts.fasta_out);
    }
  } else if (indexType == "dense") {
    PufferfishIndex pi(opts.index_dir);
    if (dump_fasta) {
      bool s = dump_index_fasta(pi, opts.fasta_out);
    }
  }


  return 0;
}
