#include <fstream>
#include <iostream>

#include "cereal/archives/json.hpp"

#include "PufferfishIndex.hpp"
#include "PufferFS.hpp"

PufferfishIndex::PufferfishIndex() {}

PufferfishIndex::PufferfishIndex(const std::string& indexDir) {
  if (!puffer::fs::DirExists(indexDir.c_str())) {
    std::cerr << "The index directory " << indexDir << " does not exist!\n";
    std::exit(1);
  }

  {
    std::ifstream infoStream(indexDir + "/info.json");
    cereal::JSONInputArchive infoArchive(infoStream);
    infoArchive(cereal::make_nvp("k", k_));
    std::cerr << "k = " << k_ << '\n';
    infoStream.close();
  }

  std::cerr << "loading mphf ... ";
  std::string hfile = indexDir + "/mphf.bin";
  std::ifstream hstream(hfile);
  hash_.reset(new boophf_t);
  hash_->load(hstream);
  hstream.close();
  std::cerr << "done\n";

  std::cerr << "loading contig boundary ... ";
  std::string bfile = indexDir + "/rank.bin";
  sdsl::load_from_file(contigBoundary_, bfile);
  contigRank_ = decltype(contigBoundary_)::rank_1_type(&contigBoundary_);
  contigSelect_ = decltype(contigBoundary_)::select_1_type(&contigBoundary_);
  std::cerr << "done\n";

  std::cerr << "loading seq ... ";
  std::string sfile = indexDir + "/seq.bin";
  sdsl::load_from_file(seq_, sfile);
  std::cerr << "done\n";

  std::cerr << "loading pos ... ";
  std::string pfile = indexDir + "/pos.bin";
  sdsl::load_from_file(pos_, pfile);
  std::cerr << "done\n";
}
