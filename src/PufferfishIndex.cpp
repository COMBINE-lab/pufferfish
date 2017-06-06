#include <fstream>
#include <iostream>

#include "cereal/archives/json.hpp"
#include "cereal/archives/binary.hpp"
#include "CLI/Timer.hpp"

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

  //std::cerr << "loading contig table ... ";
  {
    CLI::AutoTimer timer {"Loading contig table", CLI::Timer::Big};
    std::ifstream contigTableStream(indexDir + "/ctable.bin");
    cereal::BinaryInputArchive contigTableArchive(contigTableStream);
    contigTableArchive(contigTable_);
    contigTableStream.close();
  }
  //std::cerr << "done\n";

  {
    CLI::AutoTimer timer {"Loading mphf table", CLI::Timer::Big};
    std::string hfile = indexDir + "/mphf.bin";
    std::ifstream hstream(hfile);
    hash_.reset(new boophf_t);
    hash_->load(hstream);
    hstream.close();
  }

  {
    CLI::AutoTimer timer {"Loading contig boundaries", CLI::Timer::Big};
    std::string bfile = indexDir + "/rank.bin";
    sdsl::load_from_file(contigBoundary_, bfile);
    contigRank_ = decltype(contigBoundary_)::rank_1_type(&contigBoundary_);
    contigSelect_ = decltype(contigBoundary_)::select_1_type(&contigBoundary_);
  }

  {
    CLI::AutoTimer timer {"Loading sequence", CLI::Timer::Big};
    std::string sfile = indexDir + "/seq.bin";
    sdsl::load_from_file(seq_, sfile);
  }

  {
    CLI::AutoTimer timer {"Loading positions", CLI::Timer::Big};
    std::string pfile = indexDir + "/pos.bin";
    sdsl::load_from_file(pos_, pfile);
  }
}
