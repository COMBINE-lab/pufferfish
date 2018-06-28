#include <fstream>
#include <iostream>

#include "CLI/Timer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "PufferFS.hpp"
#include "PufferfishIndex.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/archives/json.hpp"

#include "jellyfish/mer_dna.hpp"

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
    infoArchive(cereal::make_nvp("num_kmers", numKmers_));
    std::cerr << "k = " << k_ << '\n';
    std::cerr << "num kmers = " << numKmers_ << '\n';
    infoStream.close();
    twok_ = 2 * k_;
  } 

  {
    CLI::AutoTimer timer{"Loading contig table", CLI::Timer::Big};
    std::ifstream contigTableStream(indexDir + "/ctable.bin");
    cereal::BinaryInputArchive contigTableArchive(contigTableStream);
    contigTableArchive(refNames_);
    contigTableArchive(contigTable_);
    contigTableArchive(contigOffsets_);
    contigTableStream.close();
  }
  numContigs_ = contigOffsets_.size()-1;

  {
    std::string rlPath = indexDir + "/reflengths.bin";
    if (puffer::fs::FileExists(rlPath.c_str())) {
      CLI::AutoTimer timer{"Loading reference lengths", CLI::Timer::Big};
      std::ifstream refLengthStream(rlPath);
      cereal::BinaryInputArchive refLengthArchive(refLengthStream);
      refLengthArchive(refLengths_);
    } else {
      refLengths_ = std::vector<uint32_t>(refNames_.size(), 1000);
    }
  }

  {
    CLI::AutoTimer timer{"Loading eq table", CLI::Timer::Big};
    std::ifstream eqTableStream(indexDir + "/eqtable.bin");
    cereal::BinaryInputArchive eqTableArchive(eqTableStream);
    eqTableArchive(eqClassIDs_);
    eqTableArchive(eqLabels_);
    eqTableStream.close();
  }

  {
    CLI::AutoTimer timer{"Loading mphf table", CLI::Timer::Big};
    std::string hfile = indexDir + "/mphf.bin";
    std::ifstream hstream(hfile);
    hash_.reset(new boophf_t);
    hash_->load(hstream);
    hstream.close();
    hash_raw_ = hash_.get();
  }

  {
    CLI::AutoTimer timer{"Loading contig boundaries", CLI::Timer::Big};
    std::string bfile = indexDir + "/rank.bin";
    sdsl::load_from_file(contigBoundary_, bfile);
    contigRank_ = decltype(contigBoundary_)::rank_1_type(&contigBoundary_);
    contigSelect_ = decltype(contigBoundary_)::select_1_type(&contigBoundary_);
  }
  /*
  selectPrecomp_.reserve(numContigs_+1);
  selectPrecomp_.push_back(0);
  for (size_t i = 1; i < numContigs_; ++i) {
    selectPrecomp_.push_back(contigSelect_(i));
  }
  selectPrecomp_.push_back(contigSelect_(numContigs_));
  */

  {
    CLI::AutoTimer timer{"Loading sequence", CLI::Timer::Big};
    std::string sfile = indexDir + "/seq.bin";
    sdsl::load_from_file(seq_, sfile);
    lastSeqPos_ = seq_.size() - k_;
  }

  {
    CLI::AutoTimer timer{"Loading positions", CLI::Timer::Big};
    std::string pfile = indexDir + "/pos.bin";
    sdsl::load_from_file(pos_, pfile);
  }

  {
    CLI::AutoTimer timer{"Loading edges", CLI::Timer::Big};
    std::string pfile = indexDir + "/edge.bin";
    sdsl::load_from_file(edge_, pfile);
  }
  /*
  {
    CLI::AutoTimer timer{"Loading edges", CLI::Timer::Big};
    std::string pfile = indexDir + "/revedge.bin";
    sdsl::load_from_file(revedge_, pfile);
  }
  */
}

/**
 * Returns a ProjectedHits object containing all of the reference loci matching this
 * provided Canonical kmer (including the oritentation of the match).  The provided
 * QueryCache argument will be used to avoid redundant rank / select operations if feasible.
 */
auto PufferfishIndex::getRefPos(CanonicalKmer& mer, util::QueryCache& qc)
    -> util::ProjectedHits {
  using IterT = std::vector<util::Position>::iterator;
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);
  if (res < numKmers_) {
    uint64_t pos = pos_[res];
    // if using quasi-dictionary idea (https://arxiv.org/pdf/1703.00667.pdf)
    /* 
    uint64_t hashbits = pos & 0xF;
    pos = pos >> 4;
    if ((km & 0xF) != hashbits) { 
        return {std::numeric_limits<uint32_t>::max(),
          std::numeric_limits<uint64_t>::max(),
          std::numeric_limits<uint32_t>::max(),
          true,
          0,
          k_,
          core::range<IterT>{}};
    }
    */
    uint64_t twopos = pos << 1;
    uint64_t fk = seq_.get_int(twopos, twok_);
    // say how the kmer fk matches mer; either
    // identity, twin (i.e. rev-comp), or no match
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      // the index of this contig
      auto rank = contigRank_(pos);
      // the reference information in the contig table
      auto contigIterRange = contigRange(rank);
      // start position of this contig
      uint64_t sp = 0;
      uint64_t contigEnd = 0;
      if (rank == qc.prevRank) {
        sp = qc.contigStart;
        contigEnd = qc.contigEnd;
      } else {
        sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
        contigEnd = contigSelect_(rank + 1);
        qc.prevRank = rank;
        qc.contigStart = sp;
        qc.contigEnd = contigEnd;
      }

      // relative offset of this k-mer in the contig
      uint32_t relPos = static_cast<uint32_t>(pos - sp);

      // start position of the next contig - start position of this one
      auto clen = static_cast<uint64_t>(contigEnd + 1 - sp);
      // auto clen =
      // cPosInfo_[rank].length();//static_cast<uint64_t>(contigSelect_(rank +
      // 1) + 1 - sp);

      // how the k-mer hits the contig (true if k-mer in fwd orientation, false
      // otherwise)
      bool hitFW = (keq == KmerMatchType::IDENTITY_MATCH);
      return {static_cast<uint32_t>(rank),
              pos,
              relPos,
              hitFW,
              static_cast<uint32_t>(clen),
              k_,
              contigIterRange};
              //core::range<IterT>{pvec.begin(), pvec.end()}};
    } else {
      return {std::numeric_limits<uint32_t>::max(),
              std::numeric_limits<uint64_t>::max(),
              std::numeric_limits<uint32_t>::max(),
              true,
              0,
              k_,
              core::range<IterT>{}};
    }
  }

  return {std::numeric_limits<uint32_t>::max(),
          std::numeric_limits<uint64_t>::max(),
          std::numeric_limits<uint32_t>::max(),
          true,
          0,
          k_,
          core::range<IterT>{}};
}

auto PufferfishIndex::getRefPos(CanonicalKmer& mer) -> util::ProjectedHits {
  using IterT = std::vector<util::Position>::iterator;
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);
  if (res < numKmers_) {
    uint64_t pos = pos_[res];
    uint64_t twopos = pos << 1;
    uint64_t fk = seq_.get_int(twopos, twok_);
    // say how the kmer fk matches mer; either
    // identity, twin (i.e. rev-comp), or no match
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      // the index of this contig
      auto rank = contigRank_(pos);
      // the reference information in the contig table
      auto contigIterRange = contigRange(rank);
      // start position of this contig
      uint64_t sp =
          (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
      uint64_t contigEnd = contigSelect_(rank + 1);

      // relative offset of this k-mer in the contig
      uint32_t relPos = static_cast<uint32_t>(pos - sp);

      // start position of the next contig - start position of this one
      auto clen = static_cast<uint64_t>(contigEnd + 1 - sp);
      // auto clen =
      // cPosInfo_[rank].length();//static_cast<uint64_t>(contigSelect_(rank +
      // 1) + 1 - sp);

      // how the k-mer hits the contig (true if k-mer in fwd orientation, false
      // otherwise)
      bool hitFW = (keq == KmerMatchType::IDENTITY_MATCH);
      return {static_cast<uint32_t>(rank),
              pos,
              relPos,
              hitFW,
              static_cast<uint32_t>(clen),
              k_,
             contigIterRange};
    } else {
      return {std::numeric_limits<uint32_t>::max(),
              std::numeric_limits<uint64_t>::max(),
              std::numeric_limits<uint32_t>::max(),
              true,
              0,
              k_,
              core::range<IterT>{}};
    }
  }

  return {std::numeric_limits<uint32_t>::max(),
          std::numeric_limits<uint64_t>::max(),
          std::numeric_limits<uint32_t>::max(),
          true,
          0,
          k_,
          core::range<IterT>{}};
}
