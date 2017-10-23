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
    // contigTableArchive(cPosInfo_);
    contigTableArchive(contigTable_);
    contigTableStream.close();
  }
  numContigs_ = contigTable_.size();

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
}

PufferfishIndex::EqClassID PufferfishIndex::getEqClassID(uint32_t contigID) {
  return eqClassIDs_[contigID];
}

const PufferfishIndex::EqClassLabel&
PufferfishIndex::getEqClassLabel(uint32_t contigID) {
  return eqLabels_[getEqClassID(contigID)];
}

uint64_t PufferfishIndex::getRawPos(CanonicalKmer& mer) {
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);
  if (res < numKmers_) {
    uint64_t pos = pos_[res];
    uint64_t fk = seq_.get_int(2 * pos, twok_);
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      return pos;
    }
  }
  return std::numeric_limits<uint64_t>::max();
}

bool PufferfishIndex::contains(CanonicalKmer& mer) {
  return isValidPos(getRawPos(mer));
}

bool PufferfishIndex::isValidPos(uint64_t pos) {
  return pos != std::numeric_limits<uint64_t>::max();
}

uint32_t PufferfishIndex::contigID(CanonicalKmer& mer) {
  auto km = mer.getCanonicalWord();
  size_t res = hash_raw_->lookup(km);
  if (res < numKmers_) {
    uint64_t pos = pos_[res];
    uint64_t fk = seq_.get_int(2 * pos, twok_);
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      auto rank = contigRank_(pos);
      return rank;
    }
  }
  return std::numeric_limits<uint32_t>::max();
}

// this function would store the part of the
// contig that is required
// if hitFW is true then we will clip
// the part towards the end of the contig
/*
void PufferfishIndex::getContigSeq(bool hitFW, uint64_t globalPos, int clipLen
,std::string& contigStr){ auto queryPos = kit->second ;

  uint64_t globalPos = phits.globalPos_ ;
  uint64_t twoGlobalPos = globalPos << 1 ;


}*/

/*
//TODO does the orientation of the contig matter while chopping sequence ?
//not needed
void PufferfishIndex::getRawSeq(util::ProjectedHits& phits,
CanonicalKmerIterator& kit,  std::string& contigStr, int readLen){
  //auto queryPos = kit->second ;

  //auto remainingContigLen = phits.contigLen_ - (phits.contigPos_ + k_) ;
  //auto remainingReadLen = readLen - (kit->second + k_ ) ;

    uint64_t globalPos = phits.globalPos_ ;
    //uint64_t twoGlobalPos = globalPos << 1 ;

    //index of the contig
    auto rank = contigRank_(globalPos) ;
    //start position of this contig
    uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) +
1;

    //clip from both side for now
    //go to right by read length
    //if contigLen < 2*readLen then
    //get back the entire contig
    uint64_t clipLen = 2*readLen ;
    uint64_t clipStart = sp ;
    if(phits.contigLen_ < clipLen){
      clipLen = phits.contigLen_ ;
    }else{
      clipStart = std::max(sp, globalPos - readLen) ;
      clipLen = std::min(readLen,static_cast<int>(globalPos - clipStart)) +
std::min(phits.contigPos_ + 101, phits.contigLen_ - phits.contigPos_) ;
    }
    uint64_t twoSp = clipStart << 1 ;
    auto numOfKmers = clipLen / k_ ;
    auto resLen  = clipLen % k_ ;
    uint64_t i{0} ;
    while(i < numOfKmers){
      uint64_t fk = seq_.get_int(twoSp + i*twok_ , twok_) ;
      CanonicalKmer mer ;
      mer.fromNum(fk) ;
      contigStr += mer.to_str() ;
      i++ ;
    }

    //TODO check boundary
    //residual kmer we might clip the wrong end
    if(resLen > 0){
      uint64_t fk = seq_.get_int(twoSp + i*twok_ , twok_) ;
      CanonicalKmer mer ;
      mer.fromNum(fk) ;
      contigStr += mer.to_str().substr(0,resLen) ;
      contigStr += "\0" ;
    }


    //return contigStr ;

}
*/

auto PufferfishIndex::getRefPos(CanonicalKmer& mer, util::QueryCache& qc)
    -> util::ProjectedHits {
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
      auto& pvec = contigTable_[rank];
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
              core::range<IterT>{pvec.begin(), pvec.end()}};
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
      auto& pvec = contigTable_[rank];
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
              core::range<IterT>{pvec.begin(), pvec.end()}};
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

uint32_t PufferfishIndex::k() { return k_; }

/**
 * Return the position list (ref_id, pos) corresponding to a contig.
 */
const std::vector<util::Position>&
PufferfishIndex::refList(uint64_t contigRank) {
  return contigTable_[contigRank];
}

const std::string& PufferfishIndex::refName(uint64_t refRank) {
  return refNames_[refRank];
}

const std::vector<std::string>& PufferfishIndex::getRefNames() {
  return refNames_;
}
