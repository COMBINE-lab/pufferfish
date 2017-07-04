#include <fstream>
#include <iostream>
#include <bitset>

#include "cereal/archives/json.hpp"
#include "cereal/archives/binary.hpp"
#include "CLI/Timer.hpp"

#include "PufferfishSparseIndex.hpp"
#include "PufferFS.hpp"

PufferfishSparseIndex::PufferfishSparseIndex() {}

PufferfishSparseIndex::PufferfishSparseIndex(const std::string& indexDir) {
  if (!puffer::fs::DirExists(indexDir.c_str())) {
    std::cerr << "The index directory " << indexDir << " does not exist!\n";
    std::exit(1);
  }

  {
    std::ifstream infoStream(indexDir + "/info.json");
    cereal::JSONInputArchive infoArchive(infoStream);
    infoArchive(cereal::make_nvp("k", k_));
    infoArchive(cereal::make_nvp("num_kmers", numKmers_));
    infoArchive(cereal::make_nvp("num_sampled_kmers", numSampledKmers_));
    infoArchive(cereal::make_nvp("extension_size", extensionSize_));
    std::cerr << "k = " << k_ << '\n';
    std::cerr << "num kmers = " << numKmers_ << '\n';
    std::cerr << "num sampled kmers = " << numSampledKmers_ << '\n';
    std::cerr << "extension size = " << extensionSize_ << '\n';
    twok_ = 2*k_;
    infoStream.close();
  }

  //std::cerr << "loading contig table ... ";
  {
    CLI::AutoTimer timer {"Loading contig table", CLI::Timer::Big};
    std::ifstream contigTableStream(indexDir + "/ctable.bin");
    cereal::BinaryInputArchive contigTableArchive(contigTableStream);
    contigTableArchive(refNames_);
    contigTableArchive(contigTable_);
    contigTableStream.close();
  }
  numContigs_ = contigTable_.size();

  {
    CLI::AutoTimer timer {"Loading eq table", CLI::Timer::Big};
    std::ifstream eqTableStream(indexDir + "/eqtable.bin");
    cereal::BinaryInputArchive eqTableArchive(eqTableStream);
    eqTableArchive(eqClassIDs_);
    eqTableArchive(eqLabels_);
    eqTableStream.close();
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
    lastSeqPos_ = seq_.size() - k_;
  }
  /*
  {
    CLI::AutoTimer timer {"Loading positions", CLI::Timer::Big};
    std::string pfile = indexDir + "/pos.bin";
    sdsl::load_from_file(pos_, pfile);
  }*/

  {
    CLI::AutoTimer timer {"Loading presence vector", CLI::Timer::Big};
    std::string bfile = indexDir + "/presence.bin";
    sdsl::load_from_file(presenceVec_, bfile);
    presenceRank_ = decltype(presenceVec_)::rank_1_type(&presenceVec_);
    presenceSelect_ = decltype(presenceVec_)::select_1_type(&presenceVec_);
  }
  {
    CLI::AutoTimer timer {"Loading canonical vector", CLI::Timer::Big};
    std::string pfile = indexDir + "/canonical.bin";
    sdsl::load_from_file(canonicalNess_, pfile);
  }
  {
    CLI::AutoTimer timer {"Loading sampled positions", CLI::Timer::Big};
    std::string pfile = indexDir + "/sample_pos.bin";
    sdsl::load_from_file(sampledPos_, pfile);
  }

  {
    CLI::AutoTimer timer {"Loading extension vector", CLI::Timer::Big};
    std::string pfile = indexDir + "/extension.bin";
    sdsl::load_from_file(auxInfo_, pfile);
    std::string pfileSize = indexDir + "/extensionSize.bin";
    sdsl::load_from_file(extSize_, pfileSize);
  }

  {
    CLI::AutoTimer timer {"Loading direction vector", CLI::Timer::Big};
    std::string pfile = indexDir + "/direction.bin";
    sdsl::load_from_file(directionVec_, pfile);
  }

}

PufferfishSparseIndex::EqClassID PufferfishSparseIndex::getEqClassID(uint32_t contigID) {
  return eqClassIDs_[contigID];
}

const PufferfishSparseIndex::EqClassLabel& PufferfishSparseIndex::getEqClassLabel(uint32_t contigID) {
  return eqLabels_[getEqClassID(contigID)];
}

//auto endContigMap() -> decltype(contigTable_.begin()) { return contigTable_.end(); }
uint64_t PufferfishSparseIndex::getRawPos(CanonicalKmer& mer)  {
  auto km = mer.getCanonicalWord();
  size_t res = hash_->lookup(km);
  uint64_t pos =
    (res < numKmers_) ? pos_[res] : std::numeric_limits<uint64_t>::max();
  if (pos <= seq_.size() - k_) {
    uint64_t fk = seq_.get_int(2 * pos, 2 * k_);
    my_mer fkm;
    fkm.word__(0) = fk;
    auto keq = mer.isEquivalent(fkm);
    if (keq != KmerMatchType::NO_MATCH) {
      //}mer.fwWord() == fkm.word(0) or mer.rcWord() == fkm.word(0)) {
      return pos;
    }
  }
  pos = std::numeric_limits<uint64_t>::max();
  return pos;
}

bool PufferfishSparseIndex::contains(CanonicalKmer& mer) {
  return isValidPos(getRawPos(mer));
}

bool PufferfishSparseIndex::isValidPos(uint64_t pos) {
  return pos != std::numeric_limits<uint64_t>::max();
}

uint32_t PufferfishSparseIndex::contigID(CanonicalKmer& mer) {
    auto km = mer.getCanonicalWord();
    size_t res = hash_->lookup(km);
    uint64_t pos =
      (res < numKmers_) ? pos_[res] : std::numeric_limits<uint64_t>::max();
    if (pos <= seq_.size() - k_) {
      uint64_t fk = seq_.get_int(2 * pos, 2 * k_);
      my_mer fkm;
      fkm.word__(0) = fk;
      auto keq = mer.isEquivalent(fkm);
      if (keq != KmerMatchType::NO_MATCH) {
        auto rank = contigRank_(pos);
        return rank;
      }
    }
    return std::numeric_limits<uint32_t>::max();
  }

auto PufferfishSparseIndex::getRefPosHelper_(CanonicalKmer& mer, uint64_t pos, util::QueryCache& qc, bool didWalk) -> util::ProjectedHits {
  using IterT = std::vector<util::Position>::iterator;
  if (pos <= lastSeqPos_) {
    uint64_t twopos = pos << 1;
    uint64_t fk = seq_.get_int(twopos, twok_);
    // say how the kmer fk matches mer; either
    // identity, twin (i.e. rev-comp), or no match
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      // the index of this contig
      auto rank = contigRank_(pos);
      // make sure that the rank vector, from the 0th through k-1st position
      // of this k-mer is all 0s
      auto rankInterval = (didWalk) ? contigBoundary_.get_int(pos, (k_-1)) : 0;
      //auto rankEnd = contigRank_(pos + k_ - 1);
      if (rankInterval > 0) {
        return {std::numeric_limits<uint32_t>::max(), true, 0, k_, core::range<IterT>{}};
      }
      // the reference information in the contig table
      auto& pvec = contigTable_[rank];
      // start position of this contig
      uint64_t sp = 0;
      uint64_t contigEnd = 0;
      if ( rank == qc.prevRank ) {
        sp = qc.contigStart;
        contigEnd = qc.contigEnd;
      } else {
        sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
        contigEnd = contigSelect_(rank + 1);
        qc.prevRank = rank;
        qc.contigStart= sp;
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
      return {relPos, hitFW, static_cast<uint32_t>(clen), k_,
              core::range<IterT>{pvec.begin(), pvec.end()}};
    } else {
      return {std::numeric_limits<uint32_t>::max(), true, 0, k_,
              core::range<IterT>{}};
    }
  }

  return {std::numeric_limits<uint32_t>::max(), true, 0, k_,
          core::range<IterT>{}};
}



auto PufferfishSparseIndex::getRefPosHelper_(CanonicalKmer& mer, uint64_t pos, bool didWalk) -> util::ProjectedHits {

  using IterT = std::vector<util::Position>::iterator;
 if (pos <= lastSeqPos_) {
   uint64_t twopos = pos << 1;
   uint64_t fk = seq_.get_int(twopos, twok_);
    // say how the kmer fk matches mer; either
    // identity, twin (i.e. rev-comp), or no match
    auto keq = mer.isEquivalent(fk);
    if (keq != KmerMatchType::NO_MATCH) {
      // the index of this contig
      auto rank = contigRank_(pos);
      // make sure that the rank vector, from the 0th through k-1st position
      // of this k-mer is all 0s
      auto rankInterval = (didWalk) ? contigBoundary_.get_int(pos, (k_-1)) : 0;
      //auto rankEnd = contigRank_(pos + k_ - 1);
      if (rankInterval > 0) {
        return {std::numeric_limits<uint32_t>::max(), true, 0, k_, core::range<IterT>{}};
      }

      // the reference information in the contig table
      auto& pvec = contigTable_[rank];
      // start position of this contig
      uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(contigSelect_(rank)) + 1;
      uint64_t contigEnd = contigSelect_(rank + 1);

      // relative offset of this k-mer in the contig
      uint32_t relPos = static_cast<uint32_t>(pos - sp);

      // start position of the next contig - start position of this one
      auto clen = static_cast<uint64_t>(contigEnd + 1 - sp);

      bool hitFW = (keq == KmerMatchType::IDENTITY_MATCH);
      // how the k-mer hits the contig (true if k-mer in fwd orientation, false
      // otherwise)
      return {relPos, hitFW, static_cast<uint32_t>(clen), k_,
              core::range<IterT>{pvec.begin(), pvec.end()}};
    } else {
      return {std::numeric_limits<uint32_t>::max(), true, 0, k_,
              core::range<IterT>{}};
    }
  }
 return {std::numeric_limits<uint32_t>::max(), true, 0, k_,
     core::range<IterT>{}};
}


auto PufferfishSparseIndex::getRefPos(CanonicalKmer mern, util::QueryCache& qc) -> util::ProjectedHits {
  using IterT = std::vector<util::Position>::iterator;
  util::ProjectedHits emptyHit{std::numeric_limits<uint32_t>::max(), true, 0, k_, core::range<IterT>{}};

  bool didWalk{false};

  auto km = mern.getCanonicalWord() ;
  CanonicalKmer mer = mern;
  if (!mer.isFwCanonical()) {mer.swap();}

  // lookup this k-mer
  size_t idx = hash_->lookup(km);

  // if the index is invalid, it's clearly not present
  if(idx >= numKmers_) {
    return emptyHit;
  }

  uint64_t pos{0} ;
	auto currRank = (idx == 0) ? 0 : presenceRank_(idx) ;

	if(presenceVec_[idx] == 1){
		pos = sampledPos_[currRank];
	}else{
    didWalk = true;
    int signedShift{0};
		int inLoop = 0 ;

    /*
    do{
			if(inLoop >= 1){
        return emptyHit;
			}
    */

			auto extensionPos = idx - currRank ;
			uint64_t extensionWord = auxInfo_[extensionPos] ;

      if (!canonicalNess_[extensionPos] and mer.isFwCanonical()) {
        mer.swap();
      }

      bool shiftFw = (directionVec_[extensionPos] == 1);
      // + 1 because we encode 1 as 00, 2 as 01, etc.
      int32_t llimit = extensionSize_ - static_cast<int32_t>(extSize_[extensionPos] + 1);

      if(shiftFw){
        for( int32_t i = extensionSize_; i > llimit; --i) {
          uint32_t ssize = 2 * (i-1);
          int currCode = static_cast<int>((extensionWord & (0x3 << ssize)) >> ssize);
          mer.shiftFw(currCode) ;
          --signedShift;
        }
      } else {
        for( int32_t i = extensionSize_; i > llimit; --i) {
          uint32_t ssize = 2 * (i-1);
          int currCode = static_cast<int>((extensionWord & (0x3 << ssize)) >> ssize);
          mer.shiftBw(currCode) ;
          ++signedShift;
        }
      }


      km = mer.getCanonicalWord() ;
			idx = hash_->lookup(km) ;

      if(idx >= numKmers_) {
        return emptyHit;
      }

			currRank = (idx == 0) ? 0 : presenceRank_(idx) ;
			inLoop++ ;

      /*
		}while(presenceVec_[idx] != 1) ;
      */

      // if we didn't find a present kmer after extension, this is a no-go
      if (presenceVec_[idx] != 1) {
        return emptyHit;
      }
    auto sampledPos = sampledPos_[currRank];
    pos = sampledPos + signedShift;
	}
	//end of sampling based pos detection
  return getRefPosHelper_(mern, pos, qc, didWalk);
}

auto PufferfishSparseIndex::getRefPos(CanonicalKmer mern) -> util::ProjectedHits {
  using IterT = std::vector<util::Position>::iterator;
  util::ProjectedHits emptyHit{std::numeric_limits<uint32_t>::max(), true, 0, k_, core::range<IterT>{}};

  bool didWalk{false};

  auto km = mern.getCanonicalWord() ;
  CanonicalKmer mer = mern;
  if (!mer.isFwCanonical()) {mer.swap();}

  // lookup this k-mer
  size_t idx = hash_->lookup(km);

  // if the index is invalid, it's clearly not present
  if(idx >= numKmers_) {
    return emptyHit;
  }

  uint64_t pos{0} ;
	auto currRank = (idx == 0) ? 0 : presenceRank_(idx) ;

	if(presenceVec_[idx] == 1){
		pos = sampledPos_[currRank];
	}else{
    didWalk = true;
    int signedShift{0};
		int inLoop = 0 ;

    /*
    do{
			if(inLoop >= 1){
        return emptyHit;
			}
    */

			auto extensionPos = idx - currRank ;
			uint64_t extensionWord = auxInfo_[extensionPos] ;

      if (!canonicalNess_[extensionPos] and mer.isFwCanonical()) {
        mer.swap();
      }

      bool shiftFw = (directionVec_[extensionPos] == 1);
      // + 1 because we encode 1 as 00, 2 as 01, etc.
      int32_t llimit = extensionSize_ - static_cast<int32_t>(extSize_[extensionPos] + 1);

      if(shiftFw){
        for( int32_t i = extensionSize_; i > llimit; --i) {
          uint32_t ssize = 2 * (i-1);
          int currCode = static_cast<int>((extensionWord & (0x3 << ssize)) >> ssize);
          mer.shiftFw(currCode) ;
          --signedShift;
        }
      } else {
        for( int32_t i = extensionSize_; i > llimit; --i) {
          uint32_t ssize = 2 * (i-1);
          int currCode = static_cast<int>((extensionWord & (0x3 << ssize)) >> ssize);
          mer.shiftBw(currCode) ;
          ++signedShift;
        }
      }


      km = mer.getCanonicalWord() ;
			idx = hash_->lookup(km) ;

      if(idx >= numKmers_) {
        return emptyHit;
      }

			currRank = (idx == 0) ? 0 : presenceRank_(idx) ;
			inLoop++ ;

      /*
		}while(presenceVec_[idx] != 1) ;
      */

      // if we didn't find a present kmer after extension, this is a no-go
      if (presenceVec_[idx] != 1) {
        return emptyHit;
      }
    auto sampledPos = sampledPos_[currRank];
    pos = sampledPos + signedShift;
	}
	//end of sampling based pos detection
  return getRefPosHelper_(mern, pos, didWalk);
  
}

uint32_t PufferfishSparseIndex::k() { return k_; }

/**
 * Return the position list (ref_id, pos) corresponding to a contig.
 */
const std::vector<util::Position>& PufferfishSparseIndex::refList(uint64_t contigRank) {
  return contigTable_[contigRank];
}

const std::string& PufferfishSparseIndex::refName(uint64_t refRank) {
  return refNames_[refRank];
}
