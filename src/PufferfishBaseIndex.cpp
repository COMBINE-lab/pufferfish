#include "PufferfishBaseIndex.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "PufferfishLossyIndex.hpp"

// from : https://www.fluentcpp.com/2017/05/19/crtp-helper/
template <typename T>
T& PufferfishBaseIndex<T>::underlying() { return static_cast<T&>(*this); }
template <typename T>
T const& PufferfishBaseIndex<T>::underlying() const { return static_cast<T const&>(*this); }

//template <typename T>
//PufferfishBaseIndex<T>::PufferfishBaseIndex() {}
//PufferfishBaseIndex::PufferfishBaseIndex(const std::string& indexDir) {}

template <typename T>
pufferfish::types::EqClassID PufferfishBaseIndex<T>::getEqClassID(uint32_t contigID) {
  return underlying().eqClassIDs_[contigID];
}

template <typename T>
const pufferfish::types::EqClassLabel&
PufferfishBaseIndex<T>::getEqClassLabel(uint32_t contigID) {
  return underlying().eqLabels_[getEqClassID(contigID)];
}

/** Seems not to be used, ignore for now **/
/*
uint64_t PufferfishBaseIndex::getRawPos(CanonicalKmer& mer) {
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
*/

/** Seems not to be used, ignore for now **/
/*
bool PufferfishBaseIndex::contains(CanonicalKmer& mer) {
  return isValidPos(getRawPos(mer));
}
*/

/** Seems not to be used, ignore for now **/
/*
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
*/

template <typename T>
bool PufferfishBaseIndex<T>::isValidPos(uint64_t pos) {
  return pos != std::numeric_limits<uint64_t>::max();
}

template <typename T>
std::string PufferfishBaseIndex<T>::getSeqStr(size_t globalPos, size_t length, bool isFw) {
    auto& seq_ = underlying().seq_;
	std::string outstr = "";
	uint64_t validLength = 0;
	uint64_t word = 0;
	uint8_t base = 0;
	while (length > 0) {
	validLength = std::min(length, (size_t)32);
	length -= validLength;
 	word = seq_.get_int(2*globalPos, 2*validLength);
	globalPos += validLength;
	if (isFw)
	  for(uint64_t i = 0; i < 2*validLength ;i+=2){
      base = (word >> i) & 0x03;
	    switch(base){
    	case 0:
        outstr += 'A';
        break ;
	    case 1:
        outstr += 'C';
	      break ;
    	case 2:
        outstr += 'G';
    	  break ;
	    case 3:
        outstr += 'T';
	      break ;
    	}
    }
	else
    for(uint64_t i = 0; i < 2*validLength ;i+=2){
      base = (word >> i) & 0x03;
      switch(base){
      case 0:
        outstr = 'T' + outstr;
        break ;
      case 1:
        outstr = 'G' + outstr;
        break ;
      case 2:
        outstr = 'C' + outstr;
        break ;
      case 3:
        outstr = 'A' + outstr;
        break ;
      }
    }
	}
  return outstr;
}


/** Seems unused, skip for now **/
/*
inline uint64_t our_read_int(uint64_t* word, uint8_t offset, const uint8_t len)
{
    uint64_t w1 = (*word)>>offset;
    if ((offset+len) > 64) { // if offset+len > 64
        return w1 |  // w1 or w2 adepted:
               ((*(word+1) & sdsl::bits::lo_set[(offset+len)&0x3F])   // set higher bits zero
                << (64-offset));  // move bits to the left
    } else {
        return w1 & sdsl::bits::lo_set[len];
    }
}
*/

/**
 * Returns a ProjectedHits object containing all of the reference loci matching this
 * provided Canonical kmer (including the oritentation of the match).  The provided
 * QueryCache argument will be used to avoid redundant rank / select operations if feasible.
 */

template <typename T>
auto PufferfishBaseIndex<T>::getRefPos(CanonicalKmer& mer, util::QueryCache& qc)
    -> util::ProjectedHits {
    return underlying().getRefPos(mer, qc);
    }

template <typename T>
auto PufferfishBaseIndex<T>::getRefPos(CanonicalKmer& mer) -> util::ProjectedHits {
    return underlying().getRefPos(mer);
}

template <typename T>
uint32_t PufferfishBaseIndex<T>::k() { return underlying().k_; }


template <typename T>
CanonicalKmer PufferfishBaseIndex<T>::getStartKmer(uint64_t rank){
  T& derived = underlying();
  auto& rankSelDict = derived.rankSelDict;
  auto& seq_ = derived.seq_;
  auto k_ = derived.k_;
  CanonicalKmer::k(k_) ;
  CanonicalKmer kb ;
  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict->select(rank)) + 1;
  uint64_t fk = seq_.get_int(2*sp, 2*k_) ;
  kb.fromNum(fk) ;
  return kb ;
}

template <typename T>
CanonicalKmer PufferfishBaseIndex<T>::getEndKmer(uint64_t rank){
  T& derived = underlying();
  auto& rankSelDict = derived.rankSelDict;
  auto& seq_ = derived.seq_;
  auto k_ = derived.k_;
  CanonicalKmer::k(k_) ;
  CanonicalKmer kb ;
  uint64_t contigEnd = rankSelDict->select(rank + 1);

  uint64_t fk = seq_.get_int(2*(contigEnd - k_ + 1), 2*k_) ;
  kb.fromNum(fk) ;
  return kb ;
}

/** Seems not to be used, so skip for now **/
/*
std::vector<CanonicalKmer> PufferfishIndex::getNextKmerOnGraph(uint64_t rank, util::Direction dir, bool isCurContigFwd){
  //get the edge vec
  std::vector<CanonicalKmer> nextKmers ;
  uint8_t edgeVec = edge_[rank] ;
  uint8_t mask = 1 ;
  std::vector<char> nuclmap = {'C','G','T','A','C','G','T','A'} ;
  std::map<char, char> cMap = {{'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}} ;

  if(dir == util::Direction::FORWARD){
    // We need to append so let's concentrate on the lower 4 bits
    auto ke = getEndKmer(rank) ;
    auto ktmp = ke ;
    for(uint8_t i=0; i < 4; ++i){
      ktmp = ke ;
      if(edgeVec & (mask << i)){
        char c = nuclmap[i] ;
        char charToAdd = (isCurContigFwd) ? c : cMap[c] ;
        ktmp.shiftFw(charToAdd) ;
        nextKmers.push_back(ktmp) ;
      }
    }
  }else{
    auto kb = getStartKmer(rank) ;
    auto ktmp = kb ;
    for(uint8_t i=4; i < 8; ++i){
      ktmp = kb ;
      if(edgeVec & (mask << i)){
        char c = nuclmap[i] ;
        char charToAdd = (isCurContigFwd) ? c : cMap[c] ;
        ktmp.shiftBw(charToAdd) ;
        nextKmers.push_back(ktmp) ;
      }
    }
  }
  return nextKmers ;
}
*/

template <typename T>
uint32_t PufferfishBaseIndex<T>::getContigLen(uint64_t rank){
  T& derived = underlying();
  auto& rankSelDict = derived.rankSelDict;
  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict->select(rank)) + 1;
  uint64_t contigEnd = rankSelDict->select(rank + 1);
  return (static_cast<uint32_t>(contigEnd - sp + 1)) ;
}

template <typename T>
uint64_t PufferfishBaseIndex<T>::getGlobalPos(uint64_t rank){
  T& derived = underlying();
  auto& rankSelDict = derived.rankSelDict;
  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict->select(rank)) + 1;
  return sp ;
}

template <typename T>
auto  PufferfishBaseIndex<T>::getContigBlock(uint64_t rank)->util::ContigBlock{
  T& derived = underlying();
  auto& rankSelDict = derived.rankSelDict;
  auto& seq_ = derived.seq_;
  auto  k_ = derived.k_;
  CanonicalKmer::k(k_) ;
  CanonicalKmer kb;
  CanonicalKmer ke;

  uint64_t sp = (rank == 0) ? 0 : static_cast<uint64_t>(rankSelDict->select(rank)) + 1;
  uint64_t contigEnd = rankSelDict->select(rank+1) ;

  uint32_t clen = static_cast<uint32_t>(contigEnd - sp + 1) ;
  uint64_t fk = seq_.get_int(2*sp, 2*k_) ;
  kb.fromNum(fk) ;

  fk = seq_.get_int(2*(contigEnd - k_ + 1), 2*k_) ;
  ke.fromNum(fk) ;

  std::string seq = getSeqStr(sp,clen) ;

  return {rank,sp,clen,seq};
}

/**
 * Return the position list (ref_id, pos) corresponding to a contig.
 */
template <typename T>
const core::range<std::vector<util::Position>::iterator>
PufferfishBaseIndex<T>::refList(uint64_t contigRank) {
  return contigRange(contigRank);
}



template <typename T>
const std::string& PufferfishBaseIndex<T>::refName(uint64_t refRank) {
  return underlying().refNames_[refRank];
}

template <typename T>
uint32_t PufferfishBaseIndex<T>::refLength(uint64_t refRank) const {
  return underlying().refLengths_[refRank];
}

template <typename T>
const std::vector<std::string>& PufferfishBaseIndex<T>::getRefNames() {
  return underlying().refNames_;
}

template <typename T>
const std::vector<uint32_t>& PufferfishBaseIndex<T>::getRefLengths() const {
  return underlying().refLengths_;
}

template <typename T>
uint8_t PufferfishBaseIndex<T>::getEdgeEntry(uint64_t contigRank) const {
    return underlying().edge_[contigRank];
}

template <typename T>
typename PufferfishBaseIndex<T>::seq_vector_t& PufferfishBaseIndex<T>::getSeq() {
  return underlying().seq_;
}

template <typename T>
typename PufferfishBaseIndex<T>::edge_vector_t& PufferfishBaseIndex<T>::getEdge() {
  return underlying().edge_;
}

template class  PufferfishBaseIndex<PufferfishIndex>;
template class  PufferfishBaseIndex<PufferfishSparseIndex>;
template class  PufferfishBaseIndex<PufferfishLossyIndex>;
