#ifndef HIT_COLLECTOR_HPP
#define HIT_COLLECTOR_HPP

#include "Util.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "jellyfish/mer_dna.hpp"
#include "edlib.h"

#include <algorithm>
#include <iostream>
#include <iterator>


#define JUMPSIZE 10

template<typename PufferfishIndexT> class MemCollector {
  using RawHitMap = std::map<uint32_t, std::vector<util::HitQueryPos>>;
public:
  MemCollector(PufferfishIndexT* pfi) : pfi_(pfi) {}


  bool hitsToMappings(PufferfishIndexT& pi,
                      int32_t readLen,
                      uint32_t k,
                      std::vector<std::pair<int, util::ProjectedHits>>& hits,
                      std::vector<util::QuasiAlignment>& qhits,
                      std::vector<std::string>& refBlocks,
                      bool verbose = false) {
    (void)verbose;
    //std::map<uint32_t, MappingInfo> mappings;
    std::map<uint32_t, util::QuasiAlignment> mappings;

    if (hits.empty()) { return false; }

    // find the smallest target set
    std::vector<std::pair<int, util::ProjectedHits>>::iterator minHit;
    uint32_t minSize{std::numeric_limits<uint32_t>::max()};
    for (auto hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
      auto& lab = pi.getEqClassLabel(hitIt->second.contigID());
      if (lab.size() < minSize) {
        minSize = lab.size();
        minHit = hitIt;
      }
    }

    // Put the targets from the smallest set into the map
    auto& lab = pi.getEqClassLabel(minHit->second.contigID());
    for (auto l : lab) {
      mappings[l] = { l, std::numeric_limits<int32_t>::max(), true, static_cast<uint32_t>(readLen) };// 0, true };
      //{l, std::numeric_limits<int32_t>::max(), true, readLen, 0, true}, // the read mapping
                      //{l, std::numeric_limits<int32_t>::max(), true, -1, 0, false}, // the mate mapping
                      //ms, -1 };
    }

    auto mapEnd = mappings.end();
    auto hEnd = hits.end();
    // a mapping must have this many hits to remain active
    size_t currThresh={0};
    for (auto hIter = hits.begin(); hIter != hEnd; ++hIter) {
      ++currThresh;
      auto& projHits = hIter->second;
      for (auto& posIt : projHits.refRange) {
        // if we will produce a mapping for this target
        auto mapIt = mappings.find(posIt.transcript_id());
        if (mapIt != mapEnd) {
          auto refPos = projHits.decodeHit(posIt);
          int offset = refPos.isFW ? -hIter->first : (hIter->first + k) - readLen;
          if (mapIt->second.active and mapIt->second.numHits == currThresh - 1) {
            //if (verbose) {std::cerr << "thresh = " << currThresh << ", nhits = " << mapIt->second.read.numHits << "\n";}
            if (mapIt->second.numHits > 0) {
              auto d = std::abs(mapIt->second.pos - (refPos.pos + offset));
              if (d > readLen / 5) { continue; }
            }
            ++(mapIt->second.numHits);
          } else if (mapIt->second.numHits < currThresh - 1) {
            mapIt->second.active = false;
            continue;
          }
          //if (verbose) { std::cerr << "hit : t = " << pi.refName(mapIt->first) << ", pos = " << refPos.pos + offset << ", fw = " << refPos.isFW << "\n";}
          if (offset < mapIt->second.pos) {
            mapIt->second.pos = refPos.pos + offset;
            mapIt->second.fwd = refPos.isFW;
          }
        }
      }
    }

    for (auto mapIt = mappings.begin(); mapIt != mappings.end(); ++mapIt) {
      if (mapIt->second.numHits != currThresh) {
        mapIt->second.active = false;
      } else {
        qhits.push_back(mapIt->second);
      }
    }

    return true;
  }

  int updateHitMap(util::ProjectedHits& phits, RawHitMap& rawHitMap, int posIn, int32_t readLen, int32_t k){
      int32_t jump{0};


      /**
       * Sigh ... we really need to deal with this right-to-left encoding deal.  It makes
       * the below the reverse of what we would normally expect
       **/
      jump = phits.contigLen_ - phits.contigPos_;
    for(auto& rpos : phits.refRange){
      auto tid = rpos.transcript_id();//pfi_->refName(rpos.transcript_id()) ;
      auto refInfo = phits.decodeHit(rpos) ;

      int offset = 0;//refInfo.isFW ? -posIn : (posIn + k) - readLen;
      rawHitMap[tid].emplace_back(static_cast<uint32_t>(posIn), refInfo.pos + offset, refInfo.isFW) ;

    }
    //int dist =  phits.contigLen_ - pfi_->k() - phits.contigPos_ ;
    return jump ;
  }



  void processRawHitMap(RawHitMap& rawHitMap, std::vector<util::QuasiAlignment>& hits, int readLen){
    for(auto& rh : rawHitMap){
      auto tid = rh.first ;
      std::sort(rh.second.begin(),
                rh.second.end(),
                [](const util::HitQueryPos& a, const util::HitQueryPos& b)-> bool {
                  return a.pos < b.pos ;
                });

      std::map<int32_t , int> hitCov ;
      auto it = rh.second.begin() ;
      while(it != rh.second.end()){
        //int numOfCov = 0 ;
        int32_t hitKey = it->pos - it->queryPos ;
        if(hitCov.count(hitKey) == 0){
          hitCov[hitKey] = 1 ;
        }else{
          hitCov[hitKey] += 1 ;
        }
        ++it ;
      }

      using pair_type = decltype(hitCov)::value_type ;
      auto maxCov = std::max_element(hitCov.begin(), hitCov.end(),
                                       [](const pair_type& p1 , const pair_type& p2) -> bool {
                                         return p1.second > p2.second ;
                                       });
      bool isFwd = rh.second[0].queryFwd ;
      int32_t hitPos = maxCov->first ;

      hits.emplace_back(tid,hitPos,isFwd,readLen) ;
    }
  }


  size_t expandHitEfficient(util::ProjectedHits& hit, pufferfish::CanonicalKmerIterator& kit) {
		auto& allContigs = pfi_->getSeq();
  		//startPos points to the next kmer in contig (which can be the left or right based on the orientation of match)
    size_t cStartPos = hit.globalPos_-hit.contigPos_-1; //next kmer in the read
    size_t cEndPos = cStartPos+hit.contigLen_; 
		size_t cCurrPos = hit.globalPos_; //start from next character if fw match
    if (hit.contigOrientation_) { //if match is fw, go to the next k-mer in the contig
      cCurrPos+=k;
    }
    kit++; //go to the next kmer in the read (In read we always move fw.
		size_t readStart = kit->second;
    pufferfish::CanonicalKmerIterator kit_end;
    bool stillMatch = true;
    CanonicalKmer cmer;
    while (stillMatch && cCurrPos < cEndPos && cCurrPos > cStartPos && kit != kit_end) { // over kmers in contig
      if (hit.contigOrientation_) {
        auto baseCnt = k < cEndPos - cCurrPos? k : cEndPos - cCurrPos;
        uint64_t fk = allContigs.get_int(2*(cCurrPos), 2*baseCnt);
        cmer.fromNum(fk);
        cCurrPos += baseCnt;
        for (size_t i = 0; i < baseCnt && kit != kit_end; i++) {
          auto readCode = (kit->first.fwWord()) >> (2*(k-1)) & 0x3;
          auto contigCode = (fk >> (2*i)) & 0x3;
          if (readCode != contigCode) {
            stillMatch = false;
            break;
          }
          hit.k_++;
          kit++;
        }
      }
      else {
        auto baseCnt = k < cCurrPos - cStartPos? k : cCurrPos - cStartPos;
        uint64_t fk = allContigs.get_int(2*(cCurrPos-baseCnt), 2*baseCnt);
        cmer.fromNum(fk);
        cCurrPos -= baseCnt;
        for (size_t i = baseCnt-1; i >= 0 && kit != kit_end; i--) {
          auto readCode = kit->first.rcWord() & 0x3;
          auto contigCode = (fk >> (2*i)) & 0x3;
          if (readCode != contigCode) {
            stillMatch = false;
            break;
          }
          hit.k_++;
          kit++;
        }
      }
    }
    return readStart;
  }


  /*
  // simple expansion 
  size_t expandHit(util::ProjectedHits& hit, pufferfish::CanonicalKmerIterator& kit) {
		auto& allContigs = pfi_->getSeq();
  		//startPos points to the next kmer in contig (which can be the left or right based on the orientation of match)
		size_t startPos = hit.globalPos_+hit.contigPos_+(hit.contigOrientation_?+1:-1);
		//next kmer in the read
		kit++;
		size_t readStart = kit->second;
		CanonicalKmer mer;
		uint64_t fk = allContigs.get_int(2*startPos, 2*pfi_->k());
		mer.fromNum(fk);
		bool boundaryReached = false;//TODO FIX IT!!
		//simple case: extend the match toward the end of the read
		while(!boundaryReached && kit->first == mer) {
			if (hit.contigOrientation_) {
				startPos++;
				boundaryReached = startPos-hit.globalPos_ >= hit.contigLen_-pfi_->k()+1;
			} else {
				startPos--;
				boundaryReached = startPos < hit.globalPos_;
			}
			//next kmer in contig
			fk = allContigs.get_int(2*startPos, 2*pfi_->k());
			mer.fromNum(fk);
			//next kmer in read
			kit++;
			//increase match length
			hit.k_++;
		}
		//ContigPos_ should be updated since now the match starting position has moved back toward the begining of the contig bc of the rc match
		if (!hit.contigOrientation_)
			hit.contigPos_ = startPos;
	//Expand the match toward the beginning of the read
		//
		//return the new read starting position
		return readStart;
  }
  */

  bool operator()(std::string& read,
                  std::vector<util::QuasiAlignment>& hits,
                  util::MateStatus mateStatus,
                  bool consistentHits,
                  std::vector<std::string>& refBlocks) {
    (void) mateStatus;
    (void) consistentHits;
    using CanonicalKmerIterator = pufferfish::CanonicalKmerIterator ;
    using ProjectedHits = util::ProjectedHits ;
    using QueryCache = util::QueryCache ;

    uint32_t readLen = static_cast<uint32_t>(read.length()) ;
    if(refBlocks.size() != readLen)
        refBlocks.resize(readLen) ;

    ProjectedHits phits ;
    RawHitMap rawHitMap ;

    std::vector<std::pair<int, util::ProjectedHits>> rawHits;

    k = pfi_->k() ;
    CanonicalKmer::k(k);
    CanonicalKmerIterator kit_end;
    CanonicalKmerIterator kit1(read) ;
    QueryCache qc ;

    //string block iterator
    decltype(refBlocks.begin()) bl ;
    //initialize
    bl = refBlocks.begin() ;
    while(kit1 != kit_end) {
      auto phits = pfi_->getRefPos(kit1->first, qc);
      if (!phits.empty()) {
          // kit1 gets updated inside expandHitEfficient function
          size_t readPos = expandHitEfficient(phits, kit1);
          //std::cerr<<"Final match length: " << phits.k_ << "\n";
          rawHits.push_back(std::make_pair(readPos, phits));
      }
      else ++kit1;
    }

    if(rawHits.size() > 0){
      //processRawHitMap(rawHitMap, hits, readLen) ;
      hitsToMappings(*pfi_, readLen, k, rawHits, hits, refBlocks);
      return true ;
    }
    return false ;

  }

private:
  PufferfishIndexT* pfi_ ;
  size_t k ;
  AlignerEngine ae_ ;
};
#endif

