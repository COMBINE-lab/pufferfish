#ifndef HIT_COLLECTOR_HPP
#define HIT_COLLECTOR_HPP

#include "Util.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "CanonicalKmer.hpp"
#include "CanonicalKmerIterator.hpp"
#include "jellyfish/mer_dna.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>


#define JUMPSIZE 10

template<typename PufferfishIndexT> class HitCollector {
  using RawHitMap = std::map<uint32_t, std::vector<util::HitQueryPos>>;
public:
  HitCollector(PufferfishIndexT* pfi) : pfi_(pfi) {}


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

      int offset = refInfo.isFW ? -posIn : (posIn + k) - readLen;
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


  bool operator()(std::string& read,
                  std::vector<util::QuasiAlignment>& hits,
                  util::MateStatus mateStatus,
                  bool consistentHits = false) {
    (void) mateStatus;
    (void) consistentHits;
    //using QuasiAlignment = util::QuasiAlignment ;
    //using MateStatus = util::MateStatus ;
    using CanonicalKmerIterator = pufferfish::CanonicalKmerIterator ;
    using ProjectedHits = util::ProjectedHits ;
    using QueryCache = util::QueryCache ;
    //using HitQueryPos = util::HitQueryPos ;
    //using RawHitMap = std::map<std::string,std::vector<HitQueryPos> > ;


    //read is same as seq in @Rob's mapper
    int32_t readLen = static_cast<int32_t>(read.length()) ;

    ProjectedHits phits ;
    RawHitMap rawHitMap ;

    //size of the kmer
    k = pfi_->k() ;
    int32_t ks = static_cast<int32_t>(k);
    CanonicalKmer::k(pfi_->k()) ;
    //end of the kmer
    CanonicalKmerIterator kit_end ;
    //start of kmer

    //print log
    //std::cerr << "\n In hit collector now with k " << k << "\n" ;


    CanonicalKmerIterator kit1(read) ;
    QueryCache qc ;
    bool lastSearch{false} ;
    bool done{false} ;
    size_t x{0} ;


    //I am rewriting my mapper to match
    //Rob's convension, the += operator in
    //my code is not working as I imagined.

    while(kit1 != kit_end and !done) {
    ++x;
    //if (x % 1000 == 0) {
    //std::cerr << "read_num = " << read_num << ", x = " << x << ", pos = " << kit1->second << "\n";
      //}
      //get position
    auto phits = pfi_->getRefPos(kit1->first, qc);
    if (!phits.empty()) {

      //@rob stores the quesry position and the phits pair in
      //a vector, I used a processed hit map
      //they are almost same with some suttle differences
      //where I used the selective alignment kind of idea to
      //refine reads.
      //hits.push_back(std::make_pair(kit1->second, phits));
      // we looked with the fw version of this k-mer
      size_t jump{0} ;
      bool queryFW = kit1->first.isFwCanonical();
      // the k-mer we queried with hit the contig in the forward orientation
      bool hitFW = phits.contigOrientation_;
      //bool ore = (queryFW == hitFW) ;



    //this line has to be removed later
    jump = phits.contigLen_ - phits.contigPos_;
      //this is just a check let's keep it here for now

      if (queryFW == hitFW) { // jump to the end of the contig
        if (jump <= 0) {
         // std::cerr << "(1) rnum = " << read_num
        std::cerr   << ", queryFW = " << queryFW
                    << ", hitFW = " << hitFW
                    << ", jump = " << jump
                    << ", pos = " << kit1->second
                    << ", x = " << x
                    << ", kmer = " << read.substr(kit1->second, k)
                    << ", read = " <<  read << "\n"; std::exit(1);
        }
      } else {
        // k-mer is RC, but is fw on contig = read is rc on contig
        jump = phits.contigLen_ - phits.contigPos_;
        if (jump <= 0) {
        //  std::cerr << "(2) rnum = " << read_num
            std::cerr<< ", queryFW = " << queryFW
                    << ", hitFW = " << hitFW
                    << ", jump = " << jump
                    << ", pos = " << kit1->second
                    << ", x = " << x << ", kmer = " << read.substr(kit1->second, k) << ", read = " <<  read << "\n"; std::exit(1);}
      }

      //update the raw hit map
      //and get the next queryPos
      jump = updateHitMap(phits, rawHitMap, kit1->second, readLen, k) ;

      // the position where we should look
      if (lastSearch or done) { done = true; continue; }
      int32_t newPos = kit1->second + jump;
      if (newPos > readLen - k) {
        lastSearch = true;
        newPos = readLen - k;
      }
      kit1.jumpTo(newPos);
      continue;
    } else {
      ++kit1;
    }
  }

    /* old implementation
     * which was also based
     * on the logic of jump
    while(kit1 != kit_end){

      std::cerr << "\nGoing over k-mer " << kit1->first.to_str() << "\n" ;
      phits = pfi_->getRefPos(kit1->first, qc) ;
      if(!phits.empty())
        break ;
      ++kit1 ;
    }

    //TODO load the sequence and
    //properly extend the read to find
    //next safe jump without SNP
    //have to change projectedHits

    int counter_check = 0 ;

    if(kit1 != kit_end){

      if (counter_check == 20){
    	  std::exit(1) ;
      }else{
    	  counter_check++ ;
      }
      std::cerr << "\n in the loop forever \n" ;
      auto queryPos = kit1->second ;
      auto dist = updateHitMap(phits, rawHitMap, queryPos) ;
      //Either I jump by 10 or end of Contig ;
      //
      int advanceBy = std::min(JUMPSIZE, dist) ;
      kit1 += advanceBy ;
      while(kit1 != kit_end and kit1.kmerIsValid()){
          if (counter_check == 20){
              std::exit(1) ;
          }else{
              counter_check++ ;
          }
    	std::cerr << "\n in the inner loop forever \n" ;
        queryPos = kit1->second ;
        phits = pfi_->getRefPos(kit1->first,qc);
        dist = updateHitMap(phits, rawHitMap, queryPos) ;
        advanceBy = std::min(JUMPSIZE, dist) ;
        kit1 += advanceBy ;
      }
    }*/


    if(rawHitMap.size() > 0){
      processRawHitMap(rawHitMap, hits, readLen) ;
      return true ;
    }

    return false ;



  }

private:
  PufferfishIndexT* pfi_ ;
  int k ;
};
#endif
