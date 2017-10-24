#include "OurGFAReader.hpp"
#include "CanonicalKmer.hpp"
#include "cereal/archives/binary.hpp"
#include "xxhash.h"
#include <algorithm>
#include <string>



enum class Direction : bool { PREPEND, APPEND } ;

uint8_t encodeEdge(char c, Direction dir){
  std::map<char,uint8_t> shift_table = {{'A',3}, {'T',2}, {'G',1}, {'C',0}};
  uint8_t val{1} ;
  if(dir == Direction::APPEND)
    return (val << shift_table[c]) ;
  else
    return (val << (shift_table[c]+4));
}



std::vector<std::pair<uint64_t, bool>>
PosFinder::explode(const stx::string_view str, const char& ch) {
  std::string next;
  std::vector<std::pair<uint64_t, bool>> result;
  // For each character in the string
  for (auto it = str.begin(); it != str.end(); it++) {
    // If we've hit the terminal character
    if (*it == '+' or *it == '-') {
      bool orientation = true;
      // If we have some characters accumulated
      // Add them to the result vector
      if (!next.empty()) {
        if (*it == '-') {
          orientation = false;
        }
        result.emplace_back(std::stoll(next), orientation);
        next.clear();
      }
    } else if (*it != ch) {
      // Accumulate the next character into the sequence
      next += *it;
    }
  }
  if (!next.empty())
    result.emplace_back(std::stoll(next),
                        true); // this case shouldn't even happen
  return result;
}

bool PosFinder::is_number(const std::string& s) {
  return !s.empty() && std::find_if(s.begin(), s.end(), [](char c) {
                         return !std::isdigit(c);
                       }) == s.end();
}

// Avoiding un-necessary stream creation + replacing strings with string view
// is a bit > than a 2x win!
// implementation from : https://marcoarena.wordpress.com/tag/string_view/
std::vector<stx::string_view> PosFinder::split(stx::string_view str,
                                               char delims) {
  std::vector<stx::string_view> ret;

  stx::string_view::size_type start = 0;
  auto pos = str.find_first_of(delims, start);
  while (pos != stx::string_view::npos) {
    if (pos != start) {
      ret.push_back(str.substr(start, pos - start));
    }
    start = pos + 1;
    pos = str.find_first_of(delims, start);
  }
  if (start < str.length()) {
    ret.push_back(str.substr(start, str.length() - start));
  }
  return ret;
}

PosFinder::PosFinder(const char* gfaFileName, size_t input_k) {
  filename_ = std::string(gfaFileName);
  std::cerr << "Reading GFA file " << gfaFileName << "\n";
  file.reset(new zstr::ifstream(gfaFileName));
  k = input_k;
}

/*
void PosFinder::scanContigLengths() {
  std::string ln;
  std::string tag, id, value;
  size_t contig_cnt{0};
  size_t contig_len{0};
  size_t ref_cnt{0};
  while(std::getline(*file, ln)) {
    char firstC = ln[0];
    if (firstC != 'S') continue;
    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = split(lnview, '\t');
    if (is_number(id)) {
      contig_len += splited[1].length();
      tag = splited[0].to_string();
      id = splited[1].to_string();
      value = splited[2].to_string();
      contigid2seq[id] = value;
      contig_cnt++;
    }
  }
}
*/

size_t PosFinder::fillContigInfoMap_() {
  std::string ln;
  std::string tag, id, value;
  size_t contig_ctr{0};
  //size_t contig_len{0};
  while (std::getline(*file, ln)) {
    char firstC = ln[0];
    if (firstC != 'S')
      continue;
    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = split(lnview, '\t');
    try {
      auto nid = std::stoull(splited[1].to_string());
      (void)nid;
      auto clen = splited[2].length();
      contigid2seq[nid] = {contig_ctr, 0, static_cast<uint32_t>(clen)};
      ++contig_ctr;
      // contigid2seq[id] = value;
      // contig_cnt++;
    } catch (std::exception& e) {
      // not numeric contig
    }
  }
  size_t total_len{0};
  for (auto& kv : contigid2seq) {
    kv.second.offset = total_len;
    total_len += kv.second.length;
  }
  return total_len;
}

void PosFinder::encodeSeq(sdsl::int_vector<2>& seqVec, size_t offset,
                          stx::string_view str) {
  for (size_t i = 0; i < str.length(); ++i) {
    auto c = my_mer::code(str[i]);
    seqVec[offset + i] = c;
  }
}

sdsl::int_vector<2>& PosFinder::getContigSeqVec() { return seqVec_; }
sdsl::int_vector<8>& PosFinder::getEdgeVec() { return edgeVec_; }



void PosFinder::parseFile() {
  size_t total_len = fillContigInfoMap_();
  file.reset(new zstr::ifstream(filename_));
  std::cerr << "total contig length = " << total_len << "\n";
  std::cerr << "packing contigs into contig vector\n";
  seqVec_ = sdsl::int_vector<2>(total_len, 0);

  std::string ln;
  std::string tag, id, value;
  size_t contig_cnt{0};
  size_t ref_cnt{0};

  k = k + 1 ;
  CanonicalKmer::k(k) ;


  // start and end kmer-hash over the contigs
  // might get deprecated later
  uint64_t maxnid{0} ;


  while (std::getline(*file, ln)) {
    char firstC = ln[0];
    if (firstC != 'S' and firstC != 'P')
      continue;
    stx::string_view lnview(ln);
    std::vector<stx::string_view> splited = split(lnview, '\t');
    tag = splited[0].to_string();
    id = splited[1].to_string();
    // value = splited[2].to_string();
    // A segment line
    if (tag == "S") {
      try {
        uint64_t nid = std::stoll(id);
        if(nid > maxnid)
          maxnid = nid ;
        encodeSeq(seqVec_, contigid2seq[nid].offset, splited[2]);

        // contigid2seq[nid] = value;
      } catch (std::exception& e) {
        // not a numeric contig id
      }
      contig_cnt++;
    }


    // A path line
    if (tag == "P") {

      auto pvalue = splited[2];
      std::vector<std::pair<uint64_t, bool>> contigVec = explode(pvalue, ',');
      //go over the contigVec
      //uint64_t kn{0} ;
      //uint64_t knn{0} ;

      // parse value and add all conitgs to contigVec

      path[ref_cnt] = contigVec;
      refMap.push_back(id);
      ref_cnt++;
      // refMap[ref_cnt] = id;
      // refIDs[id] = ref_cnt++;
    }
  }

  //Initialize edgeVec_
  //bad way, have to re-think
  edgeVec_ = sdsl::int_vector<8>(contig_cnt, 0) ;

  std::map<char, char> cMap = {{'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}} ;
  
  for(auto const& ent: path){
    const std::vector<std::pair<uint64_t, bool>>& contigs = ent.second;

    for(size_t i = 0 ; i < contigs.size() - 1 ; i++){
      auto cid = contigs[i].first ;
      bool ore = contigs[i].second ;
      size_t forder = contigid2seq[cid].fileOrder ;
      auto nextcid = contigs[i+1].first ;
      bool nextore = contigs[i+1].second ;
      // a+,b+ end kmer of a , start kmer of b
      // a+,b- end kmer of a , rc(end kmer of b)
      // a-,b+ rc(start kmer of a) , start kmer of b
      // a-,b- rc(start kmer of a) , rc(end kmer of b)


      //uint64_t kn = (!ore)? (seqVec_.get_int(2 * contigid2seq[cid].offset, 2*k)) : (seqVec_.get_int(2 * (contigid2seq[cid].offset + contigid2seq[cid].length - k), 2 * k)) ;
      uint64_t knn = (nextore)? (seqVec_.get_int(2 * contigid2seq[nextcid].offset, 2*k)) : (seqVec_.get_int(2 * (contigid2seq[nextcid].offset + contigid2seq[nextcid].length - k), 2 * k)) ;


      CanonicalKmer sk, skk ;
      //sk.fromNum(kn) ;
      skk.fromNum(knn) ;

      //validation, to be deprecated later
      std::string kmer, nkmer, k_1mer, nk_1mer;
      //kmer = (ore)?sk.to_str():sk.rcMer() ;
      nkmer = skk.to_str();

      if(!ore and !nextore){
        edgeVec_[forder] |=  encodeEdge(nkmer[0], Direction::PREPEND);
      }else if(!ore and nextore){
        edgeVec_[forder] |=  encodeEdge(cMap[nkmer[k-1]], Direction::PREPEND);
      }else if(ore and nextore){
        edgeVec_[forder] |=  encodeEdge(nkmer[k-1], Direction::APPEND);
      }else if(ore and !nextore){
        edgeVec_[forder] |=  encodeEdge(cMap[nkmer[0]], Direction::APPEND);
      } 
      
      

      /*
      if(ore){
        kmer = sk.to_str() ;
      }else{
        kmer = sk.rcMer().to_str() ;
      }

      if(nextore){
        nkmer = skk.to_str() ;
      }else{
        nkmer = skk.rcMer().to_str() ;
      }

      if( ore != nextore){
        if(!ore)
          std::reverse(kmer.begin(),kmer.end());
        else if(!nextore)
          std::reverse(nkmer.begin(),nkmer.end());
      }

      k_1mer = kmer.substr(1,k-1) ;
      nk_1mer = nkmer.substr(0,k-1) ;

      if(ore == nextore and !ore){
        k_1mer = kmer.substr(0,k-1) ;
        nk_1mer = nkmer.substr(1,k-1) ;
      }
      */

      //std::cerr << ore << "\t" << kmer << "\t" << nextore << "\t" << nkmer << "\n" ;
      //{ std::cerr << "bingo\n" ;}
      //end validation

    }

  }

  k = k - 1;

  std::cerr << " Total # of Contigs : " << contig_cnt
            << " Total # of numerical Contigs : " << contigid2seq.size()
            << "\n\n";
}

// spp::sparse_hash_map<uint64_t, std::string>& PosFinder::getContigNameMap() {
spp::sparse_hash_map<uint64_t, util::PackedContigInfo>&
PosFinder::getContigNameMap() {
  return contigid2seq;
}
spp::sparse_hash_map<std::string, std::string>& PosFinder::getContigIDMap() {
  return seq2contigid;
}

/*
spp::sparse_hash_map<uint32_t, std::string>& PosFinder::getRefIDs() {
  return refMap;
}
*/
std::vector<std::string>& PosFinder::getRefIDs() { return refMap; }

std::map<std::pair<std::string, bool>, bool, util::cmpByPair>&
PosFinder::getPathStart() {
  return pathStart;
}
std::map<std::pair<std::string, bool>, bool, util::cmpByPair>&
PosFinder::getPathEnd() {
  return pathEnd;
}

std::vector<std::pair<std::string, std::string>>& PosFinder::getNewSegments() {
  return newSegments;
}

pufg::Graph& PosFinder::getSemiCG() { return semiCG; }

/*
void PosFinder::writeFile(std::string fileName){
    std::ofstream gfa_file(fileName) ;
    for(auto& cseq : contigid2seq){
        gfa_file << "S" << "\t" << cseq.first <<"\t" << cseq.second << "\n" ;
    }
    for(auto& p : path){
        auto tid = p.first ;
        gfa_file << "P" << "\t" << tid << "\t"  ;
        auto vec = p.second ;
        for(size_t i = 0 ; i < vec.size()-1 ; i++){
            gfa_file << vec[i].first << ((vec[i].second)?"+":"-") << "," ;
        }
        gfa_file << vec[vec.size()-1].first <<
((vec[vec.size()-1].second)?"+":"-") << "\t*\n";

    }
}
*/

void PosFinder::mapContig2Pos() {
  uint64_t pos = 0;
  uint64_t accumPos;
  uint64_t currContigLength = 0;
  uint64_t total_output_lines = 0;
  for (auto const& ent : path) {
    const uint64_t& tr = ent.first;
    const std::vector<std::pair<uint64_t, bool>>& contigs = ent.second;
    accumPos = 0;
    for (size_t i = 0; i < contigs.size(); i++) {
      if (contig2pos.find(contigs[i].first) == contig2pos.end()) {
        contig2pos[contigs[i].first] = {};
        total_output_lines += 1;
      }
      if (contigid2seq.find(contigs[i].first) == contigid2seq.end()) {
        std::cerr << contigs[i].first << "\n";
      }
      pos = accumPos;
      currContigLength = contigid2seq[contigs[i].first].length;
      accumPos += currContigLength - k;
      (contig2pos[contigs[i].first])
          .push_back(util::Position(tr, pos, contigs[i].second));
    }
  }
  std::cerr << "\nTotal # of segments we have position for : "
            << total_output_lines << "\n";
}

void PosFinder::clearContigTable() {
  refMap.clear();
  contig2pos.clear();
}

// Note : We assume that odir is the name of a valid (i.e., existing) directory.
void PosFinder::serializeContigTable(const std::string& odir) {
  std::string ofile = odir + "/ctable.bin";
  std::string eqfile = odir + "/eqtable.bin";
  std::ofstream ct(ofile);
  std::ofstream et(eqfile);
  cereal::BinaryOutputArchive ar(ct);
  cereal::BinaryOutputArchive eqAr(et);
  {
    // We want to iterate over the contigs in precisely the
    // order they appear in the contig array (i.e., the iterator
    // order of contigid2seq).
    std::vector<std::string> refNames;
    refNames.reserve(refMap.size());
    for (size_t i = 0; i < refMap.size(); ++i) {
      refNames.push_back(refMap[i]);
    }
    ar(refNames);

    class VecHasher {
    public:
      size_t operator()(const std::vector<uint32_t>& vec) const {
        return XXH64(const_cast<std::vector<uint32_t>&>(vec).data(),
                     vec.size() * sizeof(decltype(vec.front())), 0);
      }
    };

    // Write out contig offsets and lengths
    /*
    {
      std::vector<util::ContigPosInfo> cpi;
      cpi.reserve(contigid2seq.size());
      for (auto& kv : contigid2seq) {
        cpi.push_back({kv.second.offset, kv.second.length});
      }
      ar(cpi);
    }
    */

    spp::sparse_hash_map<std::vector<uint32_t>, uint32_t, VecHasher> eqMap;
    std::vector<uint32_t> eqIDs;
    std::vector<std::vector<util::Position>> cpos;

    for (auto& kv : contigid2seq) {
      cpos.push_back(contig2pos[kv.first]);
      std::vector<uint32_t> tlist;
      for (auto& p : contig2pos[kv.first]) {
        tlist.push_back(p.transcript_id());
      }
      std::sort(tlist.begin(), tlist.end());
      tlist.erase(std::unique(tlist.begin(), tlist.end()), tlist.end());
      size_t eqID = eqMap.size();
      if (eqMap.contains(tlist)) {
        eqID = eqMap[tlist];
      } else {
        eqMap[tlist] = eqID;
      }
      eqIDs.push_back(eqID);
      // ar(contig2pos[kv.first]);
      /*
      auto& ent = contig2pos[kv.first];
      ct << kv.first;
      for (auto const & ent2 : ent) {
        ct << '\t' << ent2.transcript_id << ',' << ent2.pos << ',' <<
      ent2.orien;
      }
      ct << '\n';
      */
    }
    std::cerr << "there were " << eqMap.size() << " equivalence classes\n";
    eqAr(eqIDs);
    eqIDs.clear();
    eqIDs.shrink_to_fit();
    std::vector<std::vector<uint32_t>> eqLabels;
    eqLabels.reserve(eqMap.size());
    for (auto& kv : eqMap) {
      eqLabels.push_back(kv.first);
    }
    std::sort(eqLabels.begin(), eqLabels.end(),
              [&](const std::vector<uint32_t>& l1,
                  const std::vector<uint32_t>& l2) -> bool {
                return eqMap[l1] < eqMap[l2];
              });
    eqAr(eqLabels);
    ar(cpos);
  }
  /*
    ct << refIDs.size() << '\n';
    for (auto const & ent : refIDs) {
        ct << ent.first << '\t' << ent.second << '\n';
    }
    for (auto const & ent : contig2pos) {
        ct << ent.first;
        for (auto const & ent2 : ent.second) {
            ct << '\t' << ent2.transcript_id << ',' << ent2.pos << ',' <<
    ent2.orien;
        }
        ct << '\n';
    }
  */
}

void PosFinder::deserializeContigTable() {
  // TODO read the file in the same order as you've written it down.
}
/*
int main(int argc, char* argv[]) {
    //std::vector<std::string> gfa_file = {argv[1]};
    //std::vector<std::string> read_file = {argv[2]};
    std::cerr<<" Starting ... \n";
    PosFinder pf(argv[1], 30);
    pf.parseFile();

    pf.mapContig2Pos();

    std::ofstream ct("contig_table.tsv");
    for (auto const & ent : pf.contig2pos) {
        ct << ent.first;
        for (auto const & ent2 : ent.second) {
            ct << '\t' << ent2.transcript_id << ',' << ent2.pos << ',' <<
ent2.orien;
        }
        ct << '\n';
    }
    ct.close();
}
*/
