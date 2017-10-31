#ifndef GRAPH_TRAVERSER_HPP
#define GRAPH_TRAVERSER_HPP


#include "Util.hpp"

enum Task {
           SUCCESS,
           FAILURE
};


template <typename PufferfishIndexT> class GraphTraverser {
public:
  Task doBFS(size_t tid, size_t tpos, bool moveFw, util::ContigBlock& curContig, size_t startp, util::ContigBlock& endContig, size_t endp, uint32_t threshold, std::string& seq) {
    if (curContig.contigIdx_ == endContig.contigIdx_) {
      if (!endContig.isDummy()) {
        append(seq, curContig, startp, endp, moveFw);
        return Task::SUCCESS;
      }
      else { // if the end contig is a dummy
        
      }
    }

    if (!isCompatible(curContig, tid, tpos)) return Task::FAILURE;

    return Task::SUCCESS;
 }


private:
  std::vector<util::ContigBlock> fetchSuccessors(util::ContigBlock& contig, bool moveFw) { return NULL;}
  bool isCompatible(util::ContigBlock& contig, size_t tid, size_t tpos) { return false;}
  size_t distance(size_t startp, size_t endp, bool moveFw) {if(moveFw)return endp-startp; else return startp-endp;}
  void append(std::string& seq, util::ContigBlock& contig, size_t startp, size_t endp, bool moveFw) {
    if(moveFw)
      seq += contig.getsubstr(startp, endp-startp);
    else {
      seq += rc(contig.getsubstr(endp, startp-endp)); // we are always building the seq by moving forward in transcript, so we always append any substring that we construct
    }
  }
  void cutoff(std::string& seq, size_t len) {
    seq = seq.substr(0, seq.length()-len);
  }

  std::string& rc(std::string& str) {
    for (int i = 0; i < str.length()/2; i++) {
      char tmp = str[i];
      str[i] = rev(str[str.length()-1-i]);
      str[str.length()-1-i] = rev(tmp);
    }
    return str;
  }

  char rev(const char& c) {
    switch(c){
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    }
  }
};

#endif
