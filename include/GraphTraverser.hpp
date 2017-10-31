#ifndef GRAPH_TRAVERSER_HPP
#define GRAPH_TRAVERSER_HPP


#include "Util.hpp"

enum Task {
           SUCCESS,
           FAILURE
};


template <typename PufferfishIndexT> class GraphTraverser {
public:
  Task doBFS(size_t tid, size_t tpos, bool moveFw, util::ContigBlock curContig, size_t startp, util::ContigBlock endContig, size_t endp, uint32_t threshold, std::string& seq) {}


private:
  std::vector<util::ContigBlock> fetchSuccessors(util::ContigBlock& contig, bool moveFw) { return NULL;}
  bool isCompatible(util::ContigBlock& contig, size_t tid, size_t tpos) { return false;}

};

#endif
