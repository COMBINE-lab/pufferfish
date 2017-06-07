#ifndef _UTIL__H
#define _UTIL__H

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "string_view.hpp"
#include "core/range.hpp"

#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"

namespace util {

  class IndexOptions {
  public:
    uint32_t k{31};
    std::string gfa_file;
    std::string cfile;
    std::string rfile;
    std::string outdir;
  };

  class TestOptions {
  public:
  };

  class ValidateOptions {
  public:
    std::string indexDir;
    std::string refFile;
  };

  // For the time being, assume < 4B contigs
  // and that each contig is < 4B bases
  struct Position {
    //std::string transcript_id;
    uint32_t transcript_id_;
    uint32_t pos_;
    //bool orien;
    Position() {
      transcript_id_ = std::numeric_limits<decltype(transcript_id_)>::max();
      pos_ = std::numeric_limits<decltype(pos_)>::max();
    }

    Position(uint32_t tid, uint32_t tpos, bool torien) {
      transcript_id_ = tid;
      pos_ = tpos;
      setOrientation(torien);
      //orien = torien;
    }

    void setOrientation(bool orientation) {
      if (orientation) {
        pos_ |= 1 << 31;
      } else {
        pos_ &= 0x7FFFFFFF;
      }
    }

    inline uint32_t transcript_id() { return transcript_id_; }
    inline uint32_t pos() { return (pos_ & 0x7FFFFFFF); }
    inline bool orientation() { return (pos_ & 0x7FFFFFFF) >> 31; }

    template <class Archive>
    void serialize(Archive& ar) {
      ar(transcript_id_, pos_);
    }

  private:
    //uint32_t orientMask_
  };

  // Structure to hold a list of "projected" (i.e. reference) hits
  // for a k-mer
  struct ProjectedHits {
    uint32_t conigPos_;
    bool contigOrientation_; // true for fw, false for rc
    core::range<std::vector<util::Position>::iterator> refRange;
  };

  char complement(char& c);
  std::string revcomp(std::string s);
  std::vector<std::pair<std::string, bool> > explode(const stx::string_view str, const char& ch);
	bool is_number(const std::string& s);
  // Avoiding un-necessary stream creation + replacing strings with string view
  // is a bit > than a 2x win!
  // implementation from : https://marcoarena.wordpress.com/tag/string_view/
  std::vector<stx::string_view> split(stx::string_view str, char delims);
}

#endif // _UTIL__H
