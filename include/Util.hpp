#ifndef _UTIL__H
#define _UTIL__H

#include "core/range.hpp"
#include "string_view.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <vector>


#include "cereal/types/string.hpp"
#include "cereal/types/vector.hpp"
#include "jellyfish/mer_dna.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"


#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

namespace util {
using namespace std;

  static constexpr int8_t rc_table[128] = {
    78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 15
    78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 31
    78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 787
    78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 63
    78, 84, 78, 71, 78,  78,  78, 67, 78, 78, 78, 78,  78, 78, 78, 78, // 79
    78, 78,  78, 78,  65, 65, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 95
    78, 84, 78, 71, 78,  78,  78, 67, 78, 78, 78, 78,  78, 78, 78, 78, // 101
    78, 78,  78, 78,  65, 65, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78  // 127
  };


// Adapted from
// https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/8c9933a1685e0ab50c7d8b7926c9068bc0c9d7d2/src/main.c#L36
 inline void reverseRead(std::string& seq,
         std::string& qual,
         std::string& readWork,
         std::string& qualWork) {

     readWork.resize(seq.length(), 'A');
     qualWork.resize(qual.length(), 'I');
     int32_t end = seq.length()-1, start = 0;
     //readWork[end] = '\0';
     //qualWork[end] = '\0';
     while (LIKELY(start < end)) {
         readWork[start] = (char)rc_table[(int8_t)seq[end]];
         readWork[end] = (char)rc_table[(int8_t)seq[start]];
         qualWork[start] = qual[end];
         qualWork[end] = qual[start];
         ++ start;
         -- end;
     }
     // If odd # of bases, we still have to complement the middle
     if (start == end) {
         readWork[start] = (char)rc_table[(int8_t)seq[start]];
         // but don't need to mess with quality
         // qualWork[start] = qual[start];
     }
     //std::swap(seq, readWork);
     //std::swap(qual, qualWork);
 }

 // Adapted from
 // https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/8c9933a1685e0ab50c7d8b7926c9068bc0c9d7d2/src/main.c#L36
 // Don't modify the qual
 inline void reverseRead(std::string& seq,
         std::string& readWork) {

     readWork.resize(seq.length(), 'A');
     int32_t end = seq.length()-1, start = 0;
     //readWork[end] = '\0';
     //qualWork[end] = '\0';
     while (LIKELY(start < end)) {
         readWork[start] = (char)rc_table[(int8_t)seq[end]];
         readWork[end] = (char)rc_table[(int8_t)seq[start]];
         ++ start;
         -- end;
     }
     // If odd # of bases, we still have to complement the middle
     if (start == end) {
         readWork[start] = (char)rc_table[(int8_t)seq[start]];
         // but don't need to mess with quality
         // qualWork[start] = qual[start];
     }
     //std::swap(seq, readWork);
     //std::swap(qual, qualWork);
 }

 inline std::string reverseComplement(std::string& seq) {
     std::string work;
     reverseRead(seq, work);
     return work;
 }

struct cmpByPair {
  bool operator()(std::pair<uint64_t, bool> a,
                  std::pair<uint64_t, bool> b) const {
    return (a.first != b.first) ? (a.first < b.first) : (a.second < b.second);
  }
};

// We need a wraper that will provide a "type" field
template <typename> struct Void { typedef void type; };

// The "default" type T does not have a key
template <typename T, typename Sfinae = void>
struct has_key : std::false_type {};

// Use decltype to access the key_type field of the
// provided type T if it exists.  If T has a key_type trait
// then this will succeed, and we can instantiate this struct.
// If T doesn't have a key_type trait, then SFINAE means that that
// struct won't be instantiated and we'll fall back to the default
// above.
template <typename T>
struct has_key<T, typename Void<decltype(typename T::key_type{})>::type>
    : std::true_type {};

/*
Print a map-like container
 */
template <typename T, bool b> struct ContainerPrinter {
  static string str(T& container) {
    ostringstream s;
    /*
     Loop over and print out the map using the new c++11
     style for loop.  The new-style loop will work for any
     type that has begin() and end() methods that return an iterator,
     OR if there are _free_ functions begin() and end() which
     consume the type and return an iterator.
    */
    s << "[";
    for (auto& e : container) {
      s << e.first << " : " << e.second << ", ";
    }
    s.seekp(-2, ios_base::end);
    s << "]";
    return s.str();
  }
};

/*
Print a list-like container
 */
template <typename T> struct ContainerPrinter<T, false> {
  static string str(T& container) {
    ostringstream s;
    s << "[";
    for (auto& e : container) {
      s << e << ", ";
    }
    s.seekp(-2, ios_base::end);
    s << "]";
    return s.str();
  }
};

/*
Compile-time selection between list-like and map-like printing.
 */
template <typename T> string str(T& container) {
  return ContainerPrinter<T, has_key<T>::value>::str(container);
}

class IndexOptions {
public:
  uint32_t k{31};
  std::string gfa_file;
  std::string cfile;
  std::string rfile;
  std::string outdir;
  bool isSparse{false};
  uint32_t extensionSize{4};
  uint32_t sampleSize{9};
};

class TestOptions {
public:
};

class ValidateOptions {
public:
  std::string indexDir;
  std::string refFile;
};

class AlignmentOpts{
	public:

    std::string indexDir;
    std::string read1;
    std::string read2;
    std::string unmatedReads;
    uint32_t numThreads{1};
    uint32_t maxNumHits{200};
    std::string outname;
    double quasiCov{0.0};
    bool pairedEnd{false};
    bool noOutput{false};
    bool sensitive{false};
    bool strictCheck{false};
    bool fuzzy{false};
    bool consistentHits{false};
    bool quiet{false};
	bool writeOrphans{false} ;
};




//Mapped object contains all the information
//about mapping the struct is a bit changed from
//quasi mapping
enum class MateStatus : uint8_t {
        SINGLE_END = 0,
        PAIRED_END_LEFT = 1,
        PAIRED_END_RIGHT = 2,
        PAIRED_END_PAIRED = 3 };


  struct MemInfo {
    size_t tpos;
    uint32_t rpos;
    uint32_t memlen;

    MemInfo(size_t tposIn, uint32_t rposIn, uint32_t memlenIn) :
      tpos(tposIn), rpos(rposIn), memlen(memlenIn){}
  };

  struct MemCluster {
    size_t tid;
    bool ori;
    std::vector<MemInfo> mems;

    MemCluster(size_t tidIn, bool oriIn, size_t tposIn, uint32_t rposIn, uint32_t memlenIn) :
      tid(tidIn), ori(oriIn){
      mems.emplace_back(tposIn, rposIn, memlenIn);
    }
  };

struct QuasiAlignment {
  	QuasiAlignment() :
    tid(std::numeric_limits<uint32_t>::max()),
		pos(std::numeric_limits<int32_t>::max()),
		fwd(true),
		fragLen(std::numeric_limits<uint32_t>::max()),
		readLen(std::numeric_limits<uint32_t>::max()),
		isPaired(false){}

        QuasiAlignment(uint32_t tidIn, int32_t posIn,
                bool fwdIn, uint32_t readLenIn,
                uint32_t fragLenIn = 0,
                bool isPairedIn = false) :
            tid(tidIn), pos(posIn), fwd(fwdIn),
            fragLen(fragLenIn), readLen(readLenIn),
            isPaired(isPairedIn){}

        QuasiAlignment(QuasiAlignment&& other) = default;
        QuasiAlignment& operator=(QuasiAlignment&) = default;
        QuasiAlignment& operator=(QuasiAlignment&& o) = default;
        QuasiAlignment(const QuasiAlignment& o) = default;
        QuasiAlignment(QuasiAlignment& o) = default;

        // Some convenience functions to allow salmon interop
/*
#ifdef RAPMAP_SALMON_SUPPORT
        inline uint32_t transcriptID() const { return tid; }
        inline double score() { return 1.0; }
        inline uint32_t fragLength() const { return fragLen; }

        inline uint32_t fragLengthPedantic(uint32_t txpLen) const {
            if (mateStatus != rapmap::utils::MateStatus::PAIRED_END_PAIRED
                or fwd == mateIsFwd) {
                return 0;
            }
            int32_t p1 = fwd ? pos : matePos;
            p1 = (p1 < 0) ? 0 : p1;
            p1 = (p1 > txpLen) ? txpLen : p1;
            int32_t p2 = fwd ? matePos + mateLen : pos + readLen;
            p2 = (p2 < 0) ? 0 : p2;
            p2 = (p2 > txpLen) ? txpLen : p2;

            return (p1 > p2) ? p1 - p2 : p2 - p1;
        }

        inline int32_t hitPos() { return std::min(pos, matePos); }
        double logProb{HUGE_VAL};
        double logBias{HUGE_VAL};
        inline LibraryFormat libFormat() { return format; }
        LibraryFormat format;
#endif // RAPMAP_SALMON_SUPPORT
*/
        // Only 1 since the mate must have the same tid
        // we won't call *chimeric* alignments here.
        uint32_t tid;
        // Left-most position of the hit
        int32_t pos;
        // left-most position of the mate
        int32_t matePos;
        // Is the read from the forward strand
        bool fwd;
        // Is the mate from the forward strand
        bool mateIsFwd;
        // The fragment length (template length)
        // This is 0 for single-end or orphaned reads.
        uint32_t fragLen;
        // The read's length
        uint32_t readLen;
        // The mate's length
        uint32_t mateLen;
        // Is this a paired *alignment* or not
        bool isPaired;
        MateStatus mateStatus;
  bool active = true;
  uint32_t numHits = 0;
 };

// from https://github.com/cppformat/cppformat/issues/105
 class FixedBuffer : public fmt::Buffer<char> {
        public:
            FixedBuffer(char *array, std::size_t size)
                : fmt::Buffer<char>(array, size) {}

        protected:
            void grow(std::size_t size) {
                (void) size;
                throw std::runtime_error("buffer overflow");
            }
    };

class FixedWriter : public fmt::Writer {
        private:
            FixedBuffer buffer_;
        public:
            FixedWriter(char *array, std::size_t size)
                : fmt::Writer(buffer_), buffer_(array, size) {}
    };


// For the time being, assume < 4B contigs
// and that each contig is < 4B bases
struct Position {
  // std::string transcript_id;
  uint32_t transcript_id_;
  uint32_t pos_;
  // bool orien;
  Position() {
    transcript_id_ = std::numeric_limits<decltype(transcript_id_)>::max();
    pos_ = std::numeric_limits<decltype(pos_)>::max();
  }

  Position(uint32_t tid, uint32_t tpos, bool torien) {
    transcript_id_ = tid;
    pos_ = tpos;
    setOrientation(torien);
    // orien = torien;
  }

  //The most significant bit carry
  //the orientation information

  void setOrientation(bool orientation) {
    if (orientation) {
      pos_ |= 1 << 31;
    } else {
      pos_ &= 0x7FFFFFFF;
    }
  }

  inline uint32_t transcript_id() { return transcript_id_; }
  inline uint32_t pos() { return (pos_ & 0x7FFFFFFF); }
  inline bool orientation() { return (pos_ & 0x80000000); }

  template <class Archive> void serialize(Archive& ar) {
    ar(transcript_id_, pos_);
  }

private:
  // uint32_t orientMask_
};

//struct HitPos
struct HitQueryPos {
  HitQueryPos(uint32_t queryPosIn, uint32_t posIn, bool queryFwdIn) :
	queryPos(queryPosIn) , pos(posIn), queryFwd(queryFwdIn) {}

  uint32_t queryPos, pos ;
  bool queryFwd;

};

struct QueryCache {
  uint64_t prevRank;
  uint64_t contigStart;
  uint64_t contigEnd;
};

struct ContigPosInfo {
  size_t offset_;
  uint32_t length_;

  inline size_t offset() { return offset_; }
  inline uint32_t length() { return length_; }

  template <class Archive> void serialize(Archive& ar) { ar(offset_, length_); }
};

struct PackedContigInfo {
  size_t fileOrder;
  size_t offset;
  uint32_t length;
};

struct RefPos {
  uint32_t pos;
  bool isFW;
};

struct HitCounters {
        std::atomic<uint64_t> peHits{0};
        std::atomic<uint64_t> seHits{0};
        std::atomic<uint64_t> trueHits{0};
        std::atomic<uint64_t> totHits{0};
        std::atomic<uint64_t> numReads{0};
        std::atomic<uint64_t> tooManyHits{0};
        std::atomic<uint64_t> lastPrint{0};
};


// Structure to hold a list of "projected" (i.e. reference) hits
// for a k-mer
struct ProjectedHits {
  uint32_t contigIdx_;
  // The relative position of the k-mer inducing this hit on the
  // contig
  uint64_t globalPos_ ;

  uint32_t contigPos_;
  // How the k-mer inducing this hit maps to the contig
  // true for fw, false for rc
  bool contigOrientation_;
  uint32_t contigLen_;
  uint32_t k_;
  core::range<std::vector<util::Position>::iterator> refRange;

  inline bool empty() { return refRange.empty(); }

  inline uint32_t contigID() const { return contigIdx_; }
  //inline uint64_t getGlobalPos() const { return globalPos_; }
  inline RefPos decodeHit(util::Position& p) {
    // true if the contig is fowrard on the reference
    bool contigFW = p.orientation();
    // we are forward with respect to the reference if :
    // (1) contigFW and contigOrientation_
    // (2) !contigFW and !contigOrientation_
    // we are reverse complement with respect to the reference if :
    // (3) configFW and !contigOrientation_
    // (4) !configFW and contigOrientation_

    // if we're in the forward orientation, then our position is
    // just the contig offset plus or relative position
    uint32_t rpos{0};
    bool rfw{false};
    if (contigFW and contigOrientation_) {
      // kmer   :          AGC
      // contig :      ACTTAGC
      // ref    :  GCA[ACTTAGC]CA
      rpos = p.pos() + contigPos_;
      rfw = true;
    } else if (contigFW and !contigOrientation_) {
      // kmer   :          GCT
      // contig :      ACTTAGC
      // ref    :  GCA[ACTTAGC]CA
      rpos = p.pos() + contigPos_;
      rfw = false;
    } else if (!contigFW and contigOrientation_) {
      // kmer   :          AGT
      // contig :      GCTAAGT
      // ref    :  GCA[ACTTAGC]CA
      rpos = p.pos() + contigLen_ - (contigPos_ + k_);
      rfw = false;
    } else if (!contigFW and !contigOrientation_) {
      // kmer   :          ACT
      // contig :      GCTAAGT
      // ref    :  GCA[ACTTAGC]CA
      rpos = p.pos() + contigLen_ - (contigPos_ + k_);
      rfw = true;
    }
    return {rpos, rfw};
  }
};

char complement(char& c);
std::string revcomp(std::string s);
bool isRevcomp(std::string s);
std::vector<std::pair<uint64_t, bool>> explode(const stx::string_view str,
                                                  const char& ch);
bool is_number(const std::string& s);
std::vector<std::string> tokenize(const std::string& s, char delim) ;
// Avoiding un-necessary stream creation + replacing strings with string view
// is a bit > than a 2x win!
// implementation from : https://marcoarena.wordpress.com/tag/string_view/
std::vector<stx::string_view> split(stx::string_view str, char delims);

}

#endif // _UTIL__H
