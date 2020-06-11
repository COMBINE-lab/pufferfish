//
// Created by Fatemeh Almodaresi on 2019-06-26.
//

#ifndef ROOT_BINARYWRITER_H
#define ROOT_BINARYWRITER_H

#include <cstring>
#include <cerrno>
#include "compact_vector/compact_vector.hpp"

enum EntryType {
    S,P,C,L
};

class GraphdumpBinaryWriter {
private:
    std::ostream* out;
    compact::vector<uint64_t, 2> seqVec_;
    compact::vector<uint64_t, 1> rankVec_;
    uint64_t stringLengthSizeInBytes = 2;

public:
    GraphdumpBinaryWriter(): out(&std::cout) {}
    void setAllowedMaxStringLength(uint64_t allowedIn) {
        std::cerr << "allowedIn: " << allowedIn << "\n";
        if (allowedIn < (1ULL << 8)) {
            stringLengthSizeInBytes = 1;
        } else if (allowedIn < (1ULL << 16)) {
            stringLengthSizeInBytes = 2;
        } else if (allowedIn < (1ULL << 32)) {
            stringLengthSizeInBytes = 4;
        } else {
            stringLengthSizeInBytes = 8;
        }
    }
    //GraphdumpBinaryWriter(std::ostream &outin): out(outin) {}

    void setOStream(std::ostream* o) { out = o; }
    void setCapacity(uint64_t c) {
        seqVec_.reserve(c); seqVec_.clear_mem();
        rankVec_.reserve(c); rankVec_.clear_mem();
    }

    void writeAllowedMaxStringLength() { (*this)<<stringLengthSizeInBytes; }

    void flushSegments(std::string &prefix) {
        std::ofstream seqOut(prefix+"/seq.bin");
        if (seqOut.fail()) {
            std::cerr << "failed to open " << prefix+"/seq.bin" << "; error " << std::strerror(errno) << "\n";
            std::exit(1);
        }
        
        seqVec_.serialize(seqOut);
        
        std::ofstream rankOut(prefix+"/rank.bin");
        if (rankOut.fail()) {
          std::cerr << "failed to open " << prefix + "/rank.bin"
                    << "; error " << std::strerror(errno) << "\n";
          std::exit(1);
        }

        uint64_t ones=0;
        for (uint64_t idx=0; idx < rankVec_.size(); idx++) {
            if (rankVec_[idx]) {
                ones++;
            }
        }
        std::cerr << "# of ones in rank vector: " << ones << "\n";
        rankVec_.serialize(rankOut);
    }

    GraphdumpBinaryWriter& operator<<(const EntryType &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(EntryType));
        return *this;
    }

    GraphdumpBinaryWriter& operator<<(const bool &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }
    GraphdumpBinaryWriter& operator<<(const uint8_t &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }
    GraphdumpBinaryWriter& operator<<(const uint16_t &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }
    GraphdumpBinaryWriter& operator<<(const uint32_t &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }
    GraphdumpBinaryWriter& operator<<(const uint64_t &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }
    GraphdumpBinaryWriter& operator<<(const int64_t &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }


    GraphdumpBinaryWriter &operator<<(const double &inval)  {
        char *inCharPtr = const_cast<char *>(reinterpret_cast<const char *>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }
    GraphdumpBinaryWriter& operator<<(const std::string &inval)  {
        uint64_t tmp = inval.size();
        out->write(reinterpret_cast<const char *>(&tmp), stringLengthSizeInBytes);
        char* inCharPtr = const_cast<char*>(inval.c_str());
        out->write(inCharPtr, inval.size());
        return *this;
    }

    GraphdumpBinaryWriter& operator<<(const compact::vector<uint64_t, 2> &inval) {
        uint64_t w_size = inval.size();
        out->write(reinterpret_cast<char*>(&w_size), sizeof(w_size));
        uint64_t* m_mem = inval.get_words();
        out->write(reinterpret_cast<char*>(m_mem), inval.bytes());
        return *this;
    }

    void addSeq(const std::string &inval) {
        //std::cerr << seqVec_.capacity() << " " << seqVec_.size() << " ";

        for (size_t i = 0; i < inval.size(); ++i) {
            uint16_t c = 0;
            switch (inval[i]) {
                case 'C':
                    c = 1;
                    break;
                case 'G':
                    c = 2;
                    break;
                case 'T':
                    c = 3;
                    break;
                default:
                    c = 0;
            }
            seqVec_.push_back(c);
            rankVec_.push_back(0);
        }
        rankVec_[rankVec_.size()-1] = 1;
    }

};
#endif //ROOT_BINARYWRITER_H
