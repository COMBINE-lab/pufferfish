//
// Created by Fatemeh Almodaresi on 2019-06-26.
//

#ifndef ROOT_BINARYWRITER_H
#define ROOT_BINARYWRITER_H

#include "compact_vector/compact_vector.hpp"

enum EntryType {
    S,P,C,L
};

class BinaryWriter {
private:
    std::ostream* out;
    compact::vector<uint64_t, 2> seqVec_;
    compact::vector<uint64_t, 1> rankVec_;

public:
    BinaryWriter(): out(&std::cout) {}
    //BinaryWriter(std::ostream &outin): out(outin) {}

    void setOStream(std::ostream* o) { out = o; }
    void setCapacity(uint64_t c) {
        seqVec_.reserve(c);
        rankVec_.reserve(c);
    }

    void flushSegments(std::string &prefix) {
        std::ofstream seqOut(prefix+"/seq.bin");
        seqVec_.serialize(seqOut);
        std::ofstream rankOut(prefix+"/rank.bin");
        uint64_t ones=0;
        for (uint64_t idx=0; idx < rankVec_.size(); idx++) {
            if (rankVec_[idx]) {
                ones++;
            }
        }
        std::cerr << "# of ones in rank vector: " << ones << "\n";
        rankVec_.serialize(rankOut);
    }

    BinaryWriter& operator<<(const EntryType &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(EntryType));
        return *this;
    }

    BinaryWriter& operator<<(const bool &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }
    BinaryWriter& operator<<(const uint8_t &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }
    BinaryWriter& operator<<(const uint16_t &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }
    BinaryWriter& operator<<(const uint32_t &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }
    BinaryWriter& operator<<(const uint64_t &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }
    BinaryWriter& operator<<(const int64_t &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }


    BinaryWriter &operator<<(const double &inval)  {
        char *inCharPtr = const_cast<char *>(reinterpret_cast<const char *>(&inval));
        out->write(inCharPtr, sizeof(inval));
        return *this;
    }
    BinaryWriter& operator<<(const std::string &inval)  {
        if (inval.size() >= 0x100) {
            std::cerr << "ERROR!! DOESN'T SUPPORT STRING LENGTH LONGER THAN 255. String length: "
                      << inval.size() << "\n";
            std::exit(1);
        }
        auto tmp = static_cast<uint16_t>(inval.size());
        out->write(reinterpret_cast<const char *>(&tmp), sizeof(tmp));
        char* inCharPtr = const_cast<char*>(inval.c_str());
        out->write(inCharPtr, inval.size());
        return *this;
    }

    BinaryWriter& operator<<(const compact::vector<uint64_t, 2> &inval) {
        uint64_t w_size = inval.size();
        out->write(reinterpret_cast<char*>(&w_size), sizeof(w_size));
        uint64_t* m_mem = inval.get_words();
        out->write(reinterpret_cast<char*>(m_mem), inval.bytes());
        return *this;
    }

    void addSeq(const std::string &inval) {
//        std::cerr << seqVec_.capacity() << " " << seqVec_.size() << " ";

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
