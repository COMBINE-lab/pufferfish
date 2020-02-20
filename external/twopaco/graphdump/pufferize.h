//
// Created by Fatemeh Almodaresi on 8/3/18.
//

#ifndef TWOPACO_PUFFERIZE_H
#define TWOPACO_PUFFERIZE_H

#include <dnachar.h>
#include <junctionapi/junctionapi.h>

constexpr uint16_t ONE = 1;
constexpr uint16_t START = ONE << 8;
constexpr uint16_t END = ONE << 9;
constexpr uint16_t SEEN = ONE << 10;
constexpr uint16_t MERGE_LEFT = ONE << 11;
constexpr uint16_t MERGE_RIGHT = ONE << 12;
constexpr uint16_t COMPLEX = ONE << 13;

struct KmerInfo {
    //<1 bit: complex node (crop both start and end)>
    //<1 bit: merge right (crop end)>
    //<1 bit: merge left (crop start)>
    //<1 bit: set if kmer has already been seen>
    //<1 bit: set if a reference ends with the kmer in forward>
    //<1 bit: set if a reference starts with the kmer in forward>
    //<4 bits: characters that precede kmer in forward>
    //<4 bits: characters that succeed kmer in forward>
    uint16_t kinf{0};
    uint32_t len{0};

    void setPrecedingChar(bool isFw, char c) {
        if (isFw) {
            auto idx = TwoPaCo::DnaChar::MakeUpChar(c);
            kinf |= ONE << (4 + idx);
        } else {
            auto idx = TwoPaCo::DnaChar::MakeUpChar(TwoPaCo::DnaChar::ReverseChar(c));
            kinf |= ONE << (idx);
        }
    }
    void setSucceedingChar(bool isFw, char c) {
        if (isFw) {
            auto idx = TwoPaCo::DnaChar::MakeUpChar(c);
            kinf |= ONE << (idx);
        } else {
            auto idx = TwoPaCo::DnaChar::MakeUpChar(TwoPaCo::DnaChar::ReverseChar(c));
            kinf |= ONE << (4 + idx);
        }
    }

    inline void setStart() {
        kinf |= START;
    }
    inline void setEnd() {
        kinf |= END;
    }
    inline void setSeen() {
        kinf |= SEEN;
    }

    inline bool isStart() {
        return (kinf & START) != 0;
    }

    inline bool isEnd() {
        return (kinf & END) != 0;
    }

    inline bool seen() {
        return (kinf & SEEN) != 0;
    }

    inline void setCropBoth() {
        kinf |= COMPLEX;
    }

    inline bool cropBoth() {
        return  (kinf & COMPLEX) != 0;
    }

    inline void setCropStart() {
        kinf |= MERGE_LEFT;
    }

    inline bool cropStart() {
        return (kinf & MERGE_LEFT) != 0;
    }

    inline void setCropEnd() {
        kinf |= MERGE_RIGHT;
    }

    inline bool cropEnd() {
        return (kinf & MERGE_RIGHT) != 0;
    }

    uint16_t countSucceeding() {
        uint16_t cnt{0};
        for (uint16_t i = 0; i < 4; i++) {
            cnt += (kinf >> i) & ONE;
        }
        //std::cerr << kinf << "  " << cnt << "\n";
        return cnt;
    }

    uint16_t countPreceding() {
        uint16_t cnt{0};
        for (uint16_t i = 4; i < 8; i++) {
            cnt += (kinf >> i) & ONE;
        }
        //std::cerr << kinf << " " << cnt << "\n";
        return cnt;
    }
    void decideType(uint64_t k, uint64_t &approximateContigLen, uint64_t &cntr1, uint64_t &cntr2, uint64_t &cntr3, uint64_t &cntr4) {

        if (!kinf) {// if kmer doesn't exist
            return;
        }
        auto precedeCnt = countPreceding();
        auto succeedCnt = countSucceeding();
        if (precedeCnt > 1 and succeedCnt > 1) {
            setCropBoth();
            cntr1++;
        } else if (succeedCnt > 1) {
            if (isStart()) {
                setCropBoth();
                cntr2++;
            } else {
                setCropStart();
            }

        } else if (precedeCnt > 1) {
            if (isEnd()) {
                setCropBoth();
                cntr3++;
            } else {
                setCropEnd();
            }

        } else if (precedeCnt == 1 and succeedCnt == 1) {
            if (isEnd() and isStart()) {
                setCropBoth();
                cntr4++;
            } else if (isStart()) {
                setCropEnd();
            } else {
                setCropStart();
            }
        } // otherwise, we don't require to crop any nucleotides from any sides of a contig/segment
        if (cropBoth()) approximateContigLen+=(len+k);
        else approximateContigLen+=len;
    }
};

#endif //TWOPACO_PUFFERIZE_H
