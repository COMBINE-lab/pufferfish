#ifndef __MAPPINGS_H__
#define __MAPPINGS_H__

#include "spdlog/spdlog.h"
#include <fstream>
#include <iostream>
#include "taxa.h"

struct ReadInfo {
    std::string rid;
    size_t cnt;
    uint32_t len;
    std::vector<TaxaNode> mappings;
};


class Mappings {
    public:
        Mappings() {}
        void load(std::string mfileName,
                std::shared_ptr<spdlog::logger> loggerIn) {
            logger = loggerIn;
            inFile.open(mfileName, std::ios::binary);
            if (!readHeader()) {
                logger->error("Invalid header for mapping output file.");
                std::exit(1);
            }
        }

        bool hasNext() {return  inFile.is_open() && inFile.good();}
        ReadInfo& nextRead() {
            uint32_t puff_id, mcnt; 
            uint8_t rlen, lcnt, rcnt, ibeg, ilen;
                
            uint8_t readNameLength;
            // read name
            inFile.read(reinterpret_cast<char*>(&readNameLength), sizeof(readNameLength));
            char* strChar = new char[readNameLength];
            inFile.read(strChar, readNameLength);
            readInfo.rid = std::string(strChar, readNameLength);
            // mapping count
            inFile.read(reinterpret_cast<char*>(&mcnt), sizeof(uint32_t));
            inFile.read(reinterpret_cast<char*>(&rlen), sizeof(uint8_t)); // left read len
            readInfo.cnt = mcnt;
            readInfo.len = rlen;
            if (isPaired) {
                inFile.read(reinterpret_cast<char*>(&rlen), sizeof(uint8_t)); // right read len
                readInfo.len += rlen;
            }
            // std::cout << readInfo.rid << " " << readInfo.cnt << " " << readInfo.len << "\n";
            readInfo.mappings.reserve(readInfo.cnt);
            for (size_t mappingCntr = 0; mappingCntr < readInfo.cnt; mappingCntr++) {
                inFile.read(reinterpret_cast<char*>(&puff_id), sizeof(uint32_t));
                // fetch the taxon from the map
                readInfo.mappings.push_back(puff_id);
                TaxaNode& taxaPtr = readInfo.mappings.back();
                
                inFile.read(reinterpret_cast<char*>(&lcnt), sizeof(lcnt));
                if (isPaired)
                    inFile.read(reinterpret_cast<char*>(&rcnt), sizeof(rcnt));
                
                for (size_t i = 0; i < lcnt; ++i) {
                    inFile.read(reinterpret_cast<char*>(&ibeg), sizeof(ibeg));
                    inFile.read(reinterpret_cast<char*>(&ilen), sizeof(ilen));
                    taxaPtr.addInterval(ibeg, ilen, ReadEnd::LEFT);
                }
                if (isPaired)
                    for (size_t i = 0; i < rcnt; ++i) {
                        inFile.read(reinterpret_cast<char*>(&ibeg), sizeof(ibeg));
                        inFile.read(reinterpret_cast<char*>(&ilen), sizeof(ilen));
                        taxaPtr.addInterval(ibeg, ilen, ReadEnd::RIGHT);
                    }
                taxaPtr.cleanIntervals(ReadEnd::LEFT);
                taxaPtr.cleanIntervals(ReadEnd::RIGHT);
                taxaPtr.updateScore();
            }
            return readInfo;
        }

        std::string refName(size_t id) {return refNames[id];}
        size_t refLength(size_t id) {return refLengths[id];}
        
    private:
        bool readHeader() {
            size_t refCount;
            inFile.read(reinterpret_cast<char*>(&isPaired), sizeof(bool));
            inFile.read(reinterpret_cast<char*>(&refCount), sizeof(size_t));
            logger->info("Total # of References: {}", refCount);
            refLengths.reserve(refCount);
            refNames.reserve(refCount);
            uint8_t refNameSize;
            uint32_t refLen;
            
            for (size_t i = 0; i < refCount; i++) {
                inFile.read(reinterpret_cast<char*>(&refNameSize), sizeof(refNameSize));
                char* strChar = new char[refNameSize];
                inFile.read(strChar, refNameSize);
                std::string refName(strChar, refNameSize);
                inFile.read(reinterpret_cast<char*>(&refLen), sizeof(uint32_t));
                refNames.push_back(refName);
                refLengths.push_back(refLen);
                //std::cout << refName << " " << refLen << "\n";
            }
            return true;
        }
        
        std::ifstream inFile;
        bool isPaired = true;
        std::vector<uint32_t> refLengths;
        std::vector<std::string> refNames;
        std::shared_ptr<spdlog::logger> logger; 
        ReadInfo readInfo;   
};

#endif