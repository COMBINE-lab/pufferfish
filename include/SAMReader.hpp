#ifndef SAM_READER_HPP
#define SAM_READER_HPP

#include "taxa.h"
#include "spdlog/spdlog.h"
//#include "htslib/sam.h"
//#include "htslib/hts.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/BamRecord.h"
//using namespace SeqLib;


//#define SAM_GLOBAL_ARGS_INIT 


class SAMReader {

        public:
        SAMReader() {rec.init();}
        void load(std::string mfileName,
                std::shared_ptr<spdlog::logger> loggerIn) {
            logger = loggerIn;
            br.Open(mfileName);
            if (!readHeader()) {
                logger->error("Invalid header for mapping output file.");
                std::exit(1);
            }
        }

        bool nextRead(ReadInfo& rinf, bool needReadName=false) {
            rinf.mappings.clear();
            if (!hasNext) return false;
            ReadEnd re = ReadEnd::LEFT;
            std::string readName = rec.Qname();
            if (needReadName) {
                rinf.rid = readName;
            }
            rinf.mappings.push_back(rec.ChrID());
            rinf.len = rec.Length();
            rinf.cnt = 1;
            TaxaNode& taxa = rinf.mappings.back();
            taxa.setFw(!rec.ReverseFlag(), re);
            taxa.setPos(rec.Position(), re);
            taxa.addInterval(0, rec.NumMatchBases(), re);

            hasNext = br.GetNextRecord(rec);
// by default if it's paired end we assume that when we enter this function, we've already read the left pair
            if (isPaired && static_cast<uint64_t>(rec.ChrID()) == taxa.getId()) {
                re = ReadEnd::RIGHT;
                rinf.len += rec.Length();
            }
            
            while (hasNext && rec.Qname() == readName) {
                if (rec.MappedFlag()) {
                    if (re == ReadEnd::RIGHT) {
                        if (static_cast<uint64_t>(rec.ChrID()) != taxa.getId()) {
                            re = ReadEnd::LEFT;
                        }
                    } else {
                        rinf.mappings.push_back(rec.ChrID());
                        rinf.cnt++;
                        taxa = rinf.mappings.back();
                    }
                    taxa.setFw(!rec.ReverseFlag(), re);
                    taxa.setPos(rec.Position(), re);
                    taxa.addInterval(0, rec.NumMatchBases(), re);
                }                
                hasNext = br.GetNextRecord(rec);
                if (isPaired)
                    re == ReadEnd::LEFT? re = ReadEnd::RIGHT: re = ReadEnd::LEFT;
            }
            
            return true;
        }

        const std::string& refName(size_t id) {return refNames[id];}
        size_t refLength(size_t id) {return refLengths[id];}
        size_t numRefs() const { return refNames.size(); }
        bool isMappingPaired() {return isPaired;}
    private:
        bool readHeader() {
            bh = br.Header();
            bam_hdr_t *header = bh.get_();
            if (header == NULL) {
                logger->error("Failed to read SAM file header");
                std::exit(1);
            }
            // refNames
            // refLength
            std::cout << "# of targets: " << header->n_targets << "\n";
            refLengths.resize(header->n_targets);
            refNames.resize(header->n_targets);
            for (auto i = 0; i < header->n_targets; i++) {
                refNames[i] = std::string(header->target_name[i]);
                refLengths[i] = header->target_len[i];
                std::cout << refNames[i] << " " << refLengths[i] << "\n";
            }
            // isPaired
            do {
                hasNext = br.GetNextRecord(rec);
                if (!hasNext) {
                    logger->error("SAM file has no mapped records");
                    exit(1);
                }
                isPaired = rec.PairedFlag();
            } while (!rec.MappedFlag());
            std::cout   << "read name: " << rec.Qname() << " "
                        << "len: " << rec.Length() << " "
                        << "is paired: " << rec.PairedFlag() << " "
                        << "ref ID: " << rec.ChrID() << " "
                        << "ref name: " << bh.IDtoName(rec.ChrID()) << " "
                        << "pos: " << rec.Position() << " "
                        << "coverage: " << rec.NumMatchBases() << "\n";
            return true;
        }
        
        //samFile *fp;
        SeqLib::BamReader br;//(SeqLib::SAM);
        SeqLib::BamHeader bh;
        SeqLib::BamRecord rec;
    
        bool isPaired = true;
        std::vector<refLenType> refLengths;
        std::vector<std::string> refNames;
        std::shared_ptr<spdlog::logger> logger;
        //std::vector<SeqLib::BamRecord> recs;
        bool hasNext = true;

};

#endif