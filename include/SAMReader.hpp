#ifndef SAM_READER_HPP
#define SAM_READER_HPP

#include "taxa.h"
#include "spdlog/spdlog.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/BamRecord.h"
//using namespace SeqLib;

class SAMReader {

public:
    SAMReader() { rec.init(); }

    ~SAMReader() { br.Close(); }

    void load(std::string mfileName,
              std::shared_ptr<spdlog::logger> loggerIn) {
        logger = loggerIn;
        br.Open(mfileName);
        if (!readHeader()) {
            logger->error("Invalid header for mapping output file.");
            std::exit(1);
        }
    }

    bool nextRead(ReadInfo &rinf, bool needReadName = false) {
        bool hasBean = false;
        (void)hasBean;
        rinf.mappings.clear();
        TaxaNode dummy;
        TaxaNode *taxa = &dummy;
        // we should either fill up a valid mapping record
        // or return false (first line inside the while) if we go until the end of the file
        while (rinf.mappings.size() == 0) {
            // if last record stored in rec (but not used yet) was the eof record, return false
            if (!hasNext) return false;
            ReadEnd re = ReadEnd::LEFT; // initial read is always assumed to be left
            std::string readName = rec.Qname(); // assumption: in any case, readName is valid
            rinf.len = rec.Length();
            if (needReadName) {
                rinf.rid = readName;
            }
            if (rec.MappedFlag()) {
                /*if (!hasBean and bh.IDtoName(rec.ChrID()) == "NC_032111.1|kraken:taxid|67082") {
                    hasBean = true;
                    std::cout << ">" << rec.Qname() << "\n"
                              << rec.Sequence() << "\n";
                }*/
                rinf.mappings.push_back(rec.ChrID());
                rinf.cnt = 1;
                taxa = &rinf.mappings.back();
                taxa->setFw(!rec.ReverseFlag(), re);
                taxa->setPos(rec.Position(), re);
                //FIXME for now, add a fake interval, starting from pos "0" in read and end in pos "0+coverage"
                taxa->addInterval(0, rec.NumMatchBases(), re);
                /*  std::cout   << "read name:" << rec.Qname() << " | "
                     << "len:" << rec.Length() << " | "
                     << "is paired:" << rec.PairedFlag() << " | "
                     << "ref ID:" << rec.ChrID() << " | "
                     << "ref name:" << bh.IDtoName(rec.ChrID()) << " | "
                     << "coverage:" << rec.NumMatchBases() <<  " | "
                     << "pos:" << rec.Position() << " | "
                     << "ori:" << !rec.ReverseFlag() << " | "
                     << "is considered left:" << static_cast<bool>(re == ReadEnd::LEFT) << "\n"; */
            }
            hasNext = br.GetNextRecord(rec);
            // If it's paired end, we expect the next read should be right,
            // since we've already read the left pair
            // in case of left read not being mapped,
            // the read end will switch to left again in the loop
            if (isPaired) {
                re = ReadEnd::RIGHT;
                rinf.len += rec.Length();
            }

            while (hasNext && rec.Qname() == readName) {
                if (rec.MappedFlag()) {
                    if (re == ReadEnd::RIGHT) { // it's always false for single end reads
                        if (static_cast<uint64_t>(rec.ChrID()) != taxa->getId()) {
                            re = ReadEnd::LEFT;
                        }
                    }
                    // If you're on a new mapping and the previous one was mapped at least on left end
                    if (re == ReadEnd::LEFT) {
                        if (taxa->getIntervals(ReadEnd::LEFT).size() == 0) {
                            taxa->setId(rec.ChrID());
                        } else {
                            /*std::cout << "here: " << taxa->getIntervals(ReadEnd::LEFT).size()
                                      << " " << taxa->getIntervals(ReadEnd::RIGHT).size() << "\n";*/
                            // first calc stat for previous read (pair) mapping
                            taxa->cleanIntervals(ReadEnd::LEFT);
                            taxa->cleanIntervals(ReadEnd::RIGHT);
                            taxa->updateScore();

                            // then create a new mapping and go on
                            rinf.mappings.push_back(rec.ChrID());
                            rinf.cnt++;
                            taxa = &rinf.mappings.back();
                        }
                    }
                    taxa->setFw(!rec.ReverseFlag(), re);
                    taxa->setPos(rec.Position(), re);
                    taxa->addInterval(0, rec.NumMatchBases(), re);
                    /* std::cout   << "read name:" << rec.Qname() << " | "
                    << "len:" << rec.Length() << " | "
                    << "is paired:" << rec.PairedFlag() << " | "
                    << "ref ID:" << rec.ChrID() << " | "
                    << "ref name:" << bh.IDtoName(rec.ChrID()) << " | "
                    << "coverage:" << rec.NumMatchBases() <<  " | "
                    << "pos:" << rec.Position() << " | "
                    << "ori:" << !rec.ReverseFlag() << " | "
                    << "is considered left:" << static_cast<bool>(re == ReadEnd::LEFT) << "\n"; */
                    if (isPaired)
                        re == ReadEnd::LEFT ? re = ReadEnd::RIGHT : re = ReadEnd::LEFT;
                } else { // orphans always appear in the left
                    // if it was a not-mapped left-end and the right-end is gonna map,
                    // we put the orphan mapping in the left
                    re = ReadEnd::LEFT;
                }
                hasNext = br.GetNextRecord(rec);
            }

            // for the last read (pair) mapping record
            taxa->cleanIntervals(ReadEnd::LEFT);
            taxa->cleanIntervals(ReadEnd::RIGHT);
            taxa->updateScore();

        }
        return true;
    }

    const std::string refName(size_t id) {/* std::cout << id << "\n"; */return bh.IDtoName(id); }

    size_t refLength(size_t id) {/* std::cout << id << "\n"; */return bh.GetSequenceLength(id); }

    size_t numRefs() const { return bh.NumSequences(); }

    bool isMappingPaired() { return isPaired; }

private:
    bool readHeader() {
        bh = br.Header();
        logger->info("# of targets: {}", bh.NumSequences());
        // isPaired
        uint64_t cntr{0};
        do {
            hasNext = br.GetNextRecord(rec);
            std::cerr << cntr++ << " " << rec.MappedFlag() << " | "
                      << rec.Qname() << " | "
                      << rec.Length() << " | "
                      << rec.PairedFlag() << " | "
                      << rec.ChrID() << " | "
                      << rec.Position() << " | "
                      << rec.ReverseFlag() << "\n";
            if (!hasNext) {
                logger->error("SAM file has no mapped records");
                exit(1);
            }
            isPaired = rec.PairedFlag();
        } while (!rec.MappedFlag());
        //isPaired=true;
        return true;
    }

    SeqLib::BamReader br;//(SeqLib::SAM);
    SeqLib::BamHeader bh;
    SeqLib::BamRecord rec;
    bool hasNext = true;
    bool isPaired = true;
    std::shared_ptr<spdlog::logger> logger;


};

#endif
