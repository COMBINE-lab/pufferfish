
#include <cstdlib>
#include <string>

#include "spdlog/spdlog.h"

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/BamRecord.h"

int main(int argc, char* argv[]) {
    std::shared_ptr<spdlog::logger> logger = spdlog::stderr_color_mt("console");
    logger->info("input: {}", argv[1]);
    logger->info("output: {}", argv[2]);
    logger->info("filter thresh: {}", argv[3]);
    std::string mfileName = argv[1];
    std::string moutfileName = argv[2];
    float prob = std::stof(argv[3]);
    SeqLib::BamReader br;
    SeqLib::BamWriter bw;
    SeqLib::BamRecord rec;
    SeqLib::BamRecord rec2;
    br.Open(mfileName);   
    bw.Open(moutfileName); 
    SeqLib::BamHeader bh = br.Header();   
    bw.SetHeader(br.Header());  
    bw.WriteHeader();       
    logger->info("# of targets: {}", bh.NumSequences());            
    
    while (br.GetNextRecord(rec)) {        
        if ( (rec.MappedFlag() or rec.PairMappedFlag()) ) {
            
            std::string gis = bh.IDtoName(rec.ChrID());
            std::string readgis = rec.Qname();
            size_t s = gis.find_first_of('|')+1;
            size_t e = gis.find_first_of('|', s);
            gis = gis.substr(s, e-s-1);
            s = readgis.find_first_of('|')+1;
            e = readgis.find_first_of('|', s);
            readgis = readgis.substr(s, e-s-1);
    
            if ( gis != readgis and ((double) rand() / (RAND_MAX)) <= prob ) {
                if (rec.PairedFlag())
                    br.GetNextRecord(rec2);
                continue;
            }
        }
        bw.WriteRecord(rec);
        if (rec.PairedFlag()) {
            br.GetNextRecord(rec2);
            bw.WriteRecord(rec2);
        }
    }
    bw.Close();
}