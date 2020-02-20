#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
//#include <getopt.h>
#include <cassert>
#include <cmath>

#include "ntHashIterator.hpp"
//#include "Uncompress.h"

#ifdef _OPENMP
# include <omp.h>
#endif

#define PROGRAM "nthll"

static const char VERSION_MESSAGE[] =
    PROGRAM " 1.1.0 \n"
    "Written by Hamid Mohamadi.\n"
    "Copyright 2018 Hamid Mohamadi, Licensed under MIT License\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... FILE(S)...\n"
    "Estimates distinct number of k-mers in FILE(S).\n\n"
    "Acceptable file formats: fastq, fasta, sam, bam and in compressed formats gz, bz, zip, xz.\n"
    "Accepts a list of files by adding @ at the beginning of the list name.\n"
    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N	use N parallel threads [1] (N>=2 should be used when input files are >=2)\n"
    "  -k, --kmer=N	the length of kmer [64]\n"
    "      --help	display this help and exit\n"
    "      --version	output version information and exit\n"
    "\n"
    "Report bugs to https://github.com/bcgsc/ntCard/issues\n";

using namespace std;

namespace opt {
unsigned nThrd=1;
unsigned nHash=7;
unsigned kmLen=64;
unsigned nBuck=65536;
unsigned nBits=16;
unsigned sBuck=4194304;
unsigned sBits=22;
size_t totKmer=0;
bool canon=true;
bool samH=true;
}

static const char shortopts[] = "t:k:b:s:hc";

enum { OPT_HELP = 1, OPT_VERSION };
/*
static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 't' },
    { "kmer",	required_argument, NULL, 'k' },
    { "bit",	required_argument, NULL, 'b' },
    { "sit",	required_argument, NULL, 's' },
    { "hash",	required_argument, NULL, 'h' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};
*/
unsigned getftype(std::ifstream &in, std::string &samSeq) {
    std::string hseq;
    getline(in,hseq);
    if(hseq[0]=='>') {
        return 1;
    }
    if(hseq[0]=='@') {
        if( (hseq[1]=='H'&& hseq[2]=='D') ||
                (hseq[1]=='S'&& hseq[2]=='Q') ||
                (hseq[1]=='R'&& hseq[2]=='G') ||
                (hseq[1]=='P'&& hseq[2]=='G') ||
                (hseq[1]=='C'&& hseq[2]=='O') ) {
            return 2;
        }
        else
            return 0;
    }
    opt::samH=false;
    samSeq=hseq;
    return 2;
}

inline void ntComp(const uint64_t hVal, uint8_t *mVec) {
    if(hVal&(~((uint64_t)opt::nBuck-1))) {
        uint8_t run0 = __builtin_clzll(hVal&(~((uint64_t)opt::nBuck-1)));
        if(run0 > mVec[hVal&(opt::nBuck-1)]) mVec[hVal&(opt::nBuck-1)]=run0;
    }
}

inline void ntRead(const string &seq, uint8_t *mVec) {
    ntHashIterator itr(seq, opt::kmLen);
    while (itr != itr.end()) {
        ntComp((*itr),mVec);
        ++itr;
    }
}

void getEfq(uint8_t *mVec, std::ifstream &in) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = static_cast<bool>(getline(in, seq));
        good = static_cast<bool>(getline(in, hseq));
        good = static_cast<bool>(getline(in, hseq));
        if(good && seq.length()>=opt::kmLen)
            ntRead(seq, mVec);
        good = static_cast<bool>(getline(in, hseq));
    }
}

void getEfa(uint8_t *mVec, std::ifstream &in) {
    bool good = true;
    for(string seq, hseq; good;) {
        string line;
        good = static_cast<bool>(getline(in, seq));
        while(good&&seq[0]!='>') {
            line+=seq;
            good = static_cast<bool>(getline(in, seq));
        }
        if(line.length()>=opt::kmLen)
            ntRead(line, mVec);
    }
}

void getEsm(uint8_t *mVec, std::ifstream &in, const std::string &samSeq) {
    std::string samLine,seq;
    std::string s1,s2,s3,s4,s5,s6,s7,s8,s9,s11;
    if(opt::samH) {
        while(getline(in,samLine))
            if (samLine[0]!='@') break;
    }
    else
        samLine=samSeq;
    do {
        std::istringstream iss(samLine);
        iss>>s1>>s2>>s3>>s4>>s5>>s6>>s7>>s8>>s9>>seq>>s11;
        if(seq.length()>=opt::kmLen)
            ntRead(seq, mVec);
    } while(getline(in,samLine));
}

uint64_t getNumDistinctKmers(unsigned kmlen, const std::string& ifile) {
  /*
  {
    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 't':
            arg >> opt::nThrd;
            break;
        case 'b':
            arg >> opt::nBits;
            break;
        case 's':
            arg >> opt::sBits;
            break;
        case 'k':
            arg >> opt::kmLen;
            break;
        case 'c':
            opt::canon=true;
            break;
        case OPT_HELP:
            std::cerr << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cerr << VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        }
        if (optarg != NULL && !arg.eof()) {
            std::cerr << PROGRAM ": invalid option: `-"
                      << (char)c << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }
    if (argc - optind < 1) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }
    if (die) {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }
    vector<string> inFiles;
    for (int i = optind; i < argc; ++i) {
        string file(argv[i]);
        if(file[0]=='@') {
            string inName;
            ifstream inList(file.substr(1,file.length()).c_str());
            while(getline(inList,inName))
                inFiles.push_back(inName);
        }
        else
            inFiles.push_back(file);
    }
  */

  opt::kmLen = kmlen;
  std::vector<std::string> inFiles;
  inFiles.push_back(ifile);

    opt::nBuck = ((unsigned)1) << opt::nBits;
    opt::sBuck = ((unsigned)1) << opt::sBits;

    uint8_t *tVec = new uint8_t [opt::nBuck];
    for (unsigned j=0; j<opt::nBuck; j++) tVec[j]=0;


    //#ifdef _OPENMP
    //    omp_set_num_threads(opt::nThrd);
    //#endif

    //#pragma omp parallel
    {
        string fname;
        uint8_t *mVec = new uint8_t [opt::nBuck];
        for (unsigned j=0; j<opt::nBuck; j++) mVec[j]=0;

        //#pragma omp for schedule(dynamic) nowait
        for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
            fname = inFiles[inFiles.size()-file_i-1]; // verbose
            std::ifstream in(inFiles[inFiles.size()-file_i-1].c_str());
            std::string samSeq;
            unsigned ftype = getftype(in,samSeq);
            if(ftype==0)
                getEfq(mVec, in);
            else if (ftype==1)
                getEfa(mVec, in);
            else if (ftype==2)
                getEsm(mVec, in, samSeq);
            in.close();
        }
        //#pragma omp critical(vmrg)
        {
            for (unsigned j=0; j<opt::nBuck; j++)
                if(tVec[j]<mVec[j])
                    tVec[j]=mVec[j];
        }
        delete [] mVec;
    }

    double pEst = 0.0, zEst = 0.0, eEst = 0.0, alpha = 0.0;
    alpha = 1.4426/(1 + 1.079/opt::nBuck);
    if(opt::canon) alpha/=2;

    for (unsigned j=0; j<opt::nBuck; j++)
        pEst += 1.0/((uint64_t)1<<tVec[j]);
    zEst = 1.0/pEst;
    eEst = alpha * opt::nBuck * opt::nBuck * zEst;

    delete [] tVec;
    return static_cast<uint64_t>(eEst);
    //std::cout << "F0, Exp# of distnt kmers(k=" << opt::kmLen << "): " << (unsigned long long) eEst << "\n";
    //return 0;
}

  /*
int main(int argc, char** argv) {
    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 't':
            arg >> opt::nThrd;
            break;
        case 'b':
            arg >> opt::nBits;
            break;
        case 's':
            arg >> opt::sBits;
            break;
        case 'k':
            arg >> opt::kmLen;
            break;
        case 'c':
            opt::canon=true;
            break;
        case OPT_HELP:
            std::cerr << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cerr << VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        }
        if (optarg != NULL && !arg.eof()) {
            std::cerr << PROGRAM ": invalid option: `-"
                      << (char)c << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }
    if (argc - optind < 1) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }
    if (die) {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }
    vector<string> inFiles;
    for (int i = optind; i < argc; ++i) {
        string file(argv[i]);
        if(file[0]=='@') {
            string inName;
            ifstream inList(file.substr(1,file.length()).c_str());
            while(getline(inList,inName))
                inFiles.push_back(inName);
        }
        else
            inFiles.push_back(file);
    }

    opt::nBuck = ((unsigned)1) << opt::nBits;
    opt::sBuck = ((unsigned)1) << opt::sBits;

    uint8_t *tVec = new uint8_t [opt::nBuck];
    for (unsigned j=0; j<opt::nBuck; j++) tVec[j]=0;


#ifdef _OPENMP
    omp_set_num_threads(opt::nThrd);
#endif

    #pragma omp parallel
    {
        string fname;
        uint8_t *mVec = new uint8_t [opt::nBuck];
        for (unsigned j=0; j<opt::nBuck; j++) mVec[j]=0;

        #pragma omp for schedule(dynamic) nowait
        for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
            fname = inFiles[inFiles.size()-file_i-1]; // verbose
            std::ifstream in(inFiles[inFiles.size()-file_i-1].c_str());
            std::string samSeq;
            unsigned ftype = getftype(in,samSeq);
            if(ftype==0)
                getEfq(mVec, in);
            else if (ftype==1)
                getEfa(mVec, in);
            else if (ftype==2)
                getEsm(mVec, in, samSeq);
            in.close();
        }
        #pragma omp critical(vmrg)
        {
            for (unsigned j=0; j<opt::nBuck; j++)
                if(tVec[j]<mVec[j])
                    tVec[j]=mVec[j];
        }
        delete [] mVec;
    }

    double pEst = 0.0, zEst = 0.0, eEst = 0.0, alpha = 0.0;
    alpha = 1.4426/(1 + 1.079/opt::nBuck);
    if(opt::canon) alpha/=2;

    for (unsigned j=0; j<opt::nBuck; j++)
        pEst += 1.0/((uint64_t)1<<tVec[j]);
    zEst = 1.0/pEst;
    eEst = alpha * opt::nBuck * opt::nBuck * zEst;

    delete [] tVec;

    std::cout << "F0, Exp# of distnt kmers(k=" << opt::kmLen << "): " << (unsigned long long) eEst << "\n";
    return 0;
}
  */
