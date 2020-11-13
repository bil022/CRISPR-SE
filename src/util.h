#ifndef util_h
#define util_h

#include <stdio.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>
#include <cstring>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <getopt.h>
#include <inttypes.h>

using namespace std;
extern pthread_mutex_t mtx;

/* class to keep track of time */
#include <time.h>
class Timer
{
public:
    Timer() { reset(); }
    double tick(const char* msg) {
        double ret=0;
        clock_gettime(CLOCK_REALTIME, &end_);
        clock_gettime(CLOCK_REALTIME, &end_);
        ret=end_.tv_sec - beg_.tv_sec +
        (end_.tv_nsec - beg_.tv_nsec) / 1000000000.;
        fprintf(stderr, "%s: %lf\n", msg, ret);
        reset();
        return ret;
    }
    void reset() {
        clock_gettime(CLOCK_REALTIME, &beg_);
    }
    
private:
    timespec beg_, end_;
};

/* 
  class Option to parse arguments:
    There are to main functions: index & build
    index: Create indices for reference genome or user input
        option -s: simple format for input as one gRNA per line
    build: Build genome-wide gRNA candiate, can also be used to perform offtarget search with user input
        option -v: output off-targets details
*/

#define	TASK_HELP	0
#define	TASK_INDEX	1
#define	TASK_BUILD	2

class Option {
    int argc; const char** argv;
public:
    const char *ref, *query;
    string PAM;
    int task, threads, seed_weight, maxMismatch, maxOfftarget;
    bool simple, verbose;
    int parse(int argc, const char* argv[]) {
        this->argc=argc;
        this->argv=argv;
        task = TASK_HELP;
        threads=2;
        ref = NULL;
        query = NULL;
        simple=false;
        verbose=false;
        maxMismatch=0; // CREST-Seq
        PAM="NGG";
        seed_weight=0;
        maxOfftarget=1;
        
        const struct option long_options[] =
        {
            /* These options set a flag. */
            {"index", no_argument, &task, TASK_INDEX},	// -pr
            {"build", no_argument, &task, TASK_BUILD},	// -pr
            {"thread", required_argument, 0, 'p'},
            {"ref",  required_argument, 0, 'r'},
            {"simple",  required_argument, 0, 's'},
            {"maxMismatch",  required_argument, 0, 'm'},
            {"maxOfftarget",  required_argument, 0, 'n'},
            {"PAM", optional_argument, 0, 'P'},
            {"verbose",  no_argument, 0, 'v'},
            {0, 0, 0, 0}
        };
        
        int c, option_index=0;
        while ((c = getopt_long(argc, (char* const*)argv, "p:r:sq:m:n:vh", // verbose? help?
                                long_options, &option_index)) != -1) {
            switch (c) {
                case 0: break;
                case 'p': threads=atoi(optarg); assert(threads>0); break;
                case 'r': ref=optarg; break;
                case 's': simple=true; break;
                case 'q': query=optarg; break;
                case 'm': maxMismatch=atoi(optarg); assert(maxMismatch>=0); break;
                case 'n': maxOfftarget=atoi(optarg); assert(maxMismatch>=0); break;
                case 'P': PAM=optarg;
                case 'v': verbose=true; break;
                case 'h': task=TASK_HELP; break;
                default: cerr << c << "?\n"; task = TASK_HELP; break;
            }
        }
        
        if (ref==NULL) {
            task = TASK_HELP;
            cerr << "reference genome is required\n";
        }
        
        if (PAM.compare("NGG")!=0 && PAM.compare("NAG")!=0) {
            task = TASK_HELP;
            cerr << "Only NGG|NAG is supported\n";
        }

        if (maxMismatch==0) {
            seed_weight=1;
            maxMismatch=4;
        }
        
        switch (task) {
            case TASK_INDEX: // fa_file, !ref, !idx, !rep
            {
                refFile(".fa");
                refFile(".ref", false);
                refFile(".idx", false);
                refFile(".rep", false);
                refFile(".h", false);
            }
                //
                break;
            case TASK_BUILD: // idx
            {
                refFile(".ref");
                refFile(".idx");
                refFile(".mm", false);
            }
                break;
            default:
                task = TASK_HELP;;
        }
        
        if (task == TASK_HELP)
            usage();
        
        return task;
    }
    
    /*
     Common:
     case 'p': threads=atoi(optarg); assert(threads>0); break;
     case 'r': ref=optarg; break;
     case 'h': task=TASK_HELP; break;
     
     Index:
     case 's': simple=true; break;
     
     Build:
     case 'q': query=optarg; break;
     case 'm': maxMismatch=atoi(optarg); assert(maxMismatch>=0); break;
     case 'n': maxOfftarget=atoi(optarg); assert(maxMismatch>=0); break;
     case 'v': verbose=true; break;
     */
    
    void usage() {
        cerr << "\nProgram: Crispr-SE (CRISPR Search Engine)\n"
        "Contact:\tBin Li <bil022@ucsd.edu>\n"
        "\n"
        "Usage:\tCrispr-SE <command> [options]\n"
        "\n"
        "Command:\n\t--index\tcreate index from reference/query sequence in the FASTA format\n"
        "\t--build\tbuild whole genome single guide RNA (gRNA)\n"
        "\n"
        "Options:\n\t-p INT\tnumber of threads [2]\n"
        "\t-r STR\treference genome id (mm9, mm10, hg19, hg38, etc)\n"
        "\t-s\tThe FASTA format is simple format of 20-nt gRNA\n"
        "\t-q\tThe indexed query\n"
        "\t-m INT\tMax mismatch, 0 for CREST-Seq, 1+ for #mismatches, default: 0\n"
        "\t-n INT\tMax off-target, 0 for all, default: 1\n"
        "\t-v\tverbose mode, default: false\n"
        "\n"
        "Note: if not in verbose mode, the max off-target will be set to 1\n";
    }
    
    string refFile(const char* ext, bool exists=true) {
        return getFile(ref, ext, exists);
    }
    string getFile(const char* src, const char* ext, bool exists=true) {
        string file=src; file+=ext;
        ifstream ifs(file.c_str());
        if (ifs.good()!=exists) {
            if (exists)
                cerr << file << " not found\n";
            else
                cerr << file << " exists\n";
            exit(-1);
        }
        if (ifs.good()) {
            ifs.close();
        }
        return file;
    }
};

extern Option opt;

void mkPair();

/* class Util: convert gRNA into 2-bits code or reverse */
class Util {
public:
    static int64_t ngg2cid(char* seq) {
        assert(strlen(seq)==20);
        int64_t cid=0;
        for (int i=19; i>=0; i--) {
            cid<<=2;
            switch (seq[i]) {
                case 'A': break;
                case 'C': cid|=1L; break;
                case 'T': cid|=2L; break;
                case 'G': cid|=3L; break;
                default: assert(false);
            }
        }
        return cid;
    }
    
    static void cid2ngg(int64_t cid, char* ngg) {
        for (int i=0; i<20; i++, cid>>=2) {
            switch (cid&3) {
                case 0: ngg[i]='A'; break;
                case 1: ngg[i]='C'; break;
                case 2: ngg[i]='T'; break;
                default: ngg[i]='G';
            }
        }
        ngg[20]='\0';
        assert(!cid);
    }
};

#endif
