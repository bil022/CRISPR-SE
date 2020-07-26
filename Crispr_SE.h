#ifndef CREST_h
#define CREST_h

#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <unistd.h>
#include <map>
#include <set>
#include <vector>
#include <stdint.h>
#include <cstring>
#include "util.h"

using namespace std;

#define CRISPR_LEN 20
#define CRISPR_QUAL "IIIIIIIIIIIIIIIIIIII"

#define MANY 0x100000
#define MASK 0x0FFFFF

#ifdef DEBUG
#define ASSERT(c) assert(c)
#else
#define ASSERT(c)
#endif

#define OFFTGT 0x200000
#define SKIP 0x300000

#define Q_TOOMANY	10
#define Q_UNIQUE	30

//extern map<string, size_t> chrs;
//extern map<size_t, string> chrOf;
extern pthread_mutex_t mtx;

class Crispr_t {
public:
    Crispr_t() {
        reversed=0; many=0; chr=0; pos=0;
    }
    Crispr_t(bool reversed, size_t chr, uint32_t pos) {
        many=0;
        if (reversed) this->reversed=1; else this->reversed=0;
        assert(chr<32); assert(pos<256000000);
        this->chr=chr;
        this->pos=pos;
    }
    uint64_t reversed:1;
    uint64_t many:1;
    uint64_t chr:5; // max 32 chrosomes
    uint64_t pos:28; // max 256M bp
};

//extern map<int64_t, Crispr_t> db;
//extern map<int64_t, float> eff_db;

bool file_exist(string& file);
uint32_t udist(uint32_t a, uint32_t b);

class CrisprDB {
protected:
    map<string, size_t> chrs;
    map<size_t, string> chrOf;
    map<int64_t, Crispr_t> db;
    map<int64_t, float> eff_db;
public:
    CrisprDB(string ref) {
        this->ref=ref;
        rep=NULL;
    }
    ~CrisprDB() {
        if (rep)
            fclose(rep);
        rep=NULL;
    }
    //static int64_t ngg2cid(char* seq);
    //static void cid2ngg(int64_t cid, char* ngg);
    void loadRef();
    void loadSimple();
    void scanIneffctive();
protected:
    void saveIdx();
    bool addGuideRNA(char* ngg, string& chr, bool reversed, int pos, string&);
    // float doench(string&ref, int pos, bool reversed);
    
    map<string, size_t> header;
    static void revcomp(char* seq);
    static bool good(char* seq);
    virtual void found(int offset) { assert(false); }
    
    string ref;
    FILE* rep;
    inline void sam(FILE*fp, int64_t key, Crispr_t c, int, float doench);
    void scanNGG(string& chr, string& seq);
    void outputIdx(FILE* fp, size_t key, vector<uint32_t>& lst);
};

#define MAX_UINT32 0x7FFFFFFF

class CrisprSE : public vector<uint32_t> {
public:
    CrisprSE() { query_ptr=0; }
    size_t query_ptr; // offset for the query
    // read the idx in hex format
    void loadIdx(const char* src) {
        Timer tm;
        string idx_file=opt.getFile(src, ".idx");
        FILE* fp=fopen(idx_file.c_str(), "r");
        // clear();
        if (size()) {
            assert(query_ptr==0);
            query_ptr=size();
        }
        uint32_t seed, sz, distal, total=0;
        //while (fscanf(fp, "%X\t%u\n", &seed, &sz)==2) {
        while (getNum(fp, seed) && getNum(fp, sz, false)) {
            push_back(seed);
            push_back(sz);
            for (size_t i=0; i<sz; i++) {
                //int ret=fscanf(fp, "%X\n", &distal);
                //assert(ret==1);
                bool ret=getNum(fp, distal);
                assert(ret);
                push_back(distal);
                total++;
            }
        }
        fclose(fp);
        cerr << "Load " << total << " gRNAs in " << src << "\n";
        tm.tick(src);
    }
    
    // https://www.geeksforgeeks.org/fast-io-for-competitive-programming/
    bool getNum(FILE* fp, uint32_t &number, bool hex=true)
    {
        register int c;
        number = 0;
        bool first=true;
        while ((c=getc(fp))!=EOF) {
            if (isspace(c)) {
                if (first) {
                    continue;
                }
                break;
            }
            if (hex) {
                register int base=48; //
                if (c>=97) base=87;
                else if (c>=65) base=55;
                number = (number<<4) + c - base;
            } else {
                number = number *10 + c - 48;
            }
            first=false;
        }
        // cerr << "getNum: " << number << endl;
        return !first;
    }
};

class Doench {
public:
    static float doench(string&ref, int pos, bool reversed) {
        int s=0;
        string seq="";
        if (reversed) {
            if (pos>=6) {
                s=pos-6;
            } else {
                string Ns(6-pos, 'N');
                seq=Ns;
                pos=0;
            }
        } else {
            if (pos>=4) {
                s=pos-4;
            } else {
                string Ns(4-pos, 'N');
                seq=Ns;
                pos=0;
            }
        }
        
        seq+=ref.substr(s, 30-(int)seq.length());
        
        if (reversed) {
            assert(seq[3]=='C');
            assert(seq[4]=='C');
            reverse(seq.begin(), seq.end());
            string::iterator itr=seq.begin();
            while (itr!=seq.end()) {
                switch (*itr) {
                    case 'A': *itr='T'; break;
                    case 'C': *itr='G'; break;
                    case 'T': *itr='A'; break;
                    case 'G': *itr='C'; break;
                }
                itr++;
            }
        } else {
            assert(seq[25]=='G');
            assert(seq[26]=='G');
        }
        
        return eval(seq);
    }
private:
    static float eval(string& seq) {
        float intercept =  0.59763615;
        float gcHigh    = -0.1665878;
        float gcLow     = -0.2026259;
        int i, gc=0;
        for (i=4; i<24; i++) {
            char ch=seq[i];
            if (ch=='G'||ch=='C')
                gc++;
        }
        float score=intercept;
        float weight=gcLow;
        if (gc>10)
            weight=gcHigh;
        score+=abs(10-gc)*weight;
        
        const char* pch=seq.c_str();
        for (i=0; i<70; i++) {
            int pos=offs[i], step=lens[i];
            const char* s=seqs[i]; float w=weights[i];
            if (strncmp(pch+pos, s, step)==0)
                score+=w;
        }
        
        return 1/(1+exp(-score));
    }
    
    static int offs[70];
    static const char* seqs[70];
    static int lens[70];
    static float weights[70];
};

#endif
