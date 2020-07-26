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
#ifdef PCRE
#include "pcre2.h"
#else
class Pcre2 {};
#endif

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

extern map<string, size_t> chrs;
extern map<size_t, string> chrOf;

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

extern map<int64_t, Crispr_t> db;
extern map<int64_t, float> eff_db;

bool file_exist(string& file);
uint32_t udist(uint32_t a, uint32_t b);

class CrisprDB : public Pcre2 {
public:
	CrisprDB(string ref) {
		this->ref=ref;
		fa=ref+".fa";
		rep=NULL;
	}
	~CrisprDB() {
		if (rep)
			fclose(rep);
		rep=NULL;
	}
	static int64_t ngg2cid(char* seq);
	static void cid2ngg(int64_t cid, char* ngg);
	bool add(char* ngg, string& chr, bool reversed, int pos, string&);
	float doench(string&ref, int pos, bool reversed);
	void saveRef();
	void buildIdx();
	void scanIneffctive();
	
	map<string, size_t> header;
	void load();
	void saveHeader();
	static void revcomp(char* seq);
	static bool good(char* seq);
	virtual void found(int offset) { assert(false); }
private:
	string fa, ref;
	FILE* rep;
	inline void sam(FILE*fp, int64_t key, Crispr_t c, int, float doench);
	void scan(string& chr, string& seq);
	void dumpIdx(FILE* fp, size_t key, vector<uint32_t>& lst);
};

#endif
