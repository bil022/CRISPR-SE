#include "CREST.h"

pthread_mutex_t mtx;

int64_t CrisprDB::ngg2cid(char* seq) {
	assert(strlen(seq)==CRISPR_LEN);
	int64_t cid=0;
	for (int i=CRISPR_LEN-1; i>=0; i--) {
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

void CrisprDB::cid2ngg(int64_t cid, char* ngg) {
	for (int i=0; i<CRISPR_LEN; i++, cid>>=2) {
		switch (cid&3) {
			case 0: ngg[i]='A'; break;
			case 1: ngg[i]='C'; break;
			case 2: ngg[i]='T'; break;
			default: ngg[i]='G';
		}
	}
	ngg[CRISPR_LEN]='\0';
	assert(!cid);
}

class Doench {
public:
	float eval(string& seq) {
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
} doencher;

int Doench::offs[70]={1, 2, 2, 3, 4, 4, 5, 5, 6, 6, 11, 14, 14, 15, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 23, 23, 23, 24, 24, 24, 27, 27, 28, 29, 1, 4, 5, 5, 6, 11, 11, 11, 11, 12, 13, 13, 16, 18, 18, 19, 19, 20, 20, 20, 20, 21, 22, 22, 23, 23, 24, 24, 24, 26, 28};

const char* Doench::seqs[70]={"G", "A", "C", "C", "C", "G", "A", "C", "C", "G", "A", "A", "C", "A", "C", "T", "A", "G", "C", "G", "A", "C", "G", "T", "G", "T", "C", "T", "T", "C", "G", "T", "A", "C", "T", "G", "T", "C", "G", "GT", "GC", "AA", "TA", "GG", "GG", "TA", "TC", "TT", "GG", "GA", "GC", "TG", "GG", "TC", "CC", "TG", "AC", "CG", "GA", "GG", "TC", "CG", "CT", "AA", "AG", "AG", "CG", "TG", "GT", "GG"};
int Doench::lens[70]={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
float Doench::weights[70]={-0.2753771, -0.3238875, 0.17212887, -0.1006662, -0.2018029, 0.24595663, 0.03644004, 0.09837684, -0.7411813, -0.3932644, -0.466099, 0.08537695, -0.013814, 0.27262051, -0.1190226, -0.2859442, 0.09745459, -0.1755462, -0.3457955, -0.6780964, 0.22508903, -0.5077941, -0.4173736, -0.054307, 0.37989937, -0.0907126, 0.05782332, -0.5305673, -0.8770074, -0.8762358, 0.27891626, -0.4031022, -0.0773007, 0.28793562, -0.2216372, -0.6890167, 0.11787758, -0.1604453, 0.38634258, -0.6257787, 0.30004332, -0.8348362, 0.76062777, -0.4908167, -1.5169074, 0.7092612, 0.49629861, -0.5868739, -0.3345637, 0.76384993, -0.5370252, -0.7981461, -0.6668087, 0.35318325, 0.74807209, -0.3672668, 0.56820913, 0.32907207, -0.8364568, -0.7822076, -1.029693, 0.85619782, -0.4632077, -0.5794924, 0.64907554, -0.0773007, 0.28793562, -0.2216372, 0.11787758, -0.69774};

// doench score
float CrisprDB::doench(string&ref, int pos, bool reversed) {
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
	
	return doencher.eval(seq);
}

bool CrisprDB::add(char* ngg, string& chr, bool reversed, int pos, string& ref) {
	int64_t cid=ngg2cid(ngg);
	if (chrs.find(chr)==chrs.end()) {
		size_t nChr=chrs.size();
		chrs[chr]=nChr;
		chrOf[nChr]=chr;
	}
	Crispr_t crispr(reversed, chrs[chr], pos);
	float efficiency=doench(ref, pos, reversed);
	map<int64_t, Crispr_t>::iterator other_itr=db.find(cid);
	if (other_itr!=db.end()) { // may skip the palindrome
		Crispr_t& other=other_itr->second;
		if (!other.many) {
			other.many=1; // assert(db[cid].many==1);
			sam(rep, cid, other, Q_TOOMANY, efficiency);
		}
		crispr.many=1;
		sam(rep, cid, crispr, Q_TOOMANY, efficiency);
		return false;
	}
	db[cid]=crispr;
	eff_db[cid]=efficiency;
	return true;
}

void CrisprDB::sam(FILE*fp, int64_t key, Crispr_t c, int qual, float efficiency) {
	char ngg[24], rev[24];
	cid2ngg(key, ngg);
	strncpy(rev, ngg, CRISPR_LEN); rev[CRISPR_LEN]='\0';
	if (c.reversed)
		revcomp(rev);
	fprintf(fp, "%llx:%s\t%d\t%s\t%d\t%d\t%dM\t*\t0\t0\t%s\t%s\tEF:f:%f\n",
			key, ngg, c.reversed?16:0, chrOf[c.chr].c_str(), c.pos+1, qual, CRISPR_LEN, rev, CRISPR_QUAL, efficiency);
}

void CrisprDB::saveRef() {
	string file=ref+".ref";
	if (file_exist(file))
		return;
	FILE* fp=fopen(file.c_str(), "w");
	ofstream ofs(file.c_str());
	map<int64_t, Crispr_t>::iterator db_itr;
	for (db_itr=db.begin(); db_itr!=db.end(); db_itr++) {
		int64_t k=db_itr->first; Crispr_t c=db_itr->second;
		if (c.many)
			sam(fp, k, c, Q_TOOMANY, eff_db[k]);
		else
			sam(fp, k, c, Q_UNIQUE, eff_db[k]);
	}
	fclose(fp);
}

void CrisprDB::buildIdx() {
	string file=ref+".idx";
	if (file_exist(file))
		return;
	FILE* fp=fopen(file.c_str(), "w");
	map<int64_t, Crispr_t>::iterator itr=db.begin();
	vector<uint32_t> lst;
	bool head=true;
	size_t last_seed = MANY;
	while (itr!=db.end()) {
		uint32_t seed=(uint32_t)(itr->first>>20), distal=itr->first&MASK;
		if (!head && seed!=last_seed) {
			dumpIdx(fp, last_seed, lst);
		}
		head=false;
		if (itr->second.many)
			distal|=MANY;
		lst.push_back(distal);
		last_seed=seed;
		itr++;
	}
	dumpIdx(fp, last_seed, lst);
	fclose(fp);
}

extern CrisprDB crisprDB;
extern vector<uint32_t> idx;

string inf_chr, inf_seq;

// TODO: apply seed indices
bool searchDB(size_t pos, string& sgRNA, char strand) {
	int64_t cid=CrisprDB::ngg2cid((char*)sgRNA.c_str());
	uint32_t seed=(uint32_t)(cid>>20), distal=cid&MASK;
	
	uint32_t* ptr=idx.data(), tail=(uint32_t)idx.size();
	uint32_t p, p_seed, p_sz;

	for (p=0; p<tail; p+=p_sz) {
		p_seed=ptr[p++]; p_sz=ptr[p++]; ASSERT(p_seed<=MASK);
		uint32_t ds=udist(p_seed, seed); // skip if ds >= 2
		if (ds>=2)
			continue;
		for (uint32_t pt=0,pp=p; pt<p_sz; pt++,pp++) {
			uint32_t& p_distal=ptr[pp];
			uint32_t dd=udist(p_distal, distal);
			if ((ds==0&&dd<4)||(ds==1&&dd<2)) {
#ifdef DEBUG
				char buf[24];
				int64_t id=p_seed; id<<=20; id+=p_distal;
				CrisprDB::cid2ngg(id, buf);
				cerr << pos << "@" << strand << sgRNA;
				cerr << "\tds/dd:" << ds << "/" << dd;
				cerr << "\t" << buf << endl;
#endif
				return true;
			}
		}
	}

	return false;
}

void* thIneffective(void* t) {
	long tid = (long)t;
	int nTh=opt.threads, sample=opt.ineffectiveSample;
	size_t len=inf_seq.length();
	int win_sz=(int)(len/sample); if (!win_sz) win_sz++;
	for (size_t off0=CRISPR_LEN, n=0; off0<=len; off0+=win_sz, n++) {
		if (n%nTh!=tid)
			continue;
		for (int pos=0; pos<win_sz; pos++) {
			size_t off=off0+pos, start=off-CRISPR_LEN;
			string ngg=inf_seq.substr(start, CRISPR_LEN);
			if (ngg.find_first_of("N")!=string::npos)
				continue;
			
			if (off<=(len-3)) {
				string nag=inf_seq.substr(off+1, 2);
				if ((nag.compare("GG")==0) || (nag.compare("AG")==0))
					continue;
			}

			if (!searchDB(start, ngg, '+')) {
				pthread_mutex_lock(&mtx);
				cout << start << "\t+\t" << ngg << endl;
				pthread_mutex_unlock(&mtx);
				break;
			}
			
			if (start>=3) {
				string ccn=inf_seq.substr(start-3, 2);
				if ((ccn.compare("CC")==0) || (ccn.compare("CT")==0))
					continue;
			}
			
			reverse(ngg.begin(), ngg.end());
			for (size_t i=0; i<CRISPR_LEN; i++) {
				switch (ngg[i]) {
					case 'A': ngg[i]='T'; break;
					case 'C': ngg[i]='G'; break;
					case 'T': ngg[i]='A'; break;
					case 'G': ngg[i]='C'; break;
					default: assert(false);
				}
			}

			if (!searchDB(start, ngg, '-')) {
				pthread_mutex_lock(&mtx);
				cout << start << "\t-\t" << ngg << endl;
				pthread_mutex_unlock(&mtx);
				break;
			}
		}
	}
	
	pthread_exit(NULL);
}

void runIneffective() {
	int i, nTh=opt.threads;
#ifdef DEBUG
	cerr << "Scanning " << inf_chr << " " << inf_seq.length() << " bp w. " << nTh << " threads" << endl;
	assert(inf_seq.size()<(1<<28));
#endif
	
	pthread_t threads[nTh];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	for (i=0; i<nTh; i++) {
		int ret = pthread_create(&threads[i], NULL, thIneffective, (void *)(long)i );
		assert(!ret);
	}
	pthread_attr_destroy(&attr);
	void* status;
	for( i=0; i < nTh; i++ ){
		int ret = pthread_join(threads[i], &status);
		assert(!ret); assert(!status);
	}
	pthread_exit(NULL);
}

void CrisprDB::scanIneffctive() {
	string ln;
	while (cin >> ln) {
		if (ln.size() && ln[0]=='>') {
			if (inf_seq.size()) {
				runIneffective();
			}
			inf_chr=ln.substr(1);
			inf_seq.clear();
			continue;
		}
		transform(ln.begin(), ln.end(), ln.begin(), ::toupper);
		inf_seq+=ln;
	}
	if (inf_seq.size())
		runIneffective();
}

void CrisprDB::dumpIdx(FILE* fp, size_t key, vector<uint32_t>& lst) {
	if (!lst.size())
		return;
	fprintf(fp, "%zX\t%lu\n", key, lst.size());
	vector<uint32_t>::iterator itr=lst.begin();
	while (itr!=lst.end()) {
		fprintf(fp, "%X\n", *itr);
		itr++;
	}
	lst.clear();
}

void CrisprDB::load() {
	string rep_file=ref+".rep";
	rep=fopen(rep_file.c_str(), "w");

	fstream ifs(fa.c_str());
	string ln, chr, seq;
	while (getline(ifs, ln)) {
		if (ln.length() && ln[0]=='>') {
			if (seq.length()) {
				scan(chr, seq);
				header[chr]=seq.size();
				seq.clear();
			}
			seq.clear(); chr=ln.substr(1);
			continue;
		}
		transform(ln.begin(), ln.end(), ln.begin(), ::toupper);
		seq.append(ln);
	}
	scan(chr, seq);
	header[chr]=seq.size();
	seq.clear();
}

void CrisprDB::saveHeader() {
	string file=ref+".h";
	if (!file_exist(file)) {
		ofstream ofs(file.c_str());
		map<string, size_t>::iterator itr=header.begin();
		while (itr!=header.end()) {	//@SQ	SN:chr1	LN:197195432
			ofs << "@SQ\tSN:" << itr->first << "	LN:" << itr->second << endl;
			itr++;
		}
	}
}

bool CrisprDB::good(char* seq) {
	while (*seq) {
		switch (*seq) {
			case 'A': break;
			case 'C': break;
			case 'T': break;
			case 'G': break;
			default: return false;
		}
		seq++;
	}
	return true;
}

void CrisprDB::revcomp(char* seq) {
	unsigned long n=strlen(seq), i, n2=n/2;
	// reverse complementary
	for (i=0; i<n2; i++) {
		char tmp=seq[i];
		seq[i]=seq[n-i-1];
		seq[n-i-1]=tmp;
	}
	for (i=0; i<n; i++) {
		char ch=seq[i];
		switch (ch) {
			case 'A': seq[i]='T'; break;
			case 'C': seq[i]='G'; break;
			case 'T': seq[i]='A'; break;
			case 'G': seq[i]='C'; break;
			default: assert(false);
		}
	}
}

void CrisprDB::scan(string& chr, string& seq) {
	if (!chr.size())
		return;
#ifdef DEBUG
	cerr << "Scanning " << chr << " " << seq.size() << " bp" << endl;
	assert(seq.size()<(1<<28));
#endif

	unsigned long end=seq.size()-1;
	int pos=0, lastN=-1, fwds=0,revs=0;
	bool reversed;
	char ptr[24];
	ptr[CRISPR_LEN]='\0';
	for (int i=0; i<end; i++) {
		char gc=seq[i];
		if (gc=='G' && seq[i+1]==gc) {
			if (i<21) continue;
			pos=i-21;
			strncpy(ptr, seq.c_str()+pos, CRISPR_LEN);
			if (!good(ptr))
				continue;
			reversed=false;
		} else if (gc=='C' && seq[i+1]==gc) {
			pos=i+3;
			strncpy(ptr, seq.c_str()+pos, CRISPR_LEN);
			if (!good(ptr))
				continue;
			revcomp(ptr);
			reversed=true;
		} else {
			if (gc=='N')
				lastN=i;
			continue;
		}
		
		if (strlen(ptr)!=CRISPR_LEN) continue;
		// cerr << "Adding " << ptr << endl;
		if (add(ptr, chr, reversed, pos, seq)) {
			if (reversed) revs++; else fwds++;
		}
	}

#ifdef DEBUG
	cerr << fwds << " + " << revs << " -\n";
#endif
}
