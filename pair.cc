#include "CREST.h"

Timer tm;

class SAM_t {
public:
	SAM_t() { reversed=0; crispr=0; ga=0; ef=0; }
	SAM_t(bool rev, int64_t crispr, float ef) {
		reversed=0; ga=0;
		if (rev) reversed=1;
		this->crispr=crispr;
		this->ef=ef;
		
		naturalG=0;
		if ((crispr&3)==3) naturalG=1;

		int64_t tmp=crispr>>20;
		for (int i=10; i<CRISPR_LEN; i++, tmp>>=2) {
			switch (tmp&3) {
				case 0:	break;
				case 1: case 2: ga++;
				case 3: ga+=2;
			}
		}
	}
	uint64_t reversed:1;
	uint64_t crispr:40;
	uint64_t ga:5;
	uint64_t naturalG:1;
	float ef;
};

class SAM {
public:
	void read(istream& input) {
		int64_t cid, limit=((int64_t)1)<<40; int cut;
		string rid, chr, seq;
		int strand, s;
		string cigar, Is, efs; char wild; int qual, z1, z2;
		while (!input.fail()) {
			input >> hex >> cid >> dec >> rid >> strand >> chr >> s
				>> qual >> cigar >> wild >> z1 >> z2 >> seq >> Is >> efs;
			if (input.fail())
				break;
#ifdef DEBUG
			if (!(Is.find_first_not_of("I")==string::npos))
				cerr << Is << endl;
#else
			assert(Is.find_first_not_of("I")==string::npos);
#endif
			
#ifdef DEBUG
			if (cid>=limit) {
				cerr << cid << "?" << limit << endl;
			}
#else
			assert(cid<limit);
#endif
			if (strand)
				cut=s+2;
			else
				cut=s+17;
			
			if (chrs.find(chr)==chrs.end()) {
				size_t nChr=chrs.size();
				chrs[chr]=nChr;
				chrOf[nChr]=chr;
			}
			float ef=0;
			size_t pos=efs.find_last_of(':');
			if (pos!=efs.npos)
				ef=atof(efs.c_str()+pos+1);
			SAM_t sam(strand!=0, cid, ef);
			data[chr][cut]=sam;
		}
	}
	void print() {
		char buf[24];
		map<string, map<int64_t, SAM_t> >::iterator itr=data.begin();
		while (itr!=data.end()) {
			const string& chr=itr->first;
			map<int64_t, SAM_t>& sams=itr->second;
			map<int64_t, SAM_t>::iterator sam_itr=sams.begin();
			cout << chr << endl;
			while (sam_itr!=sams.end()) {
				int64_t cut=sam_itr->first;
				SAM_t sid=sam_itr->second;
				CrisprDB::cid2ngg(sid.crispr, buf);
				cout << "\t" << cut << "\t" << buf << ":" << sid.ga << "|" << sid.naturalG << endl;
				//if (cut==62566582)
				//	cout << "Here" << endl;
				sam_itr++;
			}
			itr++;
		}
	}
	map<string, map<int64_t, SAM_t> > data;
};

/*struct samPicker {	// return true if a is better
	bool operator()( const SAM_t& a, const SAM_t& b ) const {
		if (a.ga!=b.ga)
			return a.ga>b.ga;
		if (a.naturalG!=b.naturalG)
			return a.naturalG;
		return a.reversed<b.reversed;
	}
};*/

struct samPicker
{
	inline bool operator() (const pair<SAM_t, int64_t>& pa, const pair<SAM_t, int64_t>& pb)
	{
		SAM_t a=pa.first, b=pb.first;
		return a.ef>b.ef;
/*
		if (a.ga!=b.ga)
			return a.ga>b.ga;
		if (a.naturalG!=b.naturalG)
			return a.naturalG;
		return a.reversed<b.reversed;
 */
	}
};


class SamPair : public SAM {
public:
	SamPair() {
		if (opt.ref) {
			ifstream ifs(opt.ref);
			read(ifs);
		} else {
			read(cin);
		}
#ifdef DEBUG
		print();
#endif

		string q=opt.query;
		size_t colon=q.find_first_of(':');
		size_t minus=q.find_first_of('-', colon);
		if (colon==q.npos||minus==q.npos) {
			cerr << "Query format error(chr:s-e): " << q << endl;
			return;
		}

		string chr=q.substr(0, colon);
		string ss=q.substr(colon+1, minus-colon-1);
		string es=q.substr(minus+1).c_str();
		int s, e, ret;
		ret=sscanf(ss.c_str(), "%d", &s);
		assert(ret==1);
		ret=sscanf(es.c_str(), "%d", &e);
		assert(ret==1);
		query(chr, s, e);
	}

	void query(string& chr, int s, int e) {
		// cerr << chr << ":" << s << "-" << e << endl;

		map<int64_t, SAM_t>& sams=data[chr];
		
		vector<pair<SAM_t, int64_t> > window;
		vector<pair<SAM_t, int64_t> > selected;
		
		int cur_win=0, last_win=0, half_win=(opt.win_size+1)/2, step2=opt.step_size*2-1;
		map<int64_t, SAM_t>::iterator sam_itr;
		for (sam_itr=sams.begin(); sam_itr!=sams.end(); sam_itr++) {
			int64_t cut=sam_itr->first;
			if (cut<s || cut>e) {
				continue;
			}
			SAM_t sid=sam_itr->second;
			cur_win=(int)cut/half_win;
			if (cur_win!=last_win) {
				if (window.size()) {
					sort(window.begin(), window.end(), samPicker());
					selected.push_back(*window.begin());
				}
				window.clear();
			}
			window.push_back(make_pair(sid, cut));
			last_win=cur_win;
		}

		for (size_t i=step2; i<selected.size(); i+=2) {
			pair<SAM_t, int64_t>& a=selected[i-step2], &b=selected[i];
			mkPair(chr, a, b);
		}
	}
	
	void mkPair(string& chr, pair<SAM_t, int64_t>& a, pair<SAM_t, int64_t>& b) {
		char buf[24];
		
		SAM_t ca=a.first, cb=b.first;
		float ea=ca.ef, eb=cb.ef;
		int64_t cut_a=a.second, cut_b=b.second;
		CrisprDB::cid2ngg(ca.crispr, buf);
		string seq_a=buf;
		char strand_a='+';
		if (ca.reversed) {
			strand_a='-';
			CrisprDB::revcomp(buf);
		}
		string ref_a=buf;
		
		CrisprDB::cid2ngg(cb.crispr, buf);
		string seq_b=buf;
		char strand_b='+';
		if (cb.reversed) {
			strand_b='-';
			CrisprDB::revcomp(buf);
		}
		string ref_b=buf;
		
#ifdef DEBUG
		// + rev_comp function
		int pos_a=(int)cut_a-17; if (ca.reversed) pos_a=(int)cut_a-2;
		int pos_b=(int)cut_b-17; if (cb.reversed) pos_b=(int)cut_b-2;
		int gapN=(int)(pos_b-20-pos_a);
		cout << "DBG: pos_a="<<pos_a<<" pos_b="<<pos_b <<" dist="<<gapN << endl;
		
		string I40(40, 'I');
		cout << chr << ":" << cut_a << strand_a << seq_a << "|" << cut_b << strand_b << seq_b << "\t0\t" << chr << "\t" << pos_a << "\t";
		cout << 30 << "\t" << "20M" << gapN << "N20M\t" << chr << "\t";
		cout << pos_b << "\t20\t" << ref_a << ref_b << "\t" << I40 << "\t";
		cout << "EA:f:" << ea << "\tEB:f:" << eb << endl;
#else
		int64_t dist=cut_b-cut_a;
		cout << chr << "\t" << seq_a << "\t" << cut_a << "\t" << strand_a;
		cout << "\t" << seq_b << "\t" << cut_b << "\t" << strand_b;
		cout << "\t" << dist << " bp\t" << ea << "\t" << eb << "\n";
#endif
	}
};

void mkPair() {
	SamPair pair;
}
