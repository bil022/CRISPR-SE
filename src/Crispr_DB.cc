#include "Crispr_SE.h"

// load reference genome file in fasta format
void CrisprDB::loadRef() {
    string rep_file=opt.refFile(".rep", false);
    rep=fopen(rep_file.c_str(), "w");
    
    string fa_file=opt.refFile(".fa");
    fstream ifs(fa_file.c_str());
    
    // read fasta files
    string ln, chr, seq;
    while (getline(ifs, ln)) {
        if (ln.length() && ln[0]=='>') {
            if (seq.length()) {
                // Scan NGG/NAG & keep the sequence length
                scanPAM(chr, seq);
                header[chr]=seq.size();
                seq.clear();
            }
            seq.clear(); chr=ln.substr(1);
            continue;
        }
        // Conver the sequence to upper case
        transform(ln.begin(), ln.end(), ln.begin(), ::toupper);
        seq.append(ln);
    }
    // scan the last one if any
    if (seq.length()) {
        scanPAM(chr, seq);
        header[chr]=seq.size();
        seq.clear();
    }

    // save to index file
    saveIdx();
}

// user gRNA, assuming one sgRNA per line
void CrisprDB::loadSimple() {
    string rep_file=opt.refFile(".rep", false);
    rep=fopen(rep_file.c_str(), "w");
    
    string fa_file=opt.refFile(".fa");
    fstream ifs(fa_file.c_str());
    string ln, chr="NA";
    
    // scan file in simple format
    while (getline(ifs, ln)) {
        if (ln[0]=='>') {
            chr=ln.substr(1);
            continue;
        }
        size_t len=ln.length();
        // use header to trace the position of the gRNA
        header[chr]+=len;
        assert(len==20 || len==23);
        // to upper case
        transform(ln.begin(), ln.end(), ln.begin(), ::toupper);
        if (len==23) {
            assert(ln[21]=='G' && ln[22]=='G');
            ln.resize(20);
        }
        // get position
        uint32_t pos=(uint32_t)header[chr];
        // convert in 2-bits format
        int64_t cid=Util::ngg2cid((char*)ln.c_str());
        // check if it is new chromosome
        if (chrs.find(chr)==chrs.end()) {
            size_t nChr=chrs.size();
            chrs[chr]=nChr;
            chrOf[nChr]=chr;
        }
        // make a gRNA record
        Crispr_t crispr(false, chrs[chr], pos);
        // check if it is repeated
        map<int64_t, Crispr_t>::iterator other_itr=db.find(cid);
        if (other_itr!=db.end()) { // may skip the palindrome?
            // find the other repeat
            Crispr_t& other=other_itr->second;
            // if the other repeat has not been identified as repeat, mark it, save it
            if (!other.many) {
                other.many=1; // assert(db[cid].many==1);
                sam(rep, cid, other, Q_TOOMANY, -1);
            }
            // mark current one is repeat, save it
            crispr.many=1;
            sam(rep, cid, crispr, Q_TOOMANY, -1);
        }
        // keep the last record
        db[cid]=crispr;
        eff_db[cid]=-1;
    }
    
    // save vector to index file
    saveIdx();
}

// scan gRNA from long reference sequence
void CrisprDB::scanPAM(string& chr, string& seq) {
    assert(chr.size());
#ifdef DEBUG
    cerr << "Scanning " << chr << " " << seq.size() << " bp" << endl;
    assert(seq.size()<(1<<28));
#endif
    
    char PAM1='G', PAM2='C';
    if (opt.PAM.compare("NAG")==0) {
        PAM1='A'; PAM2='T';
    }

    size_t ptail=seq.size()-1;
    int pos=0, fwds=0,revs=0;
    bool reversed;
    char ptr[24];
    ptr[CRISPR_LEN]='\0';
    for (size_t i=0; i<ptail; i++) {
        char gc=seq[i], gc2=seq[i+1];
        if (((gc==gc2) || (gc==PAM1)) && gc2=='G') {
            if (i<21) continue;
            pos=i-21;
            strncpy(ptr, seq.c_str()+pos, CRISPR_LEN);
            if (!good(ptr))
                continue;
            reversed=false;
        } else if (gc=='C' && ((gc2==gc)||(gc2==PAM2))) {
            pos=i+3;
            strncpy(ptr, seq.c_str()+pos, CRISPR_LEN);
            if (!good(ptr))
                continue;
            revcomp(ptr);
            reversed=true;
        } else {
            continue;
        }
        
        if (strlen(ptr)!=CRISPR_LEN) continue;
        // cerr << "Adding " << ptr << endl;
        // add gRNA
        if (addGuideRNA(ptr, chr, reversed, pos, seq)) {
            if (reversed) revs++; else fwds++;
        }
    }
    
#ifdef DEBUG
    cerr << fwds << " + " << revs << " -\n";
#endif
}

// check if gRNA is ACTG's
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

// reverse complementary
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

// add gRNA
bool CrisprDB::addGuideRNA(char* ngg, string& chr, bool reversed, int pos, string& seq) {
    // convert the sequence into 2-bits
    int64_t cid=Util::ngg2cid(ngg);
    if (chrs.find(chr)==chrs.end()) {
        size_t nChr=chrs.size();
        chrs[chr]=nChr;
        chrOf[nChr]=chr;
    }
    Crispr_t crispr(reversed, chrs[chr], pos);
    // calculate efficiency score
    float efficiency=Doench::doench(seq, pos, reversed);
    // check if it is repeated
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
    // update record and efficiecy score
    db[cid]=crispr;
    eff_db[cid]=efficiency;
    return true;
}

// save the record in SAM format for samtools
void CrisprDB::sam(FILE*fp, int64_t key, Crispr_t c, int qual, float efficiency) {
    char ngg[24], rev[24];
    Util::cid2ngg(key, ngg);
    strncpy(rev, ngg, CRISPR_LEN); rev[CRISPR_LEN]='\0';
    if (c.reversed)
        revcomp(rev);
    fprintf(fp, "%" PRIx64 ":%s\t%d\t%s\t%d\t%d\t%dM\t*\t0\t0\t%s\t%s\tEF:f:%f\n",
            key, ngg, c.reversed?16:0, chrOf[c.chr].c_str(), (int)(c.pos+1), qual, CRISPR_LEN, rev, CRISPR_QUAL, efficiency);
}

// save the group of gRNA with same seed region
void CrisprDB::outputIdx(FILE* fp, size_t key, vector<uint32_t>& lst) {
    assert(lst.size());
    fprintf(fp, "%zX\t%lu\n", key, lst.size());
    vector<uint32_t>::iterator itr=lst.begin();
    while (itr!=lst.end()) {
        fprintf(fp, "%X\n", *itr);
        itr++;
    }
    lst.clear();
}

// save header, gRNA sequences in SAM format and the index
void CrisprDB::saveIdx() {
    // Save header file
    string h_file=opt.refFile(".h", false);
    ofstream ofs(h_file.c_str());
    map<string, size_t>::iterator header_itr=header.begin();
    while (header_itr!=header.end()) {	//@SQ	SN:chr1	LN:197195432
        ofs << "@SQ\tSN:" << header_itr->first << "\tLN:" << header_itr->second << endl;
        header_itr++;
    }
    
    // Save reference in SAM format
    string ref_file=opt.refFile(".ref", false);
    FILE* ref_fp=fopen(ref_file.c_str(), "w");
    map<int64_t, Crispr_t>::iterator db_itr;
    for (db_itr=db.begin(); db_itr!=db.end(); db_itr++) {
        int64_t k=db_itr->first; Crispr_t c=db_itr->second;
        if (c.many)
            sam(ref_fp, k, c, Q_TOOMANY, eff_db[k]);
        else
            sam(ref_fp, k, c, Q_UNIQUE, eff_db[k]);
    }
    fclose(ref_fp);
    
    // Scan db and groups of distals with same seed
    string idx_file=opt.refFile(".idx", false);
    FILE* idx_fp=fopen(idx_file.c_str(), "w");
    map<int64_t, Crispr_t>::iterator itr=db.begin();
    vector<uint32_t> lst;
    bool head=true;
    size_t last_seed = MANY;
    while (itr!=db.end()) {
        uint32_t seed=(uint32_t)(itr->first>>20), distal=itr->first&MASK;
        if (!head && seed!=last_seed) {
            // save group of distals with same seed
            outputIdx(idx_fp, last_seed, lst);
        }
        head=false;
        if (itr->second.many)
            distal|=MANY;
        lst.push_back(distal);
        last_seed=seed;
        itr++;
    }
    
    outputIdx(idx_fp, last_seed, lst);
    fclose(idx_fp);
}

// dataset used to calculate Doench score
int Doench::offs[70]={1, 2, 2, 3, 4, 4, 5, 5, 6, 6, 11, 14, 14, 15, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 23, 23, 23, 24, 24, 24, 27, 27, 28, 29, 1, 4, 5, 5, 6, 11, 11, 11, 11, 12, 13, 13, 16, 18, 18, 19, 19, 20, 20, 20, 20, 21, 22, 22, 23, 23, 24, 24, 24, 26, 28};
const char* Doench::seqs[70]={"G", "A", "C", "C", "C", "G", "A", "C", "C", "G", "A", "A", "C", "A", "C", "T", "A", "G", "C", "G", "A", "C", "G", "T", "G", "T", "C", "T", "T", "C", "G", "T", "A", "C", "T", "G", "T", "C", "G", "GT", "GC", "AA", "TA", "GG", "GG", "TA", "TC", "TT", "GG", "GA", "GC", "TG", "GG", "TC", "CC", "TG", "AC", "CG", "GA", "GG", "TC", "CG", "CT", "AA", "AG", "AG", "CG", "TG", "GT", "GG"};
int Doench::lens[70]={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
float Doench::weights[70]={-0.2753771, -0.3238875, 0.17212887, -0.1006662, -0.2018029, 0.24595663, 0.03644004, 0.09837684, -0.7411813, -0.3932644, -0.466099, 0.08537695, -0.013814, 0.27262051, -0.1190226, -0.2859442, 0.09745459, -0.1755462, -0.3457955, -0.6780964, 0.22508903, -0.5077941, -0.4173736, -0.054307, 0.37989937, -0.0907126, 0.05782332, -0.5305673, -0.8770074, -0.8762358, 0.27891626, -0.4031022, -0.0773007, 0.28793562, -0.2216372, -0.6890167, 0.11787758, -0.1604453, 0.38634258, -0.6257787, 0.30004332, -0.8348362, 0.76062777, -0.4908167, -1.5169074, 0.7092612, 0.49629861, -0.5868739, -0.3345637, 0.76384993, -0.5370252, -0.7981461, -0.6668087, 0.35318325, 0.74807209, -0.3672668, 0.56820913, 0.32907207, -0.8364568, -0.7822076, -1.029693, 0.85619782, -0.4632077, -0.5794924, 0.64907554, -0.0773007, 0.28793562, -0.2216372, 0.11787758, -0.69774};
