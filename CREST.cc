#include "CREST.h"

// Check if a file exists
bool file_exist(string& file)
{
	std::ifstream infile(file.c_str());
	return infile.good();
}

// chr1 => 0, chr2 => 1, etc
map<string, size_t> chrs;
// 0 => chr1, 1 => chr2, etc
map<size_t, string> chrOf;
// Database of sgRNA
map<int64_t, Crispr_t> db;
// Database of effective scores
map<int64_t, float> eff_db;

#define MAX_UINT32 0x7FFFFFFF
// array representation of linked list of seed + group of distal regions
vector<uint32_t> idx;

// 4^10 = 1M indices of seeds for rapid access of linked list
class Seed_map : public vector<uint32_t> {
public:
	Seed_map() {
		resize(1<<20, MAX_UINT32);
	}
	void add(uint32_t seed, uint32_t pos) {
		(*this)[seed&0xFFFFF]=pos;
	}
	void query(uint32_t seed, uint32_t* pool) {
		seed&=0xFFFFF;
		pool[0]=seed;	
		uint32_t pidx=1;
		for (uint32_t i=0, mask=3, mm1=0; i<10; i++, mask<<=2) {
			mm1=(seed&~mask) & 0xFFFFF; if (mm1!=seed) pool[pidx++]=mm1;
			mm1=((seed&~mask)|(0x55555&mask))&0xFFFFF; if (mm1!=seed) pool[pidx++]=mm1;
			mm1=((seed&~mask)|(0xAAAAA&mask))&0xFFFFF; if (mm1!=seed) pool[pidx++]=mm1;
			mm1=((seed&~mask)|(0xFFFFF&mask))&0xFFFFF; if (mm1!=seed) pool[pidx++]=mm1;
		}
		assert(pidx==31);
	}
} seed_map;

// read the idx in hex format
void loadIdx(string ref) {
	string file=ref+".idx";
	FILE* fp=fopen(file.c_str(), "r");
	idx.clear();
	uint32_t seed, sz, distal;
	while (fscanf(fp, "%X\t%u\n", &seed, &sz)==2) {
		seed_map.add(seed, idx.size());
		idx.push_back(seed);
		idx.push_back(sz);
		for (size_t i=0; i<sz; i++) {
			int ret=fscanf(fp, "%X\n", &distal);
			assert(ret==1);
			idx.push_back(distal);
		}
	}
	fclose(fp);
}

// check the number of mismatches
inline uint32_t udist(uint32_t a, uint32_t b) {
	uint32_t diff=a^b, ret=0;
	diff&=MASK;
	while (diff) {
		if (diff&3)
			ret++;
		diff>>=2;
	}
	return ret;
}

Option opt;

#ifdef SEED
#warning WITH_SEED
// Thread to find off-targets
void *runner(void*t) {
	long tid = (long)t;
	int nThreads=opt.threads;
	uint32_t* ptr=idx.data(), tail=(uint32_t)idx.size();
	uint32_t p=0, p_seed, p_sz, job=0, pool[32];
	while (p<tail) { p_seed=ptr[p++]; p_sz=ptr[p++]; ASSERT(p_seed<=MASK);
#ifdef DEBUG
		pthread_mutex_lock(&mtx);
		cout << "seed=" << p_seed << endl;
		pthread_mutex_unlock(&mtx);
#endif
		if (job++%nThreads==tid) { // NUM_THREADS
#ifdef DEBUG
			pthread_mutex_lock(&mtx);
			cout << "tid=" << tid << endl;
			pthread_mutex_unlock(&mtx);
#endif
			seed_map.query(p_seed, pool);
			// q_seed
			//uint32_t q=0, q_seed, q_sz;
			//while (q<tail) { q_seed=ptr[q++]; q_sz=ptr[q++]; ASSERT(p_seed<=MASK);
			for (int pi=0; pi<31; pi++) {
				uint32_t q_seed=pool[pi];
				uint32_t ds=udist(p_seed, q_seed); // skip if ds >= 2
#ifdef DEBUG
				pthread_mutex_lock(&mtx);
				printf("CHK_SEED tid=%d p_seed=%0X q_seed=%0X ds=%d\n", tid, p_seed, q_seed, ds); 
				pthread_mutex_unlock(&mtx);
#endif
				if (ds<2) {
#ifdef DEBUG
					pthread_mutex_lock(&mtx);
					cout << "SAME?" << endl;	
					pthread_mutex_unlock(&mtx);
#endif
					uint32_t q=seed_map[q_seed];
					if (q==MAX_UINT32) continue;
					// cerr << "q_seed=" << q_seed << "\tq=" << q << "\tptr[q]=" << ptr[q] << endl;
					assert(ptr[q++]==q_seed);
					uint32_t q_sz=ptr[q++];
					// printf("LOG: ds=%d ps: %0x qs: %0x\n", ds, p_seed, q_seed);
					for (uint32_t qt=0,qp=q; qt<q_sz; qt++,qp++) {
						uint32_t& q_distal=ptr[qp];
#ifdef DEBUG
				pthread_mutex_lock(&mtx);
				printf("CHK_tid=%d p_seed=%0X q=%0X:%0X\n", tid, p_seed, q_seed, q_distal); 
				pthread_mutex_unlock(&mtx);
#endif

						if (q_distal&SKIP) // OFFTGT or MANY
							continue;
						for (uint32_t pt=0,pp=p; pt<p_sz; pt++,pp++) {
							uint32_t& p_distal=ptr[pp];
							uint32_t dd=udist(p_distal, q_distal);
							// printf("LOG: ds=%d dd=%d p=%0X:%0X q=%0X:%0X\n", ds, dd, p_seed, p_distal, q_seed, q_distal);

							if (!(ds || dd))  // self
								continue;
#ifdef DEBUG
							pthread_mutex_lock(&mtx);
							if ((ds==0&&dd<4)||(ds==1&&dd<2)) { printf("OFFTGT: tid=%d ds=%d dd=%d p=%0X:%0X q=%0X:%0X\n", tid, ds, dd, p_seed, p_distal, q_seed, q_distal); }
							pthread_mutex_unlock(&mtx);
#endif
							if (ds==0&&dd<4) { p_distal|=OFFTGT; q_distal|=OFFTGT; }
							else if (ds==1&&dd<2) { p_distal|=OFFTGT; q_distal|=OFFTGT; }
						}
					}
				}
			//	q+=q_sz;
			}
		}
		p+=p_sz;
	}
	pthread_exit(NULL);
}
#else
#warning WITHOUT_SEED

#ifdef CREST
#define DS_DIFF(ds) (ds<2)
#define SG_DIFF(ds, dd) (((ds<<1)+dd)<4)
#else
#define DS_DIFF(ds) (ds<mm)
#define SG_DIFF(ds, dd) ((ds+dd)<mm)
#endif

// Thread to check off-target without index
void *runner(void*t) {
	long tid = (long)t;
	int nThreads=opt.threads;
	int mm=opt.maxMismatch;
	uint32_t* ptr=idx.data(), tail=(uint32_t)idx.size();
	uint32_t p=0, p_seed, p_sz, job=0;
	while (p<tail) { p_seed=ptr[p++]; p_sz=ptr[p++]; ASSERT(p_seed<=MASK);
		if (job++%nThreads==tid) { // NUM_THREADS
			uint32_t q=0, q_seed, q_sz;
			while (q<tail) { q_seed=ptr[q++]; q_sz=ptr[q++]; ASSERT(p_seed<=MASK);
				uint32_t ds=udist(p_seed, q_seed); // skip if ds >= 2
				if (DS_DIFF(ds)) { /// if (ds<2) {
					// printf("LOG: ds=%d ps: %0x qs: %0x\n", ds, p_seed, q_seed);
					for (uint32_t qt=0,qp=q; qt<q_sz; qt++,qp++) {
						uint32_t& q_distal=ptr[qp];
						if (q_distal&SKIP) // OFFTGT or MANY
							continue;
						for (uint32_t pt=0,pp=p; pt<p_sz; pt++,pp++) {
							uint32_t& p_distal=ptr[pp];
							uint32_t dd=udist(p_distal, q_distal);
							// printf("LOG: ds=%d dd=%d p=%0X:%0X q=%0X:%0X\n", ds, dd, p_seed, p_distal, q_seed, q_distal);
							if (!(ds || dd))  // self
								continue;
#ifdef DEBUG
							if (SG_DIFF(ds, dd)) { printf("OFFTGT: ds=%d dd=%d p=%0X:%0X q=%0X:%0X\n", ds, dd, p_seed, p_distal, q_seed, q_distal); }
							// if ((ds==0&&dd<4)||(ds==1&&dd<2)) { printf("OFFTGT: ds=%d dd=%d p=%0X:%0X q=%0X:%0X\n", ds, dd, p_seed, p_distal, q_seed, q_distal); }
#endif
							if (SG_DIFF(ds, dd)) { p_distal|=OFFTGT; q_distal|=OFFTGT; }
							//if (ds==0&&dd<4) { p_distal|=OFFTGT; q_distal|=OFFTGT; }
							//else if (ds==1&&dd<2) { p_distal|=OFFTGT; q_distal|=OFFTGT; }
						}
					}
				}
				q+=q_sz;
			}
		}
		p+=p_sz;
	}
	pthread_exit(NULL);
}
#endif

// setup and run threads for off-targets scanning
void mm4(string ref) {
	int ret, i;
	void *status;
	
	int nThreads=opt.threads;
	pthread_t threads[nThreads];
	pthread_attr_t attr;
	
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	for( i=0; i < nThreads; i++ ){
		ret = pthread_create(&threads[i], NULL, runner, (void *)(long)i );
		assert(!ret);
	}
	
	pthread_attr_destroy(&attr);
	for( i=0; i < nThreads; i++ ){
		ret = pthread_join(threads[i], &status);
		assert(!ret); assert(!status);
	}

	// Save the output	
	string output=ref+".mm4";
	string ref_file=ref+".ref";
	ifstream ref_fs(ref_file.c_str());
	string ln;
	FILE* fp=fopen(output.c_str(), "w");
	uint32_t* ptr=idx.data(), tail=(uint32_t)idx.size();
	uint32_t q=0, q_seed, q_sz;
	while (q<tail) { q_seed=ptr[q++]; q_sz=ptr[q++];
		for (uint32_t qt=0,qp=q; qt<q_sz; qt++,qp++) {
			uint32_t q_distal=ptr[qp];
			if (q_distal&SKIP) // NO OFFTGT, NO MANY
				continue;
			ASSERT(q_distal<=MASK);
			int64_t key=(((int64_t)q_seed)<<CRISPR_LEN)|q_distal;
			int64_t rid;
			while (getline(ref_fs, ln)) {
				ASSERT(!ref_fs.fail());
				int ret=sscanf(ln.c_str(), "%llx", &rid); ASSERT(ret==1);
				assert(key>=rid);
				if (key==rid)
					break;
			}
			fprintf(fp, "%s\n", ln.c_str());
		}
		q+=q_sz;
	}
	fclose(fp);
	
	pthread_exit(NULL);
}

void dump(vector<uint32_t>& idx, string& gRNA) {
    uint32_t mm=(uint32_t)opt.maxMismatch;
    int64_t cid=CrisprDB::ngg2cid((char*)gRNA.c_str());
    uint32_t seed=(uint32_t)(cid>>20), distal=cid&MASK;
    
    uint32_t* ptr=idx.data(), tail=(uint32_t)idx.size();
    uint32_t p, p_seed, p_sz;
    
    for (p=0; p<tail; p+=p_sz) {
        p_seed=ptr[p++]; p_sz=ptr[p++];
        uint32_t ds=udist(p_seed, seed);
        if (ds>=mm)
            continue;
        for (uint32_t pt=0,pp=p; pt<p_sz; pt++,pp++) {
            uint32_t& p_distal=ptr[pp];
            uint32_t dd=udist(p_distal, distal);
            if ((ds+dd)<mm) {
                char buf[24];
                int64_t id=p_seed; id<<=20; id+=p_distal;
                CrisprDB::cid2ngg(id, buf);
                cerr << gRNA << "\t" << p_seed << ":" << p_distal;
                cerr << "\t" << ds << "+" << dd << "\t" << buf << "\n";
            }
        }
    }
}

int main(int argc, const char * argv[]) {
	opt.parse(argc, argv);
	
	if (opt.query) {
		mkPair();
		return 0;
	}
	
    	if (!opt.ref) {
		opt.usage();
        	return -1;
    	}

	string ref=opt.ref;
	CrisprDB crisprDB(ref);

	string ref_file=ref+".ref";
	if (!file_exist(ref_file)) {
		crisprDB.load();
		crisprDB.saveHeader();
		crisprDB.saveRef();
		crisprDB.buildIdx();
	}

	if (!idx.size())
		loadIdx(ref);
	
	string mm4_file=ref+".mm4";
	if (!file_exist(mm4_file)) {
		mm4(ref);
	}

	if (opt.ineffectiveSample) {
		crisprDB.scanIneffctive();
	}

    return 0;
}
