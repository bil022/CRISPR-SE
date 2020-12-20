#include "Crispr_SE.h"

// opt to keep options parsed from arguments
Option opt;

#define DS_LESS(ds) ((ds<<sw)<mm)
#define SG_LESS(ds, dd) (((ds<<sw)+dd)<mm)

// calculate the number of 2-bits mismatches in an efficient way
inline uint32_t udist(uint32_t a, uint32_t b) {
    uint32_t i=a^b;
    i&=MASK;
    i = i | ((i >> 1) & 0x55555555);
    i = i & 0x55555555;
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    i = (i + (i >> 4)) & 0x0f0f0f0f;
    i = i + (i >> 8);
    i = i + (i >> 16);
    return i & 0x3f;
}

CrisprSE SE;
// mutex for synchronizing
pthread_mutex_t mtx;

#define BUF_SIZE 128
#define MAX_HOLD 10240

// class to loop through the blocks of index (.idx) with seeds and distals
class SE_Iter {
protected:
    long tid;
    uint32_t* ptr;
    string ss;
    int curr_hold;
    size_t curr_len;
    char buf[BUF_SIZE];
    // use a buffer for keep the record, flush to disk once it is full
    inline void append() {
        ss+= buf; ss+= "\n";
        curr_hold++;
        curr_len=ss.length();
        if (curr_hold>=MAX_HOLD) {
            flush();
        }
    }
    // flush to the disk, use mutex for synchronizing
    inline void flush() {
        pthread_mutex_lock(&mtx);
        if (ss.length()!=curr_len) {
            cerr << ss.length() << " != " << curr_len << "\n";
            assert(false);
        }
        cout << ss;
        ss.clear();
        curr_hold=0;
        curr_len=0;
        pthread_mutex_unlock(&mtx);
    }
public:
    SE_Iter(long tid=0):tid(tid){
        curr_hold=0;
        curr_len=0;
    }
    ~SE_Iter() {
        if (ss.length())
            flush();
    }
    size_t head, tail, offset, blk_size;
    uint32_t seed, distal;
    // initialize
    void init(size_t head, size_t tail) {
        ptr=SE.data();
        this->head=head;
        this->tail=tail;
        offset=head;
        blk_size=seed=distal=0;
    }
    inline void begin() {
        offset=head;
        blk_size=seed=distal=0;
    }
    // find the next block, and set value of blk_size
    inline bool hasNext() {
        if (offset<tail) {
            seed=ptr[offset++];
            blk_size=ptr[offset++];
            return true;
        }
        return false;
    }
    // loop to next block
    inline void next() { offset+=blk_size; }
    inline uint32_t distal_val(size_t off) { distal=ptr[offset+off]; return distal; }
    inline void mark(uint32_t off) { ptr[offset+off]|=OFFTGT; }
    // find the gRNA itself
    inline void dump() {
        snprintf(buf, BUF_SIZE, "[%lX]%05X:%05X", tid, seed, distal);
        append();
    }
    // find an off-target
    inline void dump(uint32_t ref_seed, uint32_t ref_distal, int seed_diff, int distal_diff) {
        //snprintf(buf, BUF_SIZE, "[%lX]%05X:%05X-%05X:%05X=%d:%d", tid, seed, distal, ref_seed, ref_distal, seed_diff, distal_diff);
        snprintf(buf, BUF_SIZE, "[%lX]%05X:%05X-%05X:%05X=%d:%d", tid, seed, distal, ref_seed, ref_distal&MASK, seed_diff, distal_diff);
        append();
    }
};

// multi-threading
void seRunner(long tid) {
    int nThreads=opt.threads;
    uint32_t mm=opt.maxMismatch, sw=opt.seed_weight;
    bool verbose=opt.verbose;
    
    SE_Iter ref_itr, query_itr(tid);
    if (!SE.query_ptr) { // reference only, query is the reference
        ref_itr.init(0, SE.size());
        query_itr.init(0, SE.size());
    } else { // reference with user input
        ref_itr.init(0, SE.query_ptr);
        query_itr.init(SE.query_ptr, SE.size());
    }
    
    uint32_t job=0; //, step=0xFFFFF/100, bar=0;
    // loop throught blocks of distals in reference
    for (query_itr.begin(); query_itr.hasNext(); query_itr.next()) {
        if (job++%nThreads!=tid) { // only work with one of N threads
            continue;
        }
        // loop through blocks of distals in user query
        for (ref_itr.begin(); ref_itr.hasNext(); ref_itr.next()) {
            uint32_t seed_diff=udist(query_itr.seed, ref_itr.seed);
            // cerr << seed_diff << "<" << mm << "?" << DS_LESS(seed_diff) << "\n";
            // there are too many mismatches in seed region, next
            if (!DS_LESS(seed_diff))
                continue;
            // loop throught each distals in reference
            for (uint32_t query_off=0; query_off<query_itr.blk_size; query_off++) {
                uint32_t query_distal=query_itr.distal_val(query_off);
                // if no need to check all (verbose) and (it is repeated or has off-targets), next
                if (!verbose&&(query_distal&SKIP))
                    continue;
                // loop through each distals in query
                for (uint32_t ref_off=0; ref_off<ref_itr.blk_size; ref_off++) {
                    uint32_t ref_distal=ref_itr.distal_val(ref_off);
                    uint32_t distal_diff=udist(query_distal, ref_distal);
                    // when same
                    if (!(seed_diff || distal_diff))  { // self
                        if (verbose) { query_itr.dump(); }
                        continue;
                    }
                    // when has less #mismatches
                    if (SG_LESS(seed_diff, distal_diff)) {
                        // mark the record
                        query_itr.mark(query_off);
                        ref_itr.mark(ref_off);
                        if (verbose) { // print out the off-targets
                            query_itr.dump(ref_itr.seed, ref_distal, seed_diff, distal_diff);
                        } else {
                            break;
                        }
                    }
                }
            }
        }
    }
}

// pthread with timer
void *runner(void*tid) {
    Timer tm;
    seRunner((long)tid);
    char buf[32];
    snprintf(buf, 32, "p%ld", (long)tid);
    tm.tick(buf);
    pthread_exit(NULL);
}

// setup and run threads for off-targets scanning
void build() {
    int ret, i;
    void *status;
    
    // the query is reference itself (whole genome) when query is empty
    string input, output;
    if (opt.query) {
        input=opt.getFile(opt.query, ".ref");
        output=opt.getFile(opt.query, ".mm", false);
    } else {
        input=opt.getFile(opt.ref, ".ref");
        output=opt.getFile(opt.ref, ".mm", false);
    }
    
    // routine for pThreads
    int nThreads=opt.threads;
    pthread_t threads[nThreads];
    pthread_attr_t attr;
    
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    // load reference index (.idx)
    SE.loadIdx(opt.ref);
    if (opt.query) { // load query index
        SE.loadIdx(opt.query);
    }
    
    // create threads
    for( i=0; i < nThreads; i++ ){
        ret = pthread_create(&threads[i], NULL, runner, (void *)(long)i );
        assert(!ret);
    }
    
    // wait for each threads to finish
    pthread_attr_destroy(&attr);
    for( i=0; i < nThreads; i++ ){
        ret = pthread_join(threads[i], &status);
        assert(!ret); assert(!status);
    }
    
    // exit if in verbose mode
    if (opt.verbose)
        return;
    
    // save the gRNA id with no off-targets
    ifstream ref_fs(input.c_str());
    string ln;
    FILE* fp=fopen(output.c_str(), "w");
    
    SE_Iter query_itr;
    if (!SE.query_ptr) { // check reference 
        query_itr.init(0, SE.size());
    } else { // check query 
        query_itr.init(SE.query_ptr, SE.size());
    }
    
    // loop through blocks
    for (query_itr.begin(); query_itr.hasNext(); query_itr.next()) {
        uint32_t curr_seed=query_itr.seed;
        size_t blk_size=query_itr.blk_size;
        // loop through each distals
        for (uint32_t off=0; off<blk_size; off++) {
            uint32_t curr_distal=query_itr.distal_val(off);
            // if repeats or has off-targets, skip
            if (curr_distal&SKIP)
                continue;
            ASSERT(curr_distal<=MASK);
            // create the key
            int64_t rid, key=(((int64_t)curr_seed)<<CRISPR_LEN)|curr_distal;
            // loop through the reference to find the line with same key
            while (getline(ref_fs, ln)) {
                ASSERT(!ref_fs.fail());
                int ret=sscanf(ln.c_str(), "%" PRIx64 "", &rid);
                assert(ret==1);
                assert(key>=rid);
                if (key==rid)
                    break;
            }
            // print out the line
            fprintf(fp, "%s\n", ln.c_str());
        }
    }
    fclose(fp);
}

// debugging only
void dump(vector<uint32_t>& idx, string& gRNA) {
    uint32_t mm=(uint32_t)opt.maxMismatch;
    int64_t cid=Util::ngg2cid((char*)gRNA.c_str());
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
                Util::cid2ngg(id, buf);
                cerr << gRNA << "\t" << p_seed << ":" << p_distal;
                cerr << "\t" << ds << "+" << dd << "\t" << buf << "\n";
            }
        }
    }
}

// main program entry
int main(int argc, const char * argv[]) {
    // tricks to speed up the IO
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    // process arguments
    int task=opt.parse(argc, argv);
    
    switch (task) {
        case TASK_INDEX: // create index for reference or query
        {
            CrisprDB crisprDB(opt.ref);
            if (opt.simple)
                crisprDB.loadSimple();
            else
                crisprDB.loadRef();
        }
            break;
        case TASK_BUILD: // perform off-target search
        {
            build();
        }
            break;
        default:
            return -1;
    }
    
    return 0;
}
