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
//#include <thread>
//#include <mutex>

using namespace std;
extern pthread_mutex_t mtx;

#include <time.h>
class Timer
{
public:
	Timer() { reset(); }
	double tick(const char* msg) {
		double ret=0;
#ifndef __MACH__
		clock_gettime(CLOCK_REALTIME, &end_);
		ret=end_.tv_sec - beg_.tv_sec +
		(end_.tv_nsec - beg_.tv_nsec) / 1000000000.;
		fprintf(stderr, "%s: %lf\n", msg, ret);
		reset();
#endif
		return ret;
	}
	void reset() {
#ifndef __MACH__
		clock_gettime(CLOCK_REALTIME, &beg_);
#endif
	}
	
private:
	timespec beg_, end_;
};

class Option {
	int argc; const char** argv;
public:
	char *ref, *query;
	int threads, win_size, step_size, ineffectiveSample, maxMismatch, all;
	void parse(int argc, const char* argv[]) {
		this->argc=argc;
		this->argv=argv;
		ref = NULL;
		query = NULL;
		threads=2;
		win_size=200;
		step_size=10;
		ineffectiveSample=0;
		maxMismatch=4;
        all=0;

		int c;
		while ((c = getopt(argc, (char*const*)argv, "p:r:q:w:s:e:m:ah")) != -1) {
			switch (c) {
				case 'p': threads=atoi(optarg); assert(threads>0); break;
				case 'r': ref=optarg; break;
				case 'q': query=optarg; break;
				case 'w': win_size=atoi(optarg); assert(win_size>0); break;
				case 's': step_size=atoi(optarg); assert(step_size>0); break;
				case 'e': ineffectiveSample=atoi(optarg); break;
				case 'm': maxMismatch=atoi(optarg); assert(maxMismatch>0); break;
                case 'a': all=1; break;
				case 'h': usage(); exit(0);
				default: cerr << "Unknown option:" << c << endl; break;
			}
		}
	}
	bool exists(const char* file) {
		ifstream ifs(file);
		return ifs.good();
	}
	void usage() {
		cerr << "To build database:\n\t" << argv[0] << " -p 4 -r <ref>" << endl;
		cerr << "To query:\n\t" << argv[0] << " -p 4 -q <chr:s-e> -w <win_size> -s <step_size>" << endl;
		cerr << "\t-p <n>: number of threads" << endl;
		cerr << "\t-m <n>: max mismatch(exclusive)" << endl;
        cerr << "\t-a: dump all off-targets with query\n";
	}
};

extern Option opt;
void mkPair();

#endif
