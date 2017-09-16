#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <sys/time.h>
#include "../src/pmc.hpp"
#include "iomanip"
//#include <gperftools/profiler.h>

const double epsilon = 0.1;
using namespace std;

double getTimeMlsec(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec *  1e-6;
}

int main(int argc, char **argv) {
    //ProfilerStart("/Users/dangph/InfluenceMaximization/pruned-monte-carlo/pruned_perf.txt");
	double start_time = getTimeMlsec();
	std::cout << std::fixed << std::setprecision(6);
	if (argc < 5) {
		cerr << "./pmc graph k R p" << endl;
		exit(1);
	}

	string file = argv[1];
	int k = atoi(argv[2]);
	int R = atoi(argv[3]);
    double p = atof(argv[4]);

	ifstream is(file.c_str());
	vector<pair<pair<int, int>, double> > es;
	int u, v;
	for (; is >> u >> v;) {
		if (u == v) {
			continue;
		}
		es.push_back(make_pair(make_pair(u, v), p));
	}
	is.close();

	InfluenceMaximizer im;
    std::cout << "\t\tTime read=" << getTimeMlsec() - start_time << "\n";
    double start_run = getTimeMlsec();
	vector<int> seeds = im.run(es, k, R, epsilon);
	for (size_t i = 0; i < seeds.size(); i++) {
		cout << i << "-th seed =\t" << seeds[i] << endl;
	}

    std::cout << "\t\tTime run=" << getTimeMlsec() - start_run << "\n";

    //ProfilerStop();
	return 0;
}
