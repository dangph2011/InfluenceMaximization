#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <sys/time.h>
#include "../src/pmc.hpp"
#include "iomanip"
#include "../SIEA/option.h"
#include "../SIEA/getMem.hpp"
//#include <gperftools/profiler.h>

//const double epsilon = 0.1;
using namespace std;

double getTimeMlsec(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec *  1e-6;
}

int main(int argc, char **argv) {
    //ProfilerStart("/Users/dangph/InfluenceMaximization/pruned-monte-carlo/pruned_perf.txt");
    srand(time(NULL));
    OptionParser op(argc, argv);
    if (!op.validCheck()){
        printf("Parameters error, please check the readme.txt file for correct format!\n");
        return -1;
    }

	double start_time = getTimeMlsec();
	std::cout << std::fixed << std::setprecision(6);
	if (argc < 5) {
		cerr << "./pmc graph k R p" << endl;
		exit(1);
	}

    char * inFile = op.getPara("-i");
    if (inFile == NULL){
        inFile = (char*)"network.bin";
    }

    char * model = op.getPara("-m");
    if (model == NULL)
        model = (char *) "IC";

    // float scale = 1;
    // char * scaledown = op.getPara("-sd");
    // if (scaledown != NULL){
    //     scale = atof(scaledown);
    // }

    char * tmp = op.getPara("-epsilon");
    double epsilon = 0.1;
    if (tmp != NULL){
        epsilon = atof(tmp);
    }

    // double delta = 1.0/n;
    // tmp = op.getPara("-delta");
    // if (tmp != NULL){
    //     delta = atof(tmp);
    // }

    int t = 2;
    tmp = op.getPara("-t");
    if (tmp != NULL){
        t = atoi(tmp);
    }

	string file = "network.txt";
    tmp = op.getPara("-f");
    if (tmp != NULL){
        file = string(tmp);
    }

    int k = 10;
    tmp = op.getPara("-k");
    if (tmp != NULL){
        k = atoi(tmp);
    }

    int R = 200;
    tmp = op.getPara("-r");
    if (tmp != NULL){
        R = atoi(tmp);
    }

    double p = 0.01;
    tmp = op.getPara("-p");
    if (tmp != NULL){
        p = atof(tmp);
    }
    // string file = string(argv[1]);
	// int k = atoi(argv[2]);
	// int R = atoi(argv[3]);
    // double p = atof(argv[4]);

	ifstream is(file.c_str());
	vector<pair<pair<int, int>, double> > es;
	int u, v;
    //ignore first line (containing the number of vertex and edges)
    is >> u >> v;
    //cout << "U=" << u << " " << v << endl;
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
	vector<int> seeds = im.run(es, k, R, epsilon, inFile, model, t);
	for (size_t i = 0; i < seeds.size(); i++) {
		cout << i << "-th seed =\t" << seeds[i] << endl;
	}

    std::cout << "\t\tTime run=" << getTimeMlsec() - start_run << "\n";
    cout << "Memory: " << getMemValue()/1024.0 << endl;
    //ProfilerStop();
	return 0;
}
