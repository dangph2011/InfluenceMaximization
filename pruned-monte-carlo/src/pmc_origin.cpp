#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <sys/time.h>
#include "pmc.hpp"
#include <unistd.h>
#include "../SIEA/hypergraph.hpp"
#include "../SIEA/xorshift.hpp"
#include "../SIEA/sfmt/SFMT.h"

using namespace std;

UI UI_MAX = 4294967295U;
ULL ULL_MAX = 18446744073709551615ULL;
vector<vector<int> > buffer;
vector<vector<int> > bufferIndex;
vector<unsigned int> beginIndex;
vector<unsigned int> endIndex;
vector<sfmt_t> sfmtSeed;
unsigned int maxBufferSize = 10000;
double b0;

double getCurrentTimeMlsec(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec *  1e-6;
}

double calculateInfluence(double epsilon, double delta, unsigned int numNodes, unsigned int numSeeds, int t, unsigned int rn){
    cout << "Step 0" << endl;
    double threshold = 1+(2+2/3)*2*log(numNodes)*(rn-numSeeds-1)*b0;
    double threshold123 = 1+(1+epsilon)*(2+2*epsilon/3)*log(2/delta)/(epsilon*epsilon)*(rn-numSeeds-1)*b0;

    long long counter = 0;
    long long totalSampling = 0;
    double observed = 0;

    double mu0 = rn-numSeeds-1;

    cout << "Step 1 " << endl;
    double eprime = sqrt(epsilon);
    if (eprime > 1.0/2){
        eprime = 1.0/2;
    }
    threshold = 1+(1+eprime)*(2+2*eprime/3)*log(6/delta)*(mu0*b0)/epsilon;
    cout << "threshold " << threshold << " " << (2+2*epsilon/3)*log(2/delta)*(numNodes-numSeeds-1) << " " << (epsilon*epsilon) << " " << 2*(1+sqrt(epsilon)) << " " << (1+2*sqrt(epsilon)) << " " << (1+log(3/2)/log(2/delta)) << endl;

    while(observed < threshold && observed < threshold123){
        for (int i = 0; i < t; ++i){
            while (beginIndex[i] != endIndex[i] && observed < threshold){
                observed += buffer[i][beginIndex[i]] + numSeeds/b0;
                counter += 1;
                beginIndex[i] = (beginIndex[i]+1)%maxBufferSize;
            }
        }

    }

    if (observed >= threshold123){
                cout << "Estimated Influence: " << observed*b0/counter << endl;
                printf("MC-1 samples: %lld\n", counter);
                return observed*b0/counter;
        }
    double mu1 = observed/counter;

    cout << "Mu 1: " << mu1 << endl;

    cout << "Step 2" << endl;
    double Delta = 0;
    double threshold2 = 2*(1+sqrt(epsilon))/(1-sqrt(epsilon))*(1+log(3/2)/log(2/delta))*(2+2*epsilon/3)*log(2/delta)*(mu0*b0)/(mu1*epsilon);
    cout << 2*(1+sqrt(epsilon))*(1+2*sqrt(epsilon))*(1+log(3/2)/log(2/delta))*(2+2*epsilon/3) << " " << log(2/delta) << " " << (mu0)/(mu1*epsilon) << endl;
    cout << "Threshold 2: " << threshold2 << endl;
    long long counter2 = 0;
    int max = 0;
    int maxind = 0;
    int beg = 0,end = 0;

    while(counter2 < threshold2 && observed < threshold123){
                for (int i = 0; i < t; ++i){
                        while ((beginIndex[i] + 1 < endIndex[i] || (beginIndex[i] > endIndex[i] && beginIndex[i] + 1 < endIndex[i] + maxBufferSize)) && counter2 < threshold2){
                                Delta += (buffer[i][beginIndex[i]] - buffer[i][beginIndex[i]+1])^2/2;
                observed += buffer[i][beginIndex[i]] + numSeeds/b0;
                observed += buffer[i][(beginIndex[i]+1)% maxBufferSize] + numSeeds/b0;
                                counter2 += 1;
                                beginIndex[i] = (beginIndex[i]+2)%maxBufferSize;
                        }
                }
        }

    if (observed >= threshold123){
                cout << "Estimated Influence: " << observed*b0/(counter+counter2*2) << endl;
                printf("MC-1 samples: %lld\n", (counter+counter2*2));
                return observed*b0/(counter+counter2*2);
        }

    cout << "Max, Maxind: " << max << " " << maxind << " " << beg << " " << end << endl;
    totalSampling += counter2*2;
    double rho = Delta/counter2+1;
    if (rho < epsilon*mu1){
        rho = epsilon*mu1;
    }
    cout << "Step 3 " << rho << endl;

    threshold = (2+2*epsilon/3)*log(2/delta)*rho*(mu0*b0)/(epsilon*epsilon*mu1*mu1);
    cout << "Threshold 3: " << threshold << " " << counter << endl;

    while(counter < threshold && observed < threshold123){
                for (int i = 0; i < t; ++i){
                        while (beginIndex[i] != endIndex[i] && counter < threshold){
                                observed += buffer[i][beginIndex[i]] + numSeeds/b0;
                                counter += 1;
                                beginIndex[i] = (beginIndex[i]+1)%maxBufferSize;
                        }
                }
        }
    totalSampling += counter;
    cout << "Estimated Influence: " << observed*b0/totalSampling << endl;
    cout << "Total Samples: " << totalSampling << endl;
    return observed*b0/totalSampling;
}

// void runEstimate(Graph g, HyperGraph hg, vector<int> neighbors, vector<double> f1, vector<double> f2, vector<int> seeds, int t,
//                     double epsilon, double delta, int n, int numSeeds, unsigned int rn, vector<bool> link, double b0){

double runEstimate(char* inFile, char* model, int t, vector<int> seeds, double epsilon){
    //srand(time(NULL));
    // OptionParser op(argc, argv);
    // if (!op.validCheck()){
    //     printf("Parameters error, please check the readme.txt file for correct format!\n");
    //     return -1;
    // }

    // char * interFile = op.getPara("-s");
    // if (interFile == NULL){
    //     interFile = (char*) "network.seeds";
    // }
    cout << "TEST=" << " " << inFile << " " << model << " " << t << " " << epsilon << endl;
    //vector<int> seeds;

    // for (int i = 1; i <= 1000; i++) {
    //     seeds.push_back(i);
    // }

    // ifstream in(interFile);
    //     int intTmp;
    //     in >> intTmp;
    //     while (!in.eof()){
    //             seeds.push_back(intTmp);
    //             in >> intTmp;
    //     }
    //     in.close();
    // for (int se:seeds) {
    //     cout << "SEED=" << se << endl;
    // }
    unsigned int numSeeds = seeds.size();
    cout << "Seeds done" << endl;

    vector<int> neighbors;
    vector<double> f1,f2;

    // char * inFile = op.getPara("-i");
    // if (inFile == NULL){
    //     inFile = (char*)"network.bin";
    // }

    // char * model = op.getPara("-m");
    // if (model == NULL)
    //     model = (char *) "IC";

    float scale = 1;
    // char * scaledown = op.getPara("-sd");
    // if (scaledown != NULL){
    //     scale = atof(scaledown);
    // }

    Graph g;
    if (strcmp(model, "LT") == 0){
        g.readGraphLT(inFile,scale);
    } else if (strcmp(model, "IC") == 0){
        g.readGraphIC(inFile,seeds,neighbors,f1,f2);
    } else {
        printf("Incorrect model option!");
        //return -1;
    }
    b0 = f1[f1.size()-1];
    cout << "Compute b0 " << f1[f1.size()-1] << " " << endl;

    if (b0 <= 0){
                cout << "Estimated Influence: " << numSeeds << endl;
                printf("MC-1 samples: %d\n", 1);
                cout << "Time: " << 0 << endl;
                //cout << "Memory: " << getMemValue()/1024.0 << endl;
                //return 1;
        }

    int n = g.getSize();
    int m = g.getEdge();

    //char * tmp = op.getPara("-epsilon");
    //double epsilon = 0.005;
    // if (tmp != NULL){
    //     epsilon = atof(tmp);
    // }

    double delta = 1.0/n;

    // tmp = op.getPara("-delta");
    // if (tmp != NULL){
    //     delta = atof(tmp);
    // }

    // int t = 2;
    // tmp = op.getPara("-t");
    // if (tmp != NULL){
    //     t = atoi(tmp);
    // }

    vector<bool> link(n+1, false);
    for (unsigned int i = 0; i < seeds.size(); ++i){
        link[seeds[i]] = true;
    }

    HyperGraph hg(n,m);
    unsigned int rn = hg.reachableNodes(seeds,n,g);

    // char * outFile = op.getPara("-o");
    // if (outFile == NULL){
    //     outFile = (char *)"output.txt";
    // }

    //calculateInfluence(epsilon,delta,n,numSeeds,t,rn);
    //hg.pollingIC(g,neighbors,f1,f2,visit,visit_mark,num_marked,sfmtSeed[id],seeds,firstphasetime);
    //runEstimate(g, hg, neighbors, f1, f2, seeds, t, epsilon, delta, n, numSeeds, rn, link);
    sfmtSeed = vector<sfmt_t>(t+1);
        for (int i = 0; i <= t; ++i){
                sfmt_init_gen_rand(&sfmtSeed[i], rand());
        }

    buffer = vector<vector<int> >(t);
    bufferIndex = vector<vector<int> >(t);
    for (int i = 0; i < t; ++i){
        buffer[i] = vector<int>(maxBufferSize);
        bufferIndex[i] = vector<int>(maxBufferSize);
    }
    beginIndex = vector<unsigned int>(t,0);
    endIndex = vector<unsigned int>(t,0);
    vector<int> maxSeed(t);

    timespec start,stop,start1,stop1;

    clock_gettime(CLOCK_REALTIME, &start);
/*  addHyperedge(g,hg,1,totalSamples,1,seeds,prob);
 *  clock_gettime(CLOCK_REALTIME, &stop);
 */
    double firstphasetime = 0;
    double timeforpolling = 0;

    bool terminated = false;
    double estimate_influence = 0;
    omp_set_num_threads(t);
        #pragma omp parallel shared(terminated)
        {

            int id = omp_get_thread_num();
            //cout << "Thread: " << id << endl;
        if (id > 0){
            //sleep(5);
            vector<bool> visit = link;
                    vector<int> visit_mark(n,0);
                    unsigned int num_marked = 0;
            unsigned int count = 0;

            while (!terminated){
//              cout << "Generate samples" << endl;
                clock_gettime(CLOCK_REALTIME, &start1);
                    hg.pollingIC(g,neighbors,f1,f2,visit,visit_mark,num_marked,sfmtSeed[id],seeds,firstphasetime);
                clock_gettime(CLOCK_REALTIME, &stop1);
                timeforpolling += ((stop1.tv_sec - start1.tv_sec) + (stop1.tv_nsec - start1.tv_nsec)/exp(9*log(10)));
//              cout << "Generate samples" << endl;
                count++;
                while (!terminated){
                    if ((endIndex[id]+1)%maxBufferSize != beginIndex[id]){
                        #pragma omp critical
                        {
                            buffer[id][endIndex[id]] = num_marked;
                            //cout << "buffer[id]=" << id << " endIndex[id]=" << endIndex[id] << " buffer[id][endIndex[id]]=" << buffer[id][endIndex[id]] << " beginIndex[id]=" << beginIndex[id] << endl;
//                          cout << "One sample: " << num_marked << " " << numSeeds << " " << b0 << endl;
                            bufferIndex[id][endIndex[id]] = count;
                            endIndex[id] = (endIndex[id]+1)%maxBufferSize;
                            count = 0;
                        }
                        break;
                    }
                    #pragma omp flush(terminated)
                }
            }
            //cout << "Thread " << id << " terminated " << endl;
        } else {
            estimate_influence = calculateInfluence(epsilon,delta,n,numSeeds,t,rn);
            cout << "terminated" << endl;
            terminated = true;
            #pragma omp flush(terminated)
        }
    }

    clock_gettime(CLOCK_REALTIME, &stop);
    cout << "Time for first step: " << firstphasetime << endl;
    cout << "Time for polling: " << timeforpolling << endl;
    cout << "Time: " << ((stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec)/exp(9*log(10))) << endl;
    return estimate_influence;
}

inline int PrunedEstimater::unique_child(const int v) {
	int outdeg = 0, child = -1;
	for (int i = at_e[v]; i < at_e[v + 1]; i++) {
		const int u = es[i];
		if (!removed[u]) {
			outdeg++;
			child = u;
		}
	}
	if (outdeg == 0) {
		return -1;
	} else if (outdeg == 1) {
		return child;
	} else {
		return -2;
	}
}

void PrunedEstimater::init(const int _n, vector<pair<int, int> > &_es,
		vector<int> &_comp) {
	flag = true;
	//n, _n: number of scc
	//n1, _comp.size() ~ number of vertex in original graph
	n = _n;
	n1 = _comp.size();

	visited.resize(n, false);

	int m = _es.size();
	vector<int> outdeg(n), indeg(n);

	for (int i = 0; i < m; i++) {
		int a = _es[i].first, b = _es[i].second;
		outdeg[a]++;
		indeg[b]++;
	}
	es.resize(m, -1);
	rs.resize(m, -1);

	at_e.resize(n + 1, 0);
	at_r.resize(n + 1, 0);

	//Prepare occupation for each vertex base on in/out degree
	at_e[0] = at_r[0] = 0;
	for (int i = 1; i <= n; i++) {
		at_e[i] = at_e[i - 1] + outdeg[i - 1];
		at_r[i] = at_r[i - 1] + indeg[i - 1];
	}

	for (int i = 0; i < m; i++) {
		int a = _es[i].first, b = _es[i].second;
		es[at_e[a]++] = b;
		rs[at_r[b]++] = a;
	}

	//Re calculate occupation for each vertex base on in/out degree because the value changed above
	at_e[0] = at_r[0] = 0;
	for (int i = 1; i <= n; i++) {
		at_e[i] = at_e[i - 1] + outdeg[i - 1];
		at_r[i] = at_r[i - 1] + indeg[i - 1];
	}

	sigmas.resize(n);
	comp = _comp;
	//
	vector<pair<int, int> > ps;
	for (int i = 0; i < n1; i++) {
		ps.push_back(make_pair(comp[i], i));
	}
	//sort by comp, in other way move all vertex in one component adjacency
	sort(ps.begin(), ps.end());
	at_p.resize(n + 1);
	//at_p: The number of vertices in each component.
	for (int i = 0; i < n1; i++) {
		//pmoc - at_p ~~~ es - at_e ~~~ rs - at_r
		pmoc.push_back(ps[i].second);
		at_p[ps[i].first + 1]++;
	}

	//The number of vertices in each node are stored consecutive in at_p array
	for (int i = 1; i <= n; i++) {
		at_p[i] += at_p[i - 1];
	}

	memo.resize(n);
	removed.resize(n);

	//Calculate the number of vertices in each component, stored in weight variable
	//Maybe need n (the number of scc) size.
	weight.resize(n1, 0);
	for (int i = 0; i < n1; i++) {
		weight[comp[i]]++;
	}

	//Looking for hub
	//double start_time = getCurrentTimeMlsec();
	first();
	//cout << "First Hub=" << getCurrentTimeMlsec() - start_time << endl;
}

int PrunedEstimater::sigma1(const int v) {
	return sigma(comp[v]);
}

int PrunedEstimater::sigma(const int v0) {
	if (memo[v0]) {
		return sigmas[v0];
	}
	memo[v0] = true;
	if (removed[v0]) {
		return sigmas[v0] = 0;
	} else {
		int child = unique_child(v0);
		//left node, return the number of vertex in scc
		if (child == -1) {
			return sigmas[v0] = weight[v0];
		} else if (child >= 0) {
			//child = 1, travel to next level
			return sigmas[v0] = sigma(child) + weight[v0];
		} else {
			//child >= 2
			int delta = 0;
			vector<int> vec;
			visited[v0] = true;
			vec.push_back(v0);
			queue<int> Q;
			Q.push(v0);
			//always false after first() function
			bool prune = ancestor[v0];

			if (prune) {
				delta += sigma(hub);
			}

			for (; !Q.empty();) {
				const int v = Q.front();
				Q.pop();
				if (removed[v]) {
					continue;
				}
				if (prune && descendant[v]) {
					continue;
				}
				delta += weight[v];
				for (int i = at_e[v]; i < at_e[v + 1]; i++) {
					const int u = es[i];
					if (removed[u]) {
						continue;
					}
					if (!visited[u]) {
						visited[u] = true;
						vec.push_back(u);
						Q.push(u);
					}
				}
			}
			//read carefully in here
			//reset vector visited for marking travelled vertex (reduce init time of vector)
			for (int i = 0; i < vec.size(); i++) {
				visited[vec[i]] = false;
			}
			return sigmas[v0] = delta;
		}
	}
}

void PrunedEstimater::first() {
	hub = 0;
	for (int i = 0; i < n; i++) {
		if ((at_e[i + 1] - at_e[i]) + (at_r[i + 1] - at_r[i])
				> (at_e[hub + 1] - at_e[hub]) + (at_r[hub + 1] - at_r[hub])) {
			hub = i;
		}
	}

    int ancestor_num = 0;
    int descendant_num = 0;

	descendant.resize(n);
	queue<int> Q;

	//mark hub descendant
	Q.push(hub);
	for (; !Q.empty();) {
		// forall v, !remove[v]
		const int v = Q.front();
		Q.pop();
		descendant[v] = true;
        descendant_num++;
		for (int i = at_e[v]; i < at_e[v + 1]; i++) {
			const int u = es[i];
			if (!descendant[u]) {
				descendant[u] = true;
                descendant_num++;
				Q.push(u);
			}
		}
	}

	ancestor.resize(n);
	//mark hub ancestor
	Q.push(hub);
	for (; !Q.empty();) {
		const int v = Q.front();
		Q.pop();
		ancestor[v] = true;
        ancestor_num++;
		for (int i = at_r[v]; i < at_r[v + 1]; i++) {
			const int u = rs[i];
			if (!ancestor[u]) {
				ancestor[u] = true;
                ancestor_num++;
				Q.push(u);
			}
		}
	}
	ancestor[hub] = false;

    //cerr << "NUMBER OF ANCESTOR=" << ancestor_num << endl;
    //cerr << "NUMBER OF DESCENDANT=" << descendant_num << endl;
    //cerr << endl;

	//Calculate sigma based on hub technique
	for (int i = 0; i < n; i++) {
		sigma(i);
	}

	//no longer use hub
	//Unmark ancestor and descendant
	ancestor.assign(n, false);
	descendant.assign(n, false);

	//push all vertex which influent to sigma to up
	for (int i = 0; i < n1; i++) {
		up.push_back(i);
	}
}

//sum: gain of v
void PrunedEstimater::update(vector<long long> &sums) {
	//not run in the first time
	//before recomputing sigma, substracting the init value
	for (int i = 0; i < (int) up.size(); i++) {
		int v = up[i];
		//flag = true in init function
		//after that always false
		if (!flag) {
			sums[v] -= sigmas[comp[v]];
		}
	}

	//recompute sigma v
	for (int i = 0; i < (int) up.size(); i++) {
		int v = up[i];
		sums[v] += sigma1(v);
	}

	flag = false;
}

void PrunedEstimater::add(int v0) {
	v0 = comp[v0];
	queue<int> Q;
	Q.push(v0);
	removed[v0] = true;
	vector<int> rm;

	//BFS to remove tree rooted at v
	//remove from vertex v0
	//push removed vecter to queue Q
	for (; !Q.empty();) {
		const int v = Q.front();
		Q.pop();
		rm.push_back(v);
		for (int i = at_e[v]; i < at_e[v + 1]; i++) {
			const int u = es[i];
			if (!removed[u]) {
				Q.push(u);
				removed[u] = true;
			}
		}
	}

	//clear marked vertex
	up.clear();

	//make BFS run like it is
	vector<int> vec;
	//ancestor from v and descendant of v, stored in up
	for (int i = 0; i < (int) rm.size(); i++) {
		const int v = rm[i];
		memo[v] = false; // for update()// not use calculated value
		//push all change vertex in to up
		for (int j = at_p[v]; j < at_p[v + 1]; j++) {
			up.push_back(pmoc[j]);
		}
		for (int j = at_r[v]; j < at_r[v + 1]; j++) {
			const int u = rs[j];
			if (!removed[u] && !visited[u]) {
				visited[u] = true;
				vec.push_back(u);
				Q.push(u);
			}
		}
	}

	//ancestor from v, stored in up
	for (; !Q.empty();) {
		const int v = Q.front();
		Q.pop();
		memo[v] = false; // for update()// not use calculated value
		for (int j = at_p[v]; j < at_p[v + 1]; j++) {
			up.push_back(pmoc[j]);
		}
		for (int i = at_r[v]; i < at_r[v + 1]; i++) {
			const int u = rs[i];
			if (!visited[u]) {
				visited[u] = true;
				vec.push_back(u);
				Q.push(u);
			}
		}
	}

	for (int i = 0; i < vec.size(); i++) {
		visited[vec[i]] = false;
	}
}

vector<int> InfluenceMaximizer::run(vector<pair<pair<int, int>, double> > &es,
		const int k, const int R, double epsilon, char * inFile, char * model, int t) {

	double start_run = getCurrentTimeMlsec();
	n = 0;
	m = es.size();
	for (int i = 0; i < (int) es.size(); i++) {
		n = max(n, max(es[i].first.first, es[i].first.second) + 1);
	}

	sort(es.begin(), es.end());

	//List outgoing node
	es1.resize(m);
	//List incoming node
	rs1.resize(m);

	//outgoing node
	at_e.resize(n + 1);
	//incoming node
	at_r.resize(n + 1);

	vector<PrunedEstimater> infs(R);
	vector<int> seeds;
    double avr_rc1 = 0;
    double avr_rc2 = 0;
    long long seed_reachability = 0;
    double _evaluate_time = 0;
    double _find_seed_time_start = getCurrentTimeMlsec();
    int iteration_time = 1;
	for (int t = 0; t < R; t++) {
		Xorshift xs = Xorshift(t);

		int mp = 0;
		//outgoing node
		at_e.assign(n + 1, 0);
		//incoming node
		at_r.assign(n + 1, 0);
		vector<pair<int, int> > ps;
		for (int i = 0; i < m; i++) {
			if (xs.gen_double() < es[i].second) {
				es1[mp++] = es[i].first.second;
				//Count the number of out going node
				at_e[es[i].first.first + 1]++;
				//Reverse: first -> second, second -> first
				ps.push_back(make_pair(es[i].first.second, es[i].first.first));
			}
		}
		at_e[0] = 0;
		sort(ps.begin(), ps.end());
		for (int i = 0; i < mp; i++) {
			rs1[i] = ps[i].second;
			//Count the number of in coming node
			at_r[ps[i].first + 1]++;
		}
		for (int i = 1; i <= n; i++) {
			at_e[i] += at_e[i - 1];
			at_r[i] += at_r[i - 1];
		}

		vector<int> comp(n);

		//nscc: Number of strongly connected component, in other way the number of vertex in DAG
		//comp: map original vertex to each scc.
		int nscc = scc(comp);

		//Generate DAG by mapping original vertex to relatively scc.
		vector<pair<int, int> > es2;
		for (int u = 0; u < n; u++) {
			int a = comp[u];
			for (int i = at_e[u]; i < at_e[u + 1]; i++) {
				int b = comp[es1[i]];
				if (a != b) {
					es2.push_back(make_pair(a, b));
				}
			}
		}

		sort(es2.begin(), es2.end());
		es2.erase(unique(es2.begin(), es2.end()), es2.end());

		infs[t].init(nscc, es2, comp);
	}

	cout << "Time before find seed=" << getCurrentTimeMlsec() - start_run << std::endl;
	vector<long long> gain(n);
	vector<int> S;
	double start_seed = getCurrentTimeMlsec();
	//Begin to find seed
	for (int t = 0; t < k; t++) {
		for (int j = 0; j < R; j++) {
			infs[j].update(gain);
		}
		int next = 0;
		for (int i = 0; i < n; i++) {
			if (gain[i] > gain[next]) {
				next = i;
			}
		}
        seed_reachability += gain[next];
        //cerr << "GAIN " << t << " " << next << "=" << gain[next] << endl;

		S.push_back(next);
		for (int j = 0; j < R; j++) {
			infs[j].add(next);
		}
		seeds.push_back(next);
	}

    //avr_rc_std = ev.evaluate(seeds, es, epsilon, avr_rc1);
    avr_rc1 = (double)seed_reachability / R;
    //avr_rc2 = runEstimate(inFile, model, t, seeds, epsilon);
    cout << "\t\t\t\tAverage rc 1 = " << avr_rc1 << endl;
    //cout << "\t\t\t\tAverage rc 2 = " << avr_rc2 << endl;
    double _evaluate_time_start = getCurrentTimeMlsec();
    _evaluate_time += getCurrentTimeMlsec() - _evaluate_time_start;

    // if (avr_rc1/avr_rc2 - 1 < epsilon/2) {
    //     cout << "Result: TRUE";
    //     //cout << "\t\t\t\tAverage reachability standardize=" << avr_rc2/n << endl;
    //     //break;
    // } else {
    //     cout << "Result: FALSE";
    //     //cout << "\t\t\t\tAverage reachability standardize=" << avr_rc2/n << endl;
    // }

    cout << "\t\t\tFind Seed Time: " << getCurrentTimeMlsec() - _find_seed_time_start - _evaluate_time << endl;
    cout << "\t\t\tEvaluated Time: " << _evaluate_time << endl;
    cout << "\t\t\tNumber of Samples: " << R << endl;
    cout << "\t\t\tIteration Time: " << iteration_time << endl;
    cout << "\t\t\tAverage reachability standardize: " << avr_rc1/n << endl;

	cout << "Time find seed only=" << getCurrentTimeMlsec() - start_seed << endl;

    cout << getCurrentTimeMlsec() - _find_seed_time_start - _evaluate_time << endl;
    cout << _evaluate_time << endl;
    cout << R << endl;
    cout << iteration_time << endl;
    cout << avr_rc1/n << endl;

	return seeds;
}

int InfluenceMaximizer::scc(vector<int> &comp) {
	vector<bool> vis(n);
	stack<pair<int, int> > S;
	vector<int> lis;
	int k = 0;
	for (int i = 0; i < n; i++) {
		S.push(make_pair(i, 0));
	}
	for (; !S.empty();) {
		int v = S.top().first, state = S.top().second;
		S.pop();
		if (state == 0) {
			if (vis[v]) {
				continue;
			}
			vis[v] = true;
			S.push(make_pair(v, 1));
			for (int i = at_e[v]; i < at_e[v + 1]; i++) {
				int u = es1[i];
				S.push(make_pair(u, 0));
			}
		} else {
			lis.push_back(v);
		}
	}
	for (int i = 0; i < n; i++) {
		S.push(make_pair(lis[i], -1));
	}
	vis.assign(n, false);
	for (; !S.empty();) {
		int v = S.top().first, arg = S.top().second;
		S.pop();
		if (vis[v]) {
			continue;
		}
		vis[v] = true;
		comp[v] = arg == -1 ? k++ : arg;
		for (int i = at_r[v]; i < at_r[v + 1]; i++) {
			int u = rs1[i];
			S.push(make_pair(u, comp[v]));
		}
	}
	return k;
}
