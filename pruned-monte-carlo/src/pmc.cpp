#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <random>
#include "pmc.hpp"

using namespace std;
#define GAP 10

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
	//flag = true;
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

	at_e[0] = at_r[0] = 0;
	for (int i = 1; i <= n; i++) {
		at_e[i] = at_e[i - 1] + outdeg[i - 1];
		at_r[i] = at_r[i - 1] + indeg[i - 1];
	}

	sigmas.resize(n);
	comp = _comp;
	vector<pair<int, int> > ps;
	for (int i = 0; i < n1; i++) {
		ps.push_back(make_pair(comp[i], i));
	}
	sort(ps.begin(), ps.end());
	at_p.resize(n + 1);
	for (int i = 0; i < n1; i++) {
		pmoc.push_back(ps[i].second);
		at_p[ps[i].first + 1]++;
	}
	for (int i = 1; i <= n; i++) {
		at_p[i] += at_p[i - 1];
	}

	memo.resize(n);
	removed.resize(n);

	weight.resize(n1, 0);
	for (int i = 0; i < n1; i++) {
		weight[comp[i]]++;
	}

	//first();
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
		if (child == -1) {
			return sigmas[v0] = weight[v0];
		} else if (child >= 0) {
			return sigmas[v0] = sigma(child) + weight[v0];
		} else {
			int delta = 0;
			vector<int> vec;
			visited[v0] = true;
			vec.push_back(v0);
			queue<int> Q;
			Q.push(v0);
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
			for (int i = 0; i < vec.size(); i++) {
				visited[vec[i]] = false;
			}
			return sigmas[v0] = delta;
		}
	}
}

void PrunedEstimater::first() {
    flag = true;
	hub = 0;
	for (int i = 0; i < n; i++) {
		if ((at_e[i + 1] - at_e[i]) + (at_r[i + 1] - at_r[i])
				> (at_e[hub + 1] - at_e[hub]) + (at_r[hub + 1] - at_r[hub])) {
			hub = i;
		}
	}

	descendant.resize(n);
	queue<int> Q;
	Q.push(hub);
	for (; !Q.empty();) {
		// forall v, !remove[v]
		const int v = Q.front();
		Q.pop();
		descendant[v] = true;
		for (int i = at_e[v]; i < at_e[v + 1]; i++) {
			const int u = es[i];
			if (!descendant[u]) {
				descendant[u] = true;
				Q.push(u);
			}
		}
	}

	ancestor.resize(n);
	Q.push(hub);
	for (; !Q.empty();) {
		const int v = Q.front();
		Q.pop();
		ancestor[v] = true;
		for (int i = at_r[v]; i < at_r[v + 1]; i++) {
			const int u = rs[i];
			if (!ancestor[u]) {
				ancestor[u] = true;
				Q.push(u);
			}
		}
	}
	ancestor[hub] = false;

    //dangph
    memo.assign(n, false);
	removed.assign(n, false);

	for (int i = 0; i < n; i++) {
		sigma(i);
	}

	ancestor.assign(n, false);
	descendant.assign(n, false);

    up.clear();
	for (int i = 0; i < n1; i++) {
		up.push_back(i);
	}
}

void PrunedEstimater::update(vector<long long> &sums) {
	for (int i = 0; i < (int) up.size(); i++) {
		int v = up[i];
		if (!flag) {
			sums[v] -= sigmas[comp[v]];
		}
	}
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

	up.clear();

	vector<int> vec;
	for (int i = 0; i < (int) rm.size(); i++) {
		const int v = rm[i];
		memo[v] = false; // for update()
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
	// reachable to removed node
	for (; !Q.empty();) {
		const int v = Q.front();
		Q.pop();
		memo[v] = false;
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
		const int k, const int R, double epsilon){
	n = 0;
	m = es.size();
	for (int i = 0; i < (int) es.size(); i++) {
		n = max(n, max(es[i].first.first, es[i].first.second) + 1);
	}

	sort(es.begin(), es.end());

	es1.resize(m);
	rs1.resize(m);
	at_e.resize(n + 1);
	at_r.resize(n + 1);
    int infs_size = 0;
    //int step = 0;
    //infs_size = GAP*pow(2,step);
	vector<PrunedEstimater> infs;
    //vector<PrunedEstimater> infs(infs_size);
	vector<int> seeds;

    //http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1);

    double avr_rc1 = 0;
    //double avr_rc2 = 0;
    double avr_rc_std = 0;
    Evaluater ev;
    ev.init(es);
    long long seed_reachability = 0;
    while (!avr_rc_std) {
        seeds.clear();
        seed_reachability = 0;
        infs_size = infs_size + 10;
        infs.resize(infs_size);
        //generate double sample size
    	for (int t = infs_size-10; t < infs_size; t++) {
    		Xorshift xs = Xorshift(t);

    		int mp = 0;
    		at_e.assign(n + 1, 0);
    		at_r.assign(n + 1, 0);
    		vector<pair<int, int> > ps;
    		for (int i = 0; i < m; i++) {
    			if (xs.gen_double() < es[i].second) {
                //if (dis(gen) < es[i].second) {
    				es1[mp++] = es[i].first.second;
    				at_e[es[i].first.first + 1]++;
    				ps.push_back(make_pair(es[i].first.second, es[i].first.first));
    			}
    		}
    		at_e[0] = 0;
    		sort(ps.begin(), ps.end());
    		for (int i = 0; i < mp; i++) {
    			rs1[i] = ps[i].second;
    			at_r[ps[i].first + 1]++;
    		}
    		for (int i = 1; i <= n; i++) {
    			at_e[i] += at_e[i - 1];
    			at_r[i] += at_r[i - 1];
    		}

    		vector<int> comp(n);

    		int nscc = scc(comp);

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

        for (int t = 0; t < infs_size; t++) {
            infs[t].first();
        }

    	vector<long long> gain(n);
    	//vector<int> S;

    	for (int t = 0; t < k; t++) {
    		for (int j = 0; j < infs_size; j++) {
    			infs[j].update(gain);
    		}
    		int next = 0;
    		for (int i = 0; i < n; i++) {
    			if (gain[i] > gain[next]) {
    				next = i;
    			}
    		}
            seed_reachability += gain[next];
    		//S.push_back(next);
    		for (int j = 0; j < infs_size; j++) {
    			infs[j].add(next);
    		}
    		seeds.push_back(next);
    	}
        cout << "\t\tReachability=" << seed_reachability << endl;
        cout << "\t\tNumber of sample = " << infs_size << endl;
        avr_rc1 = seed_reachability / infs_size;
        //double start_evaluate = getCurrentTimeMlsec();
        avr_rc_std = ev.evaluate(seeds, es, epsilon, avr_rc1);
        cout << "\t\tAverage rc 1 = " << avr_rc1 << endl;
        cout << "\t\tAverage reachability standardize=" << avr_rc_std << endl;
        //std::cout << "\t\tTime evaluate=" << getCurrentTimeMlsec() - start_evaluate << "\n";
    }
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

double Evaluater::boundStop(int n, double epsilon){
    return (1+epsilon)*(1/(epsilon*epsilon))*n*log(n);
}

void Evaluater::init(vector<pair<pair<int, int>, double> > &es) {
    n = 0;
    m = es.size();
    for (int i = 0; i < (int) es.size(); i++) {
        n = max(n, max(es[i].first.first, es[i].first.second) + 1);
    }

    sort(es.begin(), es.end());

    //List outgoing node
    es1.resize(m);

    //outgoing node
    at_e.resize(n + 1);

    int mp = 0;
    //outgoing node
    at_e.assign(n + 1, 0);

    for (int i = 0; i < m; i++) {
        es1[mp++] = es[i].first.second;
        //Count the number of out going node
        at_e[es[i].first.first + 1]++;
    }
    at_e[0] = 0;

    for (int i = 1; i <= n; i++) {
        at_e[i] += at_e[i - 1];
    }

    removed.resize(n,false);
    visited.resize(n, false);

    // vector<bool> flip_coin;
    // flip_coin.resize(m,false);
}


double Evaluater::evaluate(vector<int> seeds, vector<pair<pair<int, int>, double> > &es, double epsilon, double avr_rc1){
    long long seed_reachability = 0;
    //double start_run = getCurrentTimeMlsec();
    //double ep2 = epsilon;
    //epsilon = 0.1;
    //http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1);

    long long bs = (long long)ceil(boundStop(n, epsilon));
    long long max_sample = bs*(epsilon+2) / (2*avr_rc1);
    //cout << "\tMax number of sample=" << max_sample << endl;
    int nu_sample = 0;

    //while (seed_reachability < bs && nu_sample < max_sample) {
    while (seed_reachability < bs) {
        removed.assign(n,false);
        for (int se : seeds) {
            if (removed[se]) {
                continue;
            }
            queue<int> Q;
            Q.push(se);
            vector<int> vec;
            vec.push_back(se);
            visited[se] = true;
            //plus itself
            seed_reachability++;
            for (; !Q.empty();) {
                const int v = Q.front();
                Q.pop();
                if (removed[v]) {
                    continue;
                }

                for (int i = at_e[v]; i < at_e[v + 1]; i++) {
                    const int u = es1[i];
                    if (removed[u]) {
                        continue;
                    }
                    //check if u was removed
                    if (!visited[u]) {
                        //flip coin to confirm edge exists
                        if (dis(gen) < es[i].second) {
                        //if (xs.gen_double() < es[i].second) {
                        //if (flip_coin[i]) {
                            seed_reachability++;
                            visited[u] = true;
                            vec.push_back(u);
                            Q.push(u);
                        }
                    }
                }
            }
            //remark visited and mark removed
            for (auto u1 : vec) {
                visited[u1] = false;
                removed[u1] = true;
            }
        }
        nu_sample++;
        //cout << "Timesa=" << times << " " << seed_reachability << " " << bs << endl;
    }
    cout << "\tNeed " << nu_sample << " to evaluate model\n";
    cout << "\tMax number of sample=" << max_sample << endl;
    cout << "\tReachability=" << seed_reachability << endl;
    cout << "\tAverage rc 2 = " << (double)seed_reachability / nu_sample << endl;
    cout << "\tRATIO=" << (avr_rc1 / ((double)seed_reachability / nu_sample)) - 1 - epsilon/2 << "\n";
    //if (nu_sample < max_sample && (avr_rc1 / ((double)seed_reachability / nu_sample)) - 1 < epsilon/2) {
    if ((avr_rc1 / ((double)seed_reachability / nu_sample)) - 1 < epsilon/2) {
        return (double)seed_reachability / nu_sample / n;
    } else {
        return 0;
    }
}
