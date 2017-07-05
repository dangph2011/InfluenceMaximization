#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <sys/time.h>
#include "pmc.hpp"

using namespace std;

double getCurrentTimeMlsec(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec *  1e-6;
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
	double start_first = getCurrentTimeMlsec();
	first();
	cerr << "First and Init time" << getCurrentTimeMlsec() - start_first << endl;
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
	/*hub = 0;
	for (int i = 0; i < n; i++) {
		if ((at_e[i + 1] - at_e[i]) + (at_r[i + 1] - at_r[i])
				> (at_e[hub + 1] - at_e[hub]) + (at_r[hub + 1] - at_r[hub])) {
			hub = i;
		}
	}

	descendant.resize(n);
	queue<int> Q;

	//mark hub descendant
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
	//mark hub ancestor
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
	ancestor[hub] = false;*/

	ancestor.assign(n, false);
	descendant.assign(n, false);
	//Calculate sigma based on hub technique
	for (int i = 0; i < n; i++) {
		sigma(i);
	}

	//no longer use hub
	//Unmark ancestor and descendant
	//ancestor.assign(n, false);
	//descendant.assign(n, false);

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
		const int k, const int R) {
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

	for (int t = 0; t < R; t++) {
		Xorshift xs = Xorshift(t);

		int mp = 0;
		//outgoing node
		at_e.assign(n + 1, 0);
		//incoming node
		at_r.assign(n + 1, 0);
		vector<pair<int, int> > ps;
		for (int i = 0; i < m; i++) {
            //if (true) {
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

	vector<long long> gain(n);
	vector<int> S;

	//Begin to find seed
	for (int t = 0; t < k; t++) {
		for (int j = 0; j < R; j++) {
			infs[j].update(gain);
		}
		int next = 0;
		for (int i = 0; i < n; i++) {
            //std::cerr << "GAIN=" << i << " " << gain[i] << "\n";
			if (gain[i] > gain[next]) {
				next = i;
			}
		}
        cerr << "\n";

		S.push_back(next);
		for (int j = 0; j < R; j++) {
			infs[j].add(next);
		}
		seeds.push_back(next);
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

	//DFS
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
	//DFS reverse to find scc
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
