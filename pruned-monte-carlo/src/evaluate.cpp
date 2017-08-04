#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <random>
#include "evaluate.hpp"
#include "pmc.hpp"

double Evaluater::boundStop(int n, double epsilon){
    return (1+epsilon)*(1/(epsilon*epsilon))*n*log(n);
}

double Evaluater::init(vector<pair<pair<int, int>, double> > &es) {
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


double Evaluater::evaluate(vector<int> seeds, vector<pair<pair<int, int>, double> > &es, double epsilon){
    long long seed_reachability = 0;
    //double start_run = getCurrentTimeMlsec();

    //http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1);

    long long bs = (long long)ceil(boundStop(n, epsilon));

    int times = 0;

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
        times++;
        //cout << "Timesa=" << times << " " << seed_reachability << " " << bs << endl;
    }
    cout << "Need " << times << " to evaluate model\n";
    cout << "Reachability=" << seed_reachability << endl;
    return (double)seed_reachability / times / n;
}
