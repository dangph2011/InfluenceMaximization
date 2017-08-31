    #include <vector>
using namespace std;

//struct Node (Vertex)
typedef struct Vertex{
    int discovery_time;
    int finish_time;
    int predecessor;
    int num_child;
    int low;
    Vertex ()
        : discovery_time(0), finish_time(0), predecessor(-1), num_child(0), low(0) {}
}Vertex;

class PrunedEstimater {
private:
	int n, n1;
	vector<int> weight, comp, sigmas;
	vector<int> pmoc;
	vector<int> at_p;

	vector<int> up;
	vector<bool> memo, removed;

	vector<int> es, rs;
	vector<int> at_e, at_r;

	vector<bool> visited;
	vector<bool> articulation;
    vector<pair<int, int > > forward_edges;
    vector<pair<int, int > > cross_edges;

	int hub;
	vector<bool> descendant, ancestor;
	bool flag;

	int sigma(const int v);
	inline int unique_child(const int v);
public:
	void first();
	void init(int n, vector<pair<int, int> > &es, vector<int> &comp);
	int sigma1(const int v);
	void add(int v);
    void add_reduce(int v);
	void update(vector<long long> &sums);
	void updateGainVertex(long long &gain_v, const int v);
	void articulationPoint();
	void articulationPointDFS();
	void dfsArticulationPoint(int &u, vector<Vertex> &l_node, int &g_time);
    void forwardEdgeDFS();
    void dfsForwardEdge(int &u, vector<Vertex> &p_node, int &g_time);
    void removeEdge();
};

class InfluenceMaximizer {
private:
	int n, m; // |V|, |E|

	vector<int> es1, at_e;
	vector<int> rs1, at_r;

	int scc(vector<int> &comp);
public:
	vector<int> run(vector<pair<pair<int, int>, double> > &es, const int k,
			const int R, double ep);
};

class Evaluater{
private:
	int n, m; // |V|, |E|
	vector<int> es1, at_e;
    vector<bool> removed;
    vector<bool> visited;

    double boundStop(int n, double epsilon);
public:
    bool evaluate(vector<int> seeds, vector<pair<pair<int, int>, double> > &es, double epsilon, double avr_rc1);
    void init(vector<pair<pair<int, int>, double> > &es);
};

// Random Number Generator
class Xorshift {
public:
	Xorshift(int seed) {
		x = _(seed, 0);
		y = _(x, 1);
		z = _(y, 2);
		w = _(z, 3);
	}

	int _(int s, int i) {
		return 1812433253 * (s ^ (s >> 30)) + i + 1;
	}

	inline int gen_int() {
		unsigned int t = x ^ (x << 11);
		x = y;
		y = z;
		z = w;
		return w = w ^ (w >> 19) ^ t ^ (t >> 8);
	}

	inline int gen_int(int n) {
		return (int) (n * gen_double());
	}

	inline double gen_double() {
		unsigned int a = ((unsigned int) gen_int()) >> 5, b =
				((unsigned int) gen_int()) >> 6;
		return (a * 67108864.0 + b) * (1.0 / (1LL << 53));
	}

private:
	unsigned int x, y, z, w;
};
