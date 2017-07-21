#include <vector>
using namespace std;

class Evaluater{
private:
	int n, m; // |V|, |E|
	vector<int> es1, at_e;
    vector<bool> removed;
    vector<bool> visited;

    double boundStop(int n, double epsilon);
public:
    double evaluate(vector<int> seeds, vector<pair<pair<int, int>, double> > &es, double epsilon);
};
