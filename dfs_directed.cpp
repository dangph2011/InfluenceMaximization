#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include <stack>
#include <queue>
#include <iomanip>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <unistd.h>

enum eColor{
    WHITE,
    GRAY,
    BLACK
};

/*
 * Get current time in milisecond.
 * @return the current time in milisecond.
 */

double getCurrentTimeMlsec(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec *  1e-6;
}

//struct Node (Vertex)
typedef struct Vertex{
    eColor color;
    uint32_t discovery_time;
    uint32_t finish_time;
    int32_t predecessor;
    uint32_t num_child;
    uint32_t low;
    Vertex ()
        : color(WHITE), discovery_time(0), finish_time(0), predecessor(-1), num_child(0), low(0) {}
}Vertex;

typedef struct GraphConnected {
    int32_t graph_id;
    int32_t u_start;
    int32_t v_end;

    GraphConnected()
        : graph_id(-1), u_start(-1), v_end(-1) {}

    GraphConnected(int32_t p_graph_id, int32_t p_u_start, int32_t p_v_end) {
        graph_id = p_graph_id;
        u_start = p_u_start;
        v_end = p_v_end;
    }
}GraphConnected;

class DGraph{
    private:
        int32_t graph_id_;
        uint32_t v_num_;
        uint32_t e_num_;
        std::vector<std::vector<uint32_t> > out_edge_list_;
        std::vector<std::vector<uint32_t> > in_edge_list_;
        std::set<uint32_t> vertex_list_;
        std::vector<uint32_t> reach_;
        std::vector<GraphConnected> graph_connected_;

    public:

        DGraph(){
            graph_id_ = -1;
        }

        DGraph(uint32_t p_size, uint32_t p_graph_id) {
            out_edge_list_.resize(p_size);
            in_edge_list_.resize(p_size);
            reach_.resize(p_size, 1);
            graph_id_ = p_graph_id;
        }

        void setGraphConnected(uint32_t p_graph_id, uint32_t pu_start, uint32_t pv_end) {
            GraphConnected l_gc(p_graph_id, pu_start, pv_end);
            graph_connected_.push_back(l_gc);
        }

        std::vector<GraphConnected> getGraphConnected() {
            return graph_connected_;
        }

        uint32_t getTotalReach() {
            uint32_t total = 0;
            for (auto &it : vertex_list_) {
                total += reach_[it];
            }
            return total;
        }

        void setNumberOfVertex(uint32_t p_vnum){
            v_num_ = p_vnum;
        }

        uint32_t getNumberOfVertex() {
            return vertex_list_.size();
        }

        std::set<uint32_t> getVertexList() {
            return vertex_list_;
        }

        void setNumberOfEdge(uint32_t p_enum) {
            e_num_ = p_enum;
        }

        uint32_t getNumberOfEdge() {
            return out_edge_list_.size();
        }

        std::vector<std::vector<uint32_t> > getOutEdgeList() {
            return out_edge_list_;
        }

        std::vector<std::vector<uint32_t> > getInEdgeList() {
            return in_edge_list_;
        }

        int32_t getGraphId(){
            return graph_id_;
        }

        void setGraphid(int32_t p_id) {
            graph_id_ = p_id;
        }

        void setReach(uint32_t p_vertex, uint32_t p_value){
            if (reach_.size() < p_vertex) {
                reach_.resize(p_vertex + 1, 1);
            }
            reach_[p_vertex] = p_value;
        }

        uint32_t getReach(uint32_t p_vertex) {
            return reach_[p_vertex];
        }

        void resizeOutEdgeList(uint32_t p_size) {
            out_edge_list_.resize(p_size);
        }

        void resizeInEdgeList(uint32_t p_size) {
            in_edge_list_.resize(p_size);
        }

        void resizeReach(uint32_t p_size) {
            reach_.resize(p_size);
        }

		void addVertex(uint32_t u) {
			if (out_edge_list_.size() < (u + 1)) {
				out_edge_list_.resize(u + 1);
			}

			if (reach_.size() < (u + 1)){
				reach_.resize(u+1, 1);
			}

			vertex_list_.insert(u);
		}

		void addEdge(uint32_t u, uint32_t v) {
			if (out_edge_list_.size() < (u + 1)) {
				out_edge_list_.resize(u + 1);
				reach_.resize(u + 1, 1);
			}
			if (out_edge_list_.size() < (v + 1)) {
				out_edge_list_.resize(v + 1);
				reach_.resize(v + 1, 1);
			}

			out_edge_list_[u].push_back(v);

            if (in_edge_list_.size() < (u + 1)) {
                in_edge_list_.resize(u + 1);
                //reach_.resize(u + 1, 1);
            }
            if (in_edge_list_.size() < (v + 1)) {
                in_edge_list_.resize(v + 1);
                //reach_.resize(v + 1, 1);
            }

            in_edge_list_[v].push_back(u);
			//reach_[u] = 1;
			//reach_[v] = 1;
		}

		void delEdge(uint32_t u, uint32_t v) {
			for (auto &it : out_edge_list_[u]) {
				if (it == v) {
					out_edge_list_[u].erase(out_edge_list_[u].begin() + (&it - &out_edge_list_[u][0]));
					break;
				}
			}

            for (auto &it : in_edge_list_[v]) {
				if (it == u) {
					in_edge_list_[v].erase(in_edge_list_[v].begin() + (&it - &in_edge_list_[v][0]));
					break;
				}
			}
		}

		void sortEdges() {
			// //sort adjacency list
			for (auto &it : out_edge_list_) {
				std::sort(it.begin(), it.end());
			}

            for (auto &it : in_edge_list_) {
				std::sort(it.begin(), it.end());
			}
		}

		void sortAndRemoveDuplicateEdges() {
			// //sort adjacency list
			for (auto &it : out_edge_list_) {
				std::sort(it.begin(), it.end());
				it.erase( unique( it.begin(), it.end() ), it.end() );
			}

            for (auto &it : in_edge_list_) {
                std::sort(it.begin(), it.end());
                it.erase( unique( it.begin(), it.end() ), it.end() );
            }
		}

        //print Vertex, Edge and graph Id connected
    	void printBridgeInfo(){
            //TODO
            for (auto &it: graph_connected_) {
        		std::cout << "GraphID=" << it.graph_id << " Start=" << it.u_start << " End=" << it.v_end << std::endl;
            }
    	}

        //print graph
        void printGraph() {
    		//
    		printVertices();
        	printEdges();
        }

    	//print graph
        void printEdges() {
    		//
    		std::cout << "Out List Edges\n";
        	for (auto &it : out_edge_list_){
        		for (auto &it1 : it) {
        			std::cout << &it - &out_edge_list_[0] << " " << it1 << std::endl;
        		}
        	}

            std::cout << "In List Edges\n";
        	for (auto &it : in_edge_list_){
        		for (auto &it1 : it) {
        			std::cout << &it - &in_edge_list_[0] << " " << it1 << std::endl;
        		}
        	}
        }

    	//print graph
        void printVertices() {
    		//
    		std::cout << "List Vertices\n";
        	for (auto &it : vertex_list_){
        		std::cout << it << " " << reach_[it] << std::endl;
        	}
        }
};

bool readEdgeList(DGraph &g, std::string p_file_name){
	std::ifstream f_in;
	f_in.open(p_file_name);
	if (f_in.fail()) {
		std::cerr << "Error: Opening file\n";
		return false;
	}
	uint32_t u, v;

	std::vector<std::pair<uint32_t, uint32_t> > edge_list;

	while (f_in >> u >> v){
		edge_list.push_back(std::pair<unsigned, unsigned>(u,v));
	}

	//Store in adjacency list
	for (auto &it : edge_list) {
        g.addVertex(it.first);
        g.addVertex(it.second);

        g.addEdge(it.first, it.second);
	}

	//g.initAndSetReach();
	// //sort adjacency list
	g.sortAndRemoveDuplicateEdges();
	return true;
}

bool readMetis(DGraph &g, std::string p_file_name) {
	std::ifstream f_in;
	std::string line;
    uint32_t f_vnum;
    uint32_t f_enum;
	uint32_t v;
	uint32_t count = 1;
	f_in.open(p_file_name);
	if (f_in.fail()) {
		std::cerr << "Error: Opening file\n";
		return false;
	}
	getline(f_in, line);
	std::istringstream iss1(line);
	iss1 >> f_vnum >> f_enum;
    g.setNumberOfVertex(f_vnum);
    g.setNumberOfEdge(f_enum);
	//std::cout << "Vertex=" << v_num_ << " Edge=" << e_num_ << std::endl;
    g.resizeOutEdgeList(f_vnum);
    g.resizeInEdgeList(f_vnum);
    g.resizeReach(f_vnum);
	//std::cout << "Vertex=" << v_num_ << " Edge=" << e_num_ << std::endl;
	while (getline(f_in, line)) {
		std::istringstream iss(line);
		if (!line.empty()) {
            g.addVertex(count);
		}

		while (iss >> v) {
            g.addVertex(v);
            g.addEdge(count, v);
		}

		count++;
		//std::cout << count << "\n";
	}

	//g.initAndSetReach();
	//sort and erase duplicate
	g.sortAndRemoveDuplicateEdges();

	return true;
}

//DFS travel
bool stronglyConnectedComponent(DGraph &g, std::vector<DGraph> &p_graph_shatter, DGraph &p_scc) {
	uint32_t g_time = 0;
    std::set<uint32_t> f_vertex_list = g.getVertexList();
    std::vector<std::vector<uint32_t> > f_out_edge_list = g.getOutEdgeList();

	//init vector of Node in graph by size
	std::vector<Vertex> l_node(f_out_edge_list.size());
	std::stack<std::pair<uint32_t, uint32_t> > f_edge_component;
	std::stack<uint32_t> f_node_component;
	uint32_t l_graph_shatter_size_begin = p_graph_shatter.size();
	uint32_t l_graph_shatter_size = p_graph_shatter.size();

	std::vector<std::vector<uint32_t> > l_graph_id(f_out_edge_list.size());
    std::vector<bool> f_stack_member(f_out_edge_list.size(), false);
    std::vector<int> f_component(f_out_edge_list.size(),-1);

	std::stack<uint32_t> S;
	//uint32_t m_num_tree_node;
	bool check = true;

	for (auto &i : f_vertex_list) {
		if (l_node[i].color == WHITE) {
			//m_num_tree_node++;
			//uint32_t s = i;
			S.push(i);
			f_node_component.push(i);
			l_node[i].color = GRAY;
			l_node[i].discovery_time = ++g_time;
			l_node[i].low = l_node[i].discovery_time;
            f_stack_member[i] = true;
			while (!S.empty()) {
				uint32_t u = S.top();
				//S.pop();
                f_stack_member[u] = true;
				check = true;
				for (auto &it : f_out_edge_list[u]) {
					if (l_node[it].color == WHITE) {
						//l_node[it].discovery_time = ++g_time;
						S.push(it);
						f_node_component.push(it);
						l_node[it].color = GRAY;
						l_node[it].discovery_time = ++g_time;
						l_node[it].low = l_node[it].discovery_time;
						l_node[it].predecessor = u;
						l_node[u].num_child++;
						f_edge_component.push(std::pair<uint32_t, uint32_t>(u, it));
                        std::cout << "1 Edge: " << u << " " << it << std::endl;
						check = false;
						break;
					} else {
                        //if ()

                        if (f_stack_member[it] == true) { //check if it
    						l_node[u].low = std::min(l_node[u].low, l_node[it].discovery_time);
                            f_edge_component.push(std::pair<uint32_t, uint32_t>(u, it));
                        }
					}
				}

				if (check) {
					if (l_node[u].color == GRAY) {
						l_node[u].color = BLACK;
						l_node[u].finish_time = ++g_time;

						if (l_node[u].predecessor > -1) {
							auto l_predecessor = l_node[u].predecessor;
							l_node[l_predecessor].low = std::min(l_node[l_predecessor].low, l_node[u].low);

                            //retrieve one component.
							if (l_node[u].low == l_node[u].discovery_time) {
                                DGraph l_g(f_out_edge_list.size(), l_graph_shatter_size);

                                //add edges into graph component
                                auto l_pop_edge = f_edge_component.top();
                                //check if edge pop not bridge
                                while (l_pop_edge.first != l_predecessor || l_pop_edge.second != u) {
                                    //add edge to current graph
                                    l_g.addEdge(l_pop_edge.first, l_pop_edge.second);
                                    f_edge_component.pop();
                                    l_pop_edge = f_edge_component.top();
                                }

                                f_edge_component.pop();
                                //add vertex into graph component
                                auto l_pop_node = f_node_component.top();
                                while (l_pop_node != u) {
                                    l_g.addVertex(l_pop_node);
                                    l_g.setReach(l_pop_node, g.getReach(l_pop_node));
                                    f_component[l_pop_node] = l_graph_shatter_size;
                                    f_stack_member[l_pop_node] = false;
                                    f_node_component.pop();
                                    l_pop_node = f_node_component.top();
                                }
                                //save edge connect two strongly connected component
                                l_g.addVertex(l_pop_node);
                                l_g.setReach(l_pop_node, g.getReach(l_pop_node));
                                f_stack_member[l_pop_node] = false;
                                f_component[l_pop_node] = l_graph_shatter_size;
                                f_node_component.pop();
                                //sort and remove duplicate edges
                                l_g.sortAndRemoveDuplicateEdges();
                                p_graph_shatter.push_back(l_g);
                                //Increase version
                                l_graph_shatter_size++;
                            }
                        }
					}
					S.pop();
				}
			}

            if (l_node[i].low == l_node[i].discovery_time) {
    			DGraph l_g(f_out_edge_list.size(),l_graph_shatter_size);

    			while (!f_edge_component.empty()) {
    				auto l_pop_edge = f_edge_component.top();
    				l_g.addEdge(l_pop_edge.first, l_pop_edge.second);
    				f_edge_component.pop();
    			}

    			while (!f_node_component.empty()) {
    				auto l_pop_node = f_node_component.top();
    				l_g.addVertex(l_pop_node);
    				l_g.setReach(l_pop_node, g.getReach(l_pop_node));
                    f_stack_member[l_pop_node] = false;
                    f_component[l_pop_node] = l_graph_shatter_size;
    				f_node_component.pop();
    			}

    			l_g.sortAndRemoveDuplicateEdges();
    			p_graph_shatter.push_back(l_g);
    			l_graph_shatter_size++;
    		}
		    l_graph_shatter_size_begin = l_graph_shatter_size;
        }

		//m_num_tree_node = g_time - m_num_tree_node;

		//TODO: Update reach for each graph shattered
		// uint32_t f_total_reach_g1;
		// uint32_t f_total_reach_g2;
		// //uint32_t f_new_graph = l_graph_shatter_size - l_graph_shatter_size_begin;
        //
		// //Only recalculate reach if origin graph forms 2 graphs ( bridge exists)
		// if (l_graph_shatter_size - l_graph_shatter_size_begin > 1 ) {
		// 	for (int i = l_graph_shatter_size_begin; i < l_graph_shatter_size - 1; i++) {
		// 		f_total_reach_g1 = p_graph_shatter[i].getTotalReach();
		// 		f_total_reach_g2 = 0;
		// 		for (int j = i+1; j < l_graph_shatter_size; j++) {
		// 			f_total_reach_g2 += p_graph_shatter[j].getTotalReach();
		// 		}
        //
		// 		auto u_bridge = p_graph_shatter[i].getGraphConnected().u_start;
		// 		auto v_bridge = p_graph_shatter[i].getGraphConnected().v_end;
		// 		// std::cout << "G1 " << i << " GraphId= " << p_graph_shatter[i].graph_id_ << " " << l_graph_shatter_size - l_graph_shatter_size_begin  << " " << " GraphConnected=" << p_graph_shatter[i].graph_connected_.graph_id << " size reach=" << p_graph_shatter[i].reach_.size() << std::endl;
		// 		// std::cout << "G2 " << i << " GraphId= " << p_graph_shatter[i+1].graph_id_ << " " << l_graph_shatter_size - l_graph_shatter_size_begin << " " << " GraphConnected=" << p_graph_shatter[i+1].graph_connected_.graph_id << " size reach=" << p_graph_shatter[i+1].reach_.size() << std::endl;
		// 		// std::cout << "UB=" << u_bridge << " VB=" << v_bridge << " IREACH=" << p_graph_shatter[i].reach_.size() << " PIREACH=" << p_graph_shatter[p_graph_shatter[i].graph_connected_.graph_id].reach_.size() << std::endl;
		// 		auto u_new_reach = p_graph_shatter[i].getReach(u_bridge) + f_total_reach_g2;
		// 		auto v_new_reach = p_graph_shatter[p_graph_shatter[i].getGraphConnected().graph_id].getReach(v_bridge) + f_total_reach_g1;
        //
		// 		//Update reach bridge
		// 		p_graph_shatter[i].setReach(u_bridge, u_new_reach);
		// 		p_graph_shatter[p_graph_shatter[i].getGraphConnected().graph_id].setReach(v_bridge, v_new_reach);
        //
		// 		//Update BC score
		// 		//g_bet_cen[u_bridge] += (f_total_reach_g1 - 1) * f_total_reach_g2;
		// 		//g_bet_cen[v_bridge] += (f_total_reach_g2 - 1) * f_total_reach_g1;
		// 	}
		// }

	}

    for (auto &i : f_vertex_list) {
        for (auto &it: f_out_edge_list[i]) {
            if (f_component[i] != f_component[it]) {
                p_graph_shatter[f_component[i]].setGraphConnected(f_component[it], i, it);
                p_scc.addVertex(f_component[i]);
                p_scc.addVertex(f_component[it]);
                p_scc.addEdge(f_component[i], f_component[it]);
            }
        }
    }
    p_scc.sortAndRemoveDuplicateEdges();

    std::cout << "-----DAG start-----\n";
    p_scc.printGraph();
    std::cout << "-----DAG end-----\n";

    std::cout << "DFS traveled\n";
    for (auto &it: l_node) {
        std::cout << "\t" << &it - &l_node[0]<< " " << it.discovery_time << " " << it.finish_time  << " " << it.low << std::endl;
    }
	return true;
}

bool articulationPoints(DGraph &g, std::set<uint32_t> &p_articulation_vertices) {
	uint32_t g_time = 0;
    std::set<uint32_t> f_vertex_list = g.getVertexList();
    std::vector<std::vector<uint32_t> > f_out_edge_list = g.getOutEdgeList();
    std::vector<std::vector<uint32_t> > f_in_edge_list = g.getInEdgeList();
	//init vector of Node in graph by size
	std::vector<Vertex> l_node(f_out_edge_list.size());
	std::stack<uint32_t> S;
	//uint32_t m_num_tree_node;
	bool check_out = true;
    bool check_in = true;

	for (auto &i : f_vertex_list) {
		//l_graph_shatter_size_start = p_graph_shatter.size();
		//l_graph_shatter_size_start = p_graph_shatter.size();
		//m_num_tree_node = g_time;
		//std::vector<DGraph> l_graph;
		//GraphComponent l_graph_component;
		if (l_node[i].color == WHITE) {
			//m_num_tree_node++;
			//uint32_t s = i;
			S.push(i);
			l_node[i].color = GRAY;
			l_node[i].discovery_time = ++g_time;
			l_node[i].low = l_node[i].discovery_time;
			while (!S.empty()) {
				uint32_t u = S.top();
				//S.pop();
				check_out = true;
                check_in = true;
				for (auto &it : f_out_edge_list[u]) {
					if (l_node[it].color == WHITE) {
						//l_node[it].discovery_time = ++g_time;
						S.push(it);
						l_node[it].color = GRAY;
						l_node[it].discovery_time = ++g_time;
						l_node[it].low = l_node[it].discovery_time;
						l_node[it].predecessor = u;
						l_node[u].num_child++;
						//f_edge_component.push(std::pair<uint32_t, uint32_t>(u, it));
                        //std::cout << "1 Edge: " << u << " " << it << std::endl;
						check_out = false;
						break;
					} else if (it != l_node[u].predecessor && u != l_node[it].predecessor) {
						l_node[u].low = std::min(l_node[u].low, l_node[it].discovery_time);
					}
				}

                for (auto &it : f_in_edge_list[u]) {
					if (l_node[it].color == WHITE) {
						//l_node[it].discovery_time = ++g_time;
						S.push(it);
						l_node[it].color = GRAY;
						l_node[it].discovery_time = ++g_time;
						l_node[it].low = l_node[it].discovery_time;
						l_node[it].predecessor = u;
						l_node[u].num_child++;
						//f_edge_component.push(std::pair<uint32_t, uint32_t>(u, it));
                        //std::cout << "2 Edge: " << u << " " << it << std::endl;
						check_out = false;
						break;
					} else if (it != l_node[u].predecessor && u != l_node[it].predecessor) {
						l_node[u].low = std::min(l_node[u].low, l_node[it].discovery_time);
					}
				}

				if (check_out && check_in) {
					if (l_node[u].color == GRAY) {
						l_node[u].color = BLACK;
						l_node[u].finish_time = ++g_time;

                        if (l_node[u].predecessor == -1) {
                            if (l_node[u].num_child > 1) {
                                p_articulation_vertices.insert(u);
                            }
                        } else { //if (l_node[u].predecessor > -1) {
							auto l_predecessor = l_node[u].predecessor;
                            l_node[l_predecessor].low = std::min(l_node[l_predecessor].low, l_node[u].low);
							if (l_node[l_predecessor].predecessor > -1 && l_node[u].low >= l_node[l_predecessor].discovery_time) {
								p_articulation_vertices.insert(l_predecessor);
							}
                        }
					}
					S.pop();
				}
			}
		}
	}

    std::cout << "DFS traveled\n";
    for (auto &it: l_node) {
        std::cout << "\t" << &it - &l_node[0]<< " " << it.low << " " << it.discovery_time << " " << it.finish_time << " "<< it.predecessor << std::endl;
    }
	return true;
}

int main(){
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    DGraph g;
	std::vector<DGraph> m_graph_shatter;

	//Read edge list format
    if (!readEdgeList(g,"data/test5.txt")) {
        return 0;
    }

    g.printGraph();

    double start_time_bridge = getCurrentTimeMlsec();
    std::set<uint32_t> m_articulation_vertices;

    DGraph m_sccs;
    stronglyConnectedComponent(g, m_graph_shatter, m_sccs);
    //print bridge
    std::cout << "SIZE SHATTER=" << m_graph_shatter.size() << "\n";

    for (auto &it : m_graph_shatter) {
    	//std::cout << &it - &m_graph_shatter[0] << "\n";
    	std::cout << "DGraph ID=" << it.getGraphId() <<"\n";
    	it.printGraph();
    	it.printBridgeInfo();
    	std::cout << "Total reach=" << it.getTotalReach() << std::endl;
    	std::cout << "----------\n";
    }

    articulationPoints(m_sccs, m_articulation_vertices);

    std::cout << "Articulation:\n";
    for (auto &it : m_articulation_vertices){
    	std::cout << "\t" << it << "\n";
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Time Bridge=" << getCurrentTimeMlsec() - start_time_bridge << "\n";

    return 0;
}
