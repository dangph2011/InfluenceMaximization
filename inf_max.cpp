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
#include <random>
#include <algorithm>

uint16_t k = 200; //The size of seed
uint16_t R = 10; //The number of DAG
double pro = 0.01; //The propagation rate

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
        uint32_t n_sccs_ = 0;
        std::vector<std::vector<uint32_t> > out_edge_list_;
        std::vector<std::vector<uint32_t> > in_edge_list_;
        std::set<uint32_t> vertex_list_;
        std::vector<uint32_t> reach_;
        std::vector<GraphConnected> graph_connected_;
        std::vector<bool> articulation_points_;
        std::vector<int32_t> vetex_component_mapping_;
        std::vector<std::vector<uint32_t> > component_vertex_;
        std::vector<bool> computed_vertex_;
        std::vector<bool> removed_vertex_;
        std::vector<uint64_t> sigma_;
        std::vector<uint64_t> weight_;
        std::vector<uint32_t> change_;

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

        std::vector<int32_t> getComponent() {
            return vetex_component_mapping_;
        }

        uint32_t getNumberOfSccs() {
            return n_sccs_;
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
            if (out_edge_list_.size() <= v_num_) {
                out_edge_list_.resize(v_num_ + 1);
            }

            if (in_edge_list_.size() <= v_num_) {
                in_edge_list_.resize(v_num_ + 1);
            }
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

        bool removeForwardEdge() {
        	uint32_t g_time = 0;
            //std::set<uint32_t> f_vertex_list_cp = g.getVertexList();
            std::vector<uint32_t> f_vertex_list(vertex_list_.begin(), vertex_list_.end());
            //std::vector<std::vector<uint32_t> > f_out_edge_list = g.getOutEdgeList();
        	//init vector of Node in graph by size
        	std::vector<Vertex> l_node(out_edge_list_.size());
        	std::stack<uint32_t> S;
        	//uint32_t m_num_tree_node;
        	bool check_out = true;
            std::vector<std::pair<uint32_t, uint32_t> > f_forward_edges;
            std::vector<std::pair<uint32_t, uint32_t> > f_cross_edges;
            for (int t = 0 ; t < f_vertex_list.size() ; t++) {
                auto i = f_vertex_list[t];
        	//for (auto &i : f_vertex_list) {
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
        				for (auto &it : out_edge_list_[u]) {
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
        					} else if (l_node[it].color == BLACK && u != l_node[it].predecessor) {
                    			if (l_node[u].discovery_time < l_node[it].discovery_time) {
                    				f_forward_edges.push_back(std::pair<uint32_t, uint32_t>(u, it));
                    			} else {
                                    f_cross_edges.push_back(std::pair<uint32_t, uint32_t>(u, it));
                                }
                            }
        				}

        				if (check_out) {
        					if (l_node[u].color == GRAY) {
        						l_node[u].color = BLACK;
        						l_node[u].finish_time = ++g_time;
        					}
        					S.pop();
        				}
        			}
        		}
        	}

            //Remove Forward Edges
            for (auto &it : f_forward_edges) {
                delEdge(it.first, it.second);
        	}

            // std::cout << "DFS traveled Forward edge\n";
            // for (auto &it: l_node) {
            //     std::cout << "\t" << &it - &l_node[0]<< " " << it.low << " " << it.discovery_time << " " << it.finish_time << " "<< it.predecessor << std::endl;
            // }
            //
            // std::cout << "---FORWARD EDGES---\n";
        	// for (auto &it : f_forward_edges) {
        	// 	std::cout << "\t" << it.first << " " << it.second;
            //     std::cout << "\n";
        	// }
            //
            // std::cout << "---CROSS EDGES---\n";
        	// for (auto &it : f_cross_edges) {
        	// 	std::cout << "\t" << it.first << " " << it.second;
            //     std::cout << "\n";
        	// }

            std::cout << "---FORWARD EDGES---\n";
        	std::cout << "\t" << f_forward_edges.size() << std::endl;

            std::cout << "---CROSS EDGES---\n";
            std::cout << "\t" << f_cross_edges.size() << std::endl;

        	return true;
        }

        bool getArticulationPoints() {
        	uint32_t g_time = 0;
            //std::set<uint32_t> f_vertex_list = g.getVertexList();
            //std::vector<std::vector<uint32_t> > f_out_edge_list = g.getOutEdgeList();
            //std::vector<std::vector<uint32_t> > f_in_edge_list = g.getInEdgeList();
            articulation_points_.assign(out_edge_list_.size(), false);
        	//init vector of Node in graph by size
        	std::vector<Vertex> l_node(out_edge_list_.size());

        	std::stack<uint32_t> S;
        	//uint32_t m_num_tree_node;
        	bool check_out = true;
            bool check_in = true;

        	for (auto &i : vertex_list_) {
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

                        //travel out going edge array first
        				for (auto &it : out_edge_list_[u]) {
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
        					} else if (it != l_node[u].predecessor) {
        						l_node[u].low = std::min(l_node[u].low, l_node[it].discovery_time);
        					}
        				}

                        //travel all out going edge of node u
                        if (check_out) {
                            for (auto &it : in_edge_list_[u]) {
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
            						check_in = false;
            						break;
            					} else if (it != l_node[u].predecessor) {
            						l_node[u].low = std::min(l_node[u].low, l_node[it].discovery_time);
            					}
            				}
                        }

                        if (check_out && check_in) {
        					if (l_node[u].color == GRAY) {
        						l_node[u].color = BLACK;
        						l_node[u].finish_time = ++g_time;

                                if (l_node[u].predecessor == -1) {
                                    if (l_node[u].num_child > 1) {
                                        articulation_points_[u] = true;
                                    }
                                } else { //if (l_node[u].predecessor > -1) {
        							auto l_predecessor = l_node[u].predecessor;
                                    l_node[l_predecessor].low = std::min(l_node[l_predecessor].low, l_node[u].low);
        							if (l_node[l_predecessor].predecessor > -1 && l_node[u].low >= l_node[l_predecessor].discovery_time) {
        								articulation_points_[l_predecessor] = true;
        							}
                                }
        					}
        					S.pop();
        				}
        			}
        		}
        	}

            auto n_ar = 0;
            for (auto it : articulation_points_) {
                if (it == true) {
                    n_ar++;
                }
            }
            std::cout << "Number of Articulation = " << n_ar << std::endl;
            // for (auto &it: l_node) {
            //     std::cout << "\t" << &it - &l_node[0]<< " " << it.low << " " << it.discovery_time << " " << it.finish_time << " "<< it.predecessor << std::endl;
            // }
        	return true;
        }

        //DFS travel
        bool stronglyConnectedComponent(std::vector<DGraph> &p_graph_shatter, DGraph &p_scc) {
            n_sccs_ = 0;
        	uint32_t g_time = 0;
            vetex_component_mapping_.assign(out_edge_list_.size(), -1);
            //std::set<uint32_t> f_vertex_list = g.getVertexList();
            //std::vector<std::vector<uint32_t> > f_out_edge_list = g.getOutEdgeList();

        	//init vector of Node in graph by size
        	std::vector<Vertex> l_node(out_edge_list_.size());
        	std::stack<std::pair<uint32_t, uint32_t> > f_edge_component;
        	std::stack<uint32_t> f_node_component;

        	uint32_t l_graph_shatter_size_begin = p_graph_shatter.size();
        	uint32_t l_graph_shatter_size = p_graph_shatter.size();

        	std::vector<std::vector<uint32_t> > l_graph_id(out_edge_list_.size());
            std::vector<bool> f_stack_member(out_edge_list_.size(), false);
            //std::vector<int> f_component(out_edge_list_.size(),-1);

        	std::stack<uint32_t> S;
        	//uint32_t m_num_tree_node;
        	bool check = true;

            for (auto i = 0; i < out_edge_list_.size(); i++) {
        	//for (auto &i : vertex_list_) {
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
        				for (auto &it : out_edge_list_[u]) {
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
                                //std::cout << "1 Edge: " << u << " " << it << std::endl;
        						check = false;
        						break;
        					} else {
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
                                        DGraph l_g(out_edge_list_.size(), l_graph_shatter_size);

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
                                            l_g.setReach(l_pop_node, getReach(l_pop_node));
                                            vetex_component_mapping_[l_pop_node] = n_sccs_;
                                            f_stack_member[l_pop_node] = false;
                                            f_node_component.pop();
                                            l_pop_node = f_node_component.top();
                                        }
                                        //save edge connect two strongly connected component
                                        l_g.addVertex(l_pop_node);
                                        l_g.setReach(l_pop_node, getReach(l_pop_node));
                                        f_stack_member[l_pop_node] = false;
                                        vetex_component_mapping_[l_pop_node] = n_sccs_;
                                        f_node_component.pop();
                                        //sort and remove duplicate edges
                                        l_g.sortAndRemoveDuplicateEdges();
                                        p_graph_shatter.push_back(l_g);
                                        //Increase version
                                        l_graph_shatter_size++;
                                        n_sccs_++;
                                    }
                                }
        					}
        					S.pop();
        				}
        			}

                    if (l_node[i].low == l_node[i].discovery_time) {
            			DGraph l_g(out_edge_list_.size(),l_graph_shatter_size);

            			while (!f_edge_component.empty()) {
            				auto l_pop_edge = f_edge_component.top();
            				l_g.addEdge(l_pop_edge.first, l_pop_edge.second);
            				f_edge_component.pop();
            			}

            			while (!f_node_component.empty()) {
            				auto l_pop_node = f_node_component.top();
            				l_g.addVertex(l_pop_node);
            				l_g.setReach(l_pop_node, getReach(l_pop_node));
                            f_stack_member[l_pop_node] = false;
                            vetex_component_mapping_[l_pop_node] = n_sccs_;
            				f_node_component.pop();
            			}

            			l_g.sortAndRemoveDuplicateEdges();
            			p_graph_shatter.push_back(l_g);
            			l_graph_shatter_size++;
                        n_sccs_++;
            		}
                }
        	}

            for (auto i = 0; i < out_edge_list_.size(); i++) {
                for (auto &it: out_edge_list_[i]) {
                    if (vetex_component_mapping_[i] != vetex_component_mapping_[it]) {
                        p_graph_shatter[vetex_component_mapping_[i]].setGraphConnected(vetex_component_mapping_[it], i, it);
                        p_scc.addVertex(vetex_component_mapping_[i]);
                        p_scc.addVertex(vetex_component_mapping_[it]);
                        p_scc.addEdge(vetex_component_mapping_[i], vetex_component_mapping_[it]);
                    }
                }
            }

            p_scc.sortAndRemoveDuplicateEdges();

            //set number of strongly connected component
            //n_sccs_ = l_graph_shatter_size - l_graph_shatter_size_begin;
            // std::cout << "-----DAG start-----\n";
            // p_scc.printGraph();
            // std::cout << "-----DAG end-----\n";
            //
            // std::cout << "DFS traveled DAG\n";
            // for (auto &it: l_node) {
            //     std::cout << "\t" << &it - &l_node[0]<< " " << it.discovery_time << " " << it.finish_time  << " " << it.low << std::endl;
            // }
        	return true;
        }

        void init(uint32_t p_n, uint32_t p_nscc, std::vector<int32_t> p_component ){
            computed_vertex_.assign(p_n, false);
            removed_vertex_.assign(p_n, false);
            sigma_.assign(p_n, 0);
            weight_.assign(p_nscc, 0);
            component_vertex_.resize(p_nscc);
            vetex_component_mapping_ = p_component;
            //calculate the number of vertex in each component
            for (auto i = 0; i < p_n; i++) {
                weight_[vetex_component_mapping_[i]]++;
                component_vertex_[vetex_component_mapping_[i]].push_back(i);
            }
            // for (auto &it : weight_) {
            //     std::cout << "\t" << it << std::endl;
            // }

            //set change
            for (auto i = 0; i < p_n; i++) {
                change_.push_back(i);
            }
        }

        void initGain(std::vector<uint64_t> &gain) {
            for (auto i = 0; i < change_.size(); i++) {
               auto v = change_[i];
               gain[v] += computeSigmaVertex(v);
            }
        }

        void updateGain(std::vector<uint64_t> &gain) {
            for (auto i = 0; i < change_.size(); i++) {
               auto v = change_[i];
               gain[v] -= sigma_[vetex_component_mapping_[v]];
            }

            for (auto i = 0; i < change_.size(); i++) {
               auto v = change_[i];
               gain[v] += computeSigmaVertex(v);
            }
        }

        uint64_t computeSigmaVertex(uint32_t v) {
            return computeSigmaComponent(vetex_component_mapping_[v]);
        }

        uint64_t computeSigmaComponent(uint32_t v) {
            if (computed_vertex_[v]) {
                return sigma_[v];
            }
            computed_vertex_[v] = true;
            if (removed_vertex_[v]) {
                return sigma_[v] = 0;
            }

            std::stack<uint32_t> Q;
            std::vector<bool> checked(out_edge_list_.size(), false);
            Q.push(v);
            auto flag = true;
            uint64_t delta = 0;
            //DFS
            while (!Q.empty()) {
                auto v0 = Q.top();
                flag = true;
                for (auto &it: out_edge_list_[v0]) {
                    //std::cout << "Test ";

                    std::cout << removed_vertex_.size() << std::endl;
                    std::cout << it << " " << removed_vertex_.size() << std::endl;
                    if (removed_vertex_[it]) {
						continue;
					}
                    //std::cout << "Test \n";
                    if (!checked[it]) {
                        checked[it] = true;
                        if (articulation_points_[it]) {
                            delta += computeSigmaComponent(it);
                            continue;
                        }
                        Q.push(it);
                        flag = false;
                        break;
                    }
                }

                if (flag) {
                    delta += weight_[v0];
                    Q.pop();
                }
            }
            return sigma_[v] = delta;
        }

        void removedVertex(uint32_t v) {
            //map v to scc contain v
            v = vetex_component_mapping_[v];
            //std::cout << "Component=" << v << std::endl;
            //remove tree from root v
            //Using DFS to traveled vertex from root v
            std::vector<uint32_t> l_remove_vertex;
            std::stack<uint32_t> Q;
            Q.push(v);
            //removed_vertex_[v] = true;
            while (!Q.empty()) {
                auto v0 = Q.top();
                auto flag = true;
                for (auto &it : out_edge_list_[v0]) {
                    if (!removed_vertex_[it]) {
                        Q.push(it);
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    removed_vertex_[v0] = true;
                    l_remove_vertex.push_back(v0);
                    Q.pop();
                }
            }

            // std::cout << "Remove:" << std::endl;
            // for (auto &it : l_remove_vertex)
            // {
            //     std::cout << "\t" << it << std::endl;
            // }

            std::vector<bool> checked(removed_vertex_.size(), false);
            //Reverse DFS to set vertex that reduces gain bases on v
            change_.clear();
            for (auto &it : l_remove_vertex) {
                Q.push(it);
                while(!Q.empty()) {
                    auto v0 = Q.top();
                    auto flag = true;
                    for (auto &it1 : in_edge_list_[v0]) {
                       if (!removed_vertex_[it1] && !checked[it1]) {
                           Q.push(it1);
                           checked[it1] = true;
                           flag = false;
                           break;
                       }
                    }
                    if (flag) {
                        computed_vertex_[v0] = false;
                        for (auto v_origin : component_vertex_[v0]) {
                            change_.push_back(v_origin);
                        }
                        Q.pop();
                    }
                }
            }

            // std::cout << "Change=" << change_.size() << std::endl;
            // for (auto &it : change_)
            // {
            //     std::cout << "\t" << it << std::endl;
            // }
        }
};

bool readEdgeListPro(std::vector<std::pair<std::pair<uint32_t, uint32_t>, double > > &p_edge_list_pro, std::string p_file_name) {
	std::ifstream f_in;
	f_in.open(p_file_name);
	if (f_in.fail()) {
		std::cerr << "Error: Opening file\n";
		return false;
	}
	uint32_t u, v;
    double p;

	//std::vector<std::pair<uint32_t, uint32_t> > edge_list;

	while (f_in >> u >> v >> p){
        if (u == v)
            continue;
        p_edge_list_pro.push_back(std::make_pair(std::make_pair(u, v), p));
	}

	return true;
}

bool readEdgeList(std::string p_file_name, std::vector<std::pair<uint32_t, uint32_t> > &edge_list, uint32_t &p_max_vertex){
	std::ifstream f_in;
	f_in.open(p_file_name);
	if (f_in.fail()) {
		std::cerr << "Error: Opening file\n";
		return false;
	}
	uint32_t u, v;
    p_max_vertex = 0;
	//std::vector<std::pair<uint32_t, uint32_t> > edge_list;

	while (f_in >> u >> v){
		edge_list.push_back(std::pair<unsigned, unsigned>(u,v));
        if (p_max_vertex <= u) {
            p_max_vertex = u + 1;
        }

        if (p_max_vertex <= v) {
            p_max_vertex = v + 1;
        }
	}

	//Store in adjacency list
	// for (auto &it : edge_list) {
    //     g.addVertex(it.first);
    //     g.addVertex(it.second);
    //
    //     g.addEdge(it.first, it.second);
	// }
    //
	// //g.initAndSetReach();
	// // //sort adjacency list
	// g.sortAndRemoveDuplicateEdges();
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

int main(int argc, char **argv) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    //srand(time(NULL));
    // std::vector<std::pair<std::pair<uint32_t, uint32_t>, double> > m_edge_list_pro;
    // if (!readEdgeListPro(m_edge_list_pro, "data/test5.txt")) {
    //     return 0;
    // }
    //
    // std::vector<DGraph>
    //DGraph g;

    if (argc < 3) {
		std::cerr << "./inf_max k R" << std::endl;
		exit(1);
	}

    k = atoi(argv[1]);
	R = atoi(argv[2]);

    std::cerr << k << " " << R << std::endl;

    std::vector<std::pair<uint32_t, uint32_t> > m_edge_list;
    uint32_t n = 0;

	//Read edge list format
    if (!readEdgeList("data/wiki-Vote.txt", m_edge_list, n)) {
       return 0;
    }

    //store scc info
    std::vector<DGraph> m_graph_shatter;
    double start_time = getCurrentTimeMlsec();
    double scc_time = 0;
    //http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1);

    std::vector<uint64_t> m_gain(n, 0);
    std::vector<uint32_t> m_seeds;

    std::vector<DGraph> m_sccs(R);
    std::vector<DGraph> g(R);
    for (auto i = 0; i < R; i++) {
        //http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
        for (auto &it : m_edge_list) {
            //Use dis to transform the random unsigned int generated by gen into a double in [0, 1)
            if (dis(gen) < pro) {
                g[i].addEdge(it.first, it.second);
            }
        }
        g[i].setNumberOfVertex(n);
        double start_time_scc = getCurrentTimeMlsec();
        g[i].stronglyConnectedComponent(m_graph_shatter, m_sccs[i]);
        scc_time += getCurrentTimeMlsec() - start_time_scc;
        //mapping between original vertex and scc (DAG)
        m_sccs[i].init(n, g[i].getNumberOfSccs(), g[i].getComponent());
        m_sccs[i].removeForwardEdge();
        m_sccs[i].getArticulationPoints();

        m_sccs[i].initGain(m_gain);
    }

    auto next = 0;
    //first seeds
    for (auto i = 0; i < n; i++) {
        if (m_gain[i] > m_gain[next]) {
            next = i;
        }
    }
    m_seeds.push_back(next);

    for (auto i = 1; i < k; i++) {
        for (auto i = 0; i < R; i++) {
            m_sccs[i].removedVertex(next);
            m_sccs[i].updateGain(m_gain);
        }
        next = 0;
        for (auto i = 0; i < n; i++) {
            if (m_gain[i] > m_gain[next]) {
                next = i;
            }
        }
        m_seeds.push_back(next);
    }

    //Seed array
    std::cout << "Seeds:" << std::endl;
    for (auto &it : m_seeds)
    {
        std::cout << "\t" << it << std::endl;
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Time run=" << getCurrentTimeMlsec() - start_time << "\n";
    std::cout << "Strongly Connected Time=" << scc_time << "\n";
    return 0;
}
