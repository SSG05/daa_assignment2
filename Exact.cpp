#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <chrono>
#include <functional>

using namespace std;

struct Edge {
    int to, rev;
    long long cap, flow;
    Edge(int t, long long c, int r) : to(t), cap(c), flow(0), rev(r) {}
};

class DensestSubgraph {
private:
    int n; 
    vector<vector<int>> graph; 
    vector<int> best_subgraph;
    double best_density; 
    int h; 
    vector<vector<int>> h_minus_1_cliques; 
    vector<vector<Edge>> flow_graph;
    vector<int> level, ptr;
    int source, sink;

public:
    DensestSubgraph(int num_vertices, int clique_size)
        : n(num_vertices), h(clique_size), best_density(0) {
        graph.resize(n);
    }

    void addEdge(int u, int v) {
        graph[u].push_back(v);
        graph[v].push_back(u);
    }

    void prepareGraph() {
        for (int i = 0; i < n; i++) {
            sort(graph[i].begin(), graph[i].end());
            auto it = unique(graph[i].begin(), graph[i].end());
            graph[i].resize(distance(graph[i].begin(), it));
        }
    }

    bool hasEdge(int u, int v) {
        return binary_search(graph[u].begin(), graph[u].end(), v);
    }

    void findHMinus1Cliques() {
        auto start_time = chrono::high_resolution_clock::now();
        cout << "Finding all (h-1)-cliques..." << endl;
        vector<int> current_clique;
        
        function<void(int)> dfs = [&](int start) {
            if (current_clique.size() == h - 1) {
                h_minus_1_cliques.push_back(current_clique);
                return;
            }

            for (int v = start; v < n; v++) {
                bool can_add = true;
                for (int u : current_clique) {
                    if (!hasEdge(u, v)) {
                        can_add = false;
                        break;
                    }
                }

                if (can_add) {
                    current_clique.push_back(v);
                    dfs(v + 1);
                    current_clique.pop_back();
                }
            }
        };

        dfs(0);
        
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();

        cout << "Found " << h_minus_1_cliques.size() << " (h-1)-cliques in "
             << duration / 1000.0 << " seconds" << endl;
    }

    bool formsHClique(const vector<int>& clique, int v) {
        for (int u : clique) {
            if (!hasEdge(u, v)) {
                return false;
            }
        }
        return true;
    }

    int degree(int v) {
        return graph[v].size();
    }

    int maxDegree() {
        int max_deg = 0;
        for (int v = 0; v < n; v++) {
            max_deg = max(max_deg, (int)graph[v].size());
        }
        return max_deg;
    }

    void addFlowEdge(int from, int to, long long capacity) {
        flow_graph[from].push_back(Edge(to, capacity, flow_graph[to].size()));
        flow_graph[to].push_back(Edge(from, 0, flow_graph[from].size() - 1));
    }

    bool bfs() {
        level.assign(flow_graph.size(), -1);
        level[source] = 0;

        queue<int> q;
        q.push(source);

        while (!q.empty() && level[sink] == -1) {
            int u = q.front();
            q.pop();

            for (auto& e : flow_graph[u]) {
                if (level[e.to] == -1 && e.flow < e.cap) {
                    level[e.to] = level[u] + 1;
                    q.push(e.to);
                }
            }
        }

        return level[sink] != -1;
    }

    long long dfs(int u, long long pushed) {
        if (u == sink || pushed == 0)
            return pushed;

        for (int& cid = ptr[u]; cid < flow_graph[u].size(); cid++) {
            Edge& e = flow_graph[u][cid];

            if (level[e.to] != level[u] + 1 || e.flow >= e.cap)
                continue;

            long long new_flow = dfs(e.to, min(pushed, e.cap - e.flow));

            if (new_flow > 0) {
                e.flow += new_flow;
                flow_graph[e.to][e.rev].flow -= new_flow;
                return new_flow;
            }
        }

        return 0;
    }

    long long maxFlow() {
        long long flow = 0;

        while (bfs()) {
            ptr.assign(flow_graph.size(), 0);
            long long new_flow;

            while ((new_flow = dfs(source, numeric_limits<long long>::max())) > 0) {
                flow += new_flow;
            }
        }

        return flow;
    }

    vector<int> findMinCut() {
        vector<bool> is_in_s(flow_graph.size(), false);
        queue<int> q;
        q.push(source);
        is_in_s[source] = true;

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (auto& e : flow_graph[u]) {
                if (!is_in_s[e.to] && e.flow < e.cap) {
                    is_in_s[e.to] = true;
                    q.push(e.to);
                }
            }
        }

        vector<int> s_vertices;
        for (int v = 0; v < n; v++) {
            if (is_in_s[v + 1]) {
                s_vertices.push_back(v);
            }
        }

        return s_vertices;
    }

    void buildFlowNetwork(double alpha) {
        int num_vertices = 1 + n + h_minus_1_cliques.size() + 1; 
        flow_graph.clear();
        flow_graph.resize(num_vertices);
        source = 0;
        sink = num_vertices - 1;
        int clique_base_idx = n + 1;
        for (int v = 0; v < n; v++) {
            addFlowEdge(source, v + 1, degree(v));
        }

        for (int v = 0; v < n; v++) {
            addFlowEdge(v + 1, sink, (long long)(alpha * n));
        }

        for (int c = 0; c < h_minus_1_cliques.size(); c++) {
            int clique_idx = clique_base_idx + c;
            const vector<int>& clique = h_minus_1_cliques[c];
            for (int v : clique) {
                addFlowEdge(clique_idx, v + 1, numeric_limits<long long>::max() / 2);
            }

            for (int v = 0; v < n; v++) {
                if (find(clique.begin(), clique.end(), v) != clique.end())
                    continue;
                if (formsHClique(clique, v)) {
                    addFlowEdge(v + 1, clique_idx, 1);
                }
            }
        }
    }

    vector<int> findDensestSubgraph() {
        auto start_time = chrono::high_resolution_clock::now();
        prepareGraph();
        findHMinus1Cliques();
        auto end_time = chrono::high_resolution_clock::now();
        auto clique_duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
        cout << "Finding (h-1)-cliques took " << clique_duration / 1000.0 << " seconds" << endl;
        start_time = chrono::high_resolution_clock::now();
        double l = 0;
        double u = n * n;
        int iterations = 0;

        while (u - l > 1.0 / (n * (n - 1))) {
            auto iter_start = chrono::high_resolution_clock::now();
            iterations++;

            double alpha = (l + u) / 2;
            cout << "Testing alpha = " << alpha << " (iteration " << iterations << ")" << endl;

            buildFlowNetwork(alpha);
            maxFlow();

            vector<int> subgraph = findMinCut();

            if (subgraph.empty()) {
                u = alpha;
            } else {
                l = alpha;
                double density = calculateDensity(subgraph);
                if (density > best_density) {
                    best_density = density;
                    best_subgraph = subgraph;
                    cout << "Found better subgraph with density " << density
                         << " and " << subgraph.size() << " vertices" << endl;
                }
            }

            auto iter_end = chrono::high_resolution_clock::now();
            auto iter_duration = chrono::duration_cast<chrono::milliseconds>(iter_end - iter_start).count();
            cout << "Iteration took " << iter_duration / 1000.0 << " seconds" << endl;
        }

        end_time = chrono::high_resolution_clock::now();
        auto search_duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
        cout << "Binary search took " << search_duration / 1000.0 << " seconds ("
             << iterations << " iterations)" << endl;

        return best_subgraph;
    }
    double calculateDensity(const vector<int>& vertices) {
        if (vertices.empty()) return 0;
        if (h == 2) {
            int edge_count = 0;
            unordered_set<int> vertex_set(vertices.begin(), vertices.end());

            for (int u : vertices) {
                for (int v : graph[u]) {
                    if (vertex_set.count(v) && u < v) { 
                        edge_count++;
                    }
                }
            }
            return (double)edge_count / vertices.size();
        }
        else {
            int clique_count = countHCliques(vertices);
            return (double)clique_count / vertices.size();
        }
    }

    int countHCliques(const vector<int>& vertices) {
        if (vertices.size() < h) return 0;

        unordered_set<int> vertex_set(vertices.begin(), vertices.end());
        int count = 0;
        vector<int> current;

        function<void(int)> dfs = [&](int start) {
            if (current.size() == h) {
                count++;
                return;
            }
            for (int i = start; i < vertices.size(); i++) {
                int v = vertices[i];
                bool can_add = true;

                for (int u : current) {
                    if (!hasEdge(u, v)) {
                        can_add = false;
                        break;
                    }
                }
                if (can_add) {
                    current.push_back(v);
                    dfs(i + 1);
                    current.pop_back();
                }
            }
        };
        dfs(0);
        return count;
    }
    pair<vector<int>, double> getResult() {
        return {best_subgraph, best_density};
    }
};

int main(int argc, char* argv[]) {
    auto total_start_time = chrono::high_resolution_clock::now();
    
    if (argc < 2 || argc > 3) {
        cerr << "Usage: " << argv[0] << " <input_file> [clique_size]" << endl;
        cerr << "Default clique_size is 2 (edge density)" << endl;
        return 1;
    }
    
    string filename = argv[1];
    ifstream input_file(filename);
    
    if (!input_file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return 1;
    }
    int h = 2;
    if (argc >= 3) {
        h = stoi(argv[2]);
        if (h < 2) {
            cerr << "Clique size must be at least 2" << endl;
            return 1;
        }
    }
    vector<pair<int, int>> edges;
    int u, v;
    int max_vertex = 0;
    
    while (input_file >> u >> v) {
        edges.push_back({u, v});
        max_vertex = max(max_vertex, max(u, v));
    }

    int n = max_vertex + 1;
    
    cout << "Graph has " << n << " vertices and " << edges.size() << " edges" << endl;
    cout << "Finding densest subgraph for " << h << "-clique density..." << endl;
    
    DensestSubgraph solver(n, h);
    for (auto& edge : edges) {
        solver.addEdge(edge.first, edge.second);
    }
    vector<int> densest = solver.findDensestSubgraph();
    auto result = solver.getResult();
    vector<int> vertices = result.first;
    double density = result.second;
    cout << "\nResults:" << endl;
    if (!vertices.empty()) {
        cout << "Densest subgraph contains " << vertices.size() << " vertices:" << endl;
        sort(vertices.begin(), vertices.end());
        for (int v : vertices) {
            cout << v << " ";
        }
        cout << "\nSubgraph " << h << "-clique density: " << density << endl;
    } else {
        cout << "No valid solution found." << endl;
    }
    auto total_end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(total_end_time - total_start_time).count();
    cout << "Total execution time: " << duration / 1000.0 << " seconds" << endl;
    
    return 0;
}