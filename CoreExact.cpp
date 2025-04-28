#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <set>
#include <utility>
#include <limits>
#include <algorithm>
#include <cmath>

using namespace std;
using Edge = pair<int, int>;
using Graph = vector<vector<int>>;

Graph readGraphFromFile(const string& filename, int& numVertices);
double rho(const Graph& G, const set<int>& U, int h);
Graph coreDecomposition(const Graph& G);
Graph pruneByPsiCore(const Graph& G);
pair<set<int>, set<int>> minSTCut(const Graph& F, int source);
Graph buildFlowNetwork(const Graph& G, const set<int>& VC);
int findMaxCore(const Graph& G);
Graph inducedSubgraph(const Graph& G, const set<int>& vertices);

set<int> CoreExact(Graph& G, Graph& Psi, int h, double& bestDensity) {
    Graph coreG = coreDecomposition(G);
    Graph psiCoreG = pruneByPsiCore(coreG);
    vector<set<int>> C;
    set<int> bestSubgraph, U;
    bestDensity = rho(psiCoreG, {}, h);
    int u = findMaxCore(G);
    vector<bool> visited(psiCoreG.size(), false);
    for (size_t i = 0; i < psiCoreG.size(); ++i) {
        if (!visited[i] && !psiCoreG[i].empty()) {
            set<int> component;
            queue<int> q;
            q.push(i);
            visited[i] = true;

            while (!q.empty()) {
                int v = q.front();
                q.pop();
                component.insert(v);
                for (int neighbor : psiCoreG[v]) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        q.push(neighbor);
                    }
                }
            }
            C.push_back(component);
        }
    }
    for (auto& VC : C) {
        if (bestDensity > findMaxCore(G)) {
            set<int> newVC;
            Graph lPsiCore = pruneByPsiCore(coreG);
            for (int v : VC) {
                if (!lPsiCore[v].empty()) newVC.insert(v);
            }
            VC = newVC;
        }
        Graph F = buildFlowNetwork(G, VC);
        auto [S, T] = minSTCut(F, 0);
        if (S.empty()) continue;
        while (u - bestDensity >= 1.0 / (VC.size() * (VC.size() - 1))) {
            double alpha = (bestDensity + u) / 2.0;

            F = buildFlowNetwork(G, VC);
            tie(S, T) = minSTCut(F, 0);
            if (S.size() == 1) {
                u = alpha;
            } else {
                if (alpha > bestDensity) {
                    for (auto it = VC.begin(); it != VC.end();) {
                        if (rand() % 2) {
                            it = VC.erase(it);
                        } else {
                            ++it;
                        }
                    }
                }
                bestDensity = alpha;
                U = S;
                U.erase(0);
                if (rho(G, U, h) > rho(G, bestSubgraph, h)) {
                    bestSubgraph = U;
                }
            }
        }
    }
    return bestSubgraph;
}

Graph readGraphFromFile(const string& filename, int& numVertices) {
    ifstream infile(filename);
    string line;
    Graph G;
    int maxNode = 0;
    vector<Edge> edges;
    while (getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        int u, v;
        ss >> u >> v;
        edges.emplace_back(u, v);
        maxNode = max({maxNode, u, v});
    }
    G.resize(maxNode + 1);
    for (auto [u, v] : edges) {
        G[u].push_back(v);
        G[v].push_back(u);
    }
    numVertices = maxNode + 1;
    return G;
}
double rho(const Graph& G, const set<int>& U, int h) {
    if (U.empty()) return 0.0;
    if (h == 2) {
        int edges = 0;
        for (int u : U) {
            for (int v : G[u]) {
                if (U.count(v)) edges++;
            }
        }
        return edges / 2.0 / U.size();
    }
    else if (h == 3) {
        int triangles = 0;
        for (int u : U) {
            for (size_t i = 0; i < G[u].size(); ++i) {
                int v = G[u][i];
                if (!U.count(v)) continue;
                for (size_t j = i + 1; j < G[u].size(); ++j) {
                    int w = G[u][j];
                    if (!U.count(w)) continue;
                    for (int x : G[v]) {
                        if (x == w) {
                            triangles++;
                            break;
                        }
                    }
                }
            }
        }
        return triangles / 3.0 / U.size();
    } else {
        cout << "h > 3 not supported yet." << endl;
        exit(1);
    }
}

Graph coreDecomposition(const Graph& G) {
    int n = G.size();
    vector<int> degree(n);
    for (int u = 0; u < n; ++u) {
        degree[u] = G[u].size();
    }
    int maxDeg = *max_element(degree.begin(), degree.end());
    vector<bool> removed(n, false);
    int k = maxDeg;
    while (k >= 1) {
        bool changed = true;
        while (changed) {
            changed = false;
            for (int u = 0; u < n; ++u) {
                if (!removed[u] && degree[u] < k) {
                    removed[u] = true;
                    changed = true;
                    for (int v : G[u]) {
                        if (!removed[v]) {
                            degree[v]--;
                        }
                    }
                }
            }
        }

        int remaining = count(removed.begin(), removed.end(), false);
        if (remaining > 0) break;
        fill(removed.begin(), removed.end(), false);
        for (int u = 0; u < n; ++u) degree[u] = G[u].size();
        k--;
    }
    Graph coreG(n);
    for (int u = 0; u < n; ++u) {
        if (!removed[u]) {
            for (int v : G[u]) {
                if (!removed[v]) {
                    coreG[u].push_back(v);
                }
            }
        }
    }
    return coreG;
}

Graph pruneByPsiCore(const Graph& G) {
    return G;
}

Graph buildFlowNetwork(const Graph& G, const set<int>& VC) {
    Graph F(G.size());
    for (int v : VC) {
        for (int u : G[v]) {
            if (VC.count(u)) {
                F[v].push_back(u);
            }
        }
    }
    return F;
}

pair<set<int>, set<int>> minSTCut(const Graph& F, int source) {
    set<int> S, T;
    vector<bool> visited(F.size(), false);
    queue<int> q;
    q.push(source);
    visited[source] = true;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        S.insert(u);
        for (int v : F[u]) {
            if (!visited[v]) {
                visited[v] = true;
                q.push(v);
            }
        }
    }
    for (int i = 0; i < F.size(); ++i) {
        if (!visited[i]) {
            T.insert(i);
        }
    }
    return {S, T};
}

int findMaxCore(const Graph& G) {
    int maxDeg = 0;
    for (const auto& neighbors : G) {
        maxDeg = max(maxDeg, static_cast<int>(neighbors.size()));
    }
    return maxDeg;
}

Graph inducedSubgraph(const Graph& G, const set<int>& vertices) {
    Graph H(vertices.size());
    vector<int> map(G.size(), -1);
    int idx = 0;
    for (int v : vertices) {
        map[v] = idx++;
    }
    for (int v : vertices) {
        for (int u : G[v]) {
            if (vertices.count(u)) {
                H[map[v]].push_back(map[u]);
            }
        }
    }
    return H;
}
int main() {
    string filename = "as20000102.txt";
    int h;

    cout << "Enter the value of h (>=2): ";
    cin >> h;
    if (h < 2) {
        cout << "Invalid h value!" << endl;
        return 1;
    }

    int numVertices;
    Graph G = readGraphFromFile(filename, numVertices);
    Graph Psi(numVertices);
    double bestDensity;
    set<int> bestSubgraph = CoreExact(G, Psi, h, bestDensity);
    cout << "\nBest Subgraph Vertices:" << endl;
    for (int v : bestSubgraph) {
        cout << v << " ";
    }
    cout << endl;
    cout << "Clique Density: " << bestDensity << endl;

    return 0;
}