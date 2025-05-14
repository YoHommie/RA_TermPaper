#include <bits/stdc++.h>
using namespace std;
using namespace chrono;

struct PathInfo {
    vector<int> path;
};

int n, k, m;
vector<vector<int>> undirected_adj, dag_adj;
vector<pair<int, int>> edge_list;

void generateRandomDAG() {
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), std::mt19937(std::random_device{}()));

    vector<int> rank(n);
    for (int i = 0; i < n; ++i)
        rank[perm[i]] = i;

    dag_adj.assign(n, {});
    for (int u = 0; u < n; ++u)
        for (int v : undirected_adj[u])
            if (rank[u] < rank[v])
                dag_adj[u].push_back(v);
}

using SparseMatrix = unordered_map<int, unordered_map<int, PathInfo>>;

SparseMatrix sparseMultiply(const SparseMatrix &A, const SparseMatrix &B) {
    SparseMatrix C;
    for (auto &[u, row] : A) {
        for (auto &[mid, p1] : row) {
            if (B.count(mid)) {
                for (auto &[v, p2] : B.at(mid)) {
                    if (!C[u].count(v)) {
                        C[u][v].path = p1.path;
                        C[u][v].path.insert(C[u][v].path.end(), p2.path.begin() + 1, p2.path.end());
                    }
                }
            }
        }
    }
    return C;
}

SparseMatrix matrixExponentiate(int power) {
    SparseMatrix result;
    for (int u = 0; u < n; ++u)
        result[u][u].path = {u};

    SparseMatrix base;
    for (int u = 0; u < n; ++u)
        for (int v : dag_adj[u])
            base[u][v].path = {u, v};

    while (power) {
        if (power & 1)
            result = sparseMultiply(result, base);
        base = sparseMultiply(base, base);
        power >>= 1;
    }
    return result;
}

bool findCycleFromPaths(const SparseMatrix &A_km1) {
    for (auto &[v, u] : edge_list) {
        if (A_km1.count(u) && A_km1.at(u).count(v)) {
            auto cycle = A_km1.at(u).at(v).path;
            cycle.push_back(u); // close the cycle
            cout << "Simple cycle of length " << k << ":\n";
            for (int x : cycle)
                cout << x << " ";
            cout << endl;
            return true;
        }
    }
    return false;
}

int main() {
    srand(time(0));
    cin >> n >> m >> k;
    undirected_adj.assign(n, {});
    for (int i = 0; i < m; ++i) {
        int u, v;
        cin >> u >> v;
        undirected_adj[u].push_back(v);
        undirected_adj[v].push_back(u);
        edge_list.push_back({u, v});
        edge_list.push_back({v, u});
    }

    int maxTries = k * k * 10;
    auto start = high_resolution_clock::now();
    for (int attempt = 1; attempt <= maxTries; ++attempt) {
        generateRandomDAG();
        auto A_km1 = matrixExponentiate(k - 1);
        if (findCycleFromPaths(A_km1)) {
            auto end = high_resolution_clock::now();
            double runtime = duration<double>(end - start).count();

            cout << "Found: 1\n";
            cout << "Iterations: " << attempt << "\n";
            cout << "Runtime(s): " << fixed << setprecision(6) << runtime << "\n";
            return 0;
        }
    }

    auto end = high_resolution_clock::now();
    double runtime = duration<double>(end - start).count();

    cout << "Found: 0\n";
    cout << "PathLength: 0\n";
    cout << "Iterations: " << maxTries << "\n";
    cout << "Runtime(s): " << fixed << setprecision(6) << runtime << "\n";
    return 0;
}