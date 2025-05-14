#include <bits/stdc++.h>
using namespace std;
using namespace chrono;

int n, k;
vector<vector<int>> adj;
vector<int> color;

using BoolMatrix = vector<vector<bool>>;

BoolMatrix bool_matrix_mult(const BoolMatrix &A, const BoolMatrix &B) {
    int N = A.size();
    BoolMatrix C(N, vector<bool>(N, false));
    for (int i = 0; i < N; ++i)
        for (int k = 0; k < N; ++k)
            if (A[i][k])
                for (int j = 0; j < N; ++j)
                    C[i][j]=C[i][j]|B[k][j];
    return C;
}

BoolMatrix find_colorful_paths(int depth, vector<int> nodes, vector<int> cols) {
    if (depth == 1) {
        BoolMatrix A(n, vector<bool>(n, false));
        for (int u : nodes)
            A[u][u] = true;
        return A;
    }

    int half = cols.size() / 2;
    vector<vector<int>> color_partitions;
    for (int mask = 1; mask < (1 << cols.size()); ++mask) {
        if (__builtin_popcount(mask) != half) continue;
        vector<int> C1, C2;
        for (int i = 0; i < (int)cols.size(); ++i) {
            if (mask & (1 << i)) C1.push_back(cols[i]);
            else C2.push_back(cols[i]);
        }
        color_partitions.push_back(C1);
    }

    BoolMatrix result(n, vector<bool>(n, false));

    for (auto &C1 : color_partitions) {
        vector<int> C2;
        set_difference(cols.begin(), cols.end(), C1.begin(), C1.end(), back_inserter(C2));

        unordered_set<int> C1_set(C1.begin(), C1.end()), C2_set(C2.begin(), C2.end());

        vector<int> V1, V2;
        for (int v = 0; v < n; ++v) {
            if (C1_set.count(color[v])) V1.push_back(v);
            else if (C2_set.count(color[v])) V2.push_back(v);
        }

        auto A1 = find_colorful_paths(depth / 2, V1, C1);
        auto A2 = find_colorful_paths(depth / 2, V2, C2);

        BoolMatrix B(n, vector<bool>(n, false));
        for (int u : V1)
            for (int v : adj[u])
                if (C2_set.count(color[v]))
                    B[u][v] = true;

        auto temp = bool_matrix_mult(A1, B);
        auto combined = bool_matrix_mult(temp, A2);

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                result[i][j]=result[i][j]|combined[i][j];
    }

    return result;
}

int main() {
    srand(time(0));
    int m;
    cin >> n >> m >> k;

    adj.assign(n, {});
    color.resize(n);
    for (int i = 0; i < n; ++i)
        cin >> color[i]; // ∈ [0, k-1]

    for (int i = 0; i < m; ++i) {
        int u, v; cin >> u >> v;
        adj[u].push_back(v);
        adj[v].push_back(u); // or not if directed
    }

    vector<int> all_colors(k);
    iota(all_colors.begin(), all_colors.end(), 0);
    vector<int> all_nodes(n);
    iota(all_nodes.begin(), all_nodes.end(), 0);
    const int MAX_TRIES = 1e3;
    auto start = high_resolution_clock::now();
    for (int tries = 0; tries < MAX_TRIES; ++tries) {
        for (int i = 0; i < n; ++i) color[i]=rand()%k;
        auto paths = find_colorful_paths(k, all_nodes, all_colors);

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (paths[i][j]){
                    auto end = high_resolution_clock::now();
                    double runtime = duration<double>(end - start).count();
                    cout << "Found: 1\n";
                    cout << "Iterations: " << tries + 1 << "\n";
                    cout << "Runtime(s): " << fixed << setprecision(6) << runtime << "\n";
                    return 0;
                }
    }
    auto end = high_resolution_clock::now();
    double runtime = duration<double>(end - start).count();
    cout << "Found: 0\n";
    cout << "PathLength: 0\n";
    cout << "Iterations: " << MAX_TRIES << "\n";
    cout << "Runtime(s): " << fixed << setprecision(6) << runtime << "\n";
    return 0;
}
