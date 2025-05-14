#include <bits/stdc++.h>
using namespace std;
using namespace chrono;

int n, m, k;
vector<vector<int>> undirected_adj;
vector<vector<int>> dag;

void orientGraphAcyclically() {
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), std::mt19937(std::random_device{}()));

    vector<int> order(n);
    for (int i = 0; i < n; ++i) order[perm[i]] = i;

    dag.assign(n, {});
    for (int u = 0; u < n; ++u)
        for (int v : undirected_adj[u])
            if (order[u] < order[v])
                dag[u].push_back(v);
}

bool findLongestPathOfAtLeastK() {
    vector<int> indeg(n, 0), topo, dp(n,1);
    for (int u = 0; u < n; ++u)
        for (int v : dag[u]) indeg[v]++;

    queue<int> q;
    for (int i = 0; i < n; ++i)
        if (indeg[i] == 0) q.push(i);

    while (!q.empty()) {
        int u = q.front(); q.pop();
        topo.push_back(u);
        for (int v : dag[u])  if (--indeg[v] == 0) q.push(v);
    }
    reverse(topo.begin(),topo.end());
    for (int v:topo) for (int u:dag[v]) dp[v]=max(dp[v],dp[u]+1);
    cout<<*max_element(dp.begin(),dp.end())<<' ';
    return *max_element(dp.begin(), dp.end()) >= k;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cin >> n >> m >> k;
    undirected_adj.assign(n, {});
    for (int i = 0; i < m; ++i) {
        int u, v;
        cin >> u >> v;
        undirected_adj[u].push_back(v);
        undirected_adj[v].push_back(u);
    }

    const int MAX_TRIES = 1e5;
    auto start = high_resolution_clock::now();
    for (int tries = 0; tries < MAX_TRIES; ++tries) {
        orientGraphAcyclically();
        if (findLongestPathOfAtLeastK()) {
            auto end = high_resolution_clock::now();
            double runtime = duration<double>(end - start).count();

            cout << "Found: 1\n";
            cout << "PathLength: " << k << "\n";
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
