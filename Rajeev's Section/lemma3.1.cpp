#include <bits/stdc++.h>
using namespace std;
using namespace chrono;

int n, m, k;
vector<vector<int>> adj;
vector<int> color; // color[v] ∈ {1,...,k}

bool find_colorful_path_from_source(int s) {
    int n = adj.size();
    int total_masks = 1 << k;

    vector<vector<bool>> dp(n, vector<bool>(total_masks, false));
for (int u = 0; u < n; ++u) {
    int mask = 1 << color[u];
    dp[u][mask] = true;
}

for (int len = 1; len < k; ++len) {
    vector<vector<bool>> next_dp(n, vector<bool>(total_masks, false));
    
    for (int u = 0; u < n; ++u) {
        for (int mask = 0; mask < total_masks; ++mask) {
            if (dp[u][mask]) {
                for (int v : adj[u]) {
                    int cbit = 1 << color[v];
                    if (!(mask & cbit)) {
                        int new_mask = mask | cbit;
                        next_dp[v][new_mask] = true;
                    }
                }
            }
        }
    }
    
    dp = next_dp;
}

bool colorful_path_exists = false;
int target_mask = (1 << k) - 1;  // All k colors used

for (int u = 0; u < n; ++u) {
    if (dp[u][target_mask]) return true;
}
    return false;
}

bool find_colorful_path_from_source1(int s) {
    int n = adj.size();
    int total_masks = 1 << k;

    vector<set<set<int> > > dp(n);
    for (int u = 0; u < n; ++u) {
        set<int> s1;
        s1.insert(color[u]);
        dp[u].insert(s1);
    }

    for (int len = 1; len < k; ++len) {
    vector<set<set<int> > > next_dp(n);
    
    for (int u = 0; u < n; ++u) for(auto &i:dp[u]){
                for (int v : adj[u]) {
                    if (i.count(color[v])==0) {
                        set<int> s1=i;
                        s1.insert(color[v]);
                        next_dp[v].insert(s1);
                    }
                }
        }
        dp = next_dp;
    }
    

for (int u = 0; u < n; ++u)  for(auto &i:dp[u]) if(i.size()==k) return true;
    return false;
}

int main() {
    srand(time(0));
    cin >> n >> m >> k;
    adj.assign(n, {});
    color.resize(n);
    for (int i = 0; i < m; ++i) {
        int u, v; cin >> u >> v;
        adj[u].push_back(v);
        adj[v].push_back(u); // or not if directed
    }
    const int MAX_TRIES = 1e3;
    auto start = high_resolution_clock::now();
    for (int tries = 0; tries < MAX_TRIES; ++tries) {
        for (int i = 0; i < n; ++i) color[i]=(rand()%k); 
        for (int s = 0; s < n; ++s) {
            if (find_colorful_path_from_source1(s)) {
                auto end = high_resolution_clock::now();
                double runtime = duration<double>(end - start).count();
    
                cout << "Found: 1\n";
                cout << "Iterations: " << tries + 1 << "\n";
                cout << "Runtime(s): " << fixed << setprecision(6) << runtime << "\n";
                return 0;
            }
        }
    }
    auto end = high_resolution_clock::now();
    double runtime = duration<double>(end - start).count();

    cout << "Found: 0\n";
    cout << "PathLength: 0\n";
    cout << "Iterations: " << MAX_TRIES << "\n";
    cout << "Runtime(s): " << fixed << setprecision(6) << runtime << "\n";
}
