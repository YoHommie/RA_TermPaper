#include <iostream>
#include <vector>
#include <set>
using namespace std;

const int MAXN = 100;

int n;
vector<vector<int>> adj;
vector<vector<long long>> matPow;

vector<vector<long long>> multiply(const vector<vector<long long>> &A, const vector<vector<long long>> &B) {
    vector<vector<long long>> res(n, vector<long long>(n, 0));
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            for(int k=0; k<n; ++k)
                res[i][j] += A[i][k] * B[k][j];
    return res;
}

vector<vector<long long>> matrixExpo(vector<vector<long long>> base, int exp) {
    vector<vector<long long>> res(n, vector<long long>(n, 0));
    for(int i=0; i<n; ++i) res[i][i] = 1; // identity matrix
    while(exp) {
        if(exp & 1)
            res = multiply(res, base);
        base = multiply(base, base);
        exp >>= 1;
    }
    return res;
}

void dfs(int u, int depth, int k, vector<bool> &visited, vector<int> &path, set<vector<int>> &paths) {
    if(depth == k) {
        paths.insert(path);
        return;
    }
    for(int v=0; v<n; ++v) {
        if(adj[u][v] && !visited[v]) {
            visited[v] = true;
            path.push_back(v);
            dfs(v, depth + 1, k, visited, path, paths);
            path.pop_back();
            visited[v] = false;
        }
    }
}

int main() {
    int m, k;
    cout << "Enter number of nodes and edges: ";
    cin >> n >> m;
    adj.assign(n, vector<int>(n, 0));
    
    cout << "Enter " << m << " undirected edges (u v):" << endl;
    for(int i=0; i<m; ++i) {
        int u, v;
        cin >> u >> v;
        adj[u][v] = adj[v][u] = 1;
    }

    cout << "Enter path length k: ";
    cin >> k;

    auto adjLL = vector<vector<long long>>(adj.begin(), adj.end());
    matPow = matrixExpo(adjLL, k);
    long long totalWalks = 0;
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            totalWalks += matPow[i][j];

    cout << "\nTotal number of walks of length " << k << ": " << totalWalks << endl;

    set<vector<int>> allSimplePaths;
    for(int i=0; i<n; ++i) {
        vector<bool> visited(n, false);
        vector<int> path = {i};
        visited[i] = true;
        dfs(i, 0, k, visited, path, allSimplePaths);
    }

    cout << "\nSimple paths of length " << k << ":\n";
    for(const auto &p : allSimplePaths) {
        for(int node : p)
            cout << node << " ";
        cout << "\n";
    }

    cout << "\nTotal simple paths of length " << k << ": " << allSimplePaths.size() << endl;

    return 0;
}
