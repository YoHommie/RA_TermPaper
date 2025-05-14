#include <iostream>
#include <vector>
#include <random>
#include <cstdint>
#include <cmath>
using namespace std;

// Single-step LFSR transition (Fibonacci/Galois form)
static inline uint32_t lfsr_next(uint32_t state, uint32_t taps) {
    uint32_t lsb = state & 1u;
    state >>= 1;
    if (lsb) state ^= taps;
    return state;
}

// LFSR-based hash family
class LFSRHashFamily {
    uint32_t taps;           // tap mask for primitive polynomial
    int m;                   // register width (bits)
    int funcs;               // number of hash functions
    vector<uint32_t> seed;   // initial state per function

public:
    // taps_mask: e.g. 0x80200003 for x^32+x^22+x^2+x+1
    LFSRHashFamily(uint32_t taps_mask, int m_, int num_funcs)
      : taps(taps_mask), m(m_), funcs(num_funcs), seed(num_funcs)
    {
        mt19937_64 rnd(random_device{}());
        uniform_int_distribution<uint32_t> dist(1, (m==32 ? 0xFFFFFFFFu : (1u<<m)-1));
        for(int i = 0; i < funcs; i++)
            seed[i] = dist(rnd);
    }

    // Generate an l‐bit pseudorandom value by stepping the LFSR 'steps' times,
    // then extracting l output bits.
    uint32_t hash(int func_idx, uint32_t steps, int l) const {
        uint32_t s = seed[func_idx % funcs];
        // advance by 'steps'
        for(uint32_t i = 0; i < steps; i++)
            s = lfsr_next(s, taps);
        // extract l bits
        uint32_t h = 0;
        for(int i = 0; i < l; i++) {
            h = (h << 1) | (s & 1u);
            s = lfsr_next(s, taps);
        }
        return h;
    }

    int size() const { return funcs; }
};

// Two-stage LFSR colorer for k‐perfect hashing
class TwoStageLFSRColorer {
    int k, n;
    LFSRHashFamily H1, H2;
    int l1, l2;
public:
    TwoStageLFSRColorer(int k_, int n_)
      : k(k_), n(n_),
        H1(0x80200003u, 32, 2*k_*int(ceil(log2(n_)))),
        H2(0x80200003u, 32, k_*k_)
    {
        l1 = int(ceil(log2(k_*k_)));
        l2 = int(ceil(log2(k_)));
    }

    // Color vertex v (1..n) under function indices (i1, i2), returns in [1..k]
    int color(int v, int i1, int i2) const {
        uint32_t h1 = H1.hash(i1, v, l1) % (k*k);
        uint32_t h2 = H2.hash(i2, h1, l2) % k;
        return int(h2) + 1;
    }

    int numStage1() const { return H1.size(); }
    int numStage2() const { return H2.size(); }
};

// Derandomized color‐coding for finding a colorful path of length k
class ColorCoding {
    const vector<vector<int>>& adj;
    int k, n;
    TwoStageLFSRColorer colorer;

public:
    ColorCoding(const vector<vector<int>>& graph, int k_)
      : adj(graph), k(k_), n(int(graph.size())-1),
        colorer(k_, n) {}

    // Returns true if there's a path on k distinct‐colored vertices
    bool findColorfulPath() {
        int full = (1<<k) - 1;
        // try all two‐stage functions
        for(int i1 = 0; i1 < colorer.numStage1(); i1++) {
            for(int i2 = 0; i2 < colorer.numStage2(); i2++) {
                // assign colors
                vector<int> col(n+1);
                for(int v = 1; v <= n; v++)
                    col[v] = colorer.color(v, i1, i2);
                // dp[v][mask]: path ending at v using color‐set=mask
                vector<vector<char>> dp(n+1, vector<char>(1<<k, 0));
                // base: single‐vertex
                for(int v = 1; v <= n; v++)
                    dp[v][1<<(col[v]-1)] = 1;
                // extend
                for(int len = 1; len < k; len++) {
                    for(int v = 1; v <= n; v++) {
                        for(int mask = 0; mask <= full; mask++) {
                            if (!dp[v][mask] || __builtin_popcount(mask)!=len) continue;
                            for(int u: adj[v]) {
                                int b = 1<<(col[u]-1);
                                if (!(mask & b))
                                    dp[u][mask|b] = 1;
                            }
                        }
                    }
                }
                // check full mask
                for(int v = 1; v <= n; v++)
                    if (dp[v][full]) return true;
            }
        }
        return false;
    }
};

int main(){
    // Example graph (1-based, index 0 unused)
    vector<vector<int>> graph = {
        {}, {2,3}, {1,3,4}, {1,2,4}, {2,3,5}, {4}
    };
    int k = 4;  // look for path on 4 vertices
    ColorCoding cc(graph, k);
    bool found = cc.findColorfulPath();
    cout << "Colorful path of length " << k << (found ? " found\n" : " not found\n");
    return 0;
}
