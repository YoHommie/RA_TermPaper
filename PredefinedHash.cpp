#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <functional>
#include <algorithm>
#include <chrono>

using namespace std;

// ------------------------- StdHashFamily -------------------------
// Simulates a family of hash functions using std::hash and a seed
class StdHashFamily {
    int numFuncs;
    vector<size_t> seeds; // Different seeds for each "function"
public:
    StdHashFamily(int numFunctions) : numFuncs(numFunctions) {
        // Use simple deterministic seeds (could also use random values)
        for (int i = 0; i < numFuncs; ++i)
            seeds.push_back(i * 2654435761u); // Knuth's multiplicative hash
    }
    // Returns a hash value for x using the i-th hash "function"
    size_t hash(int x, int i) const {
        return std::hash<int>()(x ^ seeds[i % numFuncs]);
    }
    int num_functions() const { return numFuncs; }
};

// ------------------------- Two-Stage k-Perfect Hash Family -------------------------
// Constructs a two-stage perfect hash family by composing two families.
class TwoStageKPerfectHash
{
private:
    // Stage1: maps from [1, n] to [0, k^2 - 1]
    StdHashFamily stage1;
    // Stage2: maps from [0, k^2 - 1] to [0, k - 1]
    StdHashFamily stage2;
    int k; // number of colors

public:
    // The two-stage family is built with parameters:
    //   - m1 = 2 * k * ceil(log2(n_val)) functions in stage1,
    //   - m2 = k^2 functions in stage2.
    TwoStageKPerfectHash(int k_val, int n_val)
        : k(k_val),
          stage1(2 * k_val * static_cast<int>(ceil(log2(n_val)))),
          stage2(k_val * k_val)
    {
    }

    // Compose the two stages.
    // First, hash "value" with stage1 (mod k^2), then hash the intermediate value with stage2 (mod k).
    int hash(int value, int func1_idx, int func2_idx) const
    {
        int intermediate = static_cast<int>(stage1.hash(value, func1_idx) % (k * k));
        return static_cast<int>(stage2.hash(intermediate, func2_idx) % k);
    }

    int num_stage1() const { return stage1.num_functions(); }
    int num_stage2() const { return stage2.num_functions(); }
};

// ------------------------- Derandomized Color Coding with Two-Stage Hashing -------------------------
class DerandomizedColorCoding
{
private:
    vector<vector<int>> adj;           // Graph as an adjacency list (1-indexed; index 0 is unused)
    int k;                             // Target path size (number of vertices in the desired colorful path)
    TwoStageKPerfectHash twoStageHash; // Two-stage hash family for color assignment

public:
    // "graph" is assumed to have vertex indices starting at 1.
    DerandomizedColorCoding(const vector<vector<int>> &graph, int target_size)
        : adj(graph), k(target_size), twoStageHash(target_size, static_cast<int>(graph.size() - 1))
    {
    }

    // Try every combination of the two-stage hash functions.
    // For each color assignment, use dynamic programming to check for a colorful path.
    bool contains_path()
    {
        using namespace std::chrono;
        auto start = high_resolution_clock::now();

        for (int i = 0; i < twoStageHash.num_stage1(); i++)
        {
            for (int j = 0; j < twoStageHash.num_stage2(); j++)
            {
                // Build a color assignment for vertices 1 ... n.
                vector<int> colors(adj.size(), 0); // index 0 unused
                for (int v = 1; v < adj.size(); v++)
                {
                    // Add 1 so that colors are in the range 1 ... k.
                    colors[v] = twoStageHash.hash(v, i, j) + 1;
                }
                if (has_colorful_path(colors))
                {
                    auto end = high_resolution_clock::now();
                    auto duration = duration_cast<microseconds>(end - start);
                    cout << "[DerandomizedColorCoding::contains_path] Time taken: " << duration.count() << " microseconds" << endl;
                    return true;
                }
            }
        }
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(end - start);
        cout << "[DerandomizedColorCoding::contains_path] Time taken: " << duration.count() << " microseconds" << endl;
        return false;
    }

private:
    // Dynamic programming routine to check for a colorful path given vertex colors.
    // dp[v][mask] is true if there is a path ending at vertex v with colors corresponding to bitmask "mask".
    bool has_colorful_path(const vector<int> &colors)
    {
        using namespace std::chrono;
        auto start = high_resolution_clock::now();

        int n = static_cast<int>(adj.size());
        vector<vector<bool>> dp(n, vector<bool>(1 << k, false));

        // Base case: single-vertex paths.
        for (int v = 1; v < n; v++)
        {
            int mask = 1 << (colors[v] - 1);
            dp[v][mask] = true;
        }

        // Build up paths of increasing length.
        for (int len = 1; len < k; len++)
        {
            for (int v = 1; v < n; v++)
            {
                for (int mask = 0; mask < (1 << k); mask++)
                {
                    if (!dp[v][mask])
                        continue;
                    // Proceed only if the mask has exactly 'len' colors.
                    if (__builtin_popcount(mask) != len)
                        continue;

                    // Try extending the path from vertex v.
                    for (int u : adj[v])
                    {
                        int color_bit = 1 << (colors[u] - 1);
                        if ((mask & color_bit) == 0)
                        {
                            dp[u][mask | color_bit] = true;
                        }
                    }
                }
            }
        }

        // If any vertex has a path that has all k colors, return true.
        for (int v = 1; v < n; v++)
        {
            if (dp[v][(1 << k) - 1]) {
                auto end = high_resolution_clock::now();
                auto duration = duration_cast<microseconds>(end - start);
                cout << "[DerandomizedColorCoding::has_colorful_path] Time taken: " << duration.count() << " microseconds" << endl;
                return true;
            }
        }
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(end - start);
        cout << "[DerandomizedColorCoding::has_colorful_path] Time taken: " << duration.count() << " microseconds" << endl;
        return false;
    }
};

// Class to generate various test graphs.
class WrostCaseGraphForDerandomizedColorCoding
{
public:
    // Creates a complete graph with n vertices.
    // In a complete graph, every vertex is connected to every other vertex.
    // The graph uses 1-indexing; vertex 0 is unused.
    vector<vector<int>> createCompleteGraph(int n)
    {
        // Initialize graph with n+1 vertices.
        vector<vector<int>> graph(n + 1);

        // Loop over vertices 1 to n.
        for (int u = 1; u <= n; u++)
        {
            // For each vertex u, add an edge to every other vertex v.
            for (int v = 1; v <= n; v++)
            {
                // Avoid self-loops.
                if (u != v)
                {
                    graph[u].push_back(v);
                }
            }
        }
        return graph;
    }

    // Creates a layered graph with a specified number of layers and vertices per layer.
    // Each vertex in one layer is connected to all vertices in the next layer.
    // The graph is 1-indexed.
    vector<vector<int>> createLayeredGraph(int layers, int verticesPerLayer)
    {
        // Calculate total number of vertices.
        int n = layers * verticesPerLayer;
        // Initialize graph with n+1 vertices.
        vector<vector<int>> graph(n + 1);

        // Loop over layers except the last since it doesn't connect to a next layer.
        for (int layer = 0; layer < layers - 1; layer++)
        {
            // Calculate the starting index for the current and next layers.
            int startCurrent = layer * verticesPerLayer + 1;
            int startNext = (layer + 1) * verticesPerLayer + 1;

            // For each vertex in the current layer.
            for (int i = 0; i < verticesPerLayer; i++)
            {
                int u = startCurrent + i;
                // Connect with every vertex in the next layer.
                for (int j = 0; j < verticesPerLayer; j++)
                {
                    int v = startNext + j;
                    graph[u].push_back(v);
                }
            }
        }
        return graph;
    }

    // Creates an augmented chain graph.
    // Starts with a simple chain, i.e., a path from vertex 1 to vertex n.
    // Then adds extra edges (jumps) from each vertex to subsequent vertices within a specified jump range.
    // This augmentation can potentially allow skipping intermediate vertices.
    vector<vector<int>> createAugmentedChain(int n, int extraJump)
    {
        // Initialize graph with n+1 vertices.
        vector<vector<int>> graph(n + 1);

        // Build the basic chain: edge from vertex i to i+1.
        for (int i = 1; i < n; i++)
        {
            // Add the chain edge.
            graph[i].push_back(i + 1);
            // Add extra jumps from vertex i to vertices i+jump if within bounds.
            for (int jump = 2; jump <= extraJump; jump++)
            {
                if (i + jump <= n)
                {
                    graph[i].push_back(i + jump);
                }
            }
        }
        return graph;
    }

    // Generates an undirected Erdős–Rényi G(n, p) graph
    vector<vector<int>> erdos_renyi_graph(int n, double p, unsigned int seed = time(nullptr))
    {
        mt19937 gen(seed);
        uniform_real_distribution<double> dist(0.0, 1.0);
        vector<vector<int>> graph(n + 1); // 1-indexed
        // Loop over all pairs of vertices (u, v) to create edges with probability p.

        for (int u = 1; u <= n; ++u)
        {
            for (int v = u + 1; v <= n; ++v)
            {
                if (dist(gen) < p)
                {
                    graph[u].push_back(v);
                    graph[v].push_back(u);
                }
            }
        }
        return graph;
    }
};
// ------------------------- Test Graph Class -------------------------

// Utility to count edges in an adjacency list
int count_edges(const vector<vector<int>> &graph)
{
    int edges = 0;
    for (size_t i = 1; i < graph.size(); ++i)
    {
        edges += graph[i].size();
    }
    return edges;
}

// Utility to benchmark and print performance metrics
void benchmark_graph(const string &name, const vector<vector<int>> &graph, int path_length)
{
    int V = static_cast<int>(graph.size()) - 1;
    int E = count_edges(graph);
    DerandomizedColorCoding dcc(graph, path_length);

    auto start = chrono::high_resolution_clock::now();
    bool found = dcc.contains_path();
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

    cout << name << ":\n";
    cout << "  Vertices: " << V << ", Edges: " << E << ", Path length: " << path_length << endl;
    cout << "  Contains colorful path: " << boolalpha << found << endl;
    cout << "  Time taken: " << duration.count() << " microseconds\n"
         << endl;
}

int main()
{

    int skip = 10;
    for (int number_of_vertices = 10; number_of_vertices <= 500; number_of_vertices += skip)
    {
        // int target_length = int(log2(number_of_vertices)); // 4 times the log base 2 of the number of vertices
        int target_length = 8; // 4 times the log base 2 of the number of vertices

        WrostCaseGraphForDerandomizedColorCoding testGraphGen;
        vector<vector<int>> Graph = testGraphGen.createCompleteGraph(number_of_vertices);
        // vector<vector<int>> Graph = testGraphGen.createLayeredGraph(2, number_of_vertices);
        // vector<vector<int>> Graph = testGraphGen.createAugmentedChain(number_of_vertices, 3);
        // vector<vector<int>> Graph = testGraphGen.erdos_renyi_graph(number_of_vertices, 0.5);
        if (number_of_vertices < 100)
        {
            skip = 10;
        }
        else if (number_of_vertices < 500)
        {
            skip = 50;
        }
        else
        {
            skip = 100;
        }

        cout << "---------------------------------------------------" << endl;
        benchmark_graph("Complete Graph", Graph, target_length);
        // benchmark_graph("Layered Graph", Graph, target_length);
        // benchmark_graph("Augmented Chain Graph", Graph, target_length);
        // benchmark_graph("Erdos-Renyi Graph", Graph, target_length);
    }

    return 0;
}

