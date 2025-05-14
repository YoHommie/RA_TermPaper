// #include <iostream>
// #include <vector>
// #include <unordered_set>
// #include <algorithm>
// #include <cmath>

// using namespace std;

// class KPerfectHashFamily {
//     int n;          // Maximum vertex ID
//     int k;          // Target subgraph size
//     int prime;      // Prime number > n
//     vector<pair<int, int>> hash_params;  // (a, b) pairs for hash functions

// public:
//     KPerfectHashFamily(int max_vertex, int target_size) : n(max_vertex), k(target_size) {
//         prime = next_prime(n);

//         // Generate hash parameters (simplified construction)
//         for(int a = 1; a <= k; a++) {
//             for(int b = 0; b < k; b++) {
//                 hash_params.emplace_back(a, b);
//             }
//         }
//     }

//     // Get color assignment for current hash function
//     vector<int> get_colors(int func_index) const {
//         auto [a, b] = hash_params[func_index];
//         vector<int> colors(n+1);  // vertices are 1-based

//         for(int v = 1; v <= n; v++) {
//             colors[v] = ((a * v + b) % prime) % k + 1;
//         }
//         return colors;
//     }

//     size_t num_functions() const { return hash_params.size(); }

// private:
//     int next_prime(int start) {
//         while(true) {
//             start++;
//             if(is_prime(start)) return start;
//         }
//     }

//     bool is_prime(int num) {
//         if(num <= 1) return false;
//         for(int i = 2; i <= sqrt(num); i++)
//             if(num % i == 0) return false;
//         return true;
//     }
// };

// class DerandomizedColorCoding {
//     const vector<vector<int>> adj;  // adjacency list
//     const int k;                    // path length to find (k-1 edges)
//     KPerfectHashFamily hash_family;

// public:
//     DerandomizedColorCoding(const vector<vector<int>>& graph, int target_size)
//         : adj(graph), k(target_size), hash_family(graph.size()-1, target_size) {}

//     bool contains_path() {
//         // Try all hash functions in the family
//         for(size_t fi = 0; fi < hash_family.num_functions(); fi++) {
//             auto colors = hash_family.get_colors(fi);
//             if(has_colorful_path(colors)) return true;
//         }
//         return false;
//     }

// private:
//     bool has_colorful_path(const vector<int>& colors) {
//         // DP table: dp[v][mask] = exists path to v with color mask
//         vector<vector<bool>> dp(adj.size(), vector<bool>(1 << k, false));

//         // Initialize for single-node paths
//         for(int v = 1; v < adj.size(); v++) {
//             int mask = 1 << (colors[v] - 1);
//             dp[v][mask] = true;
//         }

//         // Dynamic programming for paths of increasing length
//         for(int len = 1; len < k; len++) {
//             for(int v = 1; v < adj.size(); v++) {
//                 for(int mask = 0; mask < (1 << k); mask++) {
//                     if(!dp[v][mask] || __builtin_popcount(mask) != len) continue;

//                     for(int u : adj[v]) {
//                         int color_bit = 1 << (colors[u] - 1);
//                         if(!(mask & color_bit)) {
//                             dp[u][mask | color_bit] = true;
//                         }
//                     }
//                 }
//             }
//         }

//         // Check if any vertex has a full color mask
//         for(int v = 1; v < adj.size(); v++) {
//             if(dp[v][(1 << k) - 1]) return true;
//         }
//         return false;
//     }
// };

// int main() {
//     // Example graph (1-based vertices)
//     vector<vector<int>> graph = {
//         {},          // 0 unused
//         {2, 3},      // 1
//         {1, 3, 4},   // 2
//         {1, 2, 4},   // 3
//         {2, 3, 5},   // 4
//         {4}          // 5
//     };

//     int target_length = 4; // Look for path with 4 edges (5 vertices)

//     DerandomizedColorCoding dcc(graph, target_length);
//     cout << "Contains path of length " << target_length-1 << ": "
//          << boolalpha << dcc.contains_path() << endl;

//     return 0;
// }

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>

using namespace std;

// Class representing a family of perfect hash functions for color coding.
// Each hash function maps vertex IDs to colors in a range [1, k].
class KPerfectHashFamily
{
    int n;                              // Maximum vertex ID (number of vertices in the graph)
    int k;                              // Target subgraph size (number of colors)
    int prime;                          // A prime number greater than n, used for mod operations in hash functions
    vector<pair<int, int>> hash_params; // Parameters for hash functions as (a, b) pairs

public:
    // Constructor: Initializes the hash family with max vertex and target subgraph size.
    // It calculates a prime number and generates hash functions based on (a, b) pairs.
    KPerfectHashFamily(int max_vertex, int target_size)
    {
        n = max_vertex;
        k = target_size;
        // Find the next prime number greater than n
        prime = next_prime(n);

        // Generate parameters for the hash functions.
        // For each a in [1, k] and each b in [0, k-1] we create a hash function.
        for (int a = 1; a <= k; a++)
        {
            for (int b = 0; b < k; b++)
            {
                hash_params.push_back(make_pair(a, b));
            }
        }
    }

    // Returns the color assignment for each vertex using the hash function at index 'func_index'.
    // Colors are computed as: ((a * vertex_id + b) % prime) % k + 1.
    vector<int> get_colors(int func_index) const
    {
        // Extract the current hash parameters (a, b)
        int a = hash_params[func_index].first;
        int b = hash_params[func_index].second;
        // Create a color vector for vertices, using 1-based indexing.
        vector<int> colors(n + 1);

        // Compute the color for each vertex from 1 to n.
        for (int v = 1; v <= n; v++)
        {
            colors[v] = ((a * v + b) % prime) % k + 1;
        }
        return colors;
    }

    // Returns the number of hash functions available in the family.
    int num_functions() const
    {
        return hash_params.size();
    }

private:
    // Given a starting number, finds the next prime number greater than the given start.
    int next_prime(int start)
    {
        while (true)
        {
            start++; // Increment the number
            if (is_prime(start))
                return start; // Return if the number is prime
        }
    }

    // Checks if a given number is prime.
    // Returns true if 'num' is a prime number, false otherwise.
    bool is_prime(int num)
    {
        if (num <= 1)
            return false; // Numbers <= 1 are not prime
        // Check divisibility from 2 up to sqrt(num)
        for (int i = 2; i <= sqrt(num); i++)
        {
            if (num % i == 0)
                return false;
        }
        return true;
    }
};

// Class implementing the derandomized color coding algorithm for path detection in graphs.
class DerandomizedColorCoding
{
    vector<vector<int>> adj;        // Adjacency list representation of the graph.
                                    // Indexing is 1-based, index 0 is unused.
    int k;                          // The target size; looking for paths of length k-1 (i.e., k vertices)
    KPerfectHashFamily hash_family; // An instance of the hash family used for color assignments

public:
    // Constructor: Initializes the graph and target size, and sets up the hash family.
    // The graph is assumed to be given as a vector of vectors with 0 index unused.
    DerandomizedColorCoding(const vector<vector<int>> &graph, int target_size)
        : adj(graph), k(target_size), hash_family(graph.size() - 1, target_size) {}

    // Returns true if the graph contains a path that is colorful (all vertices have distinct colors).
    // The method iterates through all hash functions in the family.
    bool contains_path()
    {
        auto start_time = chrono::high_resolution_clock::now();

        // Try each hash function produced by the hash family.
        for (int i = 0; i < hash_family.num_functions(); i++)
        {
            vector<int> colors = hash_family.get_colors(i);
            // If any hash function yields a colorful path, return true.
            if (has_colorful_path(colors))
            {
                auto end_time = chrono::high_resolution_clock::now();
                auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
                cout << "[DerandomizedColorCoding::contains_path] Time taken: " << duration << " microseconds" << endl;
                return true;
            }
        }
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
        cout << "[DerandomizedColorCoding::contains_path] Time taken: " << duration << " microseconds" << endl;
        return false;
    }

private:
    // Determines if there exists a colorful path using the given vertex color assignment.
    // Colorful means that no two vertices on the path share the same color.
    bool has_colorful_path(const vector<int> &colors)
    {
        auto start = chrono::high_resolution_clock::now();
        int n = adj.size(); // Number of vertices (including unused 0 index)
        // dp[v][mask] will be true if there is a path ending at vertex v
        // which uses a set of colors represented by the bitmask 'mask'.
        vector<vector<bool>> dp(n, vector<bool>(1 << k, false));

        // Initialization: For each vertex, set the dp value for the single color containing that vertex.
        for (int v = 1; v < n; v++)
        {
            int mask = 1 << (colors[v] - 1); // Create bitmask for the current vertex color.
            dp[v][mask] = true;              // A single vertex path exists.
        }

        // Build up the solution for paths of increasing length.
        for (int len = 1; len < k; len++)
        {
            // Consider extending paths ending at each vertex 'v'.
            for (int v = 1; v < n; v++)
            {
                // For each possible color set (bitmask).
                for (int mask = 0; mask < (1 << k); mask++)
                {
                    // If there is no path ending at vertex v with this specific set of colors, skip.
                    if (!dp[v][mask])
                        continue;

                    // Count the number of colors (bits set) in the mask.
                    int bits = 0;
                    for (int b = 0; b < k; b++)
                    {
                        if ((mask >> b) & 1)
                            bits++;
                    }
                    // If the path length (number of vertices) does not match len, skip.
                    if (bits != len)
                        continue;

                    // Try to extend the path from vertex v to each of its neighbors.
                    for (int i = 0; i < adj[v].size(); i++)
                    {
                        int u = adj[v][i];                    // Neighbor of vertex v.
                        int color_bit = 1 << (colors[u] - 1); // Bit corresponding to neighbor's color.
                        // Check if the neighbor's color is not already used in the path.
                        if ((mask & color_bit) == 0)
                        {
                            dp[u][mask | color_bit] = true; // Extend the path with the new color.
                        }
                    }
                }
            }
        }

        // Check if any vertex has a path that uses all k colors.
        // This means we have found a colorful path with k vertices.
        bool found = false;
        for (int v = 1; v < n; v++)
        {
            if (dp[v][(1 << k) - 1])
            {
                found = true;
                break;
            }
        }
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
        cout << "[DerandomizedColorCoding::has_colorful_path] Time taken: " << duration << " microseconds" << endl;
        return found;
    }
};

// // Main function to demonstrate the derandomized color coding algorithm.
// int main() {
//     // Define an example graph using an adjacency list.
//     // Index 0 is not used. Vertices are numbered from 1.
//     // Each vector inside the main graph vector represents the neighbors of that vertex.
//     vector<vector<int>> graph = {
//         {},          // 0 unused
//         {2, 3},      // Vertex 1 is connected to vertices 2 and 3.
//         {1, 3, 4},   // Vertex 2 is connected to vertices 1, 3, and 4.
//         {1, 2, 4},   // Vertex 3 is connected to vertices 1, 2, and 4.
//         {2, 3, 5},   // Vertex 4 is connected to vertices 2, 3, and 5.
//         {4}          // Vertex 5 is connected to vertex 4.
//     };

//     int target_length = 4; // Set target length for the path (number of vertices in the path).

//     // Create an instance of the DerandomizedColorCoding class with the graph and target path length.
//     DerandomizedColorCoding dcc(graph, target_length);
//     // Check if the graph contains a colorful path of the desired length and print the result.
//     cout << "Contains path of length " << (target_length - 1) << ": "
//          << boolalpha << dcc.contains_path() << endl;

//     return 0;
// }

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
        int target_length = 8; // 4 times the log base 2 of the number of vertices
        // int target_length = int(log2(number_of_vertices)); // 4 times the log base 2 of the number of vertices

        WrostCaseGraphForDerandomizedColorCoding testGraphGen;
        // vector<vector<int>> Graph = testGraphGen.createCompleteGraph(number_of_vertices);
        // vector<vector<int>> Graph = testGraphGen.createLayeredGraph(2, number_of_vertices);
        // vector<vector<int>> Graph = testGraphGen.createAugmentedChain(number_of_vertices, 3);
        vector<vector<int>> Graph = testGraphGen.erdos_renyi_graph(number_of_vertices, 0.5);
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
