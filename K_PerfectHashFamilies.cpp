#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <cmath>
#include <unordered_set>

using namespace std;

// Class for polynomial-based hash functions in GF(p)
class PolynomialHash {
private:
    vector<int> coefficients; // a_0, a_1, ..., a_(k-1)
    int prime; // Prime modulus for the finite field

public:
    // Constructor to generate a random polynomial of degree (k-1)
    PolynomialHash(int k, int p) : prime(p) {
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<int> dist(0, p - 1);
        
        // Generate k random coefficients
        coefficients.resize(k);
        for (int i = 0; i < k; i++) {
            coefficients[i] = dist(gen);
        }
        
        // Ensure the polynomial is not constant (degree at least 1)
        while (k > 1 && coefficients[k-1] == 0) {
            coefficients[k-1] = dist(gen);
        }
    }
    
    // Evaluate the polynomial at point x
    int evaluate(int x) const {
        long long result = 0;
        // Horner's method for polynomial evaluation
        for (int i = coefficients.size() - 1; i >= 0; i--) {
            result = (result * x + coefficients[i]) % prime;
        }
        return static_cast<int>(result);
    }
    
    // Hash function operator
    int operator()(int x) const {
        return evaluate(x);
    }
};

// Class for k-perfect hash family
class KPerfectHashFamily {
private:
    int k; // k-wise independence
    int n; // Universe size
    int prime; // Prime modulus for finite field
    vector<PolynomialHash> hashFunctions;

public:
    // Constructor to generate a k-perfect hash family
    KPerfectHashFamily(int k_val, int n_val, int numFunctions) : k(k_val), n(n_val) {
        // Choose a prime larger than n (for simplicity we use the next power of 2)
        prime = nextPrime(max(n, k*k));
        
        // Generate hash functions
        for (int i = 0; i < numFunctions; i++) {
            hashFunctions.emplace_back(k, prime);
        }
    }
    
    // Get a hash function from the family
    const PolynomialHash& getHashFunction(int index) const {
        return hashFunctions[index % hashFunctions.size()];
    }
    
    // Find a hash function that's perfect for a given set
    int findPerfectHashForSet(const vector<int>& set) const {
        for (size_t i = 0; i < hashFunctions.size(); i++) {
            const auto& hashFunc = hashFunctions[i];
            
            // Check if this function is perfect for the set
            unordered_set<int> hashValues;
            bool isPerfect = true;
            
            for (int value : set) {
                int hashValue = hashFunc(value);
                if (hashValues.find(hashValue) != hashValues.end()) {
                    isPerfect = false;
                    break;
                }
                hashValues.insert(hashValue);
            }
            
            if (isPerfect) {
                return i;
            }
        }
        
        // No perfect hash found, return -1
        return -1;
    }
    
    // Utility function to find next prime
    int nextPrime(int n) const {
        while (!isPrime(n)) {
            n++;
        }
        return n;
    }
    
    // Utility function to check if a number is prime
    bool isPrime(int n) const {
        if (n <= 1) return false;
        if (n <= 3) return true;
        if (n % 2 == 0 || n % 3 == 0) return false;
        
        for (int i = 5; i * i <= n; i += 6) {
            if (n % i == 0 || n % (i + 2) == 0) {
                return false;
            }
        }
        return true;
    }
    
    // Size of the family
    size_t size() const {
        return hashFunctions.size();
    }
};

// Implementation of two-stage k-perfect hash construction
class TwoStageKPerfectHash {
private:
    // KPerfectHashFamily stage1; // Maps {1,2,...,n} to {1,2,...,k²}
    // KPerfectHashFamily stage2; // Maps {1,2,...,k²} to {1,2,...,k}
    LSFRHashFamily stage1; // Maps {1,2,...,n} to {1,2,...,k²}
    LSFRHashFamily stage2; // Maps {1,2,...,k²} to {1,2,...,k}
    int k;

public:
    TwoStageKPerfectHash(int k_val, int n_val) : 
        k(k_val),
        stage1(k_val, n_val, 2 * k_val * log2(n_val)), // First stage
        stage2(k_val, k_val * k_val, k_val * k_val)         // Second stage
    {}
    
    // Hash function applying both stages
    int hash(int value, int func1_idx, int func2_idx) const {
        // First stage: hash to {0, 1, ..., k² - 1}
        int intermediate = stage1.getHashFunction(func1_idx)(value) % (k * k);
        
        // Second stage: hash to {0, 1, ..., k - 1}
        return stage2.getHashFunction(func2_idx)(intermediate) % k;
    }
    
    // Verify if a subset is hashed without collisions
    bool isKPerfectForSubset(const vector<int>& subset) const {
        if (subset.size() > k) return false; // Can't guarantee k+1 or more elements
        
        // Check if any combination of hash functions works
        for (size_t i = 0; i < stage1.num_functions(); i++) {
            for (size_t j = 0; j < stage2.num_functions(); j++) {
                unordered_set<int> hashValues;
                bool hasCollision = false;
                
                for (int value : subset) {
                    int hashValue = hash(value, i, j);
                    if (hashValues.find(hashValue) != hashValues.end()) {
                        hasCollision = true;
                        break;
                    }
                    hashValues.insert(hashValue);
                }
                
                if (!hasCollision) {
                    return true;
                }
            }
        }
        
        return false;
    }
};

// Demo function to test the implementation
void demoKPerfectHashFamily() {
    int k = 3;  // k-wise independence
    int n = 100; // Universe size
    
    cout << "Generating a " << k << "-perfect hash family for universe size " << n << endl;
    
    TwoStageKPerfectHash hashFamily(k, n);
    
    // Test with some random subsets
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(1, n);
    
    // Generate 5 random subsets of size k
    for (int t = 0; t < 5; t++) {
        vector<int> subset;
        for (int i = 0; i < k; i++) {
            subset.push_back(dist(gen));
        }
        
        cout << "Testing subset: ";
        for (int val : subset) {
            cout << val << " ";
        }
        cout << endl;
        
        bool isPerfect = hashFamily.isKPerfectForSubset(subset);
        cout << "Is k-perfect for this subset: " << (isPerfect ? "Yes" : "No") << endl;
    }
}

// Alternative implementation using Linear Feedback Shift Register (LFSR)
// class LFSRHashFamily {
// private:
//     int m; // Parameter for feedback rule length
//     vector<vector<bool>> feedbackRules; // Array of feedback rules
//     vector<vector<bool>> startSequences; // Array of start sequences

// public:
//     LFSRHashFamily(int m_val, int numFunctions) : m(m_val) {
//         random_device rd;
//         mt19937 gen(rd());
//         uniform_int_distribution<int> dist(0, 1);
        
//         // Generate feedback rules and start sequences
//         for (int i = 0; i < numFunctions; i++) {
//             // Generate a non-degenerate feedback rule
//             vector<bool> feedbackRule(m);
//             feedbackRule[0] = true; // f₀ = 1
//             for (int j = 1; j < m; j++) {
//                 feedbackRule[j] = dist(gen) == 1;
//             }
//             feedbackRules.push_back(feedbackRule);
            
//             // Generate a start sequence
//             vector<bool> startSequence(m);
//             for (int j = 0; j < m; j++) {
//                 startSequence[j] = dist(gen) == 1;
//             }
//             startSequences.push_back(startSequence);
//         }
//     }
    
//     // Generate hash value using LFSR
//     int hash(int value, int funcIdx, int n) const {
//         const auto& feedbackRule = feedbackRules[funcIdx % feedbackRules.size()];
//         const auto& startSequence = startSequences[funcIdx % startSequences.size()];
        
//         // Generate the sequence
//         vector<bool> sequence(max(n, m));
//         for (int i = 0; i < m; i++) {
//             sequence[i] = startSequence[i];
//         }
        
//         for (int i = m; i < n; i++) {
//             bool bit = false;
//             for (int j = 0; j < m; j++) {
//                 if (feedbackRule[j]) {
//                     bit ^= sequence[i - m + j];
//                 }
//             }
//             sequence[i] = bit;
//         }
        
//         // Use the value-th bit as the hash
//         return sequence[value % n] ? 1 : 0;
//     }
    
//     // Hash function that combines l bits to get a value in range [0, 2^l - 1]
//     int multibitHash(int value, int funcIdx, int n, int l) const {
//         int result = 0;
//         for (int i = 0; i < l; i++) {
//             result = (result << 1) | hash(value + i, funcIdx, n);
//         }
//         return result;
//     }
// };


// -------------------- LSFR Hash --------------------
class LSFRHash {
    private:
        uint32_t seed;
        uint32_t taps;
        int bitWidth;
    
        uint32_t step(uint32_t value) const {
            uint32_t feedback = 0;
            for (int i = 0; i < bitWidth; ++i) {
                if ((taps >> i) & 1) {
                    feedback ^= (value >> i) & 1;
                }
            }
            return ((value << 1) | feedback) & ((1u << bitWidth) - 1);
        }
    
    public:
        LSFRHash(int bits, uint32_t seed_val, uint32_t tap_mask)
            : seed(seed_val), taps(tap_mask), bitWidth(bits) {}
    
        int operator()(int x) const {
            uint32_t state = seed;
            for (int i = 0; i < x; ++i) {
                state = step(state);
            }
            return static_cast<int>(state);
        }
    };
    
    // -------------------- LSFR Hash Family --------------------
    class LSFRHashFamily {
    private:
        int k, n;
        vector<LSFRHash> hashFunctions;
        vector<uint32_t> tapMasks = {
            0b1000000000000001, 0b1100000000000011, 0b1010000000000101,
            0b1001000000000111, 0b1110000000001001
        };
    
    public:
        LSFRHashFamily(int k_val, int n_val, int numFunctions)
            : k(k_val), n(n_val) {
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<uint32_t> seed_dist(1, (1u << 16) - 1);
    
            for (int i = 0; i < numFunctions; i++) {
                uint32_t seed = seed_dist(gen);
                uint32_t tap = tapMasks[i % tapMasks.size()];
                hashFunctions.emplace_back(16, seed, tap);
            }
        }
    
        const LSFRHash &getHashFunction(int index) const {
            return hashFunctions[index % hashFunctions.size()];
        }
    
        int num_functions() const {
            return hashFunctions.size();
        }
    };
int main() {
    demoKPerfectHashFamily();
    return 0;
}

