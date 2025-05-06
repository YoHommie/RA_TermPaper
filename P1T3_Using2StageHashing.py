import random
import math

# ------------------------- Polynomial Hash -------------------------


class PolynomialHash:
    def __init__(self, k, p):
        self.prime = p
        self.coefficients = [random.randint(0, p - 1) for _ in range(k)]
        # Ensure the polynomial is not constant (for k > 1)
        while k > 1 and self.coefficients[k - 1] == 0:
            self.coefficients[k - 1] = random.randint(0, p - 1)

    def evaluate(self, x):
        result = 0
        for coeff in reversed(self.coefficients):
            result = (result * x + coeff) % self.prime
        return result

    def __call__(self, x):
        return self.evaluate(x)

# ------------------------- K-Perfect Hash Family -------------------------


class KPerfectHashFamily:
    def __init__(self, k, n, num_functions):
        self.k = k
        self.n = n
        self.prime = self.next_prime(max(n, k))
        self.hash_functions = [PolynomialHash(
            k, self.prime) for _ in range(num_functions)]

    @staticmethod
    def is_prime(num):
        if num <= 1:
            return False
        for i in range(2, int(math.sqrt(num)) + 1):
            if num % i == 0:
                return False
        return True

    @classmethod
    def next_prime(cls, start):
        while True:
            start += 1
            if cls.is_prime(start):
                return start

    def get_hash_function(self, index):
        return self.hash_functions[index % len(self.hash_functions)]

    def num_functions(self):
        return len(self.hash_functions)

# ------------------------- Two-Stage k-Perfect Hash Family -------------------------


class TwoStageKPerfectHash:
    def __init__(self, k, n):
        m1 = 2 * k * math.ceil(math.log2(n))
        m2 = k * k
        self.k = k
        self.stage1 = KPerfectHashFamily(k, n, int(m1))
        self.stage2 = KPerfectHashFamily(k, k * k, m2)

    def hash(self, value, func1_idx, func2_idx):
        intermediate = self.stage1.get_hash_function(
            func1_idx)(value) % (self.k * self.k)
        return self.stage2.get_hash_function(func2_idx)(intermediate) % self.k

    def num_stage1(self):
        return self.stage1.num_functions()

    def num_stage2(self):
        return self.stage2.num_functions()

# ------------------------- Derandomized Color Coding with Two-Stage Hashing -------------------------


class DerandomizedColorCoding:
    def __init__(self, graph, target_size):
        self.adj = graph
        self.k = target_size
        self.two_stage_hash = TwoStageKPerfectHash(target_size, len(graph) - 1)

    def contains_path(self):
        for i in range(self.two_stage_hash.num_stage1()):
            for j in range(self.two_stage_hash.num_stage2()):
                colors = [0] * len(self.adj)
                for v in range(1, len(self.adj)):
                    colors[v] = self.two_stage_hash.hash(v, i, j) + 1
                if self.has_colorful_path(colors):
                    return True
        return False

    def has_colorful_path(self, colors):
        n = len(self.adj)
        dp = [[False] * (1 << self.k) for _ in range(n)]
        # Base case: single-vertex paths
        for v in range(1, n):
            mask = 1 << (colors[v] - 1)
            dp[v][mask] = True

        # Build up paths of increasing length
        for length in range(1, self.k):
            for v in range(1, n):
                for mask in range(1 << self.k):
                    if not dp[v][mask]:
                        continue
                    if bin(mask).count('1') != length:
                        continue
                    for u in self.adj[v]:
                        color_bit = 1 << (colors[u] - 1)
                        if (mask & color_bit) == 0:
                            dp[u][mask | color_bit] = True

        # If any vertex has a path that has all k colors, return True
        for v in range(1, n):
            if dp[v][(1 << self.k) - 1]:
                return True
        return False
