"""Independent checks of the group-algebra inputs used from the earlier pass:
 (A) the six class expansions (I1)-(I6) of Prop. 2.3 in Z[S_n], 2 <= n <= 7, by explicit multiplication
     of permutations (composition right-to-left); JM elements J_k = sum_{i<k} (i k).
 (B) Lemma 5.1 (fixed-point cost)  m >= (|V| - c_V(g)) + phi  for all sequences of m <= 4 transpositions on 5 points."""
import itertools
from collections import Counter, defaultdict
from math import comb
def compose(p, q):  # (p q)(x) = p(q(x))
    return tuple(p[q[x]] for x in range(len(p)))
def transp(n, i, j):
    p = list(range(n)); p[i], p[j] = p[j], p[i]; return tuple(p)
def ctype(p):
    n = len(p); seen = [False] * n; ct = []
    for i in range(n):
        if not seen[i]:
            l = 0; j = i
            while not seen[j]: seen[j] = True; j = p[j]; l += 1
            ct.append(l)
    return tuple(sorted([c for c in ct if c > 1], reverse=True))
def mul(A, B):  # dict perm->coeff
    R = defaultdict(int)
    for p, a in A.items():
        for q, b in B.items(): R[compose(p, q)] += a * b
    return R
def add(*terms):
    R = defaultdict(int)
    for coef, A in terms:
        for p, a in A.items(): R[p] += coef * a
    return R
def as_classes(A):
    """check A is central (coeff depends on cycle type) and return dict type->coeff"""
    out = {}
    for p, a in A.items():
        if a == 0: continue
        t = ctype(p)
        if t in out: assert out[t] == a, ("not central", t)
        else: out[t] = a
    return out
def K(n, rho):
    """class sum of type rho (no 1's) in S_n"""
    R = {}
    for p in itertools.permutations(range(n)):
        if ctype(p) == tuple(rho): R[p] = 1
    return R
for n in range(2, 8):
    ident = tuple(range(n))
    J = [add(*[(1, {transp(n, i, k): 1}) for i in range(k)]) if k > 0 else {ident: 0} for k in range(n)]  # J[k] = sum_{i<k} (i k), 0-indexed
    def p(m):
        R = defaultdict(int)
        for k in range(1, n):
            P = J[k]
            for _ in range(m - 1): P = mul(P, J[k])
            for q, a in P.items(): R[q] += a
        return R
    P1, P2, P3, P4 = p(1), p(2), p(3), p(4)
    K2, K3, K22, K4, K5, K33 = (K(n, r) for r in ((2,), (3,), (2, 2), (4,), (5,), (3, 3)))
    one = {ident: 1}
    def eq(A, B):
        keys = set(A) | set(B)
        return all(A.get(k, 0) == B.get(k, 0) for k in keys)
    assert eq(P1, K2)
    assert eq(P2, add((comb(n, 2), one), (1, K3)))
    assert eq(mul(K2, K2), add((comb(n, 2), one), (3, K3), (2, K22)))
    assert eq(P3, add((1, K4), (2 * n - 3, K2)))
    assert eq(P4, add((1, K5), (3 * n - 4, K3), (4, K22), (n * (n - 1) * (4 * n - 5) // 6, one)))
    assert eq(mul(K3, K3), add((2 * comb(n, 3), one), (3 * n - 8, K3), (8, K22), (5, K5), (2, K33)))
    print(f"(A) n={n}: (I1)-(I6) hold in Z[S_n]")
# (B)
N = 5; cnt = 0
ts = [(i, j) for i in range(N) for j in range(i + 1, N)]
for m in range(0, 5):
    for seq in itertools.product(ts, repeat=m):
        g = tuple(range(N))
        for (i, j) in seq: g = compose(g, transp(N, i, j))
        V = set(x for t in seq for x in t)
        # cycles of g on V (V is g-stable)
        seen = set(); c = 0
        for x in V:
            if x not in seen:
                c += 1; y = x
                while y not in seen: seen.add(y); y = g[y]
        phi = sum(1 for x in V if g[x] == x)
        assert m >= (len(V) - c) + phi, seq
        cnt += 1
print(f"(B) Lemma 5.1 verified on all {cnt} sequences of m <= 4 transpositions on 5 points")
