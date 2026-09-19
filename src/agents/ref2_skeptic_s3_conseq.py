"""ref2_skeptic_s3_conseq.py -- referee #2: independent checks of C1 (via Schur polynomial multiplication in 3 variables,
NOT via S3), C2 (own chain counter), C3 (project d_vector / f_skew)."""
import os
import sys, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from young import f_skew, f_hook
from census import d_vector
from collections import Counter

# --- GL_3 characters as dicts weight->mult, via Weyl formula: ch = A_{lam+rho}/A_rho (exact polynomial division by decomposition)
def weyl_num(g):
    """A_g = sum_w sgn(w) x^{w g} as Counter over Z^3 exponents."""
    c = Counter()
    for s in itertools.permutations(range(3)):
        inv = sum(1 for i in range(3) for j in range(i+1,3) if s[i] > s[j])
        c[tuple(g[s[i]] for i in range(3))] += (-1)**inv
    return c
RHO = (2,1,0)
def decompose(numer):
    """Given P*A_rho (a W-alternating Laurent polynomial), return Counter of dominant weights lam with multiplicities:
    P = sum m_lam ch V_lam. Peel off the lexicographically largest monomial (which is strictly decreasing) repeatedly."""
    numer = Counter({k:v for k,v in numer.items() if v})
    res = Counter()
    while numer:
        g = max(numer)  # lex-largest exponent; must be strictly decreasing
        m = numer[g]
        assert g[0] > g[1] > g[2], (g, numer)
        lam = tuple(g[i] - RHO[i] for i in range(3))
        res[lam] += m
        for k, v in weyl_num(g).items():
            numer[k] -= m*v
            if numer[k] == 0: del numer[k]
    return res
def char(lam):
    """weights of V_lam via peeling: ch V_lam * A_rho = A_{lam+rho}. Compute ch by dividing: iterative subtraction of x^mu A_rho."""
    numer = Counter({k:v for k,v in weyl_num(tuple(lam[i]+RHO[i] for i in range(3))).items() if v})
    ch = Counter()
    Arho = weyl_num(RHO)
    while numer:
        g = max(numer); m = numer[g]
        mu = tuple(g[i]-RHO[i] for i in range(3))
        ch[mu] += m
        for k, v in Arho.items():
            kk = tuple(k[i]+mu[i] for i in range(3))
            numer[kk] -= m*v
            if numer[kk] == 0: del numer[kk]
    return ch
def nconst(lam, W):
    """number of irreducible constituents of V_lam (x) W, W given as weight multiset."""
    numer = Counter()
    for beta, m in W.items():
        for k, v in weyl_num(tuple(lam[i]+RHO[i]+beta[i] for i in range(3))).items():
            numer[k] += m*v
    return sum(decompose(numer).values())
def dual(W): return Counter({tuple(-x for x in k): v for k, v in W.items()})

bad = cnt = 0
for l1 in range(0, 6):
    for l2 in range(0, l1+1):
        for l3 in range(0, l2+1):
            lam = (l1,l2,l3)
            for m1 in range(0, 5):
                for m2 in range(0, m1+1):
                    for m3 in range(0, m2+1):
                        W = char((m1,m2,m3)); cnt += 1
                        if nconst(lam, W) != nconst(lam, dual(W)): bad += 1; print("C1 FAIL", lam, (m1,m2,m3))
print("C1 via explicit Schur-polynomial decomposition (W irreducible, lam_1<=5, mu_1<=4): %d cases, %d failures" % (cnt, bad))
# W = V^{(x) j} (tensor powers of the vector rep), lam_1 <= 6, j <= 5
V = Counter({(1,0,0):1,(0,1,0):1,(0,0,1):1})
def tensor(U, W):
    c = Counter()
    for k1,v1 in U.items():
        for k2,v2 in W.items(): c[tuple(a+b for a,b in zip(k1,k2))] += v1*v2
    return c
bad = cnt = 0
for l1 in range(0, 7):
    for l2 in range(0, l1+1):
        for l3 in range(0, l2+1):
            W = Counter({(0,0,0):1})
            for j in range(0, 6):
                cnt += 1
                if nconst((l1,l2,l3), W) != nconst((l1,l2,l3), dual(W)): bad += 1; print("C1 FAIL tensor", (l1,l2,l3), j)
                W = tensor(W, V)
print("C1 for W = V^(x)j, lam_1<=6, j<=5: %d cases, %d failures" % (cnt, bad))

# --- C2: own chain counter, 3 rows
from functools import lru_cache
@lru_cache(None)
def N(nu, j):
    if j == 0: return 1
    tot = 0
    for i in range(3):
        if i == 0 or nu[i-1] > nu[i]:
            tot += N(nu[:i] + (nu[i]+1,) + nu[i+1:], j-1)
    return tot
bad = cnt = 0
for n1 in range(0, 11):
    for n2 in range(0, n1+1):
        for n3 in range(0, n2+1):
            nu = (n1,n2,n3); nus = (n1-n3, n1-n2, 0)
            for j in range(0, 26):
                cnt += 1
                if N(nu, j) != N(nus, j): bad += 1; print("C2 FAIL", nu, j)
print("C2 chain symmetry nu_1<=10, j<=25: %d cases, %d failures" % (cnt, bad))
# is the hypothesis nu* (not e.g. nu^t) essential? show N_j((3)) vs N_j((3,3)) small values
print("N_j((3,0,0)) j<=8:", [N((3,0,0), j) for j in range(9)])
print("N_j((3,3,0)) j<=8:", [N((3,3,0), j) for j in range(9)])

# --- C3 via project code: d_j((3^K,2,2,2)) = d_j((3^K,1,1,1)) for j<=K, and reduction d_j = N_j((3)) resp N_j((3,3))
for K in range(1, 12):
    lam = (3,)*K + (2,2,2); mu = (3,)*K + (1,1,1)
    u = d_vector(lam); v = d_vector(mu)
    first = next((j for j in range(len(u)) if u[j] != v[j]), None)
    red = all(u[j] == N((3,0,0), j) and v[j] == N((3,3,0), j) for j in range(K+1))
    print("R_%d: d_j equal for j<=K: %s; first differing index %s; reduction to N_j: %s" % (K, all(u[j]==v[j] for j in range(K+1)), first, red))
# explicit skew-count reduction check for K=4, j=4 using f_skew directly (definition-level, no recursion)
K = 4; lam = (3,)*K + (2,2,2); n = sum(lam); j = 4
from young import partitions
tot = 0
for nu in partitions(n - j):
    if len(nu) <= len(lam) and all(nu[i] <= lam[i] for i in range(len(nu))):
        tot += f_skew(lam, nu)
print("K=4, j=4: sum_nu f^{lam/nu} =", tot, " N_4((3,0,0)) =", N((3,0,0), 4), " d_4 =", d_vector(lam)[4])
