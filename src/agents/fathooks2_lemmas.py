"""fathooks (angle 'fathooks', second attempt): computational verification of every lemma used in
agent_notes/fathooks.md.  Notation: lam(a,b,c,d) = ((a+c)^b, c^d), a,b,c,d >= 1, n = ab+bc+cd, corner runs ((a,b),(c,d)).
Exact integer / Fraction arithmetic throughout.

Checks (all 2-corner shapes with n <= N, default N = 40):
 (W1)  d_j(lam) = j! [t^j] G_{a,b} G_{c,d}                      for 0 <= j <= a+d          [Lemma 1]
 (W1') d_{a+d+1}(lam) = (a+d+1)! [t^{a+d+1}] G_{a,b} G_{c,d} + C(a+d, a)   (if a+d+1 <= n)  [Lemma 1]
 (W2)  d_j(lam) = j! [t^j] I(t) G_{a,d}                          for 0 <= j <= min(b,c)     [Lemma 2]
 (W2') d_{m+1}(lam) = (m+1)! [t^{m+1}] I G_{a,d} - [b=m] - [c=m],  m = min(b,c) (if m+1 <= n) [Lemma 2]
 (BOX) d_j(lam) = # j-step saturated up-chains from a^d inside the (b+d) x (a+c) box  (n <= 24) [Lemma 2]
 (READ) the reading algorithm of Theorem 3 recovers x1, x2 and N_k for all k <= x1+x2-1 from the d-vector alone
 (S1-S4) the closed formulas for the F-families and rectangles used in Theorem F (all parameter values with n <= 60)
 (THM F) Conjecture B for every 2-corner shape with >= 2 unit parameters, against ALL partitions of n <= NB (default 30),
         and against all 2-corner shapes with n <= N.
"""
import os, sys
from fractions import Fraction
from math import factorial, comb
from functools import lru_cache
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from young import partitions, f_hook, involutions, up_set, conjugate, corner_runs
from census import d_vector

DEG = 42
I = involutions

def fathook(a, b, c, d):
    return tuple([a + c] * b + [c] * d)

def fathooks_upto(N):
    for b in range(1, N + 1):
        for d in range(1, N + 1):
            for c in range(1, N + 1):
                if b * c + c * d > N: break
                for a in range(1, N + 1):
                    if a * b + b * c + c * d > N: break
                    yield (a, b, c, d)

def egf(seq):
    return [Fraction(x, factorial(j)) for j, x in enumerate(seq)] + [Fraction(0)] * (DEG + 1 - len(seq))

def mul(p, q):
    r = [Fraction(0)] * (DEG + 1)
    for i, x in enumerate(p[:DEG + 1]):
        if x == 0: continue
        for j, y in enumerate(q[:DEG + 1 - i]):
            r[i + j] += x * y
    return r

def inv(p):
    r = [Fraction(0)] * (DEG + 1); r[0] = 1 / p[0]
    for k in range(1, DEG + 1):
        r[k] = -sum(p[i] * r[k - i] for i in range(1, k + 1)) / p[0]
    return r

I_ser = egf([I(k) for k in range(DEG + 1)])
# precompute once: for each k the list of (f^rho, rho_1, l(rho)) over rho |- k
_PART = [[(f_hook(r), (r[0] if r else 0), len(r)) for r in partitions(k)] for k in range(DEG + 1)]
_D = {}
def D_ser(x):   # sum_{rho: rho_1 > x} f^rho t^|rho|/|rho|!
    if x not in _D:
        _D[x] = egf([sum(f for f, r1, l in _PART[k] if r1 > x) for k in range(DEG + 1)])
    return _D[x]
_G = {}
def G_ser(x, y):  # P_{x^y}
    key = (min(x, y), max(x, y))
    if key not in _G:
        _G[key] = egf([sum(f for f, r1, l in _PART[k] if l <= y and r1 <= x) for k in range(DEG + 1)])
    return _G[key]
def coeff(ser, j): return ser[j] * factorial(j)

def box_chains(a, b, c, d, jmax):
    rows, cols = b + d, a + c
    @lru_cache(maxsize=None)
    def cnt(rho, j):
        if j == 0: return 1
        return sum(cnt(mu, j - 1) for mu in up_set(rho) if len(mu) <= rows and mu[0] <= cols)
    return [cnt(tuple([a] * d), j) for j in range(jmax + 1)]

def read_local(dv):
    """Theorem 3 reading algorithm.  Input: d-vector of a 2-corner shape.  Output (x1, x2, {k: N_k for k <= x1+x2-1})."""
    P = egf(list(dv))
    Q = [I_ser[k] - v for k, v in enumerate(mul(P, inv(I_ser)))]       # Q = I - P/I  ==  sum_x D_x  mod t^{x1+x2+1}
    # x1 = m: first k >= 1 with Q_k != 0 is k = x1+1 ;  N_{x1} = (x1+1)! Q_{x1+1}
    k = 1
    while coeff(Q, k) == 0: k += 1
    x1 = k - 1
    N = {}
    rem = [q for q in Q]
    def peel(x):     # subtract N_x D_x
        for i in range(DEG + 1): rem[i] -= N[x] * D_ser(x)[i]
    N[x1] = int(coeff(rem, x1 + 1)); peel(x1)
    if N[x1] >= 2:
        x2 = x1
    else:
        k = x1 + 1
        while coeff(rem, k + 1) == 0:
            N[k] = 0; k += 1
        x2 = k; N[x2] = int(coeff(rem, x2 + 1)); peel(x2)
    for k in range(x2 + 1, x1 + x2):
        N[k] = int(coeff(rem, k + 1)); peel(k)
    return x1, x2, N

def main(N=40, NB=30):
    bad = 0; cnt = 0; dvs = {}
    for (a, b, c, d) in fathooks_upto(N):
        lam = fathook(a, b, c, d); n = sum(lam); cnt += 1
        assert corner_runs(lam) == [(a, b), (c, d)]
        dv = d_vector(lam); dvs.setdefault(dv, []).append((a, b, c, d))
        P1 = mul(G_ser(a, b), G_ser(c, d))
        for j in range(0, min(a + d, n) + 1):
            if dv[j] != coeff(P1, j): bad += 1; print("W1 FAIL", (a, b, c, d), j)
        if a + d + 1 <= n and dv[a + d + 1] != coeff(P1, a + d + 1) + comb(a + d, a): bad += 1; print("W1' FAIL", (a, b, c, d))
        P2 = mul(I_ser, G_ser(a, d)); m = min(b, c)
        for j in range(0, min(m, n) + 1):
            if dv[j] != coeff(P2, j): bad += 1; print("W2 FAIL", (a, b, c, d), j)
        if m + 1 <= n and dv[m + 1] != coeff(P2, m + 1) - (b == m) - (c == m): bad += 1; print("W2' FAIL", (a, b, c, d))
        if n <= 24 and tuple(box_chains(a, b, c, d, n)) != tuple(dv): bad += 1; print("BOX FAIL", (a, b, c, d))
        # READ
        xs = sorted([a, b, c, d]); x1, x2 = xs[0], xs[1]
        r1, r2, Nk = read_local(dv)
        truth = {k: xs.count(k) for k in range(x1, x1 + x2)}
        if (r1, r2) != (x1, x2) or Nk != truth: bad += 1; print("READ FAIL", (a, b, c, d), (r1, r2, Nk), truth)
    print(f"[W1,W1',W2,W2',READ: n<={N}; BOX: n<=24] {cnt} shapes, failures = {bad}")
    # Conjecture B among 2-corner shapes
    coll = 0
    for dv, L in dvs.items():
        cls = set(min(q, q[::-1]) for q in L)
        if len(cls) > 1: coll += 1; print("COLLISION", L)
    print(f"Conjecture B among 2-corner shapes, n <= {N}: non-transpose collisions = {coll}")
    # S-formulas
    def F1(c, d): return tuple([c + 1] + [c] * d)
    def F2(b, c): return tuple([c + 1] * b + [c])
    def rect(c, d): return tuple([c] * d)
    badS = 0; cntS = 0
    for c in range(1, 61):
        for d in range(1, 61):
            if c * d + c + 1 > 60: break
            cntS += 1
            x = min(c, d); n = c * d + c + 1
            dv = d_vector(F1(c, d))
            if any(dv[j] != I(j + 1) for j in range(x + 1)): badS += 1; print("S3a", c, d)
            if dv[x + 1] != I(x + 2) - (c == x) - (d == x): badS += 1; print("S3b", c, d)
            if c < d and dv[c + 2] != I(c + 3) - 2 * (c + 2) - (d == c + 1): badS += 1; print("S3c", c, d)
            if d <= c and c >= 2 and d >= 2 and dv[d + 2] != I(d + 3) - (d + 3) - 2 * (d + 2) * (c == d) - (c == d + 1): badS += 1; print("S3d", c, d)
            b = d
            if (c + 1) * b + c <= 60:
                dv = d_vector(F2(b, c)); dR = d_vector(rect(c + 1, b + 1)); y = min(b, c)
                if any(dv[j] != dR[j + 1] for j in range(len(dv))): badS += 1; print("S1", b, c)
                if any(dv[j] != I(j + 1) for j in range(y + 1)): badS += 1; print("S2a", b, c)
                if dv[y + 1] != I(y + 2) - (b == y) - (c == y): badS += 1; print("S2b", b, c)
                if dv[y + 2] != I(y + 3) - (y + 3) * (b == y) - (b == y + 1) - (y + 3) * (c == y) - (c == y + 1): badS += 1; print("S2c", b, c)
            if c * d <= 60:
                dr = d_vector(rect(c, d))
                if any(dr[j] != I(j) for j in range(x + 1)): badS += 1; print("S4a", c, d)
                if x + 1 <= c * d and dr[x + 1] != I(x + 1) - (c == x) - (d == x): badS += 1; print("S4b", c, d)
                if x + 2 <= c * d and dr[x + 2] != I(x + 2) - (x + 2) * ((c == x) + (d == x)) - (c == x + 1) - (d == x + 1): badS += 1; print("S4c", c, d)
    print(f"[S1-S4: all parameter values with n <= 60] {cntS} (c,d) pairs, failures = {badS}")
    # Theorem F against all partitions
    badB = 0; cntB = 0
    for n in range(1, NB + 1):
        groups = {}
        for lam in partitions(n): groups.setdefault(d_vector(lam), []).append(lam)
        for lam in partitions(n):
            runs = corner_runs(lam)
            if len(runs) != 2: continue
            (a, b), (c, d) = runs
            if [a, b, c, d].count(1) < 2: continue
            cntB += 1
            if any(mu != lam and mu != conjugate(lam) for mu in groups[d_vector(lam)]): badB += 1; print("THM F FAIL", lam)
    print(f"[Theorem F vs all partitions, n <= {NB}] {cntB} shapes, failures = {badB}")

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 40, int(sys.argv[2]) if len(sys.argv) > 2 else 30)
