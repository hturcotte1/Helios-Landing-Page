"""Referee #2, second pass on fathooks.md Lemma 1.  Fully independent code (no young.py / census.py).
A: d_j from the Briefing definition sum_{nu} f^{lam/nu}, with f^{lam/nu} by the Aitken/Jacobi-Trudi determinant
   |lam/nu|! det[1/(lam_i - nu_j - i + j)!]  (exact Fractions), vs the corner recursion, all partitions n<=12.
B: Lemma 1 (i),(ii) for all 2-corner shapes n<=NMAX, G_{x,y} built from hook-length f^rho over rho in the y x x box.
C: filter-level claim (every filter meeting R0 has size >= a+d+1, equality only for the hook H) for n<=NFILT,
   by enumerating subpartitions nu of lam (S = lam \ nu meets R0 iff nu_b < c).
D: boundary cases: hooks (b=c=1, a+d+1=n), window indices j with j > ab or j > cd, and sharpness at a+d+1 / a+d+2.
Usage: python3 ref2_fathooks_lemma1_v2.py [NMAX=50] [NFILT=22]
"""
import sys
from functools import lru_cache
from fractions import Fraction
from math import comb, factorial

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 50
NFILT = int(sys.argv[2]) if len(sys.argv) > 2 else 22

def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for p in range(min(n, maxp), 0, -1):
        for rest in partitions(n - p, p): yield (p,) + rest

def conj(lam): return tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()

def f_hook(lam):
    n = sum(lam)
    if n == 0: return 1
    lt = conj(lam); h = 1
    for i, r in enumerate(lam):
        for j in range(r): h *= (r - j) + (lt[j] - i) - 1
    return factorial(n) // h

def corners(lam): return [i for i in range(len(lam)) if i == len(lam) - 1 or lam[i] > lam[i + 1]]
def remove(lam, i):
    l = list(lam); l[i] -= 1
    if l[i] == 0: l.pop()
    return tuple(l)

@lru_cache(None)
def dvec(lam):
    n = sum(lam)
    if n == 0: return (1,)
    d = [0] * (n + 1); d[0] = 1
    for i in corners(lam):
        for j, v in enumerate(dvec(remove(lam, i))): d[j + 1] += v
    return tuple(d)

def det(M):
    M = [row[:] for row in M]; n = len(M); s = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return Fraction(0)
        if p != c: M[c], M[p] = M[p], M[c]; s = -s
        s *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            if f: M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return s

def f_skew_JT(lam, nu):
    """Aitken: f^{lam/nu} = |lam/nu|! det[ 1/(lam_i - nu_j - i + j)! ]_{i,j}."""
    L = len(lam); nu = tuple(nu) + (0,) * (L - len(nu)); m = sum(lam) - sum(nu)
    M = [[Fraction(1, factorial(lam[i] - nu[j] - i + j)) if lam[i] - nu[j] - i + j >= 0 else Fraction(0)
          for j in range(L)] for i in range(L)]
    v = det(M) * factorial(m)
    assert v.denominator == 1; return int(v)

def subpartitions(lam):
    def rec(i, prev, nu):
        if i == len(lam): yield tuple(x for x in nu if x > 0); return
        for r in range(min(lam[i], prev), -1, -1):
            yield from rec(i + 1, r, nu + (r,))
    yield from rec(0, lam[0] if lam else 0, ())

def fathook(a, b, c, d): return tuple([a + c] * b + [c] * d)

def twocorner(n):
    for lam in partitions(n):
        if len(corners(lam)) == 2:
            p1, p2 = sorted(set(lam), reverse=True)
            yield (p1 - p2, lam.count(p1), p2, lam.count(p2)), lam

_G = {}
def G(x, y):
    """d_j(x^y) computed independently as sum_{rho |- j, rho_1<=x, l(rho)<=y} f^rho (hook lengths)."""
    key = (min(x, y), max(x, y))
    if key not in _G:
        _G[key] = [sum(f_hook(r) for r in partitions(j) if (not r or r[0] <= x) and len(r) <= y) for j in range(x * y + 1)]
    return _G[key]

def window(a, b, c, d, j):
    A, C = G(a, b), G(c, d)
    return sum(comb(j, j1) * A[j1] * C[j - j1] for j1 in range(j + 1) if j1 < len(A) and j - j1 < len(C))

fails = 0
# ---- A ----
cnt = 0
for n in range(0, 13):
    for lam in partitions(n):
        dv = [0] * (n + 1)
        for nu in subpartitions(lam): dv[n - sum(nu)] += f_skew_JT(lam, nu)
        cnt += 1
        if tuple(dv) != dvec(lam): fails += 1; print("A FAIL", lam, dv, dvec(lam))
print(f"[A] Briefing definition (Jacobi-Trudi skew counts) vs recursion, {cnt} partitions n<=12: fails = {fails}")

# ---- B, D ----
shapes = 0; hooks = 0; beyond = 0; sharp1 = 0; sharp2 = 0; sharp2tot = 0
for n in range(3, NMAX + 1):
    for (a, b, c, d), lam in twocorner(n):
        shapes += 1; dv = dvec(lam)
        assert n == a * b + b * c + c * d and a + d + 1 <= n
        for j in range(0, a + d + 1):
            if dv[j] != window(a, b, c, d, j): fails += 1; print("B(i) FAIL", (a, b, c, d), j)
            if j > a * b or j > c * d: beyond += 1
        w1 = window(a, b, c, d, a + d + 1)
        if dv[a + d + 1] != w1 + comb(a + d, a): fails += 1; print("B(ii) FAIL", (a, b, c, d))
        if dv[a + d + 1] == w1: sharp1 += 1                    # would mean the correction is 0 -- never
        if b == 1 and c == 1:
            hooks += 1
            if not (a + d + 1 == n and dv[n] == comb(a + d, a) == f_hook(lam) and w1 == 0): fails += 1; print("D hook FAIL", lam)
        if a + d + 2 <= n:
            sharp2tot += 1
            if dv[a + d + 2] == window(a, b, c, d, a + d + 2): sharp2 += 1
print(f"[B] Lemma 1 (i),(ii) all 2-corner shapes n<={NMAX}: {shapes} shapes, fails = {fails}")
print(f"[D] hooks b=c=1 (a+d+1=n) checked: {hooks}; window indices with j>ab or j>cd exercised: {beyond}")
print(f"[D] sharpness: (i) without correction at j=a+d+1 holds for {sharp1}/{shapes} shapes; at j=a+d+2 uncorrected product holds for {sharp2}/{sharp2tot}")

# ---- C ----  filter-level: S = lam \ nu ; meets R0 iff (b,c) in S iff nu_b < c (1-indexed row b)
fc = 0; shapesC = 0
for n in range(3, NFILT + 1):
    for (a, b, c, d), lam in twocorner(n):
        shapesC += 1
        nuH = tuple(x for x in ([a + c] * (b - 1) + [c - 1] * (d + 1)) if x > 0)   # lam minus the hook H
        assert sum(lam) - sum(nuH) == a + d + 1
        seen_sizes = {}
        for nu in subpartitions(lam):
            nub = nu[b - 1] if len(nu) >= b else 0
            size = n - sum(nu)
            if nub < c:   # S meets R0
                if size < a + d + 1: fc += 1; print("C MIN FAIL", (a, b, c, d), nu)
                if size == a + d + 1 and nu != nuH: fc += 1; print("C UNIQUE FAIL", (a, b, c, d), nu)
            else:
                # S inside R1 u R2: rows 1..b of nu have length >= c, rows > b arbitrary <= c
                assert all(nu[i] >= c for i in range(b))
        # e(H) = C(a+d,a) via Jacobi-Trudi on the skew shape lam/nuH
        if f_skew_JT(lam, nuH) != comb(a + d, a): fc += 1; print("C e(H) FAIL", (a, b, c, d))
print(f"[C] filter-level minimum/uniqueness + e(H) for all 2-corner shapes n<={NFILT}: {shapesC} shapes, fails = {fc}")
print("TOTAL FAILS:", fails + fc)
