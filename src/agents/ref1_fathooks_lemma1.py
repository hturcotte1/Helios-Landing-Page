"""Referee #1, independent check of fathooks.md Lemma 1 (two-rectangle window).
No use of young.py / census.py: everything recomputed from scratch, exact integers.

Checks
 (A) d_j(a^b) = sum_{rho |- j, rho_1<=a, l(rho)<=b} f^rho   (Lemma 4.1 / fact 3, used in the proof)   n<=40
 (B) Lemma 1(i)  d_j(lam) = sum C(j,j1) d_{j1}(a^b) d_{j2}(c^d), j<=a+d              all 2-corner n<=NMAX
 (C) Lemma 1(ii) d_{a+d+1}(lam) = same + C(a+d,a)                                     all 2-corner n<=NMAX
 (D) sharpness: at j=a+d+2 the plain product formula fails for all shapes (sanity; not part of the lemma)
 (E) filter-level checks for all 2-corner shapes with n<=NFILT:
       every filter S meeting R0 has |S|>=a+d+1, with equality iff S = hook H;
       every filter with |S|<=a+d is S1 u S2 with S1,S2 filters of R1,R2; every such union is a filter;
       e(S) = C(|S|,|S1|) e(S1) e(S2) for those; e(H) = C(a+d,a);
       d_j = sum_{|S|=j} e(S) with e = number of linear extensions computed by brute force.
 (F) every partition of n<=16 with exactly 2 removable corners is ((a+c)^b,c^d) for some a,b,c,d>=1
"""
import sys
from functools import lru_cache
from math import comb, factorial
from itertools import product

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
NFILT = int(sys.argv[2]) if len(sys.argv) > 2 else 16

# ---------- own primitives ----------
def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for p in range(min(n, maxp), 0, -1):
        for rest in partitions(n - p, p): yield (p,) + rest

def conj(lam):
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()

def f_hook(lam):
    n = sum(lam)
    if n == 0: return 1
    lt = conj(lam); h = 1
    for i, r in enumerate(lam):
        for j in range(r):
            h *= (r - j) + (lt[j] - i) - 1
    assert factorial(n) % h == 0
    return factorial(n) // h

def corners(lam):
    return [i for i in range(len(lam)) if i == len(lam) - 1 or lam[i] > lam[i + 1]]

def remove(lam, i):
    l = list(lam); l[i] -= 1
    if l[i] == 0: l.pop()
    return tuple(l)

@lru_cache(maxsize=None)
def dvec(lam):
    n = sum(lam)
    if n == 0: return (1,)
    d = [0] * (n + 1); d[0] = 1
    for i in corners(lam):
        sub = dvec(remove(lam, i))
        for j, v in enumerate(sub): d[j + 1] += v
    return tuple(d)

def fathook(a, b, c, d): return tuple([a + c] * b + [c] * d)

def shapes(N):
    for b in range(1, N + 1):
        for d in range(1, N + 1):
            for c in range(1, N + 1):
                if b * c + c * d > N: break
                for a in range(1, N + 1):
                    if a * b + b * c + c * d > N: break
                    yield a, b, c, d

# ---------- (A) rectangle series ----------
@lru_cache(maxsize=None)
def g(x, y, k):   # sum over rho |- k in y x x box of f^rho
    return sum(f_hook(r) for r in partitions(k) if (not r or (r[0] <= x and len(r) <= y)))

bad = 0
for x in range(1, NMAX + 1):
    for y in range(1, NMAX // x + 1):
        dv = dvec(fathook(0, 0, x, y)[:0] + tuple([x] * y))
        for k in range(x * y + 1):
            if dv[k] != g(x, y, k): bad += 1; print("A FAIL", x, y, k)
print(f"(A) rectangles x*y<={NMAX}: failures = {bad}")

# ---------- (B),(C),(D) ----------
badB = badC = 0; cnt = 0; sharp_ok = 0; sharp_tot = 0
for a, b, c, d in shapes(NMAX):
    lam = fathook(a, b, c, d); n = sum(lam); dv = dvec(lam); cnt += 1
    def prod(j): return sum(comb(j, j1) * g(a, b, j1) * g(c, d, j - j1) for j1 in range(j + 1))
    for j in range(a + d + 1):
        if dv[j] != prod(j): badB += 1; print("B FAIL", (a, b, c, d), j, dv[j], prod(j))
    assert a + d + 1 <= n
    if dv[a + d + 1] != prod(a + d + 1) + comb(a + d, a): badC += 1; print("C FAIL", (a, b, c, d))
    if a + d + 2 <= n:
        sharp_tot += 1
        if dv[a + d + 2] != prod(a + d + 2): sharp_ok += 1
print(f"(B) Lemma1(i)  n<={NMAX}: {cnt} shapes, failures = {badB}")
print(f"(C) Lemma1(ii) n<={NMAX}: failures = {badC}")
print(f"(D) sharpness: plain product wrong at j=a+d+2 in {sharp_ok}/{sharp_tot} shapes")

# ---------- (E) filter-level ----------
def boxes(lam): return [(i, j) for i in range(1, len(lam) + 1) for j in range(1, lam[i - 1] + 1)]

def is_filter(S, lam):
    for (i, j) in S:
        if i < len(lam) and lam[i] >= j and (i + 1, j) not in S: return False
        if j < lam[i - 1] and (i, j + 1) not in S: return False
    return True

def lin_ext(S):
    """number of removal orders of S: repeatedly remove a box that has no box of S right of it or below it."""
    S = frozenset(S)
    @lru_cache(maxsize=None)
    def rec(T):
        if not T: return 1
        tot = 0
        for (i, j) in T:
            if (i + 1, j) not in T and (i, j + 1) not in T:
                tot += rec(T - {(i, j)})
        return tot
    return rec(S)

def subparts(lam):
    """all nu subset lam"""
    if not lam: yield (); return
    def rec(i, prev):
        if i == len(lam): yield (); return
        for p in range(min(lam[i], prev), -1, -1):
            if p == 0:
                yield ()
                return
            for rest in rec(i + 1, p): yield (p,) + rest
    yield from rec(0, lam[0])

badE = 0; cntE = 0
for a, b, c, d in shapes(NFILT):
    lam = fathook(a, b, c, d); n = sum(lam); cntE += 1
    B = set(boxes(lam))
    R0 = {(i, j) for (i, j) in B if i <= b and j <= c}
    R1 = {(i, j) for (i, j) in B if i <= b and j > c}
    R2 = {(i, j) for (i, j) in B if i > b and j <= c}
    assert R0 | R1 | R2 == B and len(R0) + len(R1) + len(R2) == n
    H = {(b, j) for j in range(c, a + c + 1)} | {(i, c) for i in range(b, b + d + 1)}
    assert len(H) == a + d + 1 and is_filter(H, lam)
    filters = []
    for nu in subparts(lam):
        S = B - set(boxes(nu))
        assert is_filter(S, lam)
        filters.append(frozenset(S))
    assert len(set(filters)) == len(filters)
    # d_j = sum e(S)
    dv = dvec(lam); dd = [0] * (n + 1)
    e = {S: lin_ext(S) for S in filters}
    for S in filters: dd[len(S)] += e[S]
    if tuple(dd) != dv: badE += 1; print("E d-vector FAIL", (a, b, c, d))
    # filter families
    F1 = {S for S in filters if S <= R1}; F2 = {S for S in filters if S <= R2}
    # check these are exactly the filters of the rectangles
    r1 = fathook(0, 0, a, b)[:0] + tuple([a] * b)
    if len(F1) != len(list(subparts(tuple([a] * b)))) or len(F2) != len(list(subparts(tuple([c] * d)))):
        badE += 1; print("E rect-filter count FAIL", (a, b, c, d))
    for S in filters:
        if S & R0:
            if len(S) < a + d + 1: badE += 1; print("E size FAIL", (a, b, c, d), sorted(S))
            if len(S) == a + d + 1 and S != frozenset(H): badE += 1; print("E uniqueness FAIL", (a, b, c, d), sorted(S))
        else:
            S1, S2 = S & R1, S & R2
            if S1 not in F1 or S2 not in F2: badE += 1; print("E split FAIL", (a, b, c, d))
            if e[S] != comb(len(S), len(S1)) * e[S1] * e[S2]: badE += 1; print("E product FAIL", (a, b, c, d), sorted(S))
    if e[frozenset(H)] != comb(a + d, a): badE += 1; print("E hook FAIL", (a, b, c, d))
    # every union of a filter of R1 and a filter of R2 is a filter
    fs = set(filters)
    for S1 in F1:
        for S2 in F2:
            if (S1 | S2) not in fs: badE += 1; print("E union FAIL", (a, b, c, d))
    # any other filter of size a+d+1 besides H meets R0? count
    others = [S for S in filters if len(S) == a + d + 1 and not (S <= R1 | R2)]
    if others != [frozenset(H)]: badE += 1; print("E extra filter FAIL", (a, b, c, d), others)
print(f"(E) filter-level checks n<={NFILT}: {cntE} shapes, failures = {badE}")

# ---------- (F) parametrisation ----------
badF = 0; cntF = 0
for n in range(1, 17):
    for lam in partitions(n):
        if len(corners(lam)) != 2: continue
        cntF += 1
        parts = sorted(set(lam), reverse=True)
        if len(parts) != 2: badF += 1; print("F FAIL", lam); continue
        p, q = parts; b = lam.count(p); d = lam.count(q); a = p - q; c = q
        if fathook(a, b, c, d) != lam or min(a, b, c, d) < 1: badF += 1; print("F FAIL", lam)
print(f"(F) 2-corner partitions n<=16: {cntF} shapes, parametrisation failures = {badF}")
