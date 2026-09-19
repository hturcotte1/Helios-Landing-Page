"""Referee #2: independent check of fathooks.md Lemma 1 (two-rectangle window).
Everything recomputed from scratch with exact integers (no young.py / census.py).
Usage: python3 ref2_fathooks_lemma1.py [NMAX=40] [NFILT=16]
"""
import sys
from functools import lru_cache
from math import comb, factorial
from itertools import combinations

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
NFILT = int(sys.argv[2]) if len(sys.argv) > 2 else 16

def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for p in range(min(n, maxp), 0, -1):
        for rest in partitions(n - p, p): yield (p,) + rest

def corners(lam):
    return [i for i in range(len(lam)) if i == len(lam) - 1 or lam[i] > lam[i + 1]]

def remove(lam, i):
    l = list(lam); l[i] -= 1
    if l[i] == 0: l.pop()
    return tuple(l)

@lru_cache(maxsize=None)
def dvec(lam):
    """(d_0..d_n) via d_j(lam) = sum_c d_{j-1}(lam - c)."""
    n = sum(lam)
    if n == 0: return (1,)
    d = [0] * (n + 1); d[0] = 1
    for i in corners(lam):
        for j, v in enumerate(dvec(remove(lam, i))): d[j + 1] += v
    return tuple(d)

# ---------- brute force via filters + linear extensions ----------
def boxes(lam):
    return [(i, j) for i, r in enumerate(lam) for j in range(r)]

def is_filter(S, lam):
    S = set(S)
    for (i, j) in S:
        for (i2, j2) in boxes(lam):
            if i2 >= i and j2 >= j and (i2, j2) not in S: return False
    return True

def filters(lam):
    """All order filters of lam (complements of subpartitions)."""
    out = []
    B = boxes(lam)
    # generate subpartitions nu <= lam, filter = lam \ nu
    def rec(i, prev, nu):
        if i == len(lam):
            Sn = set(nu_boxes(nu))
            out.append(frozenset(b for b in B if b not in Sn)); return
        for r in range(0, min(lam[i], prev) + 1):
            rec(i + 1, r, nu + (r,))
    def nu_boxes(nu): return [(i, j) for i, r in enumerate(nu) for j in range(r)]
    rec(0, lam[0] if lam else 0, ())
    return out

@lru_cache(maxsize=None)
def lin_ext(S):
    """number of linear extensions of (S, product order), i.e. removal orders: a box may be removed
    once every box weakly right-and-below it is gone."""
    if not S: return 1
    tot = 0
    for b in S:
        if all(not (x[0] >= b[0] and x[1] >= b[1]) for x in S if x != b):
            tot += lin_ext(S - {b})
    return tot

def dvec_brute(lam):
    n = sum(lam); d = [0] * (n + 1)
    for S in filters(lam):
        assert is_filter(S, lam)
        d[len(S)] += lin_ext(S)
    return tuple(d)

def fathook(a, b, c, d): return tuple([a + c] * b + [c] * d)

def prod_formula(a, b, c, d, j):
    A = dvec(tuple([a] * b)); C = dvec(tuple([c] * d))
    return sum(comb(j, j1) * A[j1] * C[j - j1] for j1 in range(0, j + 1) if j1 < len(A) and j - j1 < len(C))

fails = 0
# (0) sanity: recursion d-vector == brute force filter/linear-extension sum, all partitions n<=9
for n in range(0, 10):
    for lam in partitions(n):
        if dvec(lam) != dvec_brute(lam): fails += 1; print("DVEC MISMATCH", lam)
print("[0] recursion vs brute-force filter sum, all partitions n<=9: fails =", fails)

# (F) every partition of n<=NFILT with 2 removable corners has form ((a+c)^b, c^d)
cnt = 0
for n in range(1, NFILT + 1):
    for lam in partitions(n):
        if len(corners(lam)) == 2:
            parts = sorted(set(lam), reverse=True)
            p1, p2 = parts
            b = lam.count(p1); d = lam.count(p2); c = p2; a = p1 - p2
            assert fathook(a, b, c, d) == lam and min(a, b, c, d) >= 1
            cnt += 1
print("[F] 2-corner partitions n<=%d all of form ((a+c)^b,c^d):" % NFILT, cnt)

# (B),(C) Lemma 1 (i),(ii) for all 2-corner shapes n<=NMAX, plus sharpness at a+d+2
shapes = 0; sharp_ok = 0; sharp_tot = 0
for n in range(3, NMAX + 1):
    for lam in partitions(n):
        if len(corners(lam)) != 2: continue
        p1, p2 = sorted(set(lam), reverse=True)
        b = lam.count(p1); d = lam.count(p2); c = p2; a = p1 - p2
        shapes += 1
        dv = dvec(lam)
        assert a + d + 1 <= n
        for j in range(0, a + d + 1):
            if dv[j] != prod_formula(a, b, c, d, j): fails += 1; print("L1(i) FAIL", (a, b, c, d), j)
        if dv[a + d + 1] != prod_formula(a, b, c, d, a + d + 1) + comb(a + d, a): fails += 1; print("L1(ii) FAIL", (a, b, c, d))
        if a + d + 2 <= n:
            sharp_tot += 1
            if dv[a + d + 2] != prod_formula(a, b, c, d, a + d + 2) + comb(a + d, a): sharp_ok += 1
print("[B,C] Lemma 1 (i),(ii) for all 2-corner shapes n<=%d: %d shapes, fails = %d" % (NMAX, shapes, fails))
print("      sharpness: plain '+C(a+d,a)' formula at j=a+d+2 fails for %d of %d shapes (expected: not a claim)" % (sharp_ok, sharp_tot))

# (E) filter-level claims for all 2-corner shapes with n<=NFILT
fe = 0; shapesE = 0
for n in range(3, NFILT + 1):
    for lam in partitions(n):
        if len(corners(lam)) != 2: continue
        p1, p2 = sorted(set(lam), reverse=True)
        b = lam.count(p1); d = lam.count(p2); c = p2; a = p1 - p2
        shapesE += 1
        # 0-indexed: R0 = i<b, j<c ; R1 = i<b, c<=j<a+c ; R2 = b<=i<b+d, j<c
        H = frozenset([(b - 1, j) for j in range(c - 1, a + c)] + [(i, c - 1) for i in range(b - 1, b + d)])
        assert len(H) == a + d + 1 and is_filter(H, lam)
        if lin_ext(H) != comb(a + d, a): fe += 1; print("e(H) FAIL", (a, b, c, d))
        R1 = tuple([a] * b); R2 = tuple([c] * d)
        F1 = {frozenset((i, j + c) for (i, j) in S) for S in filters(R1)}
        F2 = {frozenset((i + b, j) for (i, j) in S) for S in filters(R2)}
        allF = filters(lam)
        small = set(); count_by_size = {}
        for S in allF:
            meets = any(i < b and j < c for (i, j) in S)
            if meets:
                if len(S) < a + d + 1: fe += 1; print("MIN SIZE FAIL", (a, b, c, d), S)
                if len(S) == a + d + 1 and S != H: fe += 1; print("UNIQUE HOOK FAIL", (a, b, c, d), S)
            else:
                S1 = frozenset(x for x in S if x[0] < b); S2 = frozenset(x for x in S if x[0] >= b)
                if S1 not in F1 or S2 not in F2: fe += 1; print("DECOMP FAIL", (a, b, c, d), S)
                if lin_ext(S) != comb(len(S), len(S1)) * lin_ext(S1) * lin_ext(S2): fe += 1; print("e PRODUCT FAIL", (a, b, c, d), S)
        # converse: every union S1 u S2 is a filter of lam
        setF = set(allF)
        for S1 in F1:
            for S2 in F2:
                if (S1 | S2) not in setF: fe += 1; print("UNION NOT FILTER", (a, b, c, d), S1, S2)
        # count: filters of size <= a+d are exactly the unions; size a+d+1: unions + H
        for j in range(0, a + d + 2):
            nf = sum(1 for S in allF if len(S) == j)
            nu = sum(1 for S1 in F1 for S2 in F2 if len(S1) + len(S2) == j)
            if nf != nu + (1 if j == a + d + 1 else 0): fe += 1; print("COUNT FAIL", (a, b, c, d), j, nf, nu)
print("[E] filter-level claims of the proof for all 2-corner shapes n<=%d: %d shapes, fails = %d" % (NFILT, shapesE, fe))
print("TOTAL FAILS:", fails + fe)
