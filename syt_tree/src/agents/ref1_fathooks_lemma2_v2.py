"""Referee #1, second independent check of fathooks.md Lemma 2, built differently from ref1_fathooks_lemma2.py:
 * d_j(lam) via the FILTER definition (order filters S of the diagram, |S| = j, weighted by # linear extensions
   e(S), computed by a poset DP) -- no removal recursion; also cross-checked with census.d_vector.
 * skew SYT counts f^{lam/nu} by explicit DP over subsets of the skew diagram (linear extensions), and the
   rotation claim f^{lam/nu} = f^{rho(nu)/a^d} checked box-by-box with the explicit rotation map.
 * every partition of n <= NALL with exactly two removable corners is enumerated with young.partitions and
   shown to equal lam(a,b,c,d) for a unique (a,b,c,d); Lemma 2 (i),(ii),(iii) checked on all of them.
 * (ii),(iii) additionally for all (a,b,c,d) with n <= NMAX using the filter-free chain counts.
"""
import sys
sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src')
from functools import lru_cache
from math import comb, factorial
from itertools import combinations
from young import partitions, removable_corners, corner_runs
from census import d_vector, s_up as lib_s_up

NALL = int(sys.argv[1]) if len(sys.argv) > 1 else 16
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 40

def boxes(lam):
    return [(i, j) for i, r in enumerate(lam) for j in range(r)]

def lin_ext_count(S):
    """# linear extensions of the sub-poset S (set of boxes) of the diagram order (i,j) <= (i',j') iff i<=i', j<=j'.
    We count removal orders of a filter S: remove a box only after everything weakly right-and-below is gone --
    equivalently linear extensions of the poset with (i',j') removed before (i,j). Count = linear extensions."""
    S = tuple(sorted(S)); idx = {b: k for k, b in enumerate(S)}
    below = [0] * len(S)  # bitmask of boxes weakly right/below (excluding self) inside S
    for k, (i, j) in enumerate(S):
        for l, (i2, j2) in enumerate(S):
            if l != k and i2 >= i and j2 >= j: below[k] |= 1 << l
    full = (1 << len(S)) - 1
    @lru_cache(maxsize=None)
    def rec(removed):
        if removed == full: return 1
        tot = 0
        for k in range(len(S)):
            if not (removed >> k) & 1 and (below[k] & ~removed) == 0:
                tot += rec(removed | (1 << k))
        return tot
    return rec(0)

def filters(lam):
    """all order filters (closed under moving right/down) of lam, i.e. complements of sub-partitions nu."""
    bx = boxes(lam); out = []
    for k in range(sum(lam) + 1):
        for nu in partitions(k):
            if len(nu) <= len(lam) and all(nu[i] <= lam[i] for i in range(len(nu))):
                S = [(i, j) for (i, j) in bx if j >= (nu[i] if i < len(nu) else 0)]
                out.append((nu, frozenset(S)))
    return out

def d_by_filters(lam):
    n = sum(lam); d = [0] * (n + 1)
    for nu, S in filters(lam):
        # check S really is a filter
        for (i, j) in S:
            for (i2, j2) in boxes(lam):
                if i2 >= i and j2 >= j: assert (i2, j2) in S
        d[len(S)] += lin_ext_count(S)
    return tuple(d)

def addable(lam):
    res = []
    for i in range(len(lam) + 1):
        if i == 0 or (i < len(lam) and lam[i - 1] > lam[i]) or (i == len(lam)):
            l = list(lam)
            if i == len(l): l.append(1)
            else: l[i] += 1
            res.append(tuple(l))
    return res

@lru_cache(maxsize=None)
def up_chains(rho, k, R=None, C=None):
    if k == 0: return 1
    t = 0
    for mu in addable(rho):
        if R is not None and (len(mu) > R or mu[0] > C): continue
        t += up_chains(mu, k - 1, R, C)
    return t

def involution_numbers(N):
    I = [1, 1]
    for k in range(2, N + 1): I.append(I[-1] + (k - 1) * I[-2])
    return I

def fathook(a, b, c, d): return tuple([a + c] * b + [c] * d)

def params_of(lam):
    """(a,b,c,d) with lam = ((a+c)^b, c^d), or None."""
    ds = sorted(set(lam), reverse=True)
    if len(ds) != 2: return None
    p, q = ds; b = lam.count(p); d = lam.count(q)
    return (p - q, b, q, d)

def rotation_check(a, b, c, d):
    """for every nu <= lam: explicit rotation of the skew diagram lam/nu, check it is rho/a^d with rho a partition
    in the box containing a^d, and that the skew SYT count agrees (via linear extensions of the reversed poset)."""
    lam = fathook(a, b, c, d); R, C = b + d, a + c; n = sum(lam)
    rect = tuple([a] * d); targets = set()
    for nu, S in filters(lam):
        rotS = frozenset((R - 1 - i, C - 1 - j) for (i, j) in S)
        rect_boxes = set(boxes(rect))
        union = rect_boxes | set(rotS)
        assert not (rect_boxes & set(rotS))
        # union must be a partition diagram
        rows = [0] * R
        for (i, j) in union: rows[i] += 1
        for (i, j) in union: assert j < rows[i]
        rho = tuple(x for x in rows if x > 0)
        assert all(rho[i] >= rho[i + 1] for i in range(len(rho) - 1)), (nu, rho)
        assert len(rho) <= R and (not rho or rho[0] <= C)
        assert sum(rho) == a * d + len(S)
        targets.add(rho)
        # skew SYT count of rho/a^d as linear extensions where order is reversed (fill increasing left-to-right/top-down):
        # number of fillings of rotS increasing along rows/cols = # linear extensions of rotS with (i,j) before (i',j')
        # when i<=i', j<=j'. That is lin_ext_count applied to the 180-rotated set (the DP removes maximal elements first,
        # which for the rotated set is the same count).
        f1 = lin_ext_count(S)                      # removal orders of S in lam  (= f^{lam/nu})
        f2 = lin_ext_count(rotS)                   # removal orders of rotS as a filter of rho  (= f^{rho/a^d})
        assert f1 == f2, (nu, rho, f1, f2)
    allrho = set(rho for k in range(a * d, R * C + 1) for rho in partitions(k)
                 if len(rho) <= R and rho[0] <= C and len(rho) >= d and all(rho[i] >= a for i in range(d)))
    assert targets == allrho, (a, b, c, d)

def main():
    I = involution_numbers(NMAX + 2)
    fails = 0
    # ---- part A: every two-corner partition of n <= NALL ----
    cntA = 0
    for n in range(1, NALL + 1):
        for lam in partitions(n):
            if len(removable_corners(lam)) != 2: 
                assert params_of(lam) is None; continue
            p = params_of(lam); assert p is not None
            a, b, c, d = p; assert fathook(a, b, c, d) == lam and min(p) >= 1
            assert corner_runs(lam) == [(a, b), (c, d)]
            cntA += 1
            dv = d_by_filters(lam)
            if dv != d_vector(lam): fails += 1; print("filter d-vector != library", lam)
            m = min(b, c); R, C = b + d, a + c; rect = tuple([a] * d)
            # (i)
            box = tuple(up_chains(rect, j, R, C) for j in range(n + 1))
            if box != dv: fails += 1; print("(i) FAIL", p, box, dv)
            # (ii)
            for j in range(m + 1):
                sj = up_chains(rect, j)
                ser = sum(comb(j, k) * I[k] * up_chains((), j - k, d, a) for k in range(j + 1))
                if not (dv[j] == sj == ser == lib_s_up(rect, j)): fails += 1; print("(ii) FAIL", p, j, dv[j], sj, ser)
            # (iii)
            j = m + 1; assert j <= n
            ser = sum(comb(j, k) * I[k] * up_chains((), j - k, d, a) for k in range(j + 1))
            if dv[j] != ser - (b == m) - (c == m): fails += 1; print("(iii) FAIL", p, dv[j], ser)
            if dv[j] == up_chains(rect, j): fails += 1; print("(ii) unexpectedly holds at m+1", p)
            # rotation bijection with explicit boxes
            if n <= 12: rotation_check(a, b, c, d)
    print(f"A: all {cntA} two-corner partitions of n <= {NALL}: failures so far {fails}")
    # ---- part B: (ii),(iii) for all (a,b,c,d) with n <= NMAX via chain counts and library d_vector ----
    cntB = 0
    for a in range(1, NMAX):
        for b in range(1, NMAX):
            if a * b >= NMAX: break
            for c in range(1, NMAX):
                if a * b + b * c >= NMAX: break
                for d in range(1, NMAX):
                    n = a * b + b * c + c * d
                    if n > NMAX: break
                    cntB += 1
                    lam = fathook(a, b, c, d); dv = d_vector(lam); m = min(b, c); rect = tuple([a] * d)
                    for j in range(m + 2):
                        ser = sum(comb(j, k) * I[k] * up_chains((), j - k, d, a) for k in range(j + 1))
                        want = ser if j <= m else ser - (b == m) - (c == m)
                        if dv[j] != want: fails += 1; print("B FAIL", (a, b, c, d), j, dv[j], want)
                        if j <= m and up_chains(rect, j) != ser: fails += 1; print("B s_j FAIL", (a, b, c, d), j)
                    # also (i) at j = m+1 with box-chains
                    if up_chains(rect, m + 1, b + d, a + c) != dv[m + 1]: fails += 1; print("B box FAIL", (a, b, c, d))
    print(f"B: {cntB} shapes (a,b,c,d) with n <= {NMAX}: total failures {fails}")

if __name__ == "__main__":
    main()
