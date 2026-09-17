"""Referee #1 (second pass), independent re-check of fathooks.md Lemma 1.
Self-contained: no import of young.py / census.py. Exact integer arithmetic.

 (1) d-vector by the corner recursion (own code) vs d-vector by explicit filter enumeration
     with brute-force linear-extension counting (n <= NFILT) -- validates the definition used.
 (2) Lemma 1(i)  d_j = sum_{j1+j2=j} C(j,j1) d_{j1}(a^b) d_{j2}(c^d) for 0<=j<=a+d, all 2-corner n<=NMAX.
 (3) Lemma 1(ii) d_{a+d+1} = product term + C(a+d,a), all 2-corner n<=NMAX; also a+d+1<=n always.
 (4) sharpness: product formula (without correction) fails at j=a+d+1 always; with correction fails at a+d+2
     (informational; how often).
 (5) filter-level structural claims (n<=NFILT): every filter S meeting R0 has |S| >= a+d+1, with
     equality iff S == H; every filter S with |S| <= a+d+1, S != H, is S1 u S2 with S1,S2 filters of
     R1,R2, and conversely; e(S) = C(|S|,|S1|) e(S1) e(S2); e(H) = C(a+d,a).
 (6) G_{a,b} = P_{a^b}: d_j(a^b) = sum_{rho |- j, rho_1<=a, l(rho)<=b} f^rho (used implicitly), ab<=NMAX.
 (7) every partition of n<=NPART with exactly two removable corners is ((a+c)^b, c^d), a,b,c,d>=1.
"""
import sys
from functools import lru_cache
from math import comb, factorial

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
NFILT = int(sys.argv[2]) if len(sys.argv) > 2 else 16
NPART = 20

def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0:
        yield (); return
    for p in range(min(n, maxp), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest

def conj(lam):
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()

def f_hook(lam):
    n = sum(lam)
    if n == 0: return 1
    lt = conj(lam); h = 1
    for i, r in enumerate(lam):
        for j in range(r):
            h *= (r - j) + (lt[j] - i) - 1
    q, rem = divmod(factorial(n), h)
    assert rem == 0
    return q

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
        for j, v in enumerate(sub):
            d[j + 1] += v
    return tuple(d)

def rect(a, b):  # b rows of length a
    return (a,) * b

def shape(a, b, c, d):
    return (a + c,) * b + (c,) * d

def two_corner_params(n):
    for b in range(1, n + 1):
        for d in range(1, n + 1):
            for c in range(1, n + 1):
                if b * c + c * d > n: break
                rest = n - b * c - c * d
                if rest <= 0 or rest % b: continue
                a = rest // b
                yield (a, b, c, d)

# ---------- filters and linear extensions (brute force) ----------
def boxes(lam):
    return [(i, j) for i, r in enumerate(lam) for j in range(r)]

def filters_of(lam):
    """All order filters (subsets closed under moving right/down inside lam), as frozensets."""
    bx = boxes(lam)
    # generate by complements: nu subset lam partition -> filter = lam \ nu
    res = []
    def subparts(lam):
        # all partitions nu contained in lam
        if not lam:
            yield (); return
        for first in range(lam[0], -1, -1):
            rest_lam = tuple(min(x, first) for x in lam[1:])
            for r in subparts(rest_lam):
                if first == 0:
                    yield ()
                    break
                yield (first,) + r
    seen = set()
    for nu in subparts(lam):
        nu = tuple(x for x in nu if x > 0)
        if nu in seen: continue
        seen.add(nu)
        S = frozenset(b for b in bx if b[1] >= (nu[b[0]] if b[0] < len(nu) else 0))
        res.append(S)
    return res

@lru_cache(maxsize=None)
def lin_ext(S):
    """Number of removal orders of the set S of boxes (remove a box when no box right of or below it in S)."""
    if not S: return 1
    tot = 0
    for (i, j) in S:
        if (i + 1, j) in S or (i, j + 1) in S: continue
        tot += lin_ext(S - {(i, j)})
    return tot

fails = 0
def report(msg):
    global fails
    fails += 1
    print("FAIL:", msg)

# (7)
cnt7 = 0
for n in range(1, NPART + 1):
    for lam in partitions(n):
        if len(corners(lam)) == 2:
            cnt7 += 1
            parts = sorted(set(lam), reverse=True)
            assert len(parts) == 2
            p1, p2 = parts
            b = lam.count(p1); d = lam.count(p2); c = p2; a = p1 - p2
            if shape(a, b, c, d) != lam or min(a, b, c, d) < 1:
                report(f"(7) {lam}")
print(f"(7) two-corner partitions n<={NPART}: {cnt7} all of the form ((a+c)^b,c^d)")

# (6)
cnt6 = 0
for a in range(1, NMAX + 1):
    for b in range(1, NMAX // a + 1):
        dv = dvec(rect(a, b))
        for j in range(a * b + 1):
            s = sum(f_hook(r) for r in partitions(j) if r[0] <= a and len(r) <= b) if j else 1
            if s != dv[j]:
                report(f"(6) rect {a}^{b} j={j}: {s} vs {dv[j]}")
        cnt6 += 1
print(f"(6) rectangles with ab<={NMAX}: {cnt6} checked")

# (1),(5) filter-level
cnt1 = 0
for n in range(3, NFILT + 1):
    for (a, b, c, d) in two_corner_params(n):
        lam = shape(a, b, c, d)
        dv = dvec(lam)
        F = filters_of(lam)
        assert len(set(F)) == len(F)
        # (1)
        byj = [0] * (n + 1)
        for S in F: byj[len(S)] += lin_ext(S)
        if tuple(byj) != dv:
            report(f"(1) {lam}: filter-sum {byj} vs recursion {dv}")
        # (5)
        R0 = {(i, j) for i in range(b) for j in range(c)}
        R1 = {(i, j) for i in range(b) for j in range(c, a + c)}
        R2 = {(i, j) for i in range(b, b + d) for j in range(c)}
        H = frozenset({(b - 1, j) for j in range(c - 1, a + c)} | {(i, c - 1) for i in range(b - 1, b + d)})
        F1 = set(filters_of(rect(a, b)))          # coordinates (i, j) in rectangle b x a
        F2 = set(filters_of(rect(c, d)))
        def to_R1(S): return frozenset((i, j + c) for (i, j) in S)
        def to_R2(S): return frozenset((i + b, j) for (i, j) in S)
        unions = {}
        for S1 in F1:
            for S2 in F2:
                unions[to_R1(S1) | to_R2(S2)] = (S1, S2)
        Fset = set(F)
        for U, (S1, S2) in unions.items():
            if U not in Fset:
                report(f"(5) union not a filter {lam} {sorted(U)}")
            e = lin_ext(U)
            if e != comb(len(U), len(S1)) * lin_ext(S1) * lin_ext(S2):
                report(f"(5) e(S) product rule fails {lam} {sorted(U)}")
        if H not in Fset:
            report(f"(5) H not a filter {lam}")
        if lin_ext(H) != comb(a + d, a):
            report(f"(5) e(H)={lin_ext(H)} != C({a+d},{a}) for {lam}")
        if len(H) != a + d + 1:
            report(f"(5) |H| wrong {lam}")
        for S in F:
            meets = bool(S & R0)
            if meets:
                if len(S) < a + d + 1:
                    report(f"(5) filter meeting R0 with size {len(S)} < {a+d+1}: {lam} {sorted(S)}")
                if len(S) == a + d + 1 and S != H:
                    report(f"(5) filter meeting R0 of size a+d+1 not H: {lam} {sorted(S)}")
                if S in unions:
                    report(f"(5) filter meeting R0 in unions?! {lam}")
            else:
                if S not in unions:
                    report(f"(5) filter inside R1uR2 not a union of rectangle filters: {lam} {sorted(S)}")
        cnt1 += 1
print(f"(1),(5) filter-level checks on {cnt1} two-corner shapes, n<={NFILT}")

# (2),(3),(4)
cnt2 = 0; sharp_i = 0; sharp_ii = 0; sharp_ii_tot = 0
for n in range(3, NMAX + 1):
    for (a, b, c, d) in two_corner_params(n):
        lam = shape(a, b, c, d)
        dv = dvec(lam)
        da = dvec(rect(a, b)); dc = dvec(rect(c, d))
        def prod(j):
            return sum(comb(j, j1) * da[j1] * dc[j - j1] for j1 in range(max(0, j - c * d), min(j, a * b) + 1))
        if a + d + 1 > n:
            report(f"(3) a+d+1 > n for {(a,b,c,d)}")
        for j in range(0, a + d + 1):
            if dv[j] != prod(j):
                report(f"(2) Lemma 1(i) fails {(a,b,c,d)} j={j}: {dv[j]} vs {prod(j)}")
        j = a + d + 1
        if dv[j] != prod(j) + comb(a + d, a):
            report(f"(3) Lemma 1(ii) fails {(a,b,c,d)}: {dv[j]} vs {prod(j)}+{comb(a+d,a)}")
        if dv[j] == prod(j):
            sharp_i += 1
        if a + d + 2 <= n:
            sharp_ii_tot += 1
            if dv[a + d + 2] == prod(a + d + 2) + comb(a + d, a) * 0 + 0 and False:
                pass
            # does the *product formula* also hold at a+d+2 (i.e. is the lemma's range sharp)?
            if dv[a + d + 2] == prod(a + d + 2):
                sharp_ii += 1
        cnt2 += 1
print(f"(2),(3) Lemma 1 (i),(ii) on {cnt2} two-corner shapes n<={NMAX}")
print(f"(4) shapes where product formula holds at j=a+d+1 (should be 0): {sharp_i}; "
      f"where product formula holds at a+d+2 (uncorrected): {sharp_ii}/{sharp_ii_tot}")
print("TOTAL FAILURES:", fails)
