"""ref1_skeptic_s3.py -- referee #1 independent check of skeptic2.md Section 3 (Theorem S3, C1-C3, H_k).

Part A: S3 directly (own sgn, own loops): box 16 (must give 39168 points), box 30, random large integers,
        un-normalised a (arbitrary strictly decreasing) and b (arbitrary, incl. ties, unsorted, negative).
Part B: independent enumeration of the cells of the arrangement H in R^4 (coords A,B,p,q) restricted to
        the region A>0, B>0, p>=q>=0:
        - candidate rays = primitive kernel vectors of all rank-3 triples of forms in H (check |r_i| <= 4),
        - cells = sign vectors of sums of pairwise sign-compatible subsets of <= 4 candidate rays (by Caratheodory
          every cell contains such a sum; every such sum lies in a cell), restricted to the region,
        - verify every such cell has a representative in the box {1<=A,B<=16, 0<=q<=p<=16} and S3 holds on it,
        - verify the set of sign vectors realised in the box equals the set of cells (no cell missed).
Part C: chain symmetry N_j(nu) = N_j(nu*) (own chain counter), nu_1 <= 10, j <= 24.
Part D: R_K via project d_vector (K <= 20), reduction d_j = N_j((3)) / N_j((3,3)), and H_k differing indices,
        Lemma B formula and C_1 difference (k <= 10).
"""
import os
import sys, random, itertools
from math import gcd, comb
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from census import d_vector
from young import conjugate

PERMS = list(itertools.permutations(range(3)))

def sgn(v):
    s = 1
    for i in range(3):
        for j in range(i + 1, 3):
            if v[i] == v[j]:
                return 0
            if v[i] < v[j]:
                s = -s
    return s

def sortsign(v):
    """sign of permutation sorting v decreasingly (independent formula via explicit permutation parity)."""
    if len(set(v)) < 3:
        return 0
    order = sorted(range(3), key=lambda i: -v[i])
    # parity of permutation `order`
    par = 0
    seen = [False] * 3
    for i in range(3):
        if not seen[i]:
            l = 0; j = i
            while not seen[j]:
                seen[j] = True; j = order[j]; l += 1
            par += l - 1
    return -1 if par % 2 else 1

def lhs(a, b):
    return sum(sgn(tuple(a[i] + b[s[i]] for i in range(3))) for s in PERMS)

def rhs(a, b):
    return sum(sgn(tuple(a[i] - b[s[i]] for i in range(3))) for s in PERMS)

# ---------------- Part A ----------------
def partA():
    # sanity: sgn == sortsign on random vectors
    for _ in range(2000):
        v = tuple(random.randint(-3, 3) for _ in range(3))
        assert sgn(v) == sortsign(v), v
    for BOX in (16, 30):
        cnt = bad = 0
        for A in range(1, BOX + 1):
            for B in range(1, BOX + 1):
                a = (A + B, B, 0)
                for p in range(0, BOX + 1):
                    for q in range(0, p + 1):
                        cnt += 1
                        if lhs(a, (p, q, 0)) != rhs(a, (p, q, 0)):
                            bad += 1
        print("A: box %d: %d points, %d failures" % (BOX, cnt, bad))
    random.seed(1)
    bad = 0
    for _ in range(20000):
        a = sorted((random.randint(-10**6, 10**6) for _ in range(3)), reverse=True)
        if len(set(a)) < 3:
            continue
        b = tuple(random.randint(-10**6, 10**6) for _ in range(3))
        if lhs(a, b) != rhs(a, b):
            bad += 1
    print("A: random large (un-normalised a,b): failures", bad)
    bad = 0; cnt = 0
    for a in itertools.product(range(-6, 7), repeat=3):
        if not (a[0] > a[1] > a[2]):
            continue
        for b in itertools.product(range(-4, 5), repeat=3):
            cnt += 1
            if lhs(a, b) != rhs(a, b):
                bad += 1
                if bad < 4:
                    print("   FAIL", a, b, lhs(a, b), rhs(a, b))
    print("A: exhaustive a in [-6,6]^3 strictly decreasing, b in [-4,4]^3 (ties/unsorted): %d cases, %d failures" % (cnt, bad))

# ---------------- Part B ----------------
def forms():
    """all linear forms in (A,B,p,q) whose signs appear as factors, plus A,B,p,q,p-q,A+B. Return set of primitive tuples up to sign."""
    a = ((1, 0, 0, 0), (0, 1, 0, 0), (1, 1, 0, 0))  # a1-a2 = A, a2-a3 = B, a1-a3 = A+B  (as (A,B,p,q) coefficient vectors)
    # vector representation of a = (A+B, B, 0) and b = (p,q,0) as linear maps
    avec = ((1, 1, 0, 0), (0, 1, 0, 0), (0, 0, 0, 0))
    bvec = ((0, 0, 1, 0), (0, 0, 0, 1), (0, 0, 0, 0))
    H = set()
    def norm(f):
        f = tuple(f)
        if all(c == 0 for c in f):
            return None
        g = 0
        for c in f: g = gcd(g, abs(c))
        f = tuple(c // g for c in f)
        for c in f:
            if c != 0:
                return f if c > 0 else tuple(-x for x in f)
    for s in PERMS:
        for sign in (1, -1):
            v = [tuple(avec[i][t] + sign * bvec[s[i]][t] for t in range(4)) for i in range(3)]
            for i in range(3):
                for j in range(i + 1, 3):
                    f = tuple(v[i][t] - v[j][t] for t in range(4))
                    nf = norm(f)
                    assert nf is not None
                    assert all(abs(c) <= 1 for c in f), f  # coefficients in {0,+-1} as claimed
                    H.add(nf)
    for f in ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (0, 0, 1, -1), (1, 1, 0, 0)):
        H.add(f)
    return sorted(H)

def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))

def kernel_ray(f1, f2, f3):
    M = [f1, f2, f3]
    r = []
    for k in range(4):
        cols = [c for c in range(4) if c != k]
        sub = [[M[i][c] for c in cols] for i in range(3)]
        r.append((-1) ** k * det3(sub))
    if all(x == 0 for x in r):
        return None
    g = 0
    for x in r: g = gcd(g, abs(x))
    return tuple(x // g for x in r)

def dot(f, x):
    return sum(f[i] * x[i] for i in range(4))

def signvec(H, x):
    return tuple((dot(f, x) > 0) - (dot(f, x) < 0) for f in H)

def in_region_sv(H, sv):
    idx = {f: i for i, f in enumerate(H)}
    return (sv[idx[(1, 0, 0, 0)]] == 1 and sv[idx[(0, 1, 0, 0)]] == 1
            and sv[idx[(0, 0, 1, -1)]] >= 0 and sv[idx[(0, 0, 0, 1)]] >= 0)

def eval_point(x):
    A, B, p, q = x
    a = (A + B, B, 0); b = (p, q, 0)
    return lhs(a, b), rhs(a, b)

def partB():
    H = forms()
    print("B: |H| =", len(H), H)
    rays = set()
    maxabs = 0
    for f1, f2, f3 in itertools.combinations(H, 3):
        r = kernel_ray(f1, f2, f3)
        if r is None:
            continue
        maxabs = max(maxabs, max(abs(x) for x in r))
        rays.add(r); rays.add(tuple(-x for x in r))
    rays = sorted(rays)
    print("B: candidate rays (both signs): %d, max |coord| = %d (claimed <= 4)" % (len(rays), maxabs))
    # restrict to rays in the CLOSED region A>=0,B>=0,p>=q>=0 (extreme rays of a cell in the region lie in its closure)
    rays = [r for r in rays if r[0] >= 0 and r[1] >= 0 and r[2] >= r[3] >= 0]
    print("B: candidate rays in closed region:", len(rays))
    sv = [signvec(H, r) for r in rays]
    m = len(rays)
    compat = [[all(sv[i][t] * sv[j][t] >= 0 for t in range(len(H))) for j in range(m)] for i in range(m)]
    cells = {}
    # enumerate compatible subsets of size <= 4 (cliques in compat graph)
    def rec(start, chosen, point):
        if chosen:
            s = signvec(H, point)
            if in_region_sv(H, s) and s not in cells:
                cells[s] = tuple(point)
        if len(chosen) == 4:
            return
        for i in range(start, m):
            if all(compat[i][j] for j in chosen):
                rec(i + 1, chosen + [i], tuple(point[t] + rays[i][t] for t in range(4)))
    rec(0, [], (0, 0, 0, 0))
    print("B: cells meeting the region (from ray sums):", len(cells))
    maxc = max(max(abs(c) for c in pt) for pt in cells.values())
    print("B: max |coordinate| of the ray-sum representatives:", maxc, "(claimed <= 16)")
    # sign vectors realised in the box
    box = {}
    for A in range(1, 17):
        for B in range(1, 17):
            for p in range(0, 17):
                for q in range(0, p + 1):
                    s = signvec(H, (A, B, p, q))
                    box.setdefault(s, (A, B, p, q))
    print("B: sign vectors realised in the box 16:", len(box))
    missing = [s for s in cells if s not in box]
    extra = [s for s in box if s not in cells]
    print("B: cells missing from box: %d ; box sign vectors not among cells: %d" % (len(missing), len(extra)))
    for s in missing[:5]:
        print("   MISSING cell, representative", cells[s])
    # S3 on each cell representative and constancy check: value at ray-sum rep equals value at box rep
    bad = 0
    for s, pt in cells.items():
        l1, r1 = eval_point(pt)
        if l1 != r1:
            bad += 1
        if s in box:
            l2, r2 = eval_point(box[s])
            if (l1, r1) != (l2, r2):
                print("   NOT CONSTANT on cell?", pt, box[s], (l1, r1), (l2, r2)); bad += 1
    print("B: S3 failures on cell representatives / constancy violations:", bad)
    # also a smaller box: what is the smallest box that already covers all cells?
    for BX in (4, 6, 8, 10, 12, 14):
        seen = set()
        for A in range(1, BX + 1):
            for B in range(1, BX + 1):
                for p in range(0, BX + 1):
                    for q in range(0, p + 1):
                        seen.add(signvec(H, (A, B, p, q)))
        print("B: box %d realises %d of %d cells" % (BX, len([s for s in cells if s in seen]), len(cells)))

# ---------------- Part C ----------------
def Nchains(nu, j):
    from functools import lru_cache
    @lru_cache(maxsize=None)
    def rec(k, t):
        if t == 0:
            return 1
        tot = 0
        for i in range(3):
            if i == 0 or k[i - 1] > k[i]:
                tot += rec(k[:i] + (k[i] + 1,) + k[i + 1:], t - 1)
        return tot
    return rec(tuple(nu), j)

def partC():
    bad = cnt = 0
    for n1 in range(0, 11):
        for n2 in range(0, n1 + 1):
            for n3 in range(0, n2 + 1):
                nu = (n1, n2, n3); dual = (n1 - n3, n1 - n2, 0)
                for j in range(0, 25):
                    cnt += 1
                    if Nchains(nu, j) != Nchains(dual, j):
                        bad += 1
    print("C: chain symmetry nu_1<=10, j<=24: %d cases, %d failures" % (cnt, bad))

# ---------------- Part D ----------------
def contents(lam):
    return [c - r for r, row in enumerate(lam) for c in range(row)]

def partD():
    for K in range(1, 21):
        u = d_vector((3,) * K + (2, 2, 2)); v = d_vector((3,) * K + (1, 1, 1))
        first = next(j for j in range(len(u)) if u[j] != v[j])
        ok = all(u[j] == v[j] for j in range(K + 1))
        red = all(u[j] == Nchains((3, 0, 0), j) and v[j] == Nchains((3, 3, 0), j) for j in range(K + 1))
        # transposes as in the proof
        assert conjugate((3,) * K + (2, 2, 2)) == (K + 3, K + 3, K) and conjugate((3,) * K + (1, 1, 1)) == (K + 3, K, K)
        print("D: R_%d: holds=%s first diff=%d (K+2=%d) reduction=%s" % (K, ok, first, K + 2, red))
    for k in range(1, 11):
        lam = (5,) + (3,) * k + (2, 2, 2); mu = (4, 4) + (3,) * k + (1, 1, 1); n = 3 * k + 11
        u, v = d_vector(lam), d_vector(mu)
        D = [j for j in range(n + 1) if u[j] != v[j]]
        # Lemma B formula
        x, y = d_vector((3,) * k + (2, 2, 2)), d_vector((3,) * k + (1, 1, 1))
        lb = all(u[j] - v[j] == sum(comb(j, i) * (x[j - i] - y[j - i]) for i in range(0, min(j, 2) + 1)) for j in range(k + 1))
        lb2 = all(u[j] == sum(comb(j, i) * x[j - i] for i in range(0, min(j, 2) + 1)) for j in range(k + 1))
        c1l, c1m = sum(contents(lam)), sum(contents(mu))
        c2l, c2m = sum(c * c for c in contents(lam)), sum(c * c for c in contents(mu))
        print("D: H_%d n=%d diff idx == {k+4..n-4}: %s (first %d last %d) LemmaB-diff=%s LemmaB-abs=%s f_eq=%s C2_eq=%s C1: %d %d" % (
            k, n, D == list(range(k + 4, n - 3)), D[0], D[-1], lb, lb2, u[n] == v[n], c2l == c2m, c1l, c1m))

if __name__ == "__main__":
    partA(); partB(); partC(); partD()
