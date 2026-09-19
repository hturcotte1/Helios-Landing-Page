"""ref2_skeptic_S3.py -- referee #2, independent recheck of Theorem S3 (skeptic2.md Sec.3) and its cell argument.
Exact integer arithmetic."""
import sys, random
from itertools import permutations, combinations
from math import gcd
from fractions import Fraction

PERMS = list(permutations(range(3)))

def sgn(v):
    s = 1
    for i in range(3):
        for j in range(i+1, 3):
            d = v[i]-v[j]
            if d == 0: return 0
            if d < 0: s = -s
    return s

def lhs_rhs(a, b):
    plus = sum(sgn(tuple(a[i]+b[s[i]] for i in range(3))) for s in PERMS)
    minus = sum(sgn(tuple(a[i]-b[s[i]] for i in range(3))) for s in PERMS)
    return plus, minus

# 1. general (unreduced) S3 on random a strictly decreasing, arbitrary b (no reduction assumptions)
random.seed(12345)
bad = 0
for _ in range(200000):
    R = random.choice([5, 20, 100, 10**6])
    a = sorted(random.sample(range(-R, R), 3), reverse=True)
    b = [random.randint(-R, R) for _ in range(3)]
    p, m = lhs_rhs(a, b)
    if p != m: bad += 1; print("S3 general FAIL", a, b, p, m)
print("S3 unreduced random test: 200000 samples, failures =", bad)

# 2. exhaustive check in reduced coordinates, box 16 and box 40
def box_check(BOX):
    bad = cnt = 0
    for A in range(1, BOX+1):
        for B in range(1, BOX+1):
            a = (A+B, B, 0)
            for p in range(0, BOX+1):
                for q in range(0, p+1):
                    cnt += 1
                    pl, mi = lhs_rhs(a, (p, q, 0))
                    if pl != mi: bad += 1
    return cnt, bad
for BOX in (16, 30):
    cnt, bad = box_check(BOX)
    print("S3 reduced box %d: %d points, failures = %d" % (BOX, cnt, bad))

# 3. the arrangement: forms in H (coefficients on (A,B,p,q))
forms = set()
adiff = {(0,1): (1,0,0,0), (1,2): (0,1,0,0), (0,2): (1,1,0,0)}   # a_i - a_j for a=(A+B,B,0)
bvec = [(0,0,1,0), (0,0,0,1), (0,0,0,0)]                          # b = (p,q,0)
for s in PERMS:
    for (i,j), av in adiff.items():
        bd = tuple(bvec[s[i]][t]-bvec[s[j]][t] for t in range(4))
        for sign in (1, -1):
            forms.add(tuple(av[t]+sign*bd[t] for t in range(4)))
forms |= {(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(0,0,1,-1),(1,1,0,0)}
forms = sorted(forms)
print("number of distinct forms in H:", len(forms), "all coefficients in {0,+-1}:",
      all(c in (-1,0,1) for f in forms for c in f))

def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
            + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

# extreme-ray candidates: kernels of rank-3 triples of forms
rays = set(); maxabs = 0
for tri in combinations(forms, 3):
    cof = []
    for c in range(4):
        M = [[tri[r][t] for t in range(4) if t != c] for r in range(3)]
        cof.append((-1)**c * det3(M))
    if all(x == 0 for x in cof): continue
    g = 0
    for x in cof: g = gcd(g, abs(x))
    r = tuple(x//g for x in cof)
    maxabs = max(maxabs, max(abs(x) for x in r))
    rays.add(r); rays.add(tuple(-x for x in r))
print("rank-3 triples -> %d primitive ray directions, max |entry| = %d (claim: <= 4)" % (len(rays), maxabs))

def signvec(x):
    return tuple((f[0]*x[0]+f[1]*x[1]+f[2]*x[2]+f[3]*x[3] > 0) - (f[0]*x[0]+f[1]*x[1]+f[2]*x[2]+f[3]*x[3] < 0) for f in forms)
def in_region(x):
    A,B,p,q = x
    return A > 0 and B > 0 and p >= q >= 0

cells16 = set()
for A in range(1, 17):
    for B in range(1, 17):
        for p in range(0, 17):
            for q in range(0, p+1):
                cells16.add(signvec((A,B,p,q)))
print("distinct cells (sign vectors) met by box-16 points:", len(cells16))
# cells from sums of <=4 rays lying in the region (the proof's construction), and from random far points
cells_rays = set()
raylist = sorted(rays)
for k in range(1, 5):
    for sub in combinations(raylist, k):
        y = tuple(sum(r[t] for r in sub) for t in range(4))
        if in_region(y): cells_rays.add(signvec(y))
print("cells from sums of <=4 ray directions in region:", len(cells_rays), " subset of box16 cells:", cells_rays <= cells16)
cells_big = set()
for _ in range(300000):
    R = random.choice([50, 1000, 10**7])
    x = (random.randint(1,R), random.randint(1,R), 0, 0)
    p = random.randint(0,R); q = random.randint(0,p)
    x = (x[0], x[1], p, q)
    cells_big.add(signvec(x))
print("cells from 300000 random far points:", len(cells_big), " subset of box16 cells:", cells_big <= cells16)
# also check both sides are indeed constant on cells (by sign vector) over all points seen
table = {}
bad = 0
for A in range(1, 25):
    for B in range(1, 25):
        for p in range(0, 25):
            for q in range(0, p+1):
                sv = signvec((A,B,p,q)); val = lhs_rhs((A+B,B,0),(p,q,0))
                if table.setdefault(sv, val) != val: bad += 1
print("both sides constant on cells over box 24: violations =", bad, " cells seen:", len(table))
