"""ref2_skeptic_s3_cells.py -- referee #2: independent reconstruction of the hyperplane-arrangement argument for S3.
Coordinates x = (A, B, p, q); a = (A+B, B, 0), b = (p, q, 0)."""
import itertools, math
from fractions import Fraction
PERMS = list(itertools.permutations(range(3)))
A_ = (1,0,0,0); B_ = (0,1,0,0); P_ = (0,0,1,0); Q_ = (0,0,0,1)
a_forms = [tuple(x+y for x,y in zip(A_,B_)), B_, (0,0,0,0)]   # a_i as linear forms
b_forms = [P_, Q_, (0,0,0,0)]
def sub(u,v): return tuple(x-y for x,y in zip(u,v))
def add(u,v): return tuple(x+y for x,y in zip(u,v))
def neg(u): return tuple(-x for x in u)
def ev(l, x): return sum(c*t for c,t in zip(l,x))
# forms whose signs determine each term: (a +- sigma b)_i - (a +- sigma b)_j, i<j
terms = []   # (eps, sigma) -> list of 3 forms
H = set()
for eps in (1, -1):
    for s in PERMS:
        fl = []
        for i in range(3):
            for j in range(i+1, 3):
                l = add(sub(a_forms[i], a_forms[j]), tuple(eps*c for c in sub(b_forms[s[i]], b_forms[s[j]])))
                fl.append(l)
        terms.append((eps, s, fl))
        H.update(fl)
H.update([A_, B_, P_, Q_, sub(P_, Q_), add(A_, B_)])
# normalise sign of forms (l and -l define the same hyperplane)
def canon(l):
    for c in l:
        if c != 0: return l if c > 0 else neg(l)
    return l
H = sorted(set(canon(l) for l in H if any(l)))
assert all(c in (-1,0,1) for l in H for c in l), "coefficients not in {0,+-1}"
print("|H| =", len(H), "all coefficients in {0,+-1}: True")

def sgnvec(x): return tuple((ev(l,x) > 0) - (ev(l,x) < 0) for l in H)
def sgn3(v):
    s = 1
    for i in range(3):
        for j in range(i+1,3):
            d = v[i]-v[j]
            if d == 0: return 0
            if d < 0: s = -s
    return s
def sides(x):
    A,B,p,q = x; a = (A+B,B,0); b = (p,q,0)
    L = sum(sgn3([a[i]+b[s[i]] for i in range(3)]) for s in PERMS)
    R = sum(sgn3([a[i]-b[s[i]] for i in range(3)]) for s in PERMS)
    return (L, R)
def sides_from_signs(sv):
    d = dict(zip(H, sv))
    def sg(l):
        c = canon(l); s = d[c]; return s if c == l else -s
    L = R = 0
    for eps, s, fl in terms:
        prod = 1
        for l in fl: prod *= sg(l)
        if eps == 1: L += prod
        else: R += prod
    return (L, R)

# 1. max |det| of 3x3 {0,+-1} matrices, brute force
mx = 0
def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
for ent in itertools.product((-1,0,1), repeat=9):
    m = [ent[0:3], ent[3:6], ent[6:9]]
    mx = max(mx, abs(det3(m)))
print("max |det| over all 3^9 {0,+-1} 3x3 matrices:", mx)

# 2. candidate extreme rays: kernels of rank-3 triples of forms in H (vector of signed 3x3 minors)
def minors(rows):
    v = []
    for k in range(4):
        cols = [c for c in range(4) if c != k]
        m = [[r[c] for c in cols] for r in rows]
        v.append((-1)**k * det3(m))
    return tuple(v)
rays = set(); maxc = 0
for tri in itertools.combinations(H, 3):
    v = minors(tri)
    if any(v):
        assert all(ev(l, v) == 0 for l in tri)
        g = math.gcd(*[abs(c) for c in v]); v = tuple(c//g for c in v)
        maxc = max(maxc, max(abs(c) for c in v))
        rays.add(v); rays.add(neg(v))
print("candidate primitive rays (both signs):", len(rays), " max |coordinate|:", maxc, "(proof needs <= 4)")
closed_region = lambda x: x[0] >= 0 and x[1] >= 0 and x[3] >= 0 and x[2] >= x[3]
R = sorted(r for r in rays if closed_region(r))
print("rays in closed region {A,B>=0, p>=q>=0}:", len(R), R)

# 3. rank helper
def rank(vs):
    M = [[Fraction(c) for c in v] for v in vs]; r = 0; ncol = 4
    for c in range(ncol):
        piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]
                M[i] = [x - f*y for x, y in zip(M[i], M[r])]
        r += 1
    return r
# sign vectors of sums of <= 4 linearly independent rays from R that fall in the OPEN region (A>0,B>0,p>=q>=0):
# by the proof, this superset contains a representative of every cell meeting the region.
cells_from_rays = {}
maxrep = 0
for k in range(1, 5):
    for S in itertools.combinations(R, k):
        if rank(S) < k: continue
        y = tuple(sum(r[i] for r in S) for i in range(4))
        if y[0] > 0 and y[1] > 0:
            sv = sgnvec(y)
            cells_from_rays.setdefault(sv, y); maxrep = max(maxrep, max(abs(c) for c in y))
print("distinct sign vectors from ray sums (superset of cells in region):", len(cells_from_rays), " max |coord| of reps:", maxrep)

# 4. sign vectors realised in box 16, and constancy of both sides on each sign vector
BOX = 16
box_cells = {}; inconsistent = 0; fails = 0
for A in range(1, BOX+1):
    for B in range(1, BOX+1):
        for p in range(0, BOX+1):
            for q in range(0, p+1):
                x = (A,B,p,q); sv = sgnvec(x); val = sides(x)
                if val[0] != val[1]: fails += 1
                if sv in box_cells:
                    if box_cells[sv] != val: inconsistent += 1
                else:
                    box_cells[sv] = val
                    assert sides_from_signs(sv) == val   # both sides are functions of the sign vector
print("box 16: distinct sign vectors:", len(box_cells), " S3 failures:", fails, " sign-vector/value inconsistencies:", inconsistent)
missing = [sv for sv in cells_from_rays if sv not in box_cells]
extra = [sv for sv in box_cells if sv not in cells_from_rays]
print("ray-sum cells missing from box:", len(missing), " box cells not among ray sums:", len(extra))
# every ray-sum representative also satisfies S3 (evaluate directly)
print("S3 on ray-sum representatives: failures", sum(1 for sv, y in cells_from_rays.items() if sides(y)[0] != sides(y)[1]))
# 5. smallest box that already realises all cells
for bx in range(2, 17):
    seen = set()
    for A in range(1, bx+1):
        for B in range(1, bx+1):
            for p in range(0, bx+1):
                for q in range(0, p+1):
                    seen.add(sgnvec((A,B,p,q)))
    print("  box %d realises %d of %d cells" % (bx, len(seen), len(box_cells)))
    if len(seen) == len(box_cells): break
# 6. random rational points far outside the box: sign vector must be among the box cells
import random; random.seed(7); unseen = 0
for t in range(300000):
    A = random.randint(1, 10**5); B = random.randint(1, 10**5); p = random.randint(0, 10**5); q = random.randint(0, p)
    if random.random() < 0.5:  # force coincidences to hit low-dim cells
        A, B, p, q = [random.choice([A, B, p, q, A+B, A+p, B+q, abs(A-p), abs(B-q), p-q]) for _ in range(4)]
        A = max(A, 1); B = max(B, 1); p, q = max(p, q), min(p, q)
    if sgnvec((A,B,p,q)) not in box_cells: unseen += 1
print("random far points with sign vector not realised in box 16:", unseen)
