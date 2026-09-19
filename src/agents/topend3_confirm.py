"""Exact confirmation and classification of the candidate pairs found by topend3_scan (C++, floating log f):
for each candidate pair (lam, mu) (distinct transpose classes, equal (n, C2, C1^2, C4), |dlog f| < 1e-7):
  * exact f^lam == f^mu (hook formula)   -> depth-5 collision (equal n, f, u_3, u_4, u_5)
  * h7 = C6 - 16 C1 C3 equal?            -> depth-7 (equal u_7 too)  [h7 is transpose-invariant: C1C3 even]
  * h8 = 2 C3^2 - (8n+164) C1 C3 equal?  -> depth-8
  * h9, h10 equal?                        -> depth-9, depth-10
  * r equal?  mirror type (equal multiset of |content|)?  content multiset difference.
For every pair of depth >= 7, u_6, u_7, u_8 are recomputed EXACTLY and independently as sum_{nu |- i} f^{lam/nu}
(Aitken determinants) and compared; this does not use the polynomials at all.
"""
import sys, re, time
from collections import Counter, defaultdict
from fractions import Fraction
from topend3_lib import *
from young import removable_corners
t0 = time.time()
lines = [l for l in open(sys.argv[1] if len(sys.argv) > 1 else 'topend3_scan.log') if l.startswith('CAND')]
def parse(l):
    m = re.match(r'CAND n=(\d+) \[([\d,]*)\] \[([\d,]*)\]', l)
    n = int(m.group(1)); a = tuple(int(x) for x in m.group(2).split(',')); b = tuple(int(x) for x in m.group(3).split(','))
    return n, a, b
def moments(lam):
    cs = contents(lam)
    return {k: sum(c ** k for c in cs) for k in (1, 2, 3, 4, 5, 6, 8)}
def hs(lam):
    n = sum(lam); M = moments(lam)
    C1, C2, C3, C4, C5, C6, C8 = (M[k] for k in (1, 2, 3, 4, 5, 6, 8))
    h7 = C6 - 16 * C1 * C3
    h8 = 2 * C3 * C3 - (8 * n + 164) * C1 * C3
    h9 = Fraction(C8) - 24 * C1 * C5 - (32 * n - Fraction(7720, 3)) * C1 * C3
    h10 = -24 * C1 * C2 * C3 + (12 * n * n + 1060 * n + 62664) * C1 * C3 - 360 * C1 * C5
    return (h7, h8, h9, h10)
def u_aitken(lam, i):
    return sum(f_skew_aitken(lam, nu) for nu in partitions(i))
stats = defaultdict(Counter)
byn = defaultdict(list)
for l in lines:
    n, a, b = parse(l)
    fa, fb = f_hook(a), f_hook(b)
    if fa != fb:
        stats[n]['float-only (f differs)'] += 1; continue
    Ma, Mb = moments(a), moments(b)
    assert Ma[2] == Mb[2] and Ma[1] ** 2 == Mb[1] ** 2 and Ma[4] == Mb[4]
    ha, hb = hs(a), hs(b)
    depth = 5
    if ha[0] == hb[0]: depth = 7
    if depth == 7 and ha[1] == hb[1]: depth = 8
    if depth == 8 and ha[2] == hb[2]: depth = 9
    if depth == 9 and ha[3] == hb[3]: depth = 10
    ra, rb = len(removable_corners(a)), len(removable_corners(b))
    ca, cb = Counter(contents(a)), Counter(contents(b))
    mirror = Counter(abs(c) for c in contents(a)) == Counter(abs(c) for c in contents(b))
    removed = sorted((ca - cb).elements()); added = sorted((cb - ca).elements())
    stats[n][f'depth{depth}'] += 1
    stats[n][f'depth{depth} equal r'] += (ra == rb)
    stats[n][f'depth{depth} mirror'] += mirror
    byn[n].append((depth, a, b, ra, rb, mirror, removed, added))
    extra = ''
    if depth >= 7:
        # independent exact check of u_6, u_7, u_8
        ua = [u_aitken(a, i) for i in (6, 7, 8)]; ub = [u_aitken(b, i) for i in (6, 7, 8)]
        extra = f'  [Aitken u6,u7,u8: {ua} vs {ub}]'
    print(f"n={n} depth={depth} r={ra},{rb} mirror={mirror} {a} ~ {b}  C1={Ma[1]},{Mb[1]} removed={removed} added={added}{extra}", flush=True)
print()
for n in sorted(stats):
    print(f"n={n}: " + ", ".join(f"{k}={v}" for k, v in sorted(stats[n].items())))
print(f"[{time.time()-t0:.0f}s]")
