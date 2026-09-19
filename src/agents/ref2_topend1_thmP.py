"""ref2_topend1: test Theorem P (i)-(iii) directly. For each mu |- m (m<=6) and each n<=NMAX, compute the class
expansion p_mu(J) = sum_rho' a_{mu,rho'}(n) K_rho'(n) via a_{mu,rho'}(n) = (1/n!) sum_lam f^lam chi^lam(rho') C_mu(lam)
(central element z = sum_lam omega_lam(z) e_lam with omega_lam(p_mu(J)) = C_mu(lam) by the JM theorem).
Then check: (ii) a=0 unless m(rho')<=m, and the values for n>=|rho'| are a polynomial in n of degree <= m-m(rho');
(iii) for m(rho')=m constant, zero unless l(rho')<=l(mu), and a_{mu,mu+1} = prod m_k(mu)!.
Cross-check for n<=7: brute-force group algebra product."""
import sys, time
from fractions import Fraction
from math import factorial, prod
from collections import Counter
from itertools import product
from ref2_topend1_lib import *

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 13
MMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 6

def coeffs(mu, n):
    """dict rho' (no 1's) -> a_{mu,rho'}(n) as Fraction, for all classes of S_n"""
    lams = list(partitions(n)); res = {}
    Cmu = {lam: prod(C(lam, k) for k in mu) for lam in lams}
    fl = {lam: f_hook(lam) for lam in lams}
    for rho in lams:
        red = tuple(x for x in rho if x > 1)
        s = sum(fl[lam] * chi(lam, rho) * Cmu[lam] for lam in lams)
        res[red] = Fraction(s, factorial(n))
    return res

def poly_fit_check(pts, d):
    """pts: list of (n, value); check that they lie on a polynomial of degree <= d (Newton differences)."""
    if len(pts) <= d + 1: return None  # cannot test
    xs = [Fraction(x) for x, _ in pts]; ys = [Fraction(y) for _, y in pts]
    # divided differences; degree <= d iff all (d+1)-th divided differences vanish
    dd = ys[:]
    for k in range(1, len(xs)):
        dd = [(dd[i + 1] - dd[i]) / (xs[i + k] - xs[i]) for i in range(len(dd) - 1)]
        if k == d + 1:
            return all(v == 0 for v in dd)
    return True

# brute-force group algebra for small n
def compose(p, q):  # (pq)(x) = p(q(x))
    return tuple(p[q[x]] for x in range(len(p)))
def transp(n, i, k):
    p = list(range(n)); p[i], p[k] = p[k], p[i]; return tuple(p)
def cycle_type(p):
    n = len(p); seen = [False] * n; ct = []
    for x in range(n):
        if not seen[x]:
            l = 0; y = x
            while not seen[y]: seen[y] = True; y = p[y]; l += 1
            ct.append(l)
    return tuple(sorted(ct, reverse=True))
def ga_mul(A, B):
    R = Counter()
    for a, ca in A.items():
        for b, cb in B.items():
            R[compose(a, b)] += ca * cb
    return R
def J(n, k):  # k is 1-indexed centre; points 0-indexed
    return Counter({transp(n, i, k - 1): 1 for i in range(k - 1)})
def pm(n, m):
    R = Counter()
    for k in range(1, n + 1):
        Jk = J(n, k); P = Counter({tuple(range(n)): 1})
        for _ in range(m): P = ga_mul(P, Jk)
        R.update(P)
    return R
def brute_coeffs(mu, n):
    P = Counter({tuple(range(n)): 1})
    for part in mu: P = ga_mul(P, pm(n, part))
    res = {}
    for g, c in P.items():
        red = tuple(x for x in cycle_type(g) if x > 1)
        if red in res: assert res[red] == c, ("not central!", mu, n, red)
        else: res[red] = c
    # classes absent get 0
    for rho in partitions(n):
        red = tuple(x for x in rho if x > 1)
        res.setdefault(red, 0)
    return res

t0 = time.time()
mus = [mu for m in range(1, MMAX + 1) for mu in partitions(m)]
data = {mu: {n: coeffs(mu, n) for n in range(0, NMAX + 1)} for mu in mus}
print(f"computed class expansions for {len(mus)} mu, n<={NMAX} [{time.time()-t0:.0f}s]", flush=True)
# brute cross-check n<=6 (n=7 for m<=3)
nb = 0
for mu in mus:
    for n in range(1, 7 if sum(mu) <= 4 else 6):
        bc = brute_coeffs(mu, n)
        for red, v in bc.items():
            assert data[mu][n].get(red, 0) == v, ("brute mismatch", mu, n, red, v, data[mu][n].get(red))
            nb += 1
print(f"brute-force group algebra cross-check OK ({nb} coefficients) [{time.time()-t0:.0f}s]", flush=True)

problems = []
for mu in mus:
    m = sum(mu); l = len(mu)
    rhos = set(r for n in data[mu] for r in data[mu][n])
    for rho in sorted(rhos, key=lambda r: (sum(r), r)):
        s = sum(rho); mr = s - len(rho)
        pts = [(n, data[mu][n][rho]) for n in range(s, NMAX + 1) if rho in data[mu][n]]
        vals = [v for _, v in pts]
        if mr > m:
            if any(v != 0 for v in vals): problems.append(("(ii) nonzero with m(rho')>m", mu, rho, vals))
            continue
        d = m - mr
        ok = poly_fit_check(pts, d)
        if ok is None: problems.append(("untestable (too few n)", mu, rho, len(pts), d))
        elif not ok: problems.append(("(i)/(ii) degree bound fails", mu, rho, vals))
        if mr == m:
            if len(set(vals)) != 1: problems.append(("(iii) not constant", mu, rho, vals))
            if len(rho) > l and any(v != 0 for v in vals): problems.append(("(iii) nonzero with l(rho')>l", mu, rho, vals))
    mu1 = tuple(x + 1 for x in mu)
    want = prod(factorial(c) for c in Counter(mu).values())
    got = [data[mu][n].get(mu1) for n in range(sum(mu1), NMAX + 1)]
    if any(v != want for v in got): problems.append(("(iii) a_{mu,mu+1}", mu, want, got))
print("problems:", problems if problems else "NONE", f"[{time.time()-t0:.0f}s]")
# also print the full expansions for a couple of mu at symbolic level (fit polynomials)
from sympy import interpolate, symbols, expand
nn = symbols('n')
for mu in [(1,), (2,), (1, 1), (3,), (4,), (2, 2), (3, 3), (5, 1)]:
    if mu not in data: continue
    m = sum(mu)
    out = []
    for rho in sorted(set(r for n_ in data[mu] for r in data[mu][n_]), key=lambda r: (sum(r), r)):
        s = sum(rho); mr = s - len(rho)
        if mr > m: continue
        d = m - mr
        pts = [(n_, data[mu][n_][rho]) for n_ in range(s, s + d + 1)]
        P = expand(interpolate([(x, Fraction(y)) for x, y in pts], nn)) if pts else 0
        if P != 0: out.append(f"[{P}] K_{rho}")
    print(f"p_{mu}(J) =", " + ".join(out))
