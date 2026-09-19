"""ref1_skeptic_c123.py -- referee #1 (second pass) independent check of skeptic2.md Section 3 (S3, C1, C2, C3, H_k).

All code here is written from scratch (no import of census/young/skeptic2 modules).

A. Theorem C1 tested DIRECTLY on GL_3 characters (not via S3): for lambda dominant (possibly negative entries) and W = V_mu
   (mu dominant, possibly negative), count constituents of V_lambda (x) W and V_lambda (x) W^* by (i) multiplying Weyl
   numerators A_{lambda+rho} * ch(W) as Laurent polynomials and reading off strictly-decreasing exponents, (ii) the
   Racah-Speiser signed sum with own sgn. Both counts for W and W^* must agree.
B. Chain symmetry N_j(nu) = N_j(nu*) for nu_1 <= 14, j <= 30 (own chain counter), and N_j(nu) == #const(V_nu (x) V^{(x)j}).
C. R_K for K <= 30 with own d-vector code (sub-diagram recursion) and the box-complement reduction d_j = N_j((3)), N_j((3,3)).
D. Top-end formulas (needed for 'd_{n-4} differs'): on all partitions of n <= 16,
     d_{n-3} = f*(2/3 + (C_2 - n(n-1)/2)/(n)_3),
     d_{n-4} = f*(5/12 + (C_1^2 - 3 C_2 + n(n-1))/(n)_4 + (C_2 - n(n-1)/2)/(n)_3)
   (exact rationals), so d_{n-4} is strictly increasing in C_1^2 given (n, f, C_2).
E. Theorem H(i),(ii): hook multisets equal, C_2 equal, C_1(lam) - C_1(mu) = 2, C_1(lam) = -(3k^2+9k-2)/2, for k <= 40.
F. H_k: differing indices, Lemma B formula, d_{n-4} differs, for k <= 14 (own d-vector code).
G. Independent enumeration of cells of the S3 arrangement via random+structured sampling and exact LP-free method:
   sign vectors realised by all integer points in box 24 vs box 16 (no new cells beyond box 8 expected).
"""
import os
import sys, itertools, random
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
from collections import Counter, defaultdict

OUT = []
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

PERMS = list(itertools.permutations(range(3)))
def perm_sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i+1, len(p)):
            if p[i] > p[j]: s = -s
    return s

# ---------- Laurent polynomials in 3 variables as dict exponent-tuple -> coeff ----------
def pmul(P, Q):
    R = defaultdict(int)
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            R[(e1[0]+e2[0], e1[1]+e2[1], e1[2]+e2[2])] += c1*c2
    return {e: c for e, c in R.items() if c}

def A(gamma):
    """Weyl numerator sum_w sgn(w) x^{w gamma}."""
    P = defaultdict(int)
    for p in PERMS:
        P[(gamma[p[0]], gamma[p[1]], gamma[p[2]])] += perm_sign(p)
    return {e: c for e, c in P.items() if c}

RHO = (2, 1, 0)

@lru_cache(maxsize=None)
def char_GL3(lam):
    """character of V_lam (lam dominant in Z^3) as Laurent polynomial, via semistandard tableaux / Gelfand-Tsetlin patterns
    (shift to nonnegative, use GT patterns, shift back)."""
    m = lam[2]
    l = (lam[0]-m, lam[1]-m, lam[2]-m)
    P = defaultdict(int)
    # GT patterns: top row l; middle row (a,b) with l0>=a>=l1>=b>=l2; bottom c with a>=c>=b
    for a in range(l[1], l[0]+1):
        for b in range(l[2], l[1]+1):
            for c in range(b, a+1):
                w1 = c; w2 = a + b - c; w3 = sum(l) - a - b
                P[(w1+m, w2+m, w3+m)] += 1
    return dict(P)

def decompose(P):
    """decompose a virtual character P (Laurent poly, S_3-symmetric) into irreducibles: multiply by A_rho and read
    strictly decreasing exponents. Returns dict dominant weight -> multiplicity."""
    Q = pmul(P, A(RHO))
    res = {}
    for e, c in Q.items():
        if e[0] > e[1] > e[2]:
            res[(e[0]-2, e[1]-1, e[2])] = c
    return res

def sgn3(v):
    s = 1
    for i in range(3):
        for j in range(i+1, 3):
            if v[i] == v[j]: return 0
            if v[i] < v[j]: s = -s
    return s

def dual_weight(mu):
    return (-mu[2], -mu[1], -mu[0])

def partA():
    random.seed(7)
    bad = 0; cnt = 0
    doms = [l for l in itertools.product(range(-3, 6), repeat=3) if l[0] >= l[1] >= l[2]]
    mus = [m for m in itertools.product(range(-2, 5), repeat=3) if m[0] >= m[1] >= m[2]]
    for lam in doms:
        chl = char_GL3(lam)
        for mu in mus:
            chw = char_GL3(mu)
            chwd = char_GL3(dual_weight(mu))
            # sanity: dual character = substitute x -> 1/x
            assert chwd == {(-e[0], -e[1], -e[2]): c for e, c in chw.items()}
            n1 = sum(decompose(pmul(chl, chw)).values())
            n2 = sum(decompose(pmul(chl, chwd)).values())
            # Racah-Speiser count
            a = (lam[0]+2, lam[1]+1, lam[2])
            r1 = sum(c * sgn3((a[0]+e[0], a[1]+e[1], a[2]+e[2])) for e, c in chw.items())
            r2 = sum(c * sgn3((a[0]-e[0], a[1]-e[1], a[2]-e[2])) for e, c in chw.items())
            cnt += 1
            if not (n1 == r1 and n2 == r2):
                bad += 1; log("  RS-count mismatch", lam, mu, n1, r1, n2, r2)
            if n1 != n2:
                bad += 1; log("  C1 FAILS", lam, mu, n1, n2)
            # all multiplicities nonnegative (sanity of decomposition)
            assert all(v > 0 for v in decompose(pmul(chl, chw)).values())
    log("A: C1 direct (lam in [-3,5]^3 dom, W=V_mu, mu in [-2,4]^3 dom): %d cases, %d failures" % (cnt, bad))
    # also W = tensor powers of V and V^* mixed, and W reducible
    bad = 0; cnt = 0
    V = char_GL3((1, 0, 0)); Vd = char_GL3((0, 0, -1))
    for lam in [(0,0,0),(1,0,0),(2,1,0),(3,0,0),(3,3,0),(4,1,-2),(5,2,2),(6,3,0),(2,2,0)]:
        chl = char_GL3(lam); W = {(0,0,0): 1}
        for step in range(7):
            W = pmul(W, V if step % 3 else pmul(V, Vd))
            Wd = {(-e[0], -e[1], -e[2]): c for e, c in W.items()}
            n1 = sum(decompose(pmul(chl, W)).values()); n2 = sum(decompose(pmul(chl, Wd)).values())
            cnt += 1
            if n1 != n2: bad += 1; log("  C1 FAILS (mixed tensor)", lam, step, n1, n2)
    log("A: C1 with W = mixed tensor products of V, V^*: %d cases, %d failures" % (cnt, bad))

# ---------- B: chain symmetry ----------
@lru_cache(maxsize=None)
def Nch(k, t):
    if t == 0: return 1
    tot = 0
    for i in range(3):
        if i == 0 or k[i-1] > k[i]:
            tot += Nch(k[:i] + (k[i]+1,) + k[i+1:], t-1)
    return tot

def partB():
    bad = cnt = 0
    for n1 in range(0, 15):
        for n2 in range(0, n1+1):
            for n3 in range(0, n2+1):
                nu = (n1, n2, n3); dual = (n1-n3, n1-n2, 0)
                for j in range(0, 31):
                    cnt += 1
                    if Nch(nu, j) != Nch(dual, j):
                        bad += 1; log("  chain symmetry FAILS", nu, j)
    log("B: chain symmetry nu_1 <= 14, j <= 30: %d cases, %d failures" % (cnt, bad))
    # N_j(nu) == #const(V_nu (x) V^{(x) j}) for small cases
    V = char_GL3((1, 0, 0)); bad = cnt = 0
    for nu in [(0,0,0),(1,0,0),(3,0,0),(3,3,0),(2,1,0),(4,2,1),(5,1,1),(3,2,2)]:
        chl = char_GL3(nu); W = {(0,0,0): 1}
        for j in range(0, 8):
            cnt += 1
            if sum(decompose(pmul(chl, W)).values()) != Nch(nu, j):
                bad += 1; log("  N_j != #const", nu, j)
            W = pmul(W, V)
    log("B: N_j(nu) == #const(V_nu (x) V^j): %d cases, %d failures" % (cnt, bad))

# ---------- own d-vector code ----------
def corners(lam):
    return [i for i in range(len(lam)) if lam[i] > 0 and (i == len(lam)-1 or lam[i+1] < lam[i])]

@lru_cache(maxsize=None)
def dvec(lam):
    """(d_0,...,d_n) by recursion d_j(lam) = sum_c d_{j-1}(lam - c); lam a tuple without trailing zeros."""
    n = sum(lam)
    if n == 0: return (1,)
    d = [0]*(n+1); d[0] = 1
    for i in corners(lam):
        sub = list(lam); sub[i] -= 1
        sub = tuple(x for x in sub if x > 0)
        ds = dvec(sub)
        for j in range(n):
            d[j+1] += ds[j]
    return tuple(d)

def conj(lam):
    return tuple(sum(1 for x in lam if x > c) for c in range(lam[0])) if lam else ()

def hooks(lam):
    lt = conj(lam)
    return [lam[i] - j + lt[j] - i - 1 for i in range(len(lam)) for j in range(lam[i])]

def f_hook(lam):
    n = sum(lam); p = 1
    for h in hooks(lam): p *= h
    assert factorial(n) % p == 0
    return factorial(n) // p

def contents(lam):
    return [j - i for i in range(len(lam)) for j in range(lam[i])]

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k):
            yield (k,) + rest

def partC():
    for K in range(1, 31):
        lam = (3,)*K + (2, 2, 2); mu = (3,)*K + (1, 1, 1)
        assert conj(lam) == (K+3, K+3, K) and conj(mu) == (K+3, K, K)
        u = dvec(lam); v = dvec(mu)
        ok = all(u[j] == v[j] for j in range(K+1))
        red = all(u[j] == Nch((3, 0, 0), j) and v[j] == Nch((3, 3, 0), j) for j in range(K+1))
        first = next(j for j in range(min(len(u), len(v))) if u[j] != v[j])
        log("C: R_%d holds=%s reduction=%s first_diff=%d (K+2=%d)" % (K, ok, red, first, K+2))

def partD():
    bad = cnt = 0
    for n in range(4, 17):
        for lam in partitions(n):
            d = dvec(lam); f = f_hook(lam); cs = contents(lam)
            C1 = sum(cs); C2 = sum(c*c for c in cs)
            n3 = n*(n-1)*(n-2); n4 = n3*(n-3)
            pred3 = f*(Fraction(2, 3) + Fraction(C2 - n*(n-1)//2, n3))
            pred4 = f*(Fraction(5, 12) + Fraction(C1*C1 - 3*C2 + n*(n-1), n4) + Fraction(C2 - n*(n-1)//2, n3))
            cnt += 1
            if pred3 != d[n-3] or pred4 != d[n-4] or d[n] != f or d[n-1] != f or d[n-2] != f:
                bad += 1
                if bad < 5: log("  top formula FAILS", lam, d[n-3], pred3, d[n-4], pred4)
    log("D: top-end formulas d_{n-3}, d_{n-4} on all partitions 4 <= n <= 16: %d shapes, %d failures" % (cnt, bad))

def partE():
    bad = 0
    for k in range(1, 41):
        lam = (5,) + (3,)*k + (2, 2, 2); mu = (4, 4) + (3,)*k + (1, 1, 1)
        hl, hm = Counter(hooks(lam)), Counter(hooks(mu))
        cl, cm = contents(lam), contents(mu)
        C1l, C1m = sum(cl), sum(cm); C2l, C2m = sum(c*c for c in cl), sum(c*c for c in cm)
        ok = (hl == hm) and (C2l == C2m) and (C1l - C1m == 2) and (2*C1l == -(3*k*k + 9*k - 2)) and (C1l*C1l != C1m*C1m) and (conj(lam) != mu)
        if not ok: bad += 1; log("  Thm H(i)/(ii) FAILS at k =", k, hl == hm, C2l == C2m, C1l, C1m)
    log("E: Theorem H(i),(ii) (hook multisets, C_2, C_1 = -(3k^2+9k-2)/2, C_1 diff 2, squares differ) for k <= 40: %d failures" % bad)

def partF():
    for k in range(1, 15):
        lam = (5,) + (3,)*k + (2, 2, 2); mu = (4, 4) + (3,)*k + (1, 1, 1); n = 3*k + 11
        u, v = dvec(lam), dvec(mu)
        D = [j for j in range(n+1) if u[j] != v[j]]
        x, y = dvec((3,)*k + (2, 2, 2)), dvec((3,)*k + (1, 1, 1))
        lb = all(u[j] == sum(comb(j, i) * x[j-i] for i in range(0, min(j, 2)+1)) and
                 v[j] == sum(comb(j, i) * y[j-i] for i in range(0, min(j, 2)+1)) for j in range(k+1))
        log("F: H_%d n=%d diff idx == {k+4..n-4}: %s (first %d, last %d) LemmaB-application=%s d_{n-4} differs=%s top4 equal=%s" % (
            k, n, D == list(range(k+4, n-3)), D[0], D[-1], lb, u[n-4] != v[n-4], u[n-3:] == v[n-3:]))

# ---------- G: arrangement cells ----------
def forms_H():
    avec = ((1, 1, 0, 0), (0, 1, 0, 0), (0, 0, 0, 0)); bvec = ((0, 0, 1, 0), (0, 0, 0, 1), (0, 0, 0, 0))
    H = set()
    for s in PERMS:
        for sg in (1, -1):
            v = [tuple(avec[i][t] + sg*bvec[s[i]][t] for t in range(4)) for i in range(3)]
            for i in range(3):
                for j in range(i+1, 3):
                    f = tuple(v[i][t] - v[j][t] for t in range(4))
                    assert all(abs(c) <= 1 for c in f)
                    # normalise sign
                    for c in f:
                        if c != 0:
                            if c < 0: f = tuple(-x for x in f)
                            break
                    H.add(f)
    for f in ((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(0,0,1,-1),(1,1,0,0)): H.add(f)
    return sorted(H)

def partG():
    H = forms_H()
    def sv(x): return tuple((sum(f[i]*x[i] for i in range(4)) > 0) - (sum(f[i]*x[i] for i in range(4)) < 0) for f in H)
    def S3(x):
        Aa, Bb, p, q = x; a = (Aa+Bb, Bb, 0); b = (p, q, 0)
        L = sum(sgn3((a[0]+b[s[0]], a[1]+b[s[1]], a[2]+b[s[2]])) for s in PERMS)
        R = sum(sgn3((a[0]-b[s[0]], a[1]-b[s[1]], a[2]-b[s[2]])) for s in PERMS)
        return L, R
    log("G: |H| =", len(H))
    seen = {}
    for BOX in (8, 16, 24):
        new = 0; bad = 0
        for Aa in range(1, BOX+1):
            for Bb in range(1, BOX+1):
                for p in range(0, BOX+1):
                    for q in range(0, p+1):
                        s = sv((Aa, Bb, p, q)); L, R = S3((Aa, Bb, p, q))
                        if L != R: bad += 1
                        if s not in seen: seen[s] = (L, R); new += 1
                        elif seen[s] != (L, R): log("  NOT CONSTANT on cell", (Aa, Bb, p, q))
        log("G: box %d: cells seen so far %d (new %d), S3 failures %d" % (BOX, len(seen), new, bad))
    random.seed(3); new = 0; bad = 0
    for _ in range(200000):
        x = (random.randint(1, 10**4), random.randint(1, 10**4), 0, 0)
        p = random.randint(0, 10**4); q = random.randint(0, p); x = (x[0], x[1], p, q)
        s = sv(x); L, R = S3(x)
        if L != R: bad += 1
        if s not in seen: new += 1
    log("G: 200000 random points up to 10^4: new cells %d, S3 failures %d" % (new, bad))

if __name__ == "__main__":
    partA(); partB(); partC(); partD(); partE(); partF(); partG()
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ref1_skeptic_c123.log'), 'w') as fh:
        fh.write("\n".join(OUT) + "\n")
