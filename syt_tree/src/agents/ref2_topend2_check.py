"""Referee #2 independent checks of topend.md T3.1-T3.3 (Lemma H, Lemma T, Theorem F') and T4.1 (Lemma K).
All arithmetic exact (Python int / Fraction). Own implementations of d-vector, f (hook), f_skew (chain count)."""
import sys, itertools
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src')

# ---------- own primitives ----------
def parts(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for rest in parts(n-k, k): yield (k,)+rest

def conj(l):
    return tuple(sum(1 for x in l if x > j) for j in range(l[0])) if l else ()

def rem_corners(l):
    return [i for i in range(len(l)) if i == len(l)-1 or l[i] > l[i+1]]
def add_corners(l):
    return [i for i in range(len(l)+1) if i == 0 or (i < len(l) and l[i-1] > l[i]) or (i == len(l))]
def rm(l, i):
    m = list(l); m[i] -= 1
    if m[i] == 0: m.pop()
    return tuple(m)
def ad(l, i):
    m = list(l)
    if i == len(m): m.append(1)
    else: m[i] += 1
    return tuple(m)

def hooks(l):
    lc = conj(l)
    return {(i, j): (l[i]-j-1) + (lc[j]-i-1) + 1 for i in range(len(l)) for j in range(l[i])}
def f_hookformula(l):
    n = sum(l); H = 1
    for h in hooks(l).values(): H *= h
    assert factorial(n) % H == 0
    return factorial(n)//H

@lru_cache(None)
def dvec(l):
    n = sum(l); d = [0]*(n+1); d[0] = 1
    if n == 0: return (1,)
    subs = [dvec(rm(l, i)) for i in rem_corners(l)]
    for j in range(1, n+1):
        d[j] = sum(s[j-1] for s in subs)
    return tuple(d)

@lru_cache(None)
def fskew(l, nu):
    # number of saturated chains nu -> l
    if l == nu: return 1
    tot = 0
    for i in rem_corners(l):
        m = rm(l, i)
        if len(m) >= len(nu) and all(m[k] >= nu[k] for k in range(len(nu))):
            tot += fskew(m, nu)
    return tot

# sanity vs library
from young import f_hook, f_skew
from census import d_vector
for n in range(0, 11):
    for l in parts(n):
        assert dvec(l) == tuple(d_vector(l)), l
        assert f_hookformula(l) == f_hook(l), l
        assert dvec(l)[n] == f_hookformula(l)
print("own primitives agree with library, n<=10")

def C(j, k):
    return comb(j, k) if 0 <= k <= j else 0

# ---------- Lemma H ----------
def dH(a, b, j):
    n = a+b
    if j == n: return C(n-1, a-1)
    s1 = sum(C(j, a-c) for c in range(1, a+1) if 0 <= n-j-c <= b)
    lo, hi = max(0, a+j-n), min(a-1, j)
    s2 = sum(C(j, e) for e in range(lo, hi+1))
    assert s1 == s2, (a, b, j)
    return s1

bad = 0
for n in range(1, 41):
    for a in range(1, n+1):
        b = n-a
        lam = (a,)+(1,)*b
        d = dvec(lam)
        for j in range(n+1):
            if dH(a, b, j) != d[j]:
                bad += 1; print("Lemma H FAIL", lam, j, dH(a,b,j), d[j])
print("Lemma H checked all hooks n<=40, failures:", bad)
# also check that formula with j=n (without special case) would be wrong -> confirms need for j<n clause
print("Lemma H j=n generic-sum value for (3,1,1):", sum(C(5,e) for e in range(max(0,3+5-5), min(2,5)+1)), "vs f=", f_hookformula((3,1,1)))

# ---------- Aitken (T3.0) own check ----------
def aitken(l, nu):
    ell = len(l); nu = tuple(nu)+(0,)*(ell-len(nu)); m = sum(l)-sum(nu)
    M = [[Fraction(1, factorial(l[i]-nu[j]-i+j)) if l[i]-nu[j]-i+j >= 0 else Fraction(0) for j in range(ell)] for i in range(ell)]
    # determinant by Fraction gaussian elimination
    M = [r[:] for r in M]; det = Fraction(1)
    for c in range(ell):
        p = next((r for r in range(c, ell) if M[r][c] != 0), None)
        if p is None: return 0
        if p != c: M[c], M[p] = M[p], M[c]; det = -det
        det *= M[c][c]
        for r in range(c+1, ell):
            fct = M[r][c]/M[c][c]
            for k in range(c, ell): M[r][k] -= fct*M[c][k]
    v = det*factorial(m); assert v.denominator == 1
    return int(v)
bad = 0; cnt = 0
for n in range(0, 10):
    for l in parts(n):
        for k in range(n+1):
            for nu in parts(k):
                if len(nu) <= len(l) and all(nu[i] <= l[i] for i in range(len(nu))):
                    cnt += 1
                    if aitken(l, nu) != fskew(l, nu): bad += 1; print("Aitken FAIL", l, nu)
print("Aitken checked", cnt, "pairs n<=9, failures:", bad)

# ---------- Lemma T ----------
def dT(p, q, j):
    n = p+q; s = 0
    for qp in range(0, q+1):
        pp = n-j-qp
        if qp <= pp <= p:
            s += C(j, p-pp) - C(j, p-qp+1)
    return s
bad = 0
for n in range(0, 41):
    for q in range(0, n//2+1):
        p = n-q; lam = (p, q) if q > 0 else ((p,) if p > 0 else ())
        d = dvec(lam)
        for j in range(n+1):
            if dT(p, q, j) != d[j]: bad += 1; print("Lemma T FAIL", lam, j)
        # also check the per-term formula f^{lam/nu} = C(j,p-p') - C(j,p-q'+1)
        for qp in range(0, q+1):
            for pp in range(qp, p+1):
                j = n-pp-qp; nu = (pp, qp) if qp > 0 else ((pp,) if pp > 0 else ())
                if fskew(lam, nu) != C(j, p-pp) - C(j, p-qp+1): bad += 1; print("Lemma T term FAIL", lam, nu)
print("Lemma T checked all two-row shapes n<=40 (+ per-term), failures:", bad)

# ---------- Theorem F' ----------
def check_Fprime(a, n, use_closed):
    lam = (a,)+(1,)*(n-a); mu = (2,)*(a-1)+(1,)*(n-2*a+2)
    if use_closed:
        dl = [dH(a, n-a, j) for j in range(n+1)]
        p, q = n-a+1, a-1
        dm = [dT(p, q, j) for j in range(n+1)]
    else:
        dl = dvec(lam); dm = dvec(mu)
        assert dvec(mu) == dvec(conj(mu))
    J = n-2*a+3
    ok_agree = all(dl[j] == dm[j] for j in range(0, J))
    ok_diff = (dl[J]-dm[J] == 1)
    return ok_agree, ok_diff, [dl[j]-dm[j] for j in range(n+1)]

bad = 0
for a in range(3, 9):
    for n in range(2*a-1, 41):
        ag, df, diffs = check_Fprime(a, n, False)
        if not (ag and df): bad += 1; print("F' FAIL direct", a, n, diffs)
print("Theorem F' checked directly (own d-vectors) a<=8, n<=40, failures:", bad)
bad = 0
for a in range(3, 16):
    for n in range(2*a-1, 121):
        ag, df, diffs = check_Fprime(a, n, True)
        if not (ag and df): bad += 1; print("F' FAIL closed", a, n)
print("Theorem F' checked via closed forms a<=15, n<=120, failures:", bad)
# boundary: n = 2a-2 (mu is a rectangle) -- outside the theorem's stated range
for a in range(3, 9):
    n = 2*a-2
    ag, df, diffs = check_Fprime(a, n, False)
    print("boundary n=2a-2, a=%d: agree=%s diff1=%s diffs=%s" % (a, ag, df, diffs))
# a = 2 boundary (outside range a>=3): lam=(2,1^{n-2}), mu=(2,1^{n-2}) identical
# remark for a=3: differences at d_{n-3..n}
for n in range(5, 15):
    ag, df, diffs = check_Fprime(3, n, False)
    print("a=3 n=%d diffs[n-4..n] =" % n, diffs[n-4:])

# ---------- Lemma K ----------
def contents_XY(l):
    Y = [l[i]-1-i for i in rem_corners(l)]
    X = [(l[i] if i < len(l) else 0) - i for i in add_corners(l)]
    return X, Y
bad = 0; cnt = 0
for m in range(0, 19):
    for nu in parts(m):
        X, Y = contents_XY(nu)
        assert len(X) == len(Y)+1
        fn = f_hookformula(nu); Hn = 1
        for h in hooks(nu).values(): Hn *= h
        for i in add_corners(nu):
            x = (nu[i] if i < len(nu) else 0) - i
            nu2 = ad(nu, i); cnt += 1
            lhs = Fraction(f_hookformula(nu2), fn)
            num = 1
            for y in Y: num *= abs(x-y)
            den = 1
            for xp in X:
                if xp != x: den *= abs(x-xp)
            rhs = Fraction((m+1)*num, den)
            if lhs != rhs: bad += 1; print("Lemma K FAIL", nu, x, lhs, rhs)
            # intermediate: row product R and column product Cc
            h1 = hooks(nu); h2 = hooks(nu2)
            R = Fraction(1); Cc = Fraction(1)
            for (r, c), h in h1.items():
                if r == i: R *= Fraction(h2[(r, c)], h)
                elif c == (nu[i] if i < len(nu) else 0): Cc *= Fraction(h2[(r, c)], h)
                else: assert h2[(r, c)] == h
            assert h2[(i, nu[i] if i < len(nu) else 0)] == 1
            Rp = Fraction(1)
            for xp in X:
                if xp < x: Rp *= (x-xp)
            for y in Y:
                if y < x: Rp /= (x-y)
            Cp = Fraction(1)
            for xp in X:
                if xp > x: Cp *= (xp-x)
            for y in Y:
                if y > x: Cp /= (y-x)
            if R != Rp or Cc != Cp: bad += 1; print("Lemma K intermediate FAIL", nu, x, R, Rp, Cc, Cp)
print("Lemma K checked (formula + row/column products) for all nu |- m<=18,", cnt, "addable corners, failures:", bad)
# expected count for m<=26
def pcount(n):
    return sum(1 for _ in parts(n))
print("expected #checks for m<=26:", sum(pcount(k) for k in range(1, 28)))
