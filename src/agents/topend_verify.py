"""topend: foundational checks.
 (A) MN characters (topend_lib.chi, rim-hook surgery) vs the beta-number MN in src/check_characters.py, all lam |- n <= 10, all rho;
     column orthogonality n <= 8; branching chi^lam(rho u 1^{n-i}) = sum_nu chi^nu(rho) f^{lam/nu}, n <= 12.
 (B) sigma(rho) formula vs brute force count of square roots in S_i, i <= 8.
 (C) u_i(lam) = d_{n-i} = sum_{rho |- i} sigma(rho)/z_rho chi^lam(rho u 1^{n-i}), all lam |- n <= 14, all 0<=i<=n.
 (D) group algebra Z[S_n], n <= 7: p_m(J) := sum_k J_k^m is central for m <= 6; class expansions
     p_1 = K2, p_2 = C(n,2) + K3, p_3 = K4 + (2n-3) K2, e_2(J) = K22 + K3;
     and the JM eigenvalue theorem check: sum_rho a_rho |K_rho| chi^lam(rho) = f^lam * p_m(contents) for m<=6."""
import os, sys, itertools
from fractions import Fraction
from math import comb, factorial
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from topend_lib import *
from young import partitions, f_hook, f_skew
from census import d_vector
import check_characters as cc

def check_A(nmax_cross=10, nmax_orth=8, nmax_branch=12):
    bad = 0
    for n in range(nmax_cross + 1):
        for lam in partitions(n):
            for rho in partitions(n):
                if chi(lam, rho) != cc.chi(lam, rho): bad += 1; print("MN mismatch", lam, rho)
    print(f"(A1) MN rim-hook vs beta-number MN agree for all lam, rho |- n <= {nmax_cross}: {'OK' if bad==0 else bad}")
    bad = 0
    for n in range(nmax_orth + 1):
        P = list(partitions(n))
        for rho in P:
            for tau in P:
                s = sum(chi(l, rho) * chi(l, tau) for l in P)
                if s != (z(rho) if rho == tau else 0): bad += 1
    print(f"(A2) column orthogonality n <= {nmax_orth}: {'OK' if bad==0 else bad}")
    bad = 0
    for n in range(nmax_branch + 1):
        for lam in partitions(n):
            for i in range(n + 1):
                for rho in partitions(i):
                    lhs = chi_padded(lam, rho)
                    rhs = sum(chi(nu, rho) * f_skew(lam, nu) for nu in partitions(i))
                    if lhs != rhs: bad += 1
    print(f"(A3) branching chi^lam(rho u 1^(n-i)) = sum_nu chi^nu(rho) f^(lam/nu), n <= {nmax_branch}: {'OK' if bad==0 else bad}")

def check_B(imax=8):
    bad = 0
    for i in range(imax + 1):
        perms = list(itertools.permutations(range(i)))
        from collections import Counter
        sq = Counter(tuple(p[p[t]] for t in range(i)) for p in perms)
        for rho in partitions(i):
            g = None
            # a representative of type rho
            g = []; start = 0
            for k in rho:
                g += [start + (t + 1) % k for t in range(k)]; start += k
            if sq[tuple(g)] != sigma(rho): bad += 1; print("sigma mismatch", rho, sq[tuple(g)], sigma(rho))
    print(f"(B) sigma(rho) formula vs brute force, i <= {imax}: {'OK' if bad==0 else bad}")

def check_C(nmax=14):
    bad = 0
    for n in range(nmax + 1):
        for lam in partitions(n):
            dv = d_vector(lam)
            for i in range(n + 1):
                rhs = sum(Fraction(sigma(rho), z(rho)) * chi_padded(lam, rho) for rho in partitions(i))
                if rhs != dv[n - i]: bad += 1; print("u_i mismatch", lam, i)
    print(f"(C) u_i = sum_rho sigma/z chi, all lam |- n <= {nmax}: {'OK' if bad==0 else bad}")

# ---- group algebra
def compose(p, q):  # (p q)(x) = p(q(x))
    return tuple(p[x] for x in q)
def transp(i, k, n):
    p = list(range(n)); p[i], p[k] = p[k], p[i]; return tuple(p)
def mult(A, B):
    R = {}
    for p, a in A.items():
        for q, b in B.items():
            r = compose(p, q); R[r] = R.get(r, 0) + a * b
    return {k: v for k, v in R.items() if v}
def add(A, B, c=1):
    R = dict(A)
    for q, b in B.items(): R[q] = R.get(q, 0) + c * b
    return {k: v for k, v in R.items() if v}
def JM(k, n):  # J_k, k 1-indexed: sum_{i<k} (i k), points 0..n-1
    return {transp(i, k - 1, n): 1 for i in range(k - 1)}
def power(A, m, n):
    R = {tuple(range(n)): 1}
    for _ in range(m): R = mult(R, A)
    return R
def class_expand(A, n):
    """Return dict cycletype(without ones) -> coefficient, asserting A is a class function."""
    res = {}
    for p, a in A.items():
        ct = strip_ones(cycle_type(p)); 
        if ct in res and res[ct] != a: raise ValueError("not central")
        res[ct] = a
    # zero coefficients for classes not present are fine; but check every element of a present class is present
    for p in itertools.permutations(range(n)):
        ct = strip_ones(cycle_type(p))
        if ct in res and A.get(p, 0) != res[ct]: raise ValueError("not central (missing)")
    return res

def check_D(nmax=7, mmax=6):
    for n in range(2, nmax + 1):
        e = tuple(range(n))
        Js = [JM(k, n) for k in range(1, n + 1)]
        pm = {}
        for m in range(1, mmax + 1):
            P = {}
            for J in Js: P = add(P, power(J, m, n))
            pm[m] = class_expand(P, n)  # raises if not central
        # identities
        K = lambda rho: {strip_ones(rho): 1}
        assert pm[1] == {(2,): 1}
        assert pm[2] == {(): comb(n, 2), (3,): 1} if n >= 3 else pm[2] == {(): comb(n, 2)}
        exp3 = {(4,): 1, (2,): 2 * n - 3}
        assert pm[3] == {k: v for k, v in exp3.items() if sum(k) <= n and v}, (n, pm[3])
        # e_2(J) = sum_{j<k} J_j J_k
        E2 = {}
        for j in range(n):
            for k in range(j + 1, n):
                E2 = add(E2, mult(Js[j], Js[k]))
        e2 = class_expand(E2, n)
        assert e2 == {k: v for k, v in {(2, 2): 1, (3,): 1}.items() if sum(k) <= n}, (n, e2)
        # JM eigenvalue check via characters: tr(p_m(J) | V^lam) = f^lam C_m(lam)
        for lam in partitions(n):
            for m in range(1, mmax + 1):
                tr = sum(a * class_size(rho, n) * chi_padded(lam, rho) for rho, a in pm[m].items())
                assert tr == f_hook(lam) * C(lam, m), (lam, m)
        print(f"  n={n}: p_m(J) central for m<={mmax}; expansions of p_1,p_2,p_3,e_2 as claimed; JM trace identity OK for all lam")
        print("     p_4 =", pm[4], " p_5 =", pm[5], " p_6 =", pm[6])
    print(f"(D) group-algebra checks n <= {nmax} OK")

if __name__ == "__main__":
    check_A(); check_B(); check_C(); check_D()
