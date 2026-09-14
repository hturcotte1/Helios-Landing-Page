"""topend: (A) cross-checks of the beta-number MN implementation; (B) the group-algebra identities used in the
proofs (Theorem A), checked in Z[S_n] for n <= NMAX by brute force multiplication of permutations;
(C) the JM-trace identity tr(F(J)|V^lam) = f^lam F(contents) for F = p_m, m <= 4, via characters;
(D) the 'fixed points cost' lemma  m >= m'' + f  for all sequences of m transpositions (m <= 5, points <= 6);
(E) Theorem A formulas for omega(3), omega(2,2), omega(5), omega(3,3) via MN for all lam |- n <= 16;
(F) end-to-end identities using ONLY the d-vector (Young-lattice DP) and contents, all lam |- n <= 20:
      (3u_3-2f) 2C(n,3)      = f (C_2 - C(n,2))
      (4u_4-4u_3+f) 3C(n,4)  = f (C_1^2/2 - 3C_2/2 + C(n,2))
      (5u_5-5u_4+f) 24C(n,5) = f (C_4 - (3n-10)C_2 - 2C_1^2 + n(n-1)(5n-19)/6)
      and the level-6 identity  (n)_6 u_6 / f = explicit polynomial in n, C_1^2, C_2, C_4."""
import os, sys, itertools
from fractions import Fraction
from math import comb, factorial
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from topend_mn import *
from young import partitions, f_hook, conjugate
from census import d_vector
import check_characters as cc
import topend_lib as tl

def binom(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def check_A(nmax=10):
    bad = 0
    for n in range(nmax + 1):
        P = list(partitions(n))
        for lam in P:
            for rho in P:
                a = chi(lam, rho)
                if a != cc.chi(lam, rho) or a != tl.chi(lam, rho): bad += 1
        for rho in P:            # column orthogonality
            for tau in P:
                s = sum(chi(l, rho) * chi(l, tau) for l in P)
                if s != (z(rho) if rho == tau else 0): bad += 1
        for lam in P:            # transpose: chi^{lam'}(rho) = sgn(rho) chi^lam(rho)
            for rho in P:
                if chi(conjugate(lam), rho) != (-1) ** (n - len(rho)) * chi(lam, rho): bad += 1
    print(f"(A) MN beta-numbers vs two other MN codes, orthogonality, transpose rule, n <= {nmax}: {'OK' if bad == 0 else bad}")

# ---- group algebra Z[S_n], permutations as tuples, (p*q)(x) = p(q(x)) (q applied first)
def compose(p, q): return tuple(p[x] for x in q)
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
def JM(k, n): return {transp(i, k - 1, n): 1 for i in range(k - 1)}   # J_k = sum_{i<k} (i k), k = 1..n
def power(A, m, n):
    R = {tuple(range(n)): 1}
    for _ in range(m): R = mult(R, A)
    return R
def cycle_type(p):
    n = len(p); seen = [False] * n; ct = []
    for i in range(n):
        if not seen[i]:
            l = 0; j = i
            while not seen[j]: seen[j] = True; j = p[j]; l += 1
            ct.append(l)
    return tuple(sorted(ct, reverse=True))
def class_sum(rhobar, n):
    return {p: 1 for p in itertools.permutations(range(n)) if strip_ones(cycle_type(p)) == rhobar}
def as_classes(A, n):
    """class expansion of a central element (raises if not central)."""
    res = {}
    for p in itertools.permutations(range(n)):
        ct = strip_ones(cycle_type(p)); v = A.get(p, 0)
        if ct in res and res[ct] != v: raise ValueError(("not central", ct))
        res[ct] = v
    return {k: v for k, v in res.items() if v}

def check_B(nmax=8):
    for n in range(2, nmax + 1):
        Js = [JM(k, n) for k in range(1, n + 1)]
        p = {}
        for m in range(1, 5):
            P = {}
            for J in Js: P = add(P, power(J, m, n))
            p[m] = as_classes(P, n)
        K2 = class_sum((2,), n); K3 = class_sum((3,), n)
        K2sq = as_classes(mult(K2, K2), n); K3sq = as_classes(mult(K3, K3), n)
        def expect(d): return {k: v for k, v in d.items() if v and sum(k) <= n}
        assert p[1] == expect({(2,): 1}), (n, p[1])
        assert p[2] == expect({(): comb(n, 2), (3,): 1}), (n, p[2])
        assert K2sq == expect({(): comb(n, 2), (3,): 3, (2, 2): 2}), (n, K2sq)
        assert p[3] == expect({(4,): 1, (2,): 2 * n - 3}), (n, p[3])
        assert p[4] == expect({(5,): 1, (3,): 3 * n - 4, (2, 2): 4, (): n * (n - 1) * (4 * n - 5) // 6}), (n, p[4])
        assert K3sq == expect({(): 2 * comb(n, 3), (3,): 3 * n - 8, (2, 2): 8, (5,): 5, (3, 3): 2}), (n, K3sq)
        # (C) JM trace identity via characters: sum_rho a_rho |K_rho| chi^lam(rho) = f^lam C_m(lam)
        for lam in partitions(n):
            for m in range(1, 5):
                tr = sum(a * class_size(rho, n) * chi(lam, rho) for rho, a in p[m].items())
                assert tr == f_hook(lam) * C(lam, m), (lam, m)
        print(f"  n={n}: (I1)-(I6) hold in Z[S_{n}]; JM trace identity for p_1..p_4 on every V^lam OK")
    print(f"(B),(C) group-algebra identities and JM traces, 2 <= n <= {nmax}: OK")

def check_D(mmax=5, vmax=6):
    """m >= (|V| - c_V(g)) + f for every sequence of m transpositions on {0..vmax-1}, V = points used."""
    bad = 0; total = 0
    pairs = [(i, j) for i in range(vmax) for j in range(i + 1, vmax)]
    for m in range(0, mmax + 1):
        for seq in itertools.product(pairs, repeat=m):
            V = set(x for t in seq for x in t)
            g = list(range(vmax))
            for (i, j) in reversed(seq):   # product t_1 ... t_m applied right-to-left
                g = [ (j if x == i else i if x == j else x) for x in g ]  # left-multiply by (i j)
            # cycles of g on V
            seen = set(); c = 0; f = 0
            for x in V:
                if x in seen: continue
                c += 1; y = x; l = 0
                while y not in seen: seen.add(y); y = g[y]; l += 1
                if l == 1: f += 1
            total += 1
            if m < (len(V) - c) + f: bad += 1
    print(f"(D) fixed-point-cost lemma m >= m'' + f: {total} sequences (m <= {mmax}, {vmax} points), {'OK' if bad == 0 else bad}")

def om3(n, C1, C2, C4): return C2 - binom(n, 2)
def om22(n, C1, C2, C4): return Fraction(C1 * C1, 2) - Fraction(3 * C2, 2) + binom(n, 2)
def om5(n, C1, C2, C4): return C4 - (3 * n - 10) * C2 - 2 * C1 * C1 + Fraction(n * (n - 1) * (5 * n - 19), 6)
def om33(n, C1, C2, C4):
    w3, w22, w5 = om3(n, C1, C2, C4), om22(n, C1, C2, C4), om5(n, C1, C2, C4)
    return Fraction(w3 * w3 - 2 * binom(n, 3) - (3 * n - 8) * w3 - 8 * w22 - 5 * w5, 2)

def check_E(nmax=16):
    bad = 0; cnt = 0
    for n in range(nmax + 1):
        for lam in partitions(n):
            C1, C2, C4 = C(lam, 1), C(lam, 2), C(lam, 4)
            for rho, F in (((3,), om3), ((2, 2), om22), ((5,), om5), ((3, 3), om33)):
                cnt += 1
                if omega(lam, rho) != F(n, C1, C2, C4): bad += 1; print("E mismatch", lam, rho)
    print(f"(E) Theorem A: omega(3), omega(2,2), omega(5), omega(3,3) via MN for all lam |- n <= {nmax}: {cnt} checks, {'OK' if bad == 0 else bad}")

def g6_poly(n, C1, C2, C4):
    """(n)_6 u_6 / f = sum over even types rhobar of size <= 6 of sigma I_{6-s} C(n-s, 6-s) omega(rhobar)."""
    w = {(): 1, (3,): om3(n, C1, C2, C4), (2, 2): om22(n, C1, C2, C4), (5,): om5(n, C1, C2, C4), (3, 3): om33(n, C1, C2, C4)}
    return sum(sigma(rb) * involutions(6 - sum(rb)) * binom(n - sum(rb), 6 - sum(rb)) * w[rb] for rb in w)

def check_F(nmax=20):
    bad = 0; cnt = 0
    for n in range(nmax + 1):
        for lam in partitions(n):
            dv = d_vector(lam); f = dv[n]
            C1, C2, C4 = C(lam, 1), C(lam, 2), C(lam, 4)
            u = lambda i: dv[n - i] if i <= n else 0
            cnt += 4
            if (3 * u(3) - 2 * f) * 2 * binom(n, 3) != f * om3(n, C1, C2, C4): bad += 1; print("F3", lam)
            if (4 * u(4) - 4 * u(3) + f) * 3 * binom(n, 4) != f * om22(n, C1, C2, C4): bad += 1; print("F22", lam)
            if (5 * u(5) - 5 * u(4) + f) * 24 * binom(n, 5) != f * om5(n, C1, C2, C4): bad += 1; print("F5", lam)
            if n >= 6 and Fraction(factorial(n) // factorial(n - 6) * u(6), f) != g6_poly(n, C1, C2, C4): bad += 1; print("F6", lam)
    print(f"(F) d-vector-only identities for chi(3), chi(2,2), chi(5) and the level-6 identity, all lam |- n <= {nmax}: {cnt} checks, {'OK' if bad == 0 else bad}")

if __name__ == "__main__":
    check_A(); check_B(); check_D(); check_E(); check_F()
