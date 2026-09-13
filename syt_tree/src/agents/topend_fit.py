"""topend: for every even-type class rho (no parts 1, |rho| <= 10) find the polynomial
   Omega_rho(n; C_1, C_2, ...)  with  omega_lam(K_rho) := |K_rho| chi^lam(rho u 1^{n-|rho|}) / f^lam = Omega_rho(n; C_k(lam)).
Ansatz: monomials n^j C_mu with j + |mu| + l(mu) <= |rho| and an even number of odd parts in mu (transpose symmetry).
Solve exactly (Fractions) on all lam |- n for n <= NSOLVE, then verify on all lam |- n <= NVERIFY.
Output: the polynomials, and the degree bound N_max needed by the polynomiality lemma."""
import os, sys, pickle
from fractions import Fraction
from math import comb, factorial
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from topend_lib import *
from young import partitions, f_hook

def monomials(s):
    """(j, mu) with j + |mu| + l(mu) <= s, mu with an even number of odd parts."""
    res = []
    for w in range(s + 1):            # w = |mu| + l(mu)
        for m in range((w + 1) // 2, w + 1):   # |mu| = m, l(mu) = w - m  => each part >= 1: m >= w - m
            l = w - m
            if l == 0 and m != 0: continue
            for mu in partitions(m):
                if len(mu) != l: continue
                if sum(1 for p in mu if p % 2) % 2: continue
                for j in range(s - w + 1):
                    res.append((j, mu))
    return res

def evaluate(lam, mono):
    j, mu = mono
    n = sum(lam); v = Fraction(n ** j)
    for p in mu: v *= C(lam, p)
    return v

def solve(rows, rhs):
    """Exact Gaussian elimination; returns (solution with free vars = 0, rank, consistent)."""
    m = len(rows[0]); A = [list(r) + [b] for r, b in zip(rows, rhs)]
    piv = []; r = 0
    for c in range(m):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]; A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == len(A): break
    consistent = all(all(x == 0 for x in row[:m]) is False or row[m] == 0 for row in A)
    sol = [Fraction(0)] * m
    for i, c in enumerate(piv): sol[c] = A[i][m]
    return sol, len(piv), consistent

def fit_all(NSOLVE=13, NVERIFY=16, max_size=10, verbose=True):
    classes = even_type_no_ones(max_size)
    lams = {n: list(partitions(n)) for n in range(NVERIFY + 1)}
    results = {}
    for rho in classes:
        s = sum(rho); monos = monomials(s)
        rows, rhs = [], []
        for n in range(0, NSOLVE + 1):
            for lam in lams[n]:
                rows.append([evaluate(lam, mo) for mo in monos]); rhs.append(omega(lam, rho))
        sol, rank, cons = solve(rows, rhs)
        assert cons, ("inconsistent", rho)
        # verify
        bad = 0
        for n in range(0, NVERIFY + 1):
            for lam in lams[n]:
                v = sum(a * evaluate(lam, mo) for a, mo in zip(sol, monos) if a)
                if v != omega(lam, rho): bad += 1
        # degree bound of the polynomiality lemma: for each possible class tau appearing, D_tau = max_j,mu (j + floor(|mu| - t/2)), t=|tau|
        terms = [(a, mo) for a, mo in zip(sol, monos) if a]
        Mmax = max(sum(mo[1]) for a, mo in terms)
        Nneed = 0
        for t in range(0, 2 * Mmax + 1):
            D = max(mo[0] + (sum(mo[1]) * 2 - t) // 2 for a, mo in terms)   # floor(|mu| - t/2)
            Nneed = max(Nneed, t + D)
        results[rho] = dict(terms=terms, rank=rank, nmonos=len(monos), bad=bad, Nneed=Nneed)
        if verbose:
            print(f"rho={rho}: {len(monos)} monomials, rank {rank}, verify n<={NVERIFY}: {'OK' if bad==0 else str(bad)+' FAIL'}, lemma needs n<= {Nneed}")
            print("   omega =", " + ".join(f"({a})*n^{j}*C{list(mu)}" for a, (j, mu) in terms))
    return results

if __name__ == "__main__":
    NS = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    NV = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    res = fit_all(NS, NV)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'topend_fit_results.pkl'), 'wb') as fh:
        pickle.dump({k: v['terms'] for k, v in res.items()}, fh)
