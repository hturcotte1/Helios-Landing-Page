"""topend: rigorous determination of the content polynomials Omega_rho with
      omega_lam(K_rho) = |K_rho| chi^lam(rho u 1^{n-|rho|}) / f^lam = Omega_rho(n; C_1(lam), ..., C_m(lam)),
for every even-type rho without 1's, |rho| <= 10, where m = |rho| - l(rho).
By Theorem P (report, proved), Omega_rho lies in V_m := span{ n^j C_mu : j + |mu| <= m }  (C_mu = prod C_{mu_i}).
Method: (1) evaluate omega_lam(K_rho) by Murnaghan-Nakayama on the test set T = {lam |- n : n <= NS};
        (2) check that the evaluation map V_m -> Q^T is injective (full column rank; rank computed modulo a
            large prime, which lower-bounds the rank over Q); then the unique element of V_m agreeing with omega on T
            is Omega_rho, and the identity is PROVED for all n;
        (3) solve exactly over Q (Fractions) on a pivot subset of T and re-verify on all lam |- n <= NV (>= NS).
Output: pickle {rho: [(coeff Fraction, (j, mu)), ...]} and a log."""
import os, sys, pickle, time
from fractions import Fraction
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from topend_mn import *
from young import partitions

HERE = os.path.dirname(os.path.abspath(__file__))
PRIME = (1 << 61) - 1

def monomials(m):
    """all (j, mu) with j + |mu| <= m, mu a partition (possibly empty)."""
    res = []
    for s in range(m + 1):
        for mu in partitions(s):
            for j in range(m - s + 1):
                res.append((j, mu))
    return res

def evaluate(lam, mono, Ccache):
    j, mu = mono
    v = sum(lam) ** j
    for p in mu: v *= Ccache[p]
    return v

def rank_and_pivot_rows_modp(rows):
    """Gaussian elimination mod PRIME on integer rows; returns (rank, list of row indices forming a basis of the row space)."""
    A = [[x % PRIME for x in r] for r in rows]
    m = len(A[0]); piv_rows = []; used = [False] * len(A)
    basis = []  # reduced rows
    pivcols = []
    for c in range(m):
        p = None
        for i in range(len(A)):
            if not used[i] and A[i][c] % PRIME != 0:
                p = i; break
        if p is None: continue
        used[p] = True; piv_rows.append(p)
        inv = pow(A[p][c], PRIME - 2, PRIME)
        A[p] = [(x * inv) % PRIME for x in A[p]]
        for i in range(len(A)):
            if i != p and A[i][c] % PRIME:
                f = A[i][c]; A[i] = [(x - f * y) % PRIME for x, y in zip(A[i], A[p])]
        pivcols.append(c)
    return len(piv_rows), piv_rows

def solve_exact(rows, rhs):
    """rows: square-ish integer matrix with full column rank; exact solve by Fraction elimination."""
    m = len(rows[0])
    A = [[Fraction(x) for x in r] + [Fraction(b)] for r, b in zip(rows, rhs)]
    r = 0; piv = []
    for c in range(m):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: raise ValueError("rank deficient")
        A[r], A[p] = A[p], A[r]
        inv = 1 / A[r][c]; A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    for i in range(r, len(A)):
        if A[i][m] != 0: raise ValueError("inconsistent")
    sol = [Fraction(0)] * m
    for i, c in enumerate(piv): sol[c] = A[i][m]
    return sol

def main(NS=13, NV=16, max_size=10):
    t0 = time.time()
    lams = {n: list(partitions(n)) for n in range(NV + 1)}
    Cc = {}
    for n in range(NV + 1):
        for lam in lams[n]:
            Cc[lam] = {k: C(lam, k) for k in range(0, 11)}
    T = [lam for n in range(NS + 1) for lam in lams[n]]
    results = {}
    for rho in even_types_no_ones(max_size):
        m = sum(rho) - len(rho)
        monos = monomials(m)
        rows = [[evaluate(lam, mo, Cc[lam]) for mo in monos] for lam in T]
        rank, prows = rank_and_pivot_rows_modp(rows)
        assert rank == len(monos), (rho, rank, len(monos))
        om = {lam: omega(lam, rho) for lam in T}
        rhs = [om[T[i]] for i in prows]
        sol = solve_exact([rows[i] for i in prows], rhs)
        # verify on the whole of n <= NV
        bad = 0; cnt = 0
        for n in range(NV + 1):
            for lam in lams[n]:
                v = sum(a * evaluate(lam, mo, Cc[lam]) for a, mo in zip(sol, monos) if a)
                w = om[lam] if lam in om else omega(lam, rho)
                cnt += 1
                if v != w: bad += 1
        terms = [(a, mo) for a, mo in zip(sol, monos) if a]
        results[rho] = terms
        print(f"rho={rho}: m={m}, dim V_m={len(monos)}, rank on T (n<={NS}) = {rank} (FULL), verified on {cnt} partitions n<={NV}: {'OK' if bad == 0 else str(bad) + ' FAIL'}  [{time.time()-t0:.0f}s]", flush=True)
        print("   Omega =", " + ".join(f"({a})*n^{j}*C{list(mu)}" for a, (j, mu) in terms), flush=True)
    with open(os.path.join(HERE, 'topend_omega_results.pkl'), 'wb') as fh:
        pickle.dump(results, fh)

if __name__ == "__main__":
    NS = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    NV = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    main(NS, NV)
