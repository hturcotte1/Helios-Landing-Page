"""(1) Verify the 16 polynomials Omega_rho (as printed in agent_notes/topend.md, parsed into
    topend3_omega_from_report.json) against omega_lam(K_rho) computed by my own MN implementation,
    for ALL lam |- n <= NMAX.
(2) Rank certificate: the evaluation map V_m -> Q^{partitions with n<=16} is injective for m = 2,4,6,8,
    where V_m = span{ n^j C_mu : j+|mu| <= m }.  Rank computed mod p = 2^61-1 (full rank mod p => full rank over Q).
"""
import sys, json, time
from fractions import Fraction
from math import comb
import sympy as sp
from topend3_lib import *

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 16
t0 = time.time()
d = json.load(open('topend3_omega_from_report.json'))
syms = {s: sp.Symbol(s) for s in ['n', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C8']}
polys = {}
for k, v in d.items():
    rho = tuple(int(x) for x in k.strip('()').split(',') if x.strip())
    e = sp.sympify(v, locals=syms)
    polys[rho] = sp.Poly(e, *syms.values())

def evalpoly(P, lam):
    n = sum(lam)
    vals = {'n': n}
    for k in (1, 2, 3, 4, 5, 6, 8):
        vals['C%d' % k] = C(lam, k)
    tot = Fraction(0)
    for mon, coef in P.terms():
        term = Fraction(coef.p, coef.q)
        for s, e in zip(syms, mon):
            if e:
                term *= Fraction(vals[s]) ** e
        tot += term
    return tot

# (1)
cnt = 0
for n in range(1, NMAX + 1):
    for lam in partitions(n):
        for rho, P in polys.items():
            lhs = omega(lam, rho)
            rhs = evalpoly(P, lam)
            assert lhs == rhs, (lam, rho, lhs, rhs)
            cnt += 1
    print(f"  n={n} ok [{time.time()-t0:.0f}s]", flush=True)
print(f"(1) all 16 Omega_rho verified by MN for all lam |- n <= {NMAX}: {cnt} checks OK")

# (2) rank certificate
p = (1 << 61) - 1
def monomials(m):
    mons = []
    for j in range(m + 1):
        for w in range(m + 1 - j):
            for mu in partitions(w):
                mons.append((j, mu))
    return mons
def rank_mod_p(rows):
    rows = [r[:] for r in rows]; rank = 0; ncol = len(rows[0])
    for c in range(ncol):
        piv = next((i for i in range(rank, len(rows)) if rows[i][c] % p), None)
        if piv is None: continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = pow(rows[rank][c], p - 2, p)
        rows[rank] = [(x * inv) % p for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][c] % p:
                fac = rows[i][c]
                rows[i] = [(a - fac * b) % p for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return rank
for m in (2, 4, 6, 8):
    mons = monomials(m)
    for N in (12, 14, 15, 16):
        rows = []
        for n in range(0, N + 1):
            for lam in partitions(n):
                Cv = {k: C(lam, k) for k in range(1, m + 1)}
                row = []
                for (j, mu) in mons:
                    v = n ** j
                    for part in mu: v *= Cv[part]
                    row.append(v % p)
                rows.append(row)
        rk = rank_mod_p(rows)
        print(f"(2) m={m}: dim V_m = {len(mons)}, rank of evaluation matrix on all lam |- n<={N} (mod 2^61-1) = {rk}"
              + ("  FULL RANK" if rk == len(mons) else ""), flush=True)
print(f"done [{time.time()-t0:.0f}s]")
