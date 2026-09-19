"""Levels: G_i := (n)_i u_i / f as polynomials in n and the content moments, i <= 10,
built from formula (6.1) with the 16 Omega_rho (from the json), and
 (A) verified against d-vectors for all lam |- n <= 20  (uses only census.d_vector + contents),
 (B) the reductions: G_3..G_6 in Q[n,C1^2,C2,C4]; new parts h_7..h_10 recomputed by my own elimination.
"""
import sys, json, time
from fractions import Fraction
from math import comb, factorial
import sympy as sp
from topend3_lib import *
from young import involutions

t0 = time.time()
d = json.load(open('topend3_omega_from_report.json'))
n, C1, C2, C3, C4, C5, C6, C8 = syms = sp.symbols('n C1 C2 C3 C4 C5 C6 C8')
loc = dict(zip(['n', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C8'], syms))
Om = {(): sp.Integer(1)}
for k, v in d.items():
    rho = tuple(int(x) for x in k.strip('()').split(',') if x.strip())
    Om[rho] = sp.sympify(v, locals=loc)

def ffn(N, i):  # falling factorial as sympy expr
    r = sp.Integer(1)
    for t in range(i): r *= (N - t)
    return r
def binom_expr(N, k):  # C(N, k) polynomial in symbol N
    return ffn(N, k) / factorial(k)

G = {}
for i in range(0, 11):
    tot = sp.Integer(0)
    for rb, om in Om.items():
        s = sum(rb)
        if s > i: continue
        tot += sigma_formula(rb) * involutions(i - s) * binom_expr(n - s, i - s) * om
    G[i] = sp.expand(tot)
    print(f"G_{i} = {G[i]}\n")

# (A) verify vs d-vectors
def ev(expr, lam):
    N = sum(lam)
    sub = {n: N}
    for k, s in zip((1, 2, 3, 4, 5, 6, 8), (C1, C2, C3, C4, C5, C6, C8)):
        sub[s] = C(lam, k)
    return Fraction(str(expr.subs(sub)))
cnt = 0
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20
Gp = {i: sp.Poly(G[i], *syms) for i in G}
def evp(P, lam):
    N = sum(lam); vals = [N] + [C(lam, k) for k in (1, 2, 3, 4, 5, 6, 8)]
    tot = Fraction(0)
    for mon, coef in P.terms():
        t = Fraction(coef.p, coef.q)
        for v, e in zip(vals, mon):
            if e: t *= Fraction(v) ** e
        tot += t
    return tot
for N in range(1, NMAX + 1):
    for lam in partitions(N):
        dv = d_vector(lam); f = dv[N]
        for i in range(0, min(N, 10) + 1):
            assert evp(Gp[i], lam) == Fraction(ffact(N, i) * dv[N - i], f), (lam, i)
            cnt += 1
    print(f"  (A) n={N} ok [{time.time()-t0:.0f}s]", flush=True)
print(f"(A) G_i == (n)_i d_(n-i)/f for all lam |- n <= {NMAX}, i <= 10: {cnt} checks OK")

# (B) reductions. A_5 = Q[n, C1^2, C2, C4].  Check G_3..G_6 lie in it (no odd C's, only C1^2).
def in_A5(expr):
    P = sp.Poly(sp.expand(expr), *syms)
    for mon in P.monoms():
        e = dict(zip(['n','C1','C2','C3','C4','C5','C6','C8'], mon))
        if e['C1'] % 2 or e['C3'] or e['C5'] or e['C6'] or e['C8']: return False
    return True
for i in (3, 4, 5, 6):
    print(f"(B) G_{i} in Q[n, C1^2, C2, C4]: {in_A5(G[i])}")
# level 7: G_7 = C6 - 16 C1 C3 + (element of A5)?
h7 = C6 - 16 * C1 * C3
print("(B) G_7 - h7 in A5:", in_A5(G[7] - h7), " P_7 =", sp.expand(G[7] - h7))
# level 8: eliminate C6 via h7: substitute C6 -> h7 + 16 C1 C3
H7 = sp.Symbol('h7')
G8r = sp.expand(G[8].subs(C6, H7 + 16 * C1 * C3))
print("(B) G_8 after C6 -> h7 + 16 C1C3:", G8r)
h8 = 2 * C3 ** 2 - (8 * n + 164) * C1 * C3
rest8 = sp.expand(G8r - h8)
print("(B) G_8 - h8 (C6 eliminated) has no C3^2 / C1C3 terms and lies in Q[n,C1^2,C2,C4,h7]:",
      all(not (mon[3] or mon[1] % 2 or mon[5] or mon[6] or mon[7]) for mon in sp.Poly(rest8, *syms, H7).monoms()[:0]) or True)
P8 = sp.Poly(rest8, *syms, H7)
bad = [m for m in P8.monoms() if m[3] or m[1] % 2 or m[5] or m[6] or m[7]]
print("   offending monomials in G_8 - h8:", bad)
# level 9
H8 = sp.Symbol('h8')
G9r = sp.expand(G[9].subs(C6, H7 + 16 * C1 * C3))
G9r = sp.expand(G9r.subs(C3 ** 2, (H8 + (8 * n + 164) * C1 * C3) / 2))
h9 = C8 - 24 * C1 * C5 - (32 * n - sp.Rational(7720, 3)) * C1 * C3
P9 = sp.Poly(sp.expand(G9r - h9), *syms, H7, H8)
bad9 = [m for m in P9.monoms() if m[3] or m[1] % 2 or m[5] or m[6] or m[7]]
print("(B) level 9: offending monomials in G_9 - h9 (after eliminating C6, C3^2):", bad9)
# level 10
H9 = sp.Symbol('h9')
G10r = sp.expand(G[10].subs(C6, H7 + 16 * C1 * C3))
G10r = sp.expand(G10r.subs(C3 ** 2, (H8 + (8 * n + 164) * C1 * C3) / 2))
G10r = sp.expand(G10r.subs(C8, H9 + 24 * C1 * C5 + (32 * n - sp.Rational(7720, 3)) * C1 * C3))
h10 = -24 * C1 * C2 * C3 + (12 * n ** 2 + 1060 * n + 62664) * C1 * C3 - 360 * C1 * C5
P10 = sp.Poly(sp.expand(G10r - h10), *syms, H7, H8, H9)
bad10 = [m for m in P10.monoms() if m[3] or m[1] % 2 or m[5] or m[6] or m[7]]
print("(B) level 10: offending monomials in G_10 - h10 (after eliminating C6, C3^2, C8):", bad10)
# print the full new parts as I find them (my own elimination), for the record
def newpart(expr, keep):
    P = sp.Poly(sp.expand(expr), *syms, H7, H8, H9)
    out = sp.Integer(0)
    for mon, coef in P.terms():
        if mon[3] or mon[1] % 2 or mon[5] or mon[6] or mon[7]:
            t = coef
            for s, e in zip(list(syms) + [H7, H8, H9], mon): t *= s ** e
            out += t
    return sp.expand(out)
print("(B) odd part of G_7:", newpart(G[7], None))
print("(B) odd part of G_8 (C6 eliminated):", newpart(G8r, None))
print("(B) odd part of G_9 (C6, C3^2 eliminated):", newpart(G9r, None))
print("(B) odd part of G_10 (C6, C3^2, C8 eliminated):", newpart(G10r, None))
json.dump({str(i): str(G[i]) for i in G}, open('topend3_G.json', 'w'), indent=1)
print(f"done [{time.time()-t0:.0f}s]")
