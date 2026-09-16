"""ref2_topend1: (a) parse the 16 Omega_rho polynomials from the report text (section 5.4) and verify
omega_lam(K_rho) == Omega_rho(n; C) for all lam |- n <= NV by my own MN code;
(b) recompute the rank certificate of 5.3: rank of evaluation matrix of {n^j C_mu : j+|mu|<=m} on all lam |- n<=N,
mod two primes (and exactly over Q for m<=6), for N=15,16; also check weighted degree <= m of each Omega."""
import re, sys, time
from fractions import Fraction
from sympy import symbols, sympify, Poly, Rational
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, convert_xor, implicit_multiplication
from ref2_topend1_lib import *

REPORT = '/home/user/Helios-Landing-Page/syt_tree/agent_notes/topend.md'
txt = open(REPORT).read()
blk = txt.split('### 5.4')[1].split('```')[1]
# join continuation lines
lines = []
for ln in blk.split('\n'):
    if not ln.strip(): continue
    if ln.startswith('Ω_'): lines.append(ln.strip())
    else: lines[-1] += ' ' + ln.strip()
n, C1, C2, C3, C4, C5, C6, C7, C8 = symbols('n C1 C2 C3 C4 C5 C6 C7 C8')
loc = dict(n=n, C1=C1, C2=C2, C3=C3, C4=C4, C5=C5, C6=C6, C7=C7, C8=C8)
def parse_rhs(s):
    s = s.replace(' ', '*').replace('*+*', '+').replace('*-*', '-')
    s = re.sub(r'(\d)(?=[Cn])', r'\1*', s)
    s = re.sub(r'(n)(?=[Cn])', r'\1*', s)
    s = s.replace('^', '**')
    return sympify(s, locals=loc)
table = {}
for ln in lines:
    lhs, rhs = ln.split('=', 1)
    key = lhs.strip()[2:]
    if key == '∅': rho = ()
    else: rho = tuple(int(x) for x in key.strip('()').split(','))
    table[rho] = parse_rhs(rhs.strip())
print("parsed", len(table), "polynomials:", list(table))
assert set(table) == set(even_type_no_ones(10))

NV = int(sys.argv[1]) if len(sys.argv) > 1 else 16
# weighted degree check
wt = {n: 1, C1: 1, C2: 2, C3: 3, C4: 4, C5: 5, C6: 6, C7: 7, C8: 8}
for rho, P in table.items():
    m = sum(rho) - len(rho)
    pl = Poly(P, *wt.keys())
    for mon, coef in pl.terms():
        w = sum(e * w_ for e, w_ in zip(mon, wt.values()))
        assert w <= m, (rho, mon, w, m)
        # even number of odd indices
        odd = mon[1] + mon[3] + mon[5] + mon[7]
        assert odd % 2 == 0, (rho, mon)
print("weighted degree <= m and parity OK for all 16")

# evaluate polynomials as exact functions via lambdify-free substitution: convert to dict of monomials
polys = {}
for rho, P in table.items():
    pl = Poly(P, n, C1, C2, C3, C4, C5, C6, C7, C8)
    polys[rho] = [(Fraction(int(c.p), int(c.q)), mon) for mon, c in pl.terms()]
def evalpoly(terms, vals):
    tot = Fraction(0)
    for c, mon in terms:
        v = c
        for e, x in zip(mon, vals):
            if e: v *= x ** e
        tot += v
    return tot
t0 = time.time(); cnt = 0; bad = []
for N in range(0, NV + 1):
    for lam in partitions(N):
        vals = [N] + [C(lam, k) for k in range(1, 9)]
        for rho, terms in polys.items():
            cnt += 1
            if evalpoly(terms, vals) != omega(lam, rho):
                bad.append((lam, rho))
print(f"Omega table verified on all lam |- n<={NV}: {cnt} checks, {len(bad)} failures [{time.time()-t0:.0f}s]", bad[:5])

# rank certificate
def monomials(m):
    return [(j, mu) for s in range(m + 1) for mu in partitions(s) for j in range(m - s + 1)]
def rank_modp(rows, p):
    A = [[x % p for x in r] for r in rows]; ncol = len(A[0]); r = 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(A)) if A[i][c]), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p); A[r] = [(x * inv) % p for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]:
                f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        r += 1
    return r
def rank_exact(rows):
    A = [[Fraction(x) for x in r] for r in rows]; ncol = len(A[0]); r = 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        inv = 1 / A[r][c]; A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
    return r
for m in (2, 4, 6, 8):
    monos = monomials(m)
    for N in (15, 16):
        T = [lam for k in range(N + 1) for lam in partitions(k)]
        rows = []
        for lam in T:
            Cc = {k: C(lam, k) for k in range(0, 9)}
            rows.append([sum(lam) ** j * eval_prod(Cc, mu) if False else (sum(lam) ** j) * __import__('math').prod(Cc[k] for k in mu) for j, mu in monos])
        r1 = rank_modp(rows, (1 << 61) - 1); r2 = rank_modp(rows, 1000000007)
        r3 = rank_exact(rows) if m <= 6 else None
        print(f"m={m} dim={len(monos)} N={N} |T|={len(T)} rank mod 2^61-1 = {r1}, mod 1e9+7 = {r2}, exact = {r3}", flush=True)
