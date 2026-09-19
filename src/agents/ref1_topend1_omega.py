"""Referee check of topend.md section 5.3/5.4: independent MN characters, verify the 16 Omega polynomials
(parsed from the report text) for all partitions of n <= NMAX, and independently certify the rank claim."""
import os
import sys, re, random
from fractions import Fraction
from functools import lru_cache
from math import factorial
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import sympy
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 18

def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0:
        yield (); return
    for k in range(min(n, maxp), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

def hook_f(lam):
    n = sum(lam)
    if n == 0: return 1
    lamc = [sum(1 for p in lam if p > j) for j in range(lam[0])]
    h = 1
    for i, p in enumerate(lam):
        for j in range(p):
            h *= (p - j - 1) + (lamc[j] - i - 1) + 1
    return factorial(n) // h

def beta(lam):
    l = len(lam)
    return tuple(sorted(lam[i] + (l - 1 - i) for i in range(l)))

def lam_from_beta(b):
    b = sorted(b)
    l = len(b)
    lam = [b[i] - i for i in range(l)]
    return tuple(sorted([x for x in lam if x > 0], reverse=True))

@lru_cache(maxsize=None)
def chi_beta(b, rho):
    """character of shape with beta-set b at cycle type rho (tuple, parts >= 2, decreasing); remaining boxes are 1-cycles."""
    if not rho:
        return hook_f(lam_from_beta(b))
    r = rho[0]; rest = rho[1:]
    total = 0
    bs = set(b)
    for x in b:
        y = x - r
        if y >= 0 and y not in bs:
            sign = (-1) ** sum(1 for z in b if y < z < x)
            nb = tuple(sorted((bs - {x}) | {y}))
            total += sign * chi_beta(nb, rest)
    return total

def chi(lam, rho):
    return chi_beta(beta(lam), tuple(rho))

def zee(rho):
    from collections import Counter
    z = 1
    for k, m in Counter(rho).items():
        z *= k ** m * factorial(m)
    return z

def falling(n, k):
    r = 1
    for i in range(k): r *= (n - i)
    return r

def contents(lam):
    return [j - i for i, p in enumerate(lam) for j in range(p)]

def Ck(lam, k):
    return sum(c ** k for c in contents(lam))

# --- parse report table ---
txt = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'agent_notes/topend.md')).read()
start = txt.index('### 5.4'); block = txt[start:]
block = block[block.index('```') + 3:]; block = block[:block.index('```')]
lines = block.strip().split('\n')
entries = {}
cur = None
for ln in lines:
    if ln.startswith('Ω_'):
        lhs, rhs = ln.split('=', 1)
        key = lhs.strip()[2:]
        if key == '∅': rho = ()
        else: rho = tuple(int(x) for x in key.strip('()').split(','))
        cur = rho; entries[rho] = rhs.strip()
    else:
        entries[cur] += ' ' + ln.strip()
syms = {f'C{k}': sympy.Symbol(f'C{k}') for k in range(1, 11)}; syms['n'] = sympy.Symbol('n')
tr = standard_transformations + (implicit_multiplication_application, convert_xor)
polys = {}
for rho, s in entries.items():
    e = parse_expr(re.sub(r'C(\d+)(?=[Cn])', r'C\1*', s), local_dict=syms, transformations=tr)
    polys[rho] = sympy.Poly(sympy.expand(e), *[syms['n']] + [syms[f'C{k}'] for k in range(1, 11)])
print('parsed', len(polys), 'polynomials:', list(polys))

def m_of(rho): return sum(rho) - len(rho)

# weighted degree check
for rho, P in polys.items():
    m = m_of(rho); maxw = 0
    for mon, c in P.terms():
        w = mon[0] + sum((k + 1) * e for k, e in enumerate(mon[1:]))
        maxw = max(maxw, w)
        idx_used = [k + 1 for k, e in enumerate(mon[1:]) if e]
        assert max(idx_used + [0]) <= m, (rho, mon)
    print(f'rho={rho} m={m} max weight={maxw} terms={len(P.terms())}')
    assert maxw <= m

def eval_poly(P, n, C):
    tot = Fraction(0)
    for mon, c in P.terms():
        v = Fraction(int(c.p), int(c.q)) * n ** mon[0]
        for k, e in enumerate(mon[1:]):
            if e: v *= C[k + 1] ** e
        tot += v
    return tot

# --- verify on all partitions n <= NMAX ---
bad = 0; checks = 0
for n in range(0, NMAX + 1):
    for lam in partitions(n):
        f = hook_f(lam)
        C = {k: Ck(lam, k) for k in range(1, 11)}
        for rho, P in polys.items():
            s = sum(rho)
            if s > n:
                omega = Fraction(0)
            else:
                omega = Fraction(falling(n, s) * chi(lam, rho), zee(rho) * f)
            val = eval_poly(P, n, C)
            checks += 1
            if val != omega:
                bad += 1
                if bad < 20: print('MISMATCH', n, lam, rho, omega, val)
    print(f'n={n}: cumulative checks={checks} bad={bad}', flush=True)
print('TOTAL mismatches:', bad)

# --- rank certificate: evaluation matrix of V_m on partitions n <= N ---
def monomials(m):
    mons = []
    for k in range(0, m + 1):
        for mu in partitions(k):
            for j in range(0, m + 1 - k):
                mons.append((j, mu))
    return mons

def rank_mod(rows, p):
    M = [[x % p for x in r] for r in rows]
    rank = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = None
        for r in range(rank, len(M)):
            if M[r][c]: piv = r; break
        if piv is None: continue
        M[rank], M[piv] = M[piv], M[rank]
        inv = pow(M[rank][c], p - 2, p)
        M[rank] = [(x * inv) % p for x in M[rank]]
        for r in range(len(M)):
            if r != rank and M[r][c]:
                fac = M[r][c]
                M[r] = [(a - fac * b) % p for a, b in zip(M[r], M[rank])]
        rank += 1
    return rank

for m in (2, 4, 6, 8):
    mons = monomials(m)
    for N in (15, 16):
        rows = []
        for n in range(0, N + 1):
            for lam in partitions(n):
                C = {k: Ck(lam, k) for k in range(1, m + 1)}
                row = []
                for j, mu in mons:
                    v = n ** j
                    for t in mu: v *= C[t]
                    row.append(v)
                rows.append(row)
        rk = [rank_mod(rows, p) for p in (2**61 - 1, 1000000007)]
        print(f'm={m}: dim V_m={len(mons)}, N={N}, rows={len(rows)}, rank mod (2^61-1, 1e9+7) = {rk}', flush=True)
