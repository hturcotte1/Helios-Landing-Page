"""Referee check of Theorem P (topend.md sec 5) directly in Z[S_n], n <= NMAX (default 8), no characters used.
(1) p_k(J) computed in the group algebra, centrality checked, expanded in class sums.
(2) For all mu with |mu| <= MMAX: p_mu(J) expansion coefficients a_{mu,rho'}(n) checked against Theorem P (ii),(iii)
    and the polynomial-degree bound via finite differences where enough n are available.
(3) K_rho(n) = Omega_rho(n; p(J)) verified in Z[S_n] for the 16 polynomials of 5.4 (parsed from the report)."""
import os
import sys, itertools, re
from fractions import Fraction
from collections import Counter, defaultdict
from math import factorial
import sympy
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
MMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 6

def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0:
        yield (); return
    for k in range(min(n, maxp), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

def compose(p, q):  # (pq)(x) = p(q(x))
    return tuple(p[q[x]] for x in range(len(p)))

def inverse(p):
    r = [0] * len(p)
    for i, x in enumerate(p): r[x] = i
    return tuple(r)

def ctype(p):
    n = len(p); seen = [False] * n; cyc = []
    for i in range(n):
        if not seen[i]:
            l = 0; j = i
            while not seen[j]:
                seen[j] = True; j = p[j]; l += 1
            cyc.append(l)
    return tuple(sorted(cyc, reverse=True))

def strip1(t): return tuple(x for x in t if x > 1)
def m_of(rho): return sum(rho) - len(rho)

# --- parse polynomials from report ---
txt = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'agent_notes/topend.md')).read()
start = txt.index('### 5.4'); block = txt[start:]
block = block[block.index('```') + 3:]; block = block[:block.index('```')]
entries = {}; cur = None
for ln in block.strip().split('\n'):
    if ln.startswith('Ω_'):
        lhs, rhs = ln.split('=', 1); key = lhs.strip()[2:]
        rho = () if key == '∅' else tuple(int(x) for x in key.strip('()').split(','))
        cur = rho; entries[rho] = rhs.strip()
    else: entries[cur] += ' ' + ln.strip()
syms = {f'C{k}': sympy.Symbol(f'C{k}') for k in range(1, 11)}; syms['n'] = sympy.Symbol('n')
tr = standard_transformations + (implicit_multiplication_application, convert_xor)
gens = [syms['n']] + [syms[f'C{k}'] for k in range(1, 11)]
polys = {rho: sympy.Poly(sympy.expand(parse_expr(re.sub(r'C(\d+)(?=[Cn])', r'C\1*', s), local_dict=syms, transformations=tr)), *gens) for rho, s in entries.items()}

grand_ok = True
for n in range(1, NMAX + 1):
    classes = list(partitions(n)); cidx = {c: i for i, c in enumerate(classes)}; nc = len(classes)
    perms = list(itertools.permutations(range(n)))
    bytype = defaultdict(list)
    for p in perms: bytype[ctype(p)].append(p)
    rep = {c: bytype[c][0] for c in classes}
    # structure table T[b][c] = Counter over a of #{h in K_b : type(g_c h^-1) = a}
    T = [[None] * nc for _ in range(nc)]
    for b, cb in enumerate(classes):
        Kb_inv = [inverse(h) for h in bytype[cb]]
        for c, cc in enumerate(classes):
            g = rep[cc]; cnt = [0] * nc
            for hi in Kb_inv: cnt[cidx[ctype(compose(g, hi))]] += 1
            T[b][c] = cnt
    def mul(z, w):  # class vectors (lists of Fractions/ints)
        out = [Fraction(0)] * nc
        for b in range(nc):
            if w[b] == 0: continue
            for c in range(nc):
                s = 0
                Tb = T[b][c]
                for a in range(nc):
                    if z[a] and Tb[a]: s += z[a] * Tb[a]
                out[c] += w[b] * s
        return out
    def scal(x): v = [Fraction(0)] * nc; v[cidx[tuple([1] * n)]] = Fraction(x); return v
    # sanity: K_2 * K_2 = C(n,2) + 3 K_3 + 2 K_22 (I3)
    def Kvec(rho):
        v = [Fraction(0)] * nc
        if sum(rho) <= n: v[cidx[tuple(sorted(list(rho) + [1] * (n - sum(rho)), reverse=True))]] = Fraction(1)
        return v
    ident = scal(1)
    # --- p_k(J) in the group algebra ---
    kmax = max(MMAX, 8)
    pk = {}
    for k in range(1, kmax + 1):
        acc = defaultdict(int)
        for i in range(n):  # J_{i+1} = sum_{j<i} (j i)
            trans = []
            for j in range(i):
                t = list(range(n)); t[j], t[i] = t[i], t[j]; trans.append(tuple(t))
            cur = {tuple(range(n)): 1}
            for _ in range(k):
                new = defaultdict(int)
                for x, cx in cur.items():
                    for t in trans: new[compose(x, t)] += cx
                cur = new
            for x, cx in cur.items(): acc[x] += cx
        # centrality check and conversion
        v = [None] * nc
        for x, cx in acc.items():
            a = cidx[ctype(x)]
            if v[a] is None: v[a] = cx
            elif v[a] != cx: print('NOT CENTRAL', n, k); grand_ok = False
        for a in range(nc):
            if v[a] is None:
                # class not present -> coefficient 0, but verify no member present (acc default)
                v[a] = 0
        # verify all members of each class have the same coefficient (including 0)
        for a, c in enumerate(classes):
            for x in bytype[c]:
                if acc.get(x, 0) != v[a]: print('NOT CENTRAL(0)', n, k); grand_ok = False; break
        pk[k] = [Fraction(x) for x in v]
    # sanity identities (I2),(I3)
    assert pk[2] == [a + b for a, b in zip(scal(n * (n - 1) // 2), Kvec((3,)))]
    K2sq = mul(pk[1], pk[1])
    tgt = [a + 3 * b + 2 * c for a, b, c in zip(scal(n * (n - 1) // 2), Kvec((3,)), Kvec((2, 2)))]
    assert K2sq == tgt, (n, K2sq, tgt)
    # --- (3) K_rho = Omega_rho(n; p(J)) ---
    powcache = {}
    def ppow(k, e):
        if (k, e) not in powcache:
            v = ident
            for _ in range(e): v = mul(v, pk[k])
            powcache[(k, e)] = v
        return powcache[(k, e)]
    nbad = 0
    for rho, P in polys.items():
        tot = [Fraction(0)] * nc
        for mon, c in P.terms():
            v = scal(Fraction(int(c.p), int(c.q)) * n ** mon[0])
            for k, e in enumerate(mon[1:]):
                if e: v = mul(v, ppow(k + 1, e))
            tot = [a + b for a, b in zip(tot, v)]
        if tot != Kvec(rho):
            nbad += 1; grand_ok = False
            print(f'  n={n} rho={rho}: Omega(n;p(J)) != K_rho.  diff nonzero at', [(classes[i], tot[i] - Kvec(rho)[i]) for i in range(nc) if tot[i] != Kvec(rho)[i]][:5])
    print(f'n={n}: {len(polys)} identities K_rho = Omega_rho(n;p(J)) in Q[S_n]: {len(polys)-nbad} ok, {nbad} bad', flush=True)
    # --- (2) p_mu(J) expansions, store for Theorem P checks ---
    if 'store' not in globals(): store = {}
    for m in range(1, MMAX + 1):
        for mu in partitions(m):
            v = ident
            for t in mu: v = mul(v, pk[t])
            for a, c in enumerate(classes):
                store[(mu, strip1(c), n)] = v[a]

# --- Theorem P checks over stored data ---
def mult_fact(mu):
    r = 1
    for k, e in Counter(mu).items(): r *= factorial(e)
    return r
viol = 0; tested_deg = 0; tested_iii = 0
mus = [mu for m in range(1, MMAX + 1) for mu in partitions(m)]
rhos = sorted({r for (mu, r, n) in store})
for mu in mus:
    m = sum(mu); l = len(mu)
    for rho in rhos:
        s = sum(rho)
        vals = [(n, store[(mu, rho, n)]) for n in range(max(s, 1), NMAX + 1) if (mu, rho, n) in store]
        if not vals: continue
        if m_of(rho) > m:
            if any(v != 0 for _, v in vals): print('VIOLATION (ii) nonzero:', mu, rho, vals); viol += 1
            continue
        if m_of(rho) == m:
            tested_iii += 1
            vs = {v for _, v in vals}
            if len(vs) != 1: print('VIOLATION (iii) not constant:', mu, rho, vals); viol += 1
            v0 = vals[0][1]
            if len(rho) > l and v0 != 0: print('VIOLATION (iii) l(rho)>l(mu) but nonzero:', mu, rho, v0); viol += 1
            mup1 = tuple(x + 1 for x in mu)
            if rho == mup1 and v0 != mult_fact(mu): print('VIOLATION (iii) leading coeff:', mu, rho, v0, mult_fact(mu)); viol += 1
        d = m - m_of(rho)
        # finite-difference degree test: need >= d+2 consecutive points
        if len(vals) >= d + 2:
            tested_deg += 1
            seq = [v for _, v in vals]
            for _ in range(d + 1):
                seq = [seq[i + 1] - seq[i] for i in range(len(seq) - 1)]
            if any(x != 0 for x in seq): print('VIOLATION (ii) degree bound:', mu, rho, vals); viol += 1
print(f'Theorem P checks: {len(mus)} mu, {len(rhos)} rho; (iii) cases tested {tested_iii}; degree-bound tests {tested_deg}; violations {viol}')
print('GRAND OK' if grand_ok and viol == 0 else 'PROBLEMS FOUND')
