"""ref2_topend1: brute-force Lemma 5.1 (m >= |V| - c_V(g) + phi) and Lemma 5.2 (equality => phi=0, forest without
multiple edges, each component's product is a single cycle on its vertex set) over all sequences of m transpositions
on N points, for (N,m) in a list."""
import sys
from itertools import product, combinations
def run(N, M):
    trans = list(combinations(range(N), 2))
    cnt = 0; viol51 = 0; viol52 = 0
    for m in range(0, M + 1):
        for seq in product(trans, repeat=m):
            cnt += 1
            g = list(range(N))
            # product t_1 ... t_m composed right to left: apply t_m first
            for (a, b) in reversed(seq):
                g[a], g[b] = g[b], g[a]
            # careful: this computes g = t_1 t_2 ... t_m as functions? Let's define via composition explicitly instead
            # rebuild: g(x) = t_1(t_2(...t_m(x)))
            g = list(range(N))
            for (a, b) in seq[::-1]:
                # apply transposition (a b) after current g: newg = t o g
                g = [b if y == a else a if y == b else y for y in g]
            V = set(x for t in seq for x in t)
            # cycles on V
            seen = set(); c = 0; phi = 0
            for x in V:
                if x in seen: continue
                c += 1; y = x; L = 0
                while y not in seen: seen.add(y); y = g[y]; L += 1
                if L == 1: phi += 1
            if m < len(V) - c + phi: viol51 += 1
            if m == len(V) - c:
                # lemma 5.2 claims
                ok = phi == 0 and len(set(seq)) == m
                # components
                par = {x: x for x in V}
                def f(x):
                    while par[x] != x: par[x] = par[par[x]]; x = par[x]
                    return x
                for a, b in seq: par[f(a)] = f(b)
                comps = {}
                for x in V: comps.setdefault(f(x), set()).add(x)
                for comp in comps.values():
                    e = sum(1 for a, b in seq if f(a) == f(next(iter(comp))))
                    if e != len(comp) - 1: ok = False
                    # product of component's transpositions in order
                    h = list(range(N))
                    for (a, b) in [t for t in seq if f(t[0]) in [f(next(iter(comp)))]][::-1]:
                        h = [b if y == a else a if y == b else y for y in h]
                    # single cycle on comp
                    y = next(iter(comp)); L = 0; s = set()
                    while y not in s: s.add(y); y = h[y]; L += 1
                    if s != comp: ok = False
                if not ok: viol52 += 1
    print(f"N={N} m<={M}: {cnt} sequences, Lemma5.1 violations={viol51}, Lemma5.2 violations={viol52}", flush=True)
for N, M in [(4, 6), (5, 5), (6, 5), (7, 4)]:
    run(N, M)
