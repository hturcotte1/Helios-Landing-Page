"""Brute-force Lemma 5.1 and 5.2 of topend.md: all sequences of m transpositions on N points."""
import itertools, sys
from collections import Counter
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
MMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 5
trans = list(itertools.combinations(range(N), 2))
def compose(p, q): return tuple(p[q[x]] for x in range(N))
def tperm(t):
    p = list(range(N)); p[t[0]], p[t[1]] = p[t[1]], p[t[0]]; return tuple(p)
def cycles_on(g, V):
    seen = set(); cyc = []
    for x in V:
        if x not in seen:
            c = []; j = x
            while j not in seen: seen.add(j); c.append(j); j = g[j]
            cyc.append(c)
    return cyc
cnt = 0; eq = 0; bad = 0
for m in range(0, MMAX + 1):
    for seq in itertools.product(trans, repeat=m):
        g = tuple(range(N))
        for t in seq: g = compose(g, tperm(t))
        V = set(x for t in seq for x in t)
        cyc = cycles_on(g, V)
        c = len(cyc); phi = sum(1 for x in V if g[x] == x)
        cnt += 1
        if m < len(V) - c + phi: bad += 1; print('5.1 FAILS', seq)
        if m == len(V) - c:
            eq += 1
            # Lemma 5.2: phi = 0, forest w/o multiple edges, each component's product a single cycle on its vertex set
            if phi != 0: bad += 1; print('5.2 phi', seq)
            edges = Counter(seq)
            if any(v > 1 for v in edges.values()): bad += 1; print('5.2 multi-edge', seq)
            # components via union-find
            par = {x: x for x in V}
            def find(x):
                while par[x] != x: x = par[x]
                return x
            for a, b in seq: par[find(a)] = find(b)
            comps = {}
            for x in V: comps.setdefault(find(x), set()).add(x)
            ne = Counter(find(a) for a, b in seq)
            for r, Vi in comps.items():
                if ne[r] != len(Vi) - 1: bad += 1; print('5.2 not tree', seq)
                gi = tuple(range(N))
                for t in seq:
                    if find(t[0]) == r: gi = compose(gi, tperm(t))
                cy = cycles_on(gi, Vi)
                if len(cy) != 1 or set(cy[0]) != Vi: bad += 1; print('5.2 not single cycle', seq)
    print(f'm={m} done, sequences so far {cnt}, equality cases {eq}, bad {bad}', flush=True)
print('N', N, 'MMAX', MMAX, 'total', cnt, 'equality', eq, 'bad', bad)
