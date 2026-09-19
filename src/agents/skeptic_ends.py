"""skeptic_ends.py -- both-ends agreement grid. For n <= NMAX and a grid of (a, s):
count pairs of transpose classes agreeing on d_0..d_{a-1} AND on d_{n-s+1}..d_n; list the pairs with
the largest a among those with s >= S0 (S0 = 4: f and u_3 agree; S0 = 5: f,u_3,u_4 agree)."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from skeptic_dvec import levels, transpose
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
A = [3, 4, 5, 6, 8, 10, 12, 14]
S = [3, 4, 5, 7, 8]
for n, D in levels(NMAX):
    if n < 12: continue
    reps = sorted(lam for lam in D if lam <= transpose(lam))
    V = [tuple(D[lam]) for lam in reps]
    line = []
    for s in S:
        for a in A:
            if a + s > n + 1: line.append(f"s{s}a{a}:-"); continue
            g = {}
            for v in V:
                k = (v[:a], v[n - s + 1:]); g[k] = g.get(k, 0) + 1
            line.append(f"s{s}a{a}:{sum(c*(c-1)//2 for c in g.values())}")
    print(f"n={n} " + " ".join(line), flush=True)
    for s0 in (4, 5):
        # largest a with some pair agreeing on top s0 and bottom a: group by top s0, within group longest common prefix
        best = (-1, None)
        g = {}
        for i, v in enumerate(V): g.setdefault(v[n - s0 + 1:], []).append(i)
        for L in g.values():
            if len(L) < 2: continue
            L.sort(key=lambda i: V[i])
            for x in range(len(L) - 1):
                i, j = L[x], L[x + 1]
                k = 0
                while k <= n and V[i][k] == V[j][k]: k += 1
                if k > best[0]: best = (k, (i, j))
        if best[1]:
            i, j = best[1]
            diff = [q for q in range(n + 1) if V[i][q] != V[j][q]]
            print(f"   top>={s0}: max prefix a={best[0]}  {reps[i]} vs {reps[j]}  ham={len(diff)} diff={diff}", flush=True)
