"""Test: among ALL pairs of partitions of n (n <= NMAX) with the same multiset of |content| and equal f,
is C_1(lam)^2 == C_1(mu)^2 always?  Also record the ratio f^mu/f^lam structure for single moves."""
import sys
from collections import defaultdict, Counter
from topend3_lib import *
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
tot_pairs = 0; bad = []
for n in range(2, NMAX + 1):
    B = defaultdict(list)
    for lam in partitions(n):
        key = tuple(sorted(Counter(abs(c) for c in contents(lam)).items()))
        B[(key, f_hook(lam))].append(lam)
    for (key, f), L in B.items():
        if len(L) < 2: continue
        # distinct transpose classes
        classes = {}
        for lam in L:
            classes[min(lam, conjugate(lam))] = C(lam, 1) ** 2
        vals = list(classes.items())
        for a in range(len(vals)):
            for b in range(a + 1, len(vals)):
                tot_pairs += 1
                if vals[a][1] != vals[b][1]:
                    bad.append((n, vals[a], vals[b]))
    print(f"n={n}: cumulative pairs (same |c|-multiset, same f, distinct classes) = {tot_pairs}, with different C_1^2: {len(bad)}", flush=True)
for b in bad[:30]: print("  different C_1^2:", b)
