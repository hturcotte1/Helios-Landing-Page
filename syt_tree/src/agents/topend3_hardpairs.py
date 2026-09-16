"""Task (4), own recomputation: all pairs of transpose classes of partitions of n that agree on
d_0..d_{n-k} for k = 4, 5, 6 (i.e. differ only in the top k-1 values u_0..u_{k-2}), n <= NMAX."""
import sys, time
from collections import defaultdict
from topend3_lib import *
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 32
t0 = time.time()
for n in range(4, NMAX + 1):
    classes = {}
    for lam in partitions(n):
        lt = conjugate(lam)
        key = min(lam, lt)
        if key in classes: continue
        classes[key] = d_vector(lam)
    for k in (4, 5, 6):
        buckets = defaultdict(list)
        for lam, dv in classes.items():
            buckets[dv[:n - k + 1]].append(lam)
        groups = [g for g in buckets.values() if len(g) >= 2]
        if groups:
            for g in groups:
                print(f"n={n} agree on d_0..d_(n-{k}): {g}   d-vectors differ at top: " +
                      "; ".join(str(classes[l][n - k + 1:]) for l in g))
    print(f"n={n} done [{time.time()-t0:.0f}s]", flush=True)
