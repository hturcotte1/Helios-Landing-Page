"""K_n = smallest K such that (d_0,...,d_K) already separates all transpose classes of partitions of n
(equivalently, by Theorem 2.3, the depth to which the up-census must be known)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from young import partitions, conjugate
from census import d_vector
for n in range(3, int(sys.argv[1]) + 1 if len(sys.argv) > 1 else 46):
    reps = {}
    for lam in partitions(n):
        key = min(lam, conjugate(lam))
        if key not in reps: reps[key] = d_vector(lam)
    ncl = len(reps)
    K = None
    for j in range(0, n + 1):
        if len({dv[:j + 1] for dv in reps.values()}) == ncl:
            K = j; break
    print(n, ncl, K, "n-K =", n - K if K is not None else None, flush=True)
