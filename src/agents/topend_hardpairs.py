"""topend task (4): pairs of transpose-classes of partitions of n <= N whose d-vectors agree on d_0..d_{n-k}
(i.e. differ only in the top k values u_0..u_{k-1}), for k = 4 (the question asked) and also k = 5, 6, 7, 8
to see the families.  Uses the memoised d_vector (Young-lattice recursion)."""
import os, sys, time
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from young import partitions, conjugate, corner_runs
from census import d_vector

def main(N=30, KS=(4, 5, 6, 7, 8)):
    for n in range(4, N + 1):
        reps = {}
        for lam in partitions(n):
            key = min(lam, conjugate(lam))
            if key not in reps: reps[key] = d_vector(lam)
        for k in KS:
            if n - k < 0: continue
            groups = defaultdict(list)
            for lam, dv in reps.items(): groups[dv[:n - k + 1]].append(lam)
            bad = [g for g in groups.values() if len(g) > 1]
            if bad:
                for g in bad:
                    tops = {lam: reps[lam][n - k + 1:] for lam in g}
                    print(f"n={n} agree on d_0..d_(n-{k}): {g}   tops (u_{k-1}..u_0): { {l: tuple(reversed(t)) for l, t in tops.items()} }", flush=True)
        if n % 5 == 0: print(f"  ... n={n} done", flush=True)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    main(N)
