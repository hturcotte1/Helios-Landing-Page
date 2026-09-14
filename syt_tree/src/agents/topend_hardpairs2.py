"""topend task (4), extended: pairs of transpose classes of n <= N agreeing on d_0..d_{n-k}, k in KS;
level-by-level DP on full d-vectors."""
import os, sys
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from young import partitions, conjugate, removable_corners, remove_box

def main(N=40, KS=(4, 5, 6)):
    prev = {(): (1,)}
    for n in range(1, N + 1):
        cur = {}
        for lam in partitions(n):
            d = [0] * (n + 1); d[0] = 1
            for c in removable_corners(lam):
                dm = prev[remove_box(lam, c)]
                for j in range(1, n + 1): d[j] += dm[j - 1]
            cur[lam] = tuple(d)
        prev = cur
        reps = {}
        for lam in cur:
            key = min(lam, conjugate(lam))
            if key not in reps: reps[key] = cur[lam]
        for k in KS:
            if n - k < 0: continue
            groups = defaultdict(list)
            for lam, dv in reps.items(): groups[dv[:n - k + 1]].append(lam)
            for g in groups.values():
                if len(g) > 1:
                    print(f"n={n} k={k}: {g}", flush=True)
        print(f"  n={n} done", flush=True)

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 40)
