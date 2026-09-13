"""Which small tuples of d-vector-derived invariants already separate non-transpose shapes?
Invariants: n, r=d_1, f=d_n, u3=d_{n-3}, C2 = sum of squared contents (recoverable from f and u3)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from young import partitions, conjugate, boxes, corner_runs
from census import d_vector

def C2(lam):
    return sum((j - i) ** 2 for (i, j) in boxes(lam))

def test(max_n):
    print("n : #classes | collisions for (f) | (f,C2) | (r,f) | (r,f,C2) | (n,r,m,f) | (d_1,d_2,d_3,f)")
    for n in range(4, max_n + 1):
        reps = {}
        for lam in partitions(n):
            key = min(lam, conjugate(lam))
            reps[key] = lam
        stats = {"f": {}, "fC2": {}, "rf": {}, "rfC2": {}, "nrmf": {}, "d123f": {}}
        for key, lam in reps.items():
            dv = d_vector(lam); f = dv[n]; r = dv[1]; c2 = C2(lam)
            m = min(min(a, b) for a, b in corner_runs(lam))
            stats["f"].setdefault(f, []).append(key)
            stats["fC2"].setdefault((f, c2), []).append(key)
            stats["rf"].setdefault((r, f), []).append(key)
            stats["rfC2"].setdefault((r, f, c2), []).append(key)
            stats["nrmf"].setdefault((r, m, f), []).append(key)
            stats["d123f"].setdefault((dv[1], dv[2], dv[3], f), []).append(key)
        counts = {k: sum(len(v) - 1 for v in s.values() if len(v) > 1) for k, s in stats.items()}
        print(n, len(reps), counts)
        if n <= 12:
            for k, v in stats["fC2"].items():
                if len(v) > 1: print("   (f,C2) collision", k, v)

if __name__ == "__main__":
    test(int(sys.argv[1]) if len(sys.argv) > 1 else 30)
