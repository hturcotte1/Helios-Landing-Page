"""topend: collisions of top keys up to N=58, keys (f,u3,u4,u5), (r,f,u3,u4,u5), (f,u3,u4,u5,u7), (r,f,u3,..,u7)."""
import os, sys, time
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from young import partitions, conjugate, corner_runs, removable_corners, remove_box
def main(N=58, I=7):
    prev = {(): tuple([1] + [0] * I)}
    t0 = time.time()
    for n in range(1, N + 1):
        cur = {}
        for lam in partitions(n):
            subs = [prev[remove_box(lam, c)] for c in removable_corners(lam)]
            cur[lam] = tuple(1 if i == n else (0 if i > n else sum(s[i] for s in subs)) for i in range(I + 1))
        prev = cur
        if n < 40: continue
        reps = {}
        for lam in cur:
            key = min(lam, conjugate(lam))
            if key not in reps: reps[key] = (len(removable_corners(lam)), cur[lam])
        out = []
        for name, mk in (("(f,u3,u4,u5)", lambda r, u: (u[2], u[3], u[4], u[5])), ("(r,f,u3,u4,u5)", lambda r, u: (r, u[2], u[3], u[4], u[5])),
                         ("(f,u3,u4,u5,u7)", lambda r, u: (u[2], u[3], u[4], u[5], u[7])), ("(r,f,u3,u4,u5,u7)", lambda r, u: (r, u[2], u[3], u[4], u[5], u[7]))):
            groups = defaultdict(list)
            for lam, (r, u) in reps.items(): groups[mk(r, u)].append(lam)
            bad = [g for g in groups.values() if len(g) > 1]
            out.append(f"{name}: {len(bad)}")
            for g in bad: print(f"n={n} key{name} collision: " + " ~ ".join(f"{lam} runs={corner_runs(lam)}" for lam in g), flush=True)
        print(f"n={n}: collisions " + ", ".join(out) + f"  [{time.time()-t0:.0f}s]", flush=True)
if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 58)
