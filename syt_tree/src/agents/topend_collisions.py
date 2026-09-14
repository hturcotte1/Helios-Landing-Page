"""topend task (3), detail: list every collision (pair of transpose classes of the same n with equal key) for
   key3 = (r, f, u_3)  [equivalently (n, r, f, C_2)]   and   key4 = (r, f, u_3, u_4)  [(n, r, f, C_2, C_1^2)],
and also without r; and check that key5 = (f, u_3, u_4, u_5) [(n, f, C_2, C_1^2, C_4)] has no collision, n <= N.
Level-by-level DP (only two levels of Young's lattice kept in memory)."""
import os, sys, time
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from young import partitions, conjugate, corner_runs, removable_corners, remove_box

def main(N=50, I=5):
    prev = {(): tuple([1] + [0] * I)}
    for n in range(1, N + 1):
        cur = {}
        for lam in partitions(n):
            subs = [prev[remove_box(lam, c)] for c in removable_corners(lam)]
            cur[lam] = tuple(1 if i == n else (0 if i > n else sum(s[i] for s in subs)) for i in range(I + 1))
        prev = cur
        if n < 4: continue
        reps = {}
        for lam in cur:
            key = min(lam, conjugate(lam))
            if key not in reps: reps[key] = (len(removable_corners(lam)), cur[lam])
        for name, mk in (("(r,f,u3)", lambda r, u: (r, u[2], u[3])), ("(f,u3)", lambda r, u: (u[2], u[3])),
                         ("(r,f,u3,u4)", lambda r, u: (r, u[2], u[3], u[4])), ("(f,u3,u4)", lambda r, u: (u[2], u[3], u[4])),
                         ("(f,u3,u4,u5)", lambda r, u: (u[2], u[3], u[4], u[5]))):
            groups = defaultdict(list)
            for lam, (r, u) in reps.items(): groups[mk(r, u)].append(lam)
            bad = [g for g in groups.values() if len(g) > 1]
            if name in ("(r,f,u3,u4)", "(f,u3,u4)", "(f,u3,u4,u5)") or n <= 24:
                for g in bad:
                    print(f"n={n} key{name} collision: " + " ~ ".join(f"{lam} runs={corner_runs(lam)}" for lam in g), flush=True)
            if name == "(f,u3,u4,u5)":
                print(f"n={n}: key (f,u3,u4,u5) collisions: {len(bad)}", flush=True)

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 50)
