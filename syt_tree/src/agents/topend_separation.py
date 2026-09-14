"""topend task (3): which top data separate partitions up to transpose?
For all partitions of n <= N compute r and the top vector (u_0..u_I) by the Young-lattice DP (no characters),
then for each n and each depth i report the number of collisions (pairs of transpose-classes with equal key)
for the keys  K_i = (r, f, u_3, ..., u_i)  and  K'_i = (f, u_3, ..., u_i)  (no r) and  K''_i = (u_3..u_i)/f ratios only.
Also: minimal i such that K_i separates all transpose classes of n ('top separating depth')."""
import os, sys, time
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from topend_mn import u_table, r_of
from young import partitions, conjugate

def main(N=40, I=10):
    t0 = time.time()
    U = u_table(N, I)
    print(f"u-table built for n <= {N}, i <= {I} in {time.time()-t0:.0f}s", flush=True)
    for n in range(4, N + 1):
        reps = {}
        for lam in partitions(n):
            key = min(lam, conjugate(lam))
            if key not in reps: reps[key] = (r_of(lam), U[lam])
        ncl = len(reps)
        line = f"n={n:2d} classes={ncl:6d}"
        depths = {}
        for variant, mk in (("with r", lambda r, u, i: (r,) + u[2:i + 1]), ("no r", lambda r, u, i: u[2:i + 1])):
            first = None; coll = []
            for i in range(3, I + 1):
                groups = defaultdict(list)
                for lam, (r, u) in reps.items(): groups[mk(r, u, i)].append(lam)
                ncoll = sum(len(g) - 1 for g in groups.values() if len(g) > 1)
                coll.append(ncoll)
                if ncoll == 0 and first is None: first = i
            line += f" | {variant}: collisions at i=3..{I}: {coll} sep.depth={first}"
            depths[variant] = (first, coll)
        print(line, flush=True)
        # show the surviving collisions at i = I, if any
        groups = defaultdict(list)
        for lam, (r, u) in reps.items(): groups[(r,) + u[2:I + 1]].append(lam)
        for g in groups.values():
            if len(g) > 1: print("    unresolved at i=%d (with r): %s" % (I, g), flush=True)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    main(N)
