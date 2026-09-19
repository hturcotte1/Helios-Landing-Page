"""topend: level witnesses via the u-table (n <= 40): pairs with equal u_3/f..u_{i-1}/f but different u_i/f."""
import os, sys
from fractions import Fraction
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from topend_mn import u_table
from young import partitions, conjugate
N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
U = u_table(N, 10)
found = {}
for n in range(10, N + 1):
    reps = {}
    for lam in partitions(n):
        key = min(lam, conjugate(lam))
        if key not in reps:
            u = U[lam]; reps[key] = tuple(Fraction(u[i], u[0]) for i in range(11))
    for i in range(7, 11):
        if i in found: continue
        groups = defaultdict(list)
        for lam, rt in reps.items(): groups[rt[3:i]].append(lam)
        for g in groups.values():
            if len(g) > 1 and len({reps[l][i] for l in g}) > 1:
                found[i] = n; print(f"level {i}: first witness at n={n}: {g[:2]} u_{i}/f = {[str(reps[l][i]) for l in g[:2]]}", flush=True); break
    # also: count of groups with equal ratios u_3..u_6 (i.e. equal (n,C_2,C_1^2,C_4)) and their sizes
print("found:", found)
