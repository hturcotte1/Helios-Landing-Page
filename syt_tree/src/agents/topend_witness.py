"""topend: function-level independence of the levels.  For i = 6..10 find lam, mu |- n (n <= N), not transposes,
with equal (u_3/f, ..., u_{i-1}/f) (equivalently equal g_3..g_{i-1}) but different u_i/f.  For i = 6 none can exist
(Theorem C); for i = 7..10 we exhibit the smallest n."""
import os, sys
from fractions import Fraction
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from young import partitions, conjugate
from census import d_vector
N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
found = {i: None for i in range(6, 11)}
for n in range(10, N + 1):
    reps = {}
    for lam in partitions(n):
        key = min(lam, conjugate(lam))
        if key not in reps:
            dv = d_vector(lam); reps[key] = tuple(Fraction(dv[n - i], dv[n]) for i in range(11))
    for i in range(6, 11):
        if found[i] is not None: continue
        groups = defaultdict(list)
        for lam, ratios in reps.items(): groups[ratios[3:i]].append(lam)
        for g in groups.values():
            if len(g) > 1 and len({reps[l][i] for l in g}) > 1:
                found[i] = (n, g[:2], [reps[l][i] for l in g[:2]])
                print(f"level {i}: n={n}: {g[:2]} have equal u_3/f..u_{i-1}/f but u_{i}/f = {[str(reps[l][i]) for l in g[:2]]}", flush=True)
                break
    if all(found[i] is not None for i in range(7, 11)): break
print("level 6 witness (should be None):", found[6])
