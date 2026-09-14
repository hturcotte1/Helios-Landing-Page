"""topend: does (n, C_2) separate class C = two-row + two-column + hooks up to transpose?  And (n, C_2, f)?  n <= N."""
import os, sys
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from topend_mn import C
from young import f_hook, conjugate
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
for n in range(4, N + 1):
    shapes = set()
    for b in range(1, n // 2 + 1): shapes.add((n - b, b))           # two-row (incl. hooks (n-1,1))
    for d in range(1, n - 1): shapes.add((n - d,) + (1,) * d)        # hooks
    reps = {}
    for lam in shapes:
        key = min(lam, conjugate(lam)); reps[key] = lam
    g2 = defaultdict(list); g2f = defaultdict(list)
    for lam in reps.values():
        g2[C(lam, 2)].append(lam); g2f[(C(lam, 2), f_hook(lam))].append(lam)
    c2 = [g for g in g2.values() if len(g) > 1]; c2f = [g for g in g2f.values() if len(g) > 1]
    if c2 or c2f:
        print(f"n={n}: (n,C_2) collisions {c2[:3]} ; (n,C_2,f) collisions {c2f[:3]}")
print("done to", N)
