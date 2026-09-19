"""Small-n companion of structure_recover_fast.py: for n <= 20, every ambiguous recovery is exactly {lam, lam^t}."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from structure_recover_fast import rec3, rec4, shapes, Cpow
from young import conjugate
for L, rec, N in ((3, lambda lam, n: rec3(n, Cpow(lam,1)**2, Cpow(lam,2)), 20),
                  (4, lambda lam, n: rec4(n, Cpow(lam,1)**2, Cpow(lam,2), Cpow(lam,4)), 20)):
    bad = []; amb = []
    for n in range(1, N+1):
        if L == 4 and n == 8: continue
        for lam in shapes(n, L):
            r = rec(lam, n)
            if r != [lam]:
                amb.append(n)
                if sorted(r) != sorted({lam, conjugate(lam)}): bad.append((lam, r))
    print(f"L={L}, n<={N} (n=8 skipped for L=4): ambiguous cases occur for n in {sorted(set(amb))}; cases not equal to {{lam, lam^t}}: {bad}")
