"""topend: analyse specific pairs (the depth-5 collisions at n = 49, 50): r, f, contents data, u_0..u_10 and the
first indices (from the bottom and from the top) where the d-vectors differ."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from topend_mn import C, r_of
from young import f_hook, corner_runs
from census import d_vector
PAIRS = [((9, 8, 8, 6, 3, 3, 3, 3, 3, 1, 1, 1), (11, 7, 5, 5, 4, 4, 4, 4, 3, 1, 1)),
         ((10, 10, 5, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2), (12, 7, 5, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1)),
         ((10, 10, 8, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2), (12, 11, 11, 3, 3, 2, 2, 1, 1, 1, 1, 1, 1)),
         ((12, 10, 8, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1), (12, 11, 8, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1))]
for lam, mu in PAIRS:
    n = sum(lam); assert sum(mu) == n
    dl, dm = d_vector(lam), d_vector(mu)
    print(f"n={n}: {lam} runs={corner_runs(lam)} r={r_of(lam)}  vs  {mu} runs={corner_runs(mu)} r={r_of(mu)}")
    print("   f:", f_hook(lam), f_hook(mu), " C_2:", C(lam, 2), C(mu, 2), " C_1:", C(lam, 1), C(mu, 1), " C_4:", C(lam, 4), C(mu, 4))
    print("   C_6-16C_1C_3:", C(lam, 6) - 16 * C(lam, 1) * C(lam, 3), C(mu, 6) - 16 * C(mu, 1) * C(mu, 3), " C_3^2:", C(lam, 3) ** 2, C(mu, 3) ** 2)
    print("   u_0..u_10 lam:", [dl[n - i] for i in range(11)])
    print("   u_0..u_10 mu :", [dm[n - i] for i in range(11)])
    diff = [j for j in range(n + 1) if dl[j] != dm[j]]
    print(f"   d-vectors differ exactly at j in {diff[:6]}...{diff[-6:]} (count {len(diff)}); first from bottom j={diff[0]}, first from top i=n-j={n - diff[-1]}")
