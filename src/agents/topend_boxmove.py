"""topend: 'mirror box moves'.  If mu is obtained from lam by removing a corner box of content -c and adding a box of
content +c (c != 0), then C_k(mu) = C_k(lam) for all even k and C_1(mu) = C_1(lam) + 2c... so C_1^2 is preserved iff
C_1(lam) = -c.  We list all such pairs (lam, mu) with equal f^lam = f^mu and equal C_1^2, n <= N: they are
indistinguishable by (n, f, C_2, C_4, C_6, ..., C_1^2) but distinguished by C_1 C_3 (hence by u_7)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from topend_mn import C
from young import partitions, f_hook, removable_corners, addable_corners, remove_box, add_box, conjugate, corner_runs
N = int(sys.argv[1]) if len(sys.argv) > 1 else 50
cnt = 0
for n in range(4, N + 1):
    for lam in partitions(n):
        C1 = C(lam, 1)
        for i in removable_corners(lam):
            c_rem = (lam[i] - 1) - i
            if c_rem != -C1 or c_rem == 0: continue     # need C_1(lam) = -c_rem so that C_1(mu) = +c_rem... wait: C_1(mu) = C1 - c_rem + (-c_rem) = C1 - 2 c_rem; C_1(mu)^2 = C1^2 iff C1 = c_rem.
        # (handled below properly)
        for i in removable_corners(lam):
            c_rem = (lam[i] - 1) - i
            if c_rem == 0: continue
            nu = remove_box(lam, i)
            for a in addable_corners(nu):
                c_add = (nu[a] if a < len(nu) else 0) - a
                if c_add != -c_rem: continue
                mu = add_box(nu, a)
                if mu == lam or mu == conjugate(lam) or mu < lam: continue
                if C(mu, 1) ** 2 != C1 ** 2: continue
                if f_hook(mu) != f_hook(lam): continue
                assert C(mu, 2) == C(lam, 2) and C(mu, 4) == C(lam, 4) and C(mu, 6) == C(lam, 6)
                cnt += 1
                print(f"n={n}: {lam} runs={corner_runs(lam)} -> {mu} runs={corner_runs(mu)}: box content {c_rem} -> {c_add}; f={f_hook(lam)}, C_1: {C1} -> {C(mu,1)}, C_1C_3: {C1*C(lam,3)} vs {C(mu,1)*C(mu,3)}", flush=True)
print("total pairs:", cnt)
