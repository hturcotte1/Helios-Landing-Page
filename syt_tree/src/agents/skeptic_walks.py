"""skeptic_walks.py -- 3-row chain symmetry.
N_j(nu) = number of j-step up-chains nu = kappa^0 < kappa^1 < ... < kappa^j in Young's lattice with all l(kappa^i) <= 3
        = sum_{kappa, l(kappa)<=3, |kappa/nu|=j} f^{kappa/nu}.
It depends only on (u,v) = (nu_1-nu_2+1, nu_2-nu_3+1).  Test: N_j(u,v) == N_j(v,u) for u,v <= 8, j <= 24.
Also confirm: d_j((3^K,2,2,2)) = N_j(1,4)?? / N_j(4,1) and d_j((3^K,1,1,1)) = N_j(4,1)?? for j <= K (box-complement + transpose).
"""
import sys
from functools import lru_cache
sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src/agents')
from skeptic_dvec import removable, remove_box
@lru_cache(maxsize=None)
def N(j, u, v):
    if j == 0: return 1
    tot = N(j-1, u+1, v)                       # add box to row 1
    if u >= 2: tot += N(j-1, u-1, v+1)         # add box to row 2
    if v >= 2: tot += N(j-1, u, v-1)           # add box to row 3
    return tot
bad = 0
for u in range(1, 9):
    for v in range(1, 9):
        for j in range(0, 25):
            if N(j, u, v) != N(j, v, u): bad += 1
print("symmetry N_j(u,v)=N_j(v,u) failures:", bad, " (u,v<=8, j<=24)")
memo = {(): (1,)}
def dvec(lam):
    x = memo.get(lam)
    if x is not None: return x
    n = sum(lam); vec = [1] + [0]*n
    for i in removable(lam):
        sub = dvec(remove_box(lam, i))
        for j in range(n): vec[j+1] += sub[j]
    memo[lam] = tuple(vec); return memo[lam]
for K in range(1, 13):
    a, b = dvec((3,)*K + (2,2,2)), dvec((3,)*K + (1,1,1))
    ok1 = all(a[j] == N(j, 1, 4) for j in range(K + 1))
    ok2 = all(b[j] == N(j, 4, 1) for j in range(K + 1))
    ok1b = all(a[j] == N(j, 4, 1) for j in range(K + 1))
    print(f"K={K}: d_j(3^K,2,2,2)==N_j(1,4) for j<=K: {ok1}; ==N_j(4,1): {ok1b}; d_j(3^K,1,1,1)==N_j(4,1): {ok2}")
print("N_j(1,4) for j<=10:", [N(j,1,4) for j in range(11)])
print("N_j(4,1) for j<=10:", [N(j,4,1) for j in range(11)])
