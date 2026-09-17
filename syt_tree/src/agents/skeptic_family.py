"""skeptic_family.py -- linear both-ends families found by skeptic_ends.py.
H_k:  lam(k) = (5,3^k,2,2,2)   mu(k) = (4,4,3^k,1,1,1)     n = 3k+11
H'_k: lam(k) = (6,4,3^k,2^4,1) mu(k) = (5,4,4,3^k,2,1^4)   n = 3k+19
H''_k: lam(k) = (6,4^k,3,3,2)  mu(k) = (5,5,4^k,2,1,1)     n = 4k+14
For k = 1..KMAX: exact d-vectors, first/last differing index, Hamming distance, f, C_1, C_2, C_1^2, C_4,
and the hook-length multisets (to confirm the hand proof that f agrees).
Also the reduced pairs (3^K,2,2,2) vs (3^K,1,1,1) and (3^K,2^4,1) vs (3^K,2,1^4): first differing index vs K.
"""
import sys
from math import factorial
sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src/agents')
from skeptic_dvec import transpose, removable, remove_box
memo = {(): (1,)}
def dvec(lam):
    v = memo.get(lam)
    if v is not None: return v
    n = sum(lam); vec = [1] + [0]*n
    for i in removable(lam):
        sub = dvec(remove_box(lam, i))
        for j in range(n): vec[j+1] += sub[j]
    memo[lam] = tuple(vec); return memo[lam]
def hooks(lam):
    lt = transpose(lam)
    return sorted(lam[i] - j + lt[j] - i - 1 for i in range(len(lam)) for j in range(lam[i]))
def cont(lam):
    cs = [j - i for i, p in enumerate(lam) for j in range(p)]
    return sum(cs), sum(c*c for c in cs), sum(c**4 for c in cs), sum(c**3 for c in cs)
fams = {
 'H':   lambda k: ((5,) + (3,)*k + (2,2,2), (4,4) + (3,)*k + (1,1,1)),
 "H'":  lambda k: ((6,4) + (3,)*k + (2,2,2,2,1), (5,4,4) + (3,)*k + (2,1,1,1,1)),
 "H''": lambda k: ((6,) + (4,)*k + (3,3,2), (5,5) + (4,)*k + (2,1,1)),
}
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
for name, F in fams.items():
    print("family", name)
    for k in range(1, KMAX + 1):
        lam, mu = F(k); n = sum(lam); assert sum(mu) == n
        assert transpose(lam) != mu and lam != mu
        hl, hm = hooks(lam), hooks(mu)
        dl, dm = dvec(lam), dvec(mu)
        diff = [j for j in range(n + 1) if dl[j] != dm[j]]
        cl, cm = cont(lam), cont(mu)
        print(f"  k={k:2d} n={n:3d} hooks_equal={hl==hm} f_equal={dl[n]==dm[n]} C1={cl[0]},{cm[0]} C2eq={cl[1]==cm[1]} "
              f"C4eq={cl[2]==cm[2]} C3={cl[3]},{cm[3]} first_diff={diff[0]} last_diff={diff[-1]} ham={len(diff)} "
              f"agree_bottom=d_0..d_{diff[0]-1} agree_top=d_{diff[-1]+1}..d_{n}", flush=True)
    memo.clear()
print("reduced pairs")
for K in range(1, 16):
    for (b1, b2) in (((2,2,2), (1,1,1)), ((2,2,2,2,1), (2,1,1,1,1)), ((3,3,2),(2,1,1))):
        c = 3 if len(b1) != 3 or b1[0] < 3 else 4
        lam, mu = (c,)*K + b1, (c,)*K + b2
        dl, dm = dvec(lam), dvec(mu); n = sum(lam)
        diff = [j for j in range(n+1) if dl[j] != dm[j]]
        print(f"  K={K:2d} {c}^K|{b1} vs {c}^K|{b2}: n={n} first_diff={diff[0] if diff else None} ham={len(diff)}", flush=True)
