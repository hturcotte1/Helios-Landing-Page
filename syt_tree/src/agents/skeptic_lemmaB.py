"""skeptic_lemmaB.py -- Lemma B (body insertion) and the reduced pairs.
Lemma B: lam = alpha|(c^k)|beta with alpha_last >= c >= beta_1. alpha' = alpha/(c^{l(alpha)}) (skew diagram).
For j <= k:  d_j(lam) = sum_i C(j,i) D_i(alpha') d_{j-i}((c^k)|beta),
where D_i(alpha') = number of ways to remove i boxes one at a time from the skew diagram alpha'.
Test: all alpha with |alpha| <= 8, beta with |beta| <= 5, c in 1..4, k in 1..5 (exact).
Reduced pairs: (3^K,2,2,2) vs (3^K,1,1,1), (3^K,2,2,2,2,1) vs (3^K,2,1,1,1,1), (4^K,3,3,2) vs (4^K,2,1,1): first differing index.
Family H'' for k >= 2.
"""
import sys
from math import comb
sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src/agents')
from skeptic_dvec import gen_partitions, transpose, removable, remove_box
memo = {(): (1,)}
def dvec(lam):
    v = memo.get(lam)
    if v is not None: return v
    n = sum(lam); vec = [1] + [0]*n
    for i in removable(lam):
        sub = dvec(remove_box(lam, i))
        for j in range(n): vec[j+1] += sub[j]
    memo[lam] = tuple(vec); return memo[lam]
smemo = {}
def dskew(alpha, c):
    """down-census of the skew diagram alpha/(c^{l(alpha)}): remove corners of alpha lying in columns > c"""
    key = (alpha, c)
    if key in smemo: return smemo[key]
    size = sum(max(0, p - c) for p in alpha)
    vec = [1] + [0]*size
    for i in removable(alpha):
        if alpha[i] > c:
            sub = dskew(remove_box(alpha, i), c)
            for j in range(size): vec[j+1] += sub[j]
    smemo[key] = tuple(vec); return smemo[key]

tested = fails = 0
for A in range(0, 9):
    for alpha in gen_partitions(A):
        for B in range(0, 6):
            for beta in gen_partitions(B):
                for c in range(1, 5):
                    if alpha and alpha[-1] < c: continue
                    if beta and beta[0] > c: continue
                    for k in range(1, 6):
                        lam = alpha + (c,)*k + beta
                        body = (c,)*k + beta
                        dl, da, db = dvec(lam), dskew(alpha, c), dvec(body)
                        for j in range(k + 1):
                            pred = sum(comb(j, i) * (da[i] if i < len(da) else 0) * db[j - i] for i in range(j + 1))
                            tested += 1
                            if pred != dl[j]:
                                fails += 1; print("LEMMA B FAILS", alpha, c, k, beta, j, pred, dl[j])
print(f"Lemma B: {tested} instances tested, failures = {fails}", flush=True)
memo.clear()
print("reduced pairs (first differing index; K = number of body rows)")
for K in range(1, 19):
    for c, b1, b2 in ((3, (2,2,2), (1,1,1)), (3, (2,2,2,2,1), (2,1,1,1,1)), (4, (3,3,2), (2,1,1))):
        lam, mu = (c,)*K + b1, (c,)*K + b2
        dl, dm = dvec(lam), dvec(mu); n = sum(lam)
        diff = [j for j in range(n+1) if dl[j] != dm[j]]
        print(f"  K={K:2d} {c}^K|{b1} vs {c}^K|{b2}: n={n} first_diff={diff[0]} (=K+{diff[0]-K}) last_diff={diff[-1]} ham={len(diff)}", flush=True)
memo.clear()
print("family H''")
for k in range(2, 11):
    lam, mu = (6,) + (4,)*k + (3,3,2), (5,5) + (4,)*k + (2,1,1); n = sum(lam)
    dl, dm = dvec(lam), dvec(mu)
    diff = [j for j in range(n+1) if dl[j] != dm[j]]
    print(f"  k={k} n={n} f_equal={dl[n]==dm[n]} first_diff={diff[0]} last_diff={diff[-1]} ham={len(diff)}", flush=True)
