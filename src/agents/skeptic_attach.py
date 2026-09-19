"""skeptic_attach.py -- the attachment family  kappa|nu  vs  kappa|nu^t.
Notation: kappa|nu = rows of kappa followed by rows of nu (a partition when kappa_last >= nu_1).
Lemma A (tested here): for j <= kappa_last - nu_1,  d_j(kappa|nu) = sum_i C(j,i) d_i(kappa) d_{j-i}(nu).
Corollary: d_j(kappa|nu) = d_j(kappa|nu^t) for j <= kappa_last - max(nu_1, l(nu)).
Test 1: Lemma A and the corollary for all kappa, nu with |kappa| <= 10, |nu| <= 6 (exact d-vectors).
Test 2: f-equality search: f^{kappa|nu} = f^{kappa|nu^t}, via hook products, for kappa = (N), (N,N2), (N,N,N3 ..)
        with nu non-symmetric, |nu| <= 12, and N up to 400.  Records C_1(nu) for each solution.
"""
import os
import sys
from math import comb, prod
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
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

def hookprod(lam):
    lt = transpose(lam)
    return prod(lam[i] - j + lt[j] - i - 1 for i in range(len(lam)) for j in range(lam[i]))

def c1(nu):
    return sum(j - i for i, p in enumerate(nu) for j in range(p))

# Test 1
bad = 0; tested = 0
for K in range(1, 11):
    for kappa in gen_partitions(K):
        for M in range(1, 7):
            for nu in gen_partitions(M):
                if kappa[-1] < nu[0]: continue
                lam = kappa + nu
                dl, dk, dn = dvec(lam), dvec(kappa), dvec(nu)
                J = kappa[-1] - nu[0]
                for j in range(J + 1):
                    pred = sum(comb(j, i) * (dk[i] if i <= K else 0) * (dn[j-i] if j-i <= M else 0) for i in range(j+1))
                    tested += 1
                    if pred != dl[j]: bad += 1; print("LEMMA A FAILS", kappa, nu, j, pred, dl[j])
                nut = transpose(nu)
                if kappa[-1] >= nut[0]:
                    dm = dvec(kappa + nut)
                    J2 = kappa[-1] - max(nu[0], len(nu))
                    for j in range(J2 + 1):
                        if dl[j] != dm[j]: bad += 1; print("COROLLARY FAILS", kappa, nu, j)
                    # sharpness: does it fail at J2+1 when nu != nu^t ?
print(f"Test 1: {tested} instances of Lemma A checked, failures={bad}", flush=True)

# sharpness statistics of the corollary: first differing index vs J2+1
sharp = {}
for K in range(2, 11):
    for kappa in gen_partitions(K):
        for M in range(2, 7):
            for nu in gen_partitions(M):
                nut = transpose(nu)
                if nu >= nut or kappa[-1] < max(nu[0], len(nu)): continue
                dl, dm = dvec(kappa + nu), dvec(kappa + nut)
                J2 = kappa[-1] - max(nu[0], len(nu))
                first = next((j for j in range(K + M + 1) if dl[j] != dm[j]), None)
                sharp[first - (J2 + 1) if first is not None else 'equal'] = sharp.get(first - (J2 + 1) if first is not None else 'equal', 0) + 1
print("Corollary sharpness: histogram of (first differing j) - (J2+1):", sharp, flush=True)

# Test 2: f-equality
sols = []
for M in range(2, 13):
    for nu in gen_partitions(M):
        nut = transpose(nu)
        if nu >= nut: continue
        w = max(nu[0], len(nu))
        for shape in ('N', 'NN2', 'NNN'):
            for N in range(w, 401):
                if shape == 'N': kappas = [(N,)]
                elif shape == 'NN2': kappas = [(N, N2) for N2 in range(w, N + 1)] if N <= 60 else []
                else: kappas = [(N, N, N)] if N <= 400 else []
                for kappa in kappas:
                    if hookprod(kappa + nu) == hookprod(kappa + nut):
                        sols.append((kappa, nu, c1(nu)))
                        print("f-EQUAL:", kappa, "|", nu, " vs |", nut, " C1(nu)=", c1(nu), " n=", sum(kappa)+M, flush=True)
print("total f-equal solutions:", len(sols), flush=True)
