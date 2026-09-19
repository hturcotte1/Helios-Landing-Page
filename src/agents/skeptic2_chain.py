"""skeptic2_chain.py -- proof-by-finite-check of the sign identity S3 for GL_3 and its consequences.

S3(A,B;p,q):  for a = (A+B, B, 0) [any strictly decreasing a in Z^3 has this form up to a common shift], b = (p,q,0),
    sum_{sigma in S_3} sgn(a + sigma b) = sum_{sigma in S_3} sgn(a - sigma b),
where sgn(v) = product over i<j of sgn(v_i - v_j)  (= sign of the permutation sorting v decreasingly; 0 if a tie).
Each factor sgn(v_i - v_j) is the sign of a linear form in (A,B,p,q) with coefficients in {0,+1,-1}; hence both sides are
constant on the open cells of the central hyperplane arrangement H in R^4 spanned by all these forms together with
A, B, p, q, p-q, A+B.  Every cell of H meeting the region {A>0, B>0, p>=q>=0} contains an integer point with all coordinates
<= 16 (extreme rays of a closed cell are kernels of rank-3 {0,+-1}-matrices, so have integer entries of absolute value
<= 4 (maximal 3x3 {0,+-1}-determinant); by Caratheodory a point of the cell is a positive combination of <= 4 linearly
independent extreme rays, and the sum of those rays is again in the cell).  So checking all 1<=A,B<=16, 0<=q<=p<=16 proves S3.

Consequences checked numerically as well:
  * N_j(nu) = N_j(nu*) for GL_3 (number of chains nu -> kappa of length j through partitions with <= 3 rows),
    nu* = (nu_1 - nu_3, nu_1 - nu_2, 0), for all nu with nu_1 <= 8, j <= 20;
  * R_K:  d_j((3^K,2,2,2)) = d_j((3^K,1,1,1)) for j <= K, K <= 22 (own d-vector code);
  * family H_k = {(5,3^k,2,2,2), (4,4,3^k,1,1,1)}: differing indices exactly {k+4,...,n-4}, k <= 12.
"""
import sys, os
from itertools import permutations
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from skeptic2_targeted import dvec

def sgn(v):
    s = 1
    for i in range(len(v)):
        for j in range(i + 1, len(v)):
            d = v[i] - v[j]
            if d == 0:
                return 0
            if d < 0:
                s = -s
    return s

def S3_check(BOX=16):
    bad = 0; cnt = 0
    perms = list(permutations(range(3)))
    for A in range(1, BOX + 1):
        for B in range(1, BOX + 1):
            a = (A + B, B, 0)
            for p in range(0, BOX + 1):
                for q in range(0, p + 1):
                    b = (p, q, 0)
                    plus = sum(sgn((a[0] + b[s[0]], a[1] + b[s[1]], a[2] + b[s[2]])) for s in perms)
                    minus = sum(sgn((a[0] - b[s[0]], a[1] - b[s[1]], a[2] - b[s[2]])) for s in perms)
                    cnt += 1
                    if plus != minus:
                        bad += 1
                        if bad < 5:
                            print("  S3 FAILS at A,B,p,q =", A, B, p, q, plus, minus)
    print("S3 finite check: box %d, %d parameter points, %d failures" % (BOX, cnt, bad), flush=True)
    return bad == 0

def N(nu, j, rows=3):
    """number of chains nu = k0 c k1 c ... c kj, one box at a time, all with <= rows rows."""
    nu = tuple(nu) + (0,) * (rows - len(nu))
    memo = {}
    def rec(k, t):
        if t == 0:
            return 1
        key = (k, t)
        if key in memo:
            return memo[key]
        tot = 0
        for i in range(rows):
            if i == 0 or k[i - 1] > k[i]:
                k2 = k[:i] + (k[i] + 1,) + k[i + 1:]
                tot += rec(k2, t - 1)
        memo[key] = tot
        return tot
    return rec(nu, j)

if __name__ == "__main__":
    ok = S3_check(16)
    print("S3 PROVED by finite check:", ok)
    # chain symmetry
    bad = 0; cnt = 0
    for n1 in range(0, 9):
        for n2 in range(0, n1 + 1):
            for n3 in range(0, n2 + 1):
                nu = (n1, n2, n3); dual = (n1 - n3, n1 - n2, 0)
                for j in range(0, 21):
                    cnt += 1
                    if N(nu, j) != N(dual, j):
                        bad += 1
                        print("  chain symmetry FAILS", nu, dual, j)
    print("chain symmetry N_j(nu) = N_j(nu*): %d cases (nu_1<=8, j<=20), %d failures" % (cnt, bad), flush=True)
    # R_K
    bad = 0
    for K in range(1, 23):
        u = dvec((3,) * K + (2, 2, 2)); v = dvec((3,) * K + (1, 1, 1))
        first = next(j for j in range(len(u)) if u[j] != v[j])
        assert all(u[j] == v[j] for j in range(K + 1))
        # also check the reduction d_j = N_j((3)) resp. N_j((3,3)) for j <= K
        assert all(u[j] == N((3, 0, 0), j) and v[j] == N((3, 3, 0), j) for j in range(K + 1))
        print("R_%d holds; first differing index %d (= K+2: %s); reduction to N_j verified" % (K, first, first == K + 2))
    # family H_k
    for k in range(1, 13):
        lam = (5,) + (3,) * k + (2, 2, 2); mu = (4, 4) + (3,) * k + (1, 1, 1)
        n = 3 * k + 11
        u, v = dvec(lam), dvec(mu)
        D = [j for j in range(n + 1) if u[j] != v[j]]
        print("H_%d n=%d: differing indices == {k+4..n-4}: %s  (first %d last %d ham %d)" % (
            k, n, D == list(range(k + 4, n - 3)), D[0], D[-1], len(D)))
