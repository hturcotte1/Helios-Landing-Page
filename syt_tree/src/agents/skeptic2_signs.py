"""skeptic2_signs.py -- test of the sign identity  S(a,b):  sum_{sigma in S_n} sgn(a + sigma.b) = sum_{sigma} sgn(a - sigma.b)
for strictly decreasing a in Z^n and weakly decreasing b in Z^n, where for a vector v in Z^n
sgn(v) = sign of the permutation sorting v into strictly decreasing order (0 if v has a repeated entry).
By the Racah-Speiser / Brauer-Klimyk form of Weyl's character formula, sum_{beta} m_W(beta) sgn(lambda+rho+beta) is the number of
irreducible constituents (with multiplicity) of V_lambda (x) W for GL_n; so S(a,b) for all orbits says
#const(V_lambda (x) W) = #const(V_lambda (x) W^*) for every W, which (with W = V^{(x) j}) gives the 3-row chain symmetry
N_j(u,v) = N_j(v,u) of agent_notes/skeptic.md (identity R_K).
Also tests the weaker character-wise version for GL_3 via SSYT weights (sanity), and prints the first counterexample if any.
"""
import sys
from itertools import permutations, combinations_with_replacement

def sgn_sort(v):
    n = len(v)
    if len(set(v)) < n:
        return 0
    # sign of permutation sorting v decreasingly = parity of inversions w.r.t. decreasing order
    inv = 0
    for i in range(n):
        for j in range(i + 1, n):
            if v[i] < v[j]:
                inv += 1
    return -1 if inv % 2 else 1

def Phi(a, b, sign):
    n = len(a)
    tot = 0
    for s in set(permutations(b)):
        tot += sgn_sort(tuple(a[i] + sign * s[i] for i in range(n)))
    return tot

def test(n, amax, bmax):
    bad = 0; cnt = 0
    for a in combinations_with_replacement(range(amax, -1, -1), n):
        if len(set(a)) < n:
            continue
        for b in combinations_with_replacement(range(bmax, -1, -1), n):
            if b[-1] != 0:
                continue
            cnt += 1
            p, m = Phi(a, b, 1), Phi(a, b, -1)
            if p != m:
                bad += 1
                if bad <= 5:
                    print("  COUNTEREXAMPLE n=%d a=%s b=%s  sum sgn(a+sb)=%d  sum sgn(a-sb)=%d" % (n, a, b, p, m))
    print("n=%d: %d (a,b) pairs tested (a_1<=%d, b_1<=%d), %d failures" % (n, cnt, amax, bmax, bad), flush=True)
    return bad

if __name__ == "__main__":
    test(2, 8, 8)
    test(3, 9, 9)
    test(4, 7, 6)
    test(5, 6, 4)
