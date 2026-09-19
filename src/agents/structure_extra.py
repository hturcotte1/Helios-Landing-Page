"""structure angle: extra checks.
(8) derangement transform: delta_m(lam) := sum_{tau in Der_m} chi^lam(tau^2 cup 1^{n-m}) satisfies
    (n-j)! d_j(lam) = sum_m C(n-j, m) delta_m(lam), i.e. i! u_i = sum_m C(i,m) delta_m; checked n <= 8
    with delta computed from characters (MN) by brute force over Der_m.
(9) Psi_lam(-1) = (-1)^n, Psi_lam(1) = 1 for n <= 10 (Psi from the d-vector).
(10) values n! P_lam(-1) for hooks and two-row shapes, n <= 12 (table only).
"""
import os, sys, itertools
from fractions import Fraction
from math import factorial, comb
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from young import partitions, involutions
from census import d_vector
from check_characters import chi
from structure_verify import cycle_type, w_from_d

def check_delta(N):
    bad = 0
    for n in range(1, N+1):
        for lam in partitions(n):
            lam = tuple(lam); d = d_vector(lam)
            delta = []
            for m in range(n+1):
                tot = 0
                for p in itertools.permutations(range(m)):
                    if any(p[i] == i for i in range(m)): continue
                    sq = tuple(p[p[i]] for i in range(m))
                    ct = tuple(sorted(list(cycle_type(sq)) + [1]*(n-m), reverse=True)) if m > 0 else tuple([1]*n)
                    tot += chi(lam, ct)
                delta.append(tot)
            for j in range(n+1):
                i = n-j
                if factorial(i)*d[j] != sum(comb(i, m)*delta[m] for m in range(i+1)):
                    bad += 1; print("DELTA MISMATCH", lam, j)
    print(f"(8) derangement transform i! u_i = sum_m C(i,m) delta_m, all n <= {N}: mismatches = {bad}")

def check_psi_special(N):
    bad = 0
    for n in range(1, N+1):
        for lam in partitions(n):
            w = w_from_d(tuple(lam))
            v1 = sum(Fraction(involutions(m))*w[m] for m in range(n+1))
            vm1 = sum(Fraction((-1)**m*involutions(m))*w[m] for m in range(n+1))
            if v1 != 1 or vm1 != (-1)**n:
                bad += 1; print("PSI special value", lam, v1, vm1)
    print(f"(9) Psi_lam(1) = 1 and Psi_lam(-1) = (-1)^n for all n <= {N}: violations = {bad}")

def table(N):
    print("(10) n! P_lam(-1) for hooks (n-k,1^k) and two-row (n-k,k):")
    for n in range(2, N+1):
        d = lambda lam: sum((-1)**j*d_vector(lam)[j]*factorial(n)//factorial(j) for j in range(n+1))
        hooks = [d(tuple([n-k]+[1]*k)) for k in range(0, n)]
        rows = [d(tuple(x for x in (n-k, k) if x)) for k in range(0, n//2+1)]
        print(n, "hooks:", hooks, " two-row:", rows)

if __name__ == "__main__":
    check_delta(7)
    check_psi_special(10)
    table(12)
