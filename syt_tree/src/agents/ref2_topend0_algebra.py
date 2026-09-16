"""Referee #2, independent brute-force check of Proposition 2.3 (I1)-(I6) in Z[S_n], n = 1..NMAX.
Permutations as tuples p with p[x] = image of x; (p*q)(x) = p(q(x)) (q applied first), as in the report.
Every identity is checked as a full equality of dictionaries perm -> coefficient, so any class the author
might have missed would show up as a mismatch."""
import sys, itertools
from math import comb
from collections import Counter

def compose(p, q): return tuple(p[x] for x in q)
def transp(i, k, n):
    p = list(range(n)); p[i], p[k] = p[k], p[i]; return tuple(p)
def mult(A, B):
    R = Counter()
    for p, a in A.items():
        for q, b in B.items():
            R[compose(p, q)] += a * b
    return {k: v for k, v in R.items() if v}
def add(*terms):
    R = Counter()
    for c, A in terms:
        for q, b in A.items(): R[q] += c * b
    return {k: v for k, v in R.items() if v}
def ident(n): return {tuple(range(n)): 1}
def J(k, n):  # J_k = sum_{i<k} (i k), 1-indexed k; internally 0-indexed positions i-1, k-1
    return {transp(i, k - 1, n): 1 for i in range(k - 1)}
def power(A, m, n):
    R = ident(n)
    for _ in range(m): R = mult(R, A)
    return R
def cycle_type(p):
    n = len(p); seen = [False] * n; ct = []
    for x in range(n):
        if not seen[x]:
            l = 0; y = x
            while not seen[y]: seen[y] = True; y = p[y]; l += 1
            if l > 1: ct.append(l)
    return tuple(sorted(ct, reverse=True))
def K(rho, n):
    return {p: 1 for p in itertools.permutations(range(n)) if cycle_type(p) == tuple(rho)}
def p_m(m, n):
    return add(*[(1, power(J(k, n), m, n)) for k in range(1, n + 1)]) if n >= 1 else {}

def check(n):
    K2, K3, K22, K4, K5, K33 = (K(r, n) for r in [(2,), (3,), (2, 2), (4,), (5,), (3, 3)])
    one = ident(n)
    P1, P2, P3, P4 = (p_m(m, n) for m in (1, 2, 3, 4))
    ok = {}
    ok['I1'] = P1 == K2
    ok['I2'] = P2 == add((comb(n, 2), one), (1, K3))
    ok['I3'] = mult(K2, K2) == add((comb(n, 2), one), (3, K3), (2, K22))
    ok['I4'] = P3 == add((1, K4), (2 * n - 3, K2))
    ok['I5'] = P4 == add((1, K5), (3 * n - 4, K3), (4, K22), (n * (n - 1) * (4 * n - 5) // 6, one))
    ok['I6'] = mult(K3, K3) == add((2 * comb(n, 3), one), (3 * n - 8, K3), (8, K22), (5, K5), (2, K33))
    # also: e_2(J) = K3 + K22
    e2 = add((1, mult(P1, P1)), (-1, P2)); e2 = {k: v for k, v in ((k, v // 2) for k, v in e2.items())}
    ok['e2'] = e2 == add((1, K3), (1, K22))
    # what classes actually occur in K3^2 (report on the class list)
    types = sorted(set(cycle_type(p) for p in mult(K3, K3)))
    return ok, types

if __name__ == '__main__':
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    allok = True
    for n in range(1, NMAX + 1):
        ok, types = check(n)
        allok &= all(ok.values())
        print(n, ok, 'classes in K3^2:', types, flush=True)
    print('ALL OK' if allok else 'FAILURE')
