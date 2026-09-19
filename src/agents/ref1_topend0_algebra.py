"""Referee #1, independent check of Proposition 2.3 (I1)-(I6) in Z[S_n], n <= NMAX, by brute-force
permutation multiplication.  Permutations are tuples p with p[x] = image of x; (p*q)(x) = p(q(x)).
Also: check which cycle types occur in K_3^2 (the report claims (4,2) is odd -- it is even)."""
import sys
from itertools import permutations
from math import comb
from collections import Counter

def compose(p, q): return tuple(p[x] for x in q)
def cycle_type(p):
    n = len(p); seen = [False]*n; ct = []
    for s in range(n):
        if not seen[s]:
            l = 0; x = s
            while not seen[x]:
                seen[x] = True; x = p[x]; l += 1
            ct.append(l)
    return tuple(sorted(ct, reverse=True))
def reduced(ct): return tuple(c for c in ct if c > 1)
def mult(A, B):
    R = Counter()
    for p, a in A.items():
        for q, b in B.items():
            R[compose(p, q)] += a*b
    return R
def add(A, B, c=1):
    R = Counter(A)
    for q, b in B.items(): R[q] += c*b
    return R
def transp(i, k, n):
    p = list(range(n)); p[i], p[k] = p[k], p[i]; return tuple(p)
def J(k, n):  # J_k = sum_{i<k} (i k), 1-indexed k; 0-indexed points i=0..k-2, k-1
    return Counter({transp(i, k-1, n): 1 for i in range(k-1)})
def ident(n): return Counter({tuple(range(n)): 1})
def power(A, m, n):
    R = ident(n)
    for _ in range(m): R = mult(R, A)
    return R
def classsum(rho, n):
    K = Counter()
    for p in permutations(range(n)):
        if reduced(cycle_type(p)) == rho: K[p] = 1
    return K
def clean(A): return Counter({k: v for k, v in A.items() if v})

def check(n):
    e = ident(n)
    Ks = {rho: classsum(rho, n) for rho in [(), (2,), (3,), (2,2), (4,), (5,), (3,3)]}
    K = lambda rho: Ks[rho]
    p = {m: Counter() for m in range(1, 5)}
    for k in range(1, n+1):
        Jk = J(k, n)
        for m in range(1, 5):
            p[m] = add(p[m], power(Jk, m, n))
    lhs_rhs = {
      'I1': (p[1], K((2,))),
      'I2': (p[2], add(Counter({tuple(range(n)): comb(n,2)}), K((3,)))),
      'I3': (mult(K((2,)), K((2,))), add(add(Counter({tuple(range(n)): comb(n,2)}), K((3,)), 3), K((2,2)), 2)),
      'I4': (p[3], add(K((4,)), K((2,)), 2*n-3)),
      'I5': (p[4], add(add(add(K((5,)), K((3,)), 3*n-4), K((2,2)), 4), Counter({tuple(range(n)): n*(n-1)*(4*n-5)//6}))),
      'I6': (mult(K((3,)), K((3,))), add(add(add(add(Counter({tuple(range(n)): 2*comb(n,3)}), K((3,)), 3*n-8), K((2,2)), 8), K((5,)), 5), K((3,3)), 2)),
    }
    ok = True
    for name, (L, R) in lhs_rhs.items():
        L, R = clean(L), clean(R)
        if L != R:
            ok = False
            diff = clean(add(L, R, -1))
            print(f"  n={n} {name} FAILS; differing terms by type:", Counter(reduced(cycle_type(g)) for g in diff))
    # cycle types occurring in K_3^2, with their (constant) coefficients
    K33 = clean(mult(K((3,)), K((3,))))
    types = {}
    for g, c in K33.items():
        types.setdefault(reduced(cycle_type(g)), set()).add(c)
    print(f"  n={n}: (I1)-(I6) {'OK' if ok else 'FAIL'};  K_3^2 support types/coeffs: {types}")
    return ok

if __name__ == '__main__':
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    allok = True
    for n in range(1, NMAX+1):
        allok &= check(n)
    # sign of (4,2): the report says '(4,2) is odd'
    p = (1,2,3,0,5,4)  # 4-cycle * 2-cycle on 6 points
    def sign(p):
        return (-1)**(len(p) - len(cycle_type(p)))
    print("sign of a (4,2)-permutation:", sign(p), " (so (4,2) is", "EVEN" if sign(p)==1 else "ODD", ")")
    print("ALL OK" if allok else "SOME FAIL")
