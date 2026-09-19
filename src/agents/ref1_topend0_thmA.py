"""Referee #1: independent Murnaghan-Nakayama (rim-hook removal on partitions, my own code) and
check of Theorem A (A1)-(A4) for all lam |- n <= NMAX, including n < |rho| boundary cases;
symbolic check of the (A4) expansion; cross-check against the d-vector (u_i) via census.d_vector."""
import os
import sys, os
from fractions import Fraction as Fr
from functools import lru_cache
from math import comb
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from young import partitions, f_hook
from census import d_vector

def rim_hooks(lam, h):
    """All ways to remove a border strip of size h from lam: yields (new partition, height)."""
    lam = list(lam); L = len(lam)
    # beta-numbers: b_i = lam_i + (L-1-i)
    beta = [lam[i] + (L-1-i) for i in range(L)]
    bset = set(beta)
    for i, b in enumerate(beta):
        nb = b - h
        if nb >= 0 and nb not in bset:
            height = sum(1 for c in beta if nb < c < b)  # number of beta's jumped over
            newbeta = sorted([c for c in beta if c != b] + [nb], reverse=True)
            mu = tuple(newbeta[j] - (L-1-j) for j in range(L))
            mu = tuple(x for x in mu if x > 0)
            yield mu, height

@lru_cache(maxsize=None)
def chi(lam, rho):
    """chi^lam(rho), rho a partition of |lam| (tuple, weakly decreasing)."""
    if sum(lam) == 0: return 1
    h = rho[0]; rest = rho[1:]
    return sum((-1)**ht * chi(mu, rest) for mu, ht in rim_hooks(lam, h))

def contents(lam):
    return [j - i for i, r in enumerate(lam) for j in range(r)]
def Cm(lam, m): return sum(c**m for c in contents(lam))
def pad(rho, n): return tuple(sorted(list(rho) + [1]*(n - sum(rho)), reverse=True))
def falling(n, k):
    r = 1
    for i in range(k): r *= (n - i)
    return r
def z(rho):
    from collections import Counter
    from math import factorial
    r = 1
    for k, m in Counter(rho).items(): r *= k**m * factorial(m)
    return r
def omega(lam, rho):
    n = sum(lam)
    if n < sum(rho): return Fr(0)
    size = Fr(falling(n, sum(rho)), z(rho))
    return size * chi(lam, pad(rho, n)) / f_hook(lam)

def thmA(lam):
    n = sum(lam); C1, C2, C4 = Cm(lam,1), Cm(lam,2), Cm(lam,4)
    A1 = Fr(C2 - comb(n,2))
    A2 = Fr(C1**2, 2) - Fr(3*C2, 2) + comb(n,2)
    A3 = Fr(C4 - (3*n-10)*C2 - 2*C1**2) + Fr(n*(n-1)*(5*n-19), 6)
    A4 = Fr(1,2)*(A1**2 - 2*comb(n,3) - (3*n-8)*A1 - 8*A2 - 5*A3)
    A4b = (Fr(C2**2,2) - Fr(5,2)*C4 + 3*C1**2 - Fr(n**2*C2,2) + Fr(13*n*C2,2) - 15*C2
           + Fr(n**4,8) - Fr(7*n**3,4) + Fr(47*n**2,8) - Fr(17*n,4))
    return A1, A2, A3, A4, A4b

def main(NMAX):
    # sanity of my MN: orthogonality of columns and f^lam = chi(1^n) for n <= 8
    from young import partitions as P
    for n in range(0, 9):
        parts = list(P(n))
        for lam in parts:
            assert chi(tuple(lam), tuple([1]*n)) == f_hook(lam), (lam,)
        for rho in parts:
            for tau in parts:
                s = sum(chi(tuple(l), tuple(rho))*chi(tuple(l), tuple(tau)) for l in parts)
                assert s == (z(rho) if rho == tau else 0)
    print("MN sanity (dimension, column orthogonality) n<=8 OK")
    bad = 0; cnt = 0; worst = []
    for n in range(0, NMAX+1):
        for lam in partitions(n):
            lam = tuple(lam)
            A1, A2, A3, A4, A4b = thmA(lam)
            o = [omega(lam, r) for r in [(3,), (2,2), (5,), (3,3)]]
            for name, x, y in zip(['A1','A2','A3','A4'], o, [A1,A2,A3,A4]):
                cnt += 1
                if x != y:
                    bad += 1
                    if len(worst) < 10: worst.append((n, lam, name, x, y))
            cnt += 1
            if A4 != A4b:
                bad += 1
                if len(worst) < 10: worst.append((n, lam, 'A4expansion', A4, A4b))
        print(f"n={n} done, cumulative checks={cnt}, bad={bad}", flush=True)
    print("Theorem A vs independent MN: checks", cnt, "bad", bad, worst)
    # d-vector cross-check: chi(3)=3u3-2f, chi(22)=4u4-4u3+f, chi(5)=5u5-5u4+f, chi(33)=(9/2)u6-(9/2)u5+(3/2)u3-f/2,
    # then (A1)-(A4) via omega = |K| chi / f, for n <= 18 (no characters used)
    bad2 = 0; cnt2 = 0
    for n in range(0, 19):
        for lam in partitions(n):
            lam = tuple(lam); d = d_vector(lam); f = d[n]
            u = lambda i: d[n-i] if i <= n else 0
            c3 = 3*u(3) - 2*f; c22 = 4*u(4) - 4*u(3) + f; c5 = 5*u(5) - 5*u(4) + f
            c33 = Fr(9,2)*u(6) - Fr(9,2)*u(5) + Fr(3,2)*u(3) - Fr(f,2)
            A1, A2, A3, A4, _ = thmA(lam)
            om = [Fr(2*comb(n,3)*c3, f), Fr(3*comb(n,4)*c22, f), Fr(24*comb(n,5)*c5, f), Fr(falling(n,6),18)*c33/f]
            for x, y in zip(om, [A1,A2,A3,A4]):
                cnt2 += 1
                if x != y: bad2 += 1
    print("Theorem A vs d-vector (Theorem B inversions), n<=18: checks", cnt2, "bad", bad2)

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 16)
