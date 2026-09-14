"""topend (second pass): self-contained exact tools.
 * chi(lam, rho): Murnaghan-Nakayama via beta-numbers (abacus), independent of topend_lib.chi (rim surgery)
   and of src/check_characters.py (different padding logic; cross-checked in topend_algebra.py).
 * contents / content power sums C_k(lam); n, r.
 * z_rho, sigma(rho) (square-root count), even-type test, class sizes, central character omega.
 * u-vector (u_0..u_I) for all partitions of all n <= N by a top-down DP that never touches characters:
      u_i(lam) = sum_{corners c} u_i(lam - c)  (|lam| > i),  u_i(lam) = 1 (|lam| = i).
Exact integer / Fraction arithmetic only."""
import os, sys
from fractions import Fraction
from functools import lru_cache
from math import factorial, comb
from collections import Counter
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from young import partitions, conjugate, boxes, f_hook, f_skew, involutions, removable_corners, remove_box

# ---------------------------------------------------------------- Murnaghan-Nakayama (beta numbers)
def _strip_removals(beta, k):
    """beta: strictly decreasing tuple of beta-numbers (length L, all >= 0).  Yield (newbeta, height) for
    every way to move one bead from x to x-k >= 0 with x-k unoccupied; height = # beads strictly between."""
    S = set(beta)
    for x in beta:
        y = x - k
        if y < 0 or y in S:
            continue
        h = sum(1 for b in beta if y < b < x)
        nb = tuple(sorted((set(beta) - {x}) | {y}, reverse=True))
        yield nb, h

@lru_cache(maxsize=None)
def _chi_beta(beta, rho):
    """rho: tuple of cycle lengths (any order, no zeros).  beta: beta-set with |beta| = sum(rho) + C(L,2)."""
    if not rho:
        return 1
    k, rest = rho[0], rho[1:]
    return sum((-1) ** h * _chi_beta(nb, rest) for nb, h in _strip_removals(beta, k))

def chi(lam, rho):
    """chi^lam(rho), rho a tuple of cycle lengths summing to |lam| (padded with 1's or not: we pad)."""
    n = sum(lam)
    s = sum(rho)
    if s > n:
        raise ValueError((lam, rho))
    big = tuple(sorted([p for p in rho if p > 1], reverse=True))
    rho = big + (1,) * (n - sum(big))
    L = max(len(lam), 1)
    beta = tuple(lam[i] + (L - 1 - i) for i in range(L)) if lam else (0,)
    return _chi_beta(beta, rho)

# ---------------------------------------------------------------- contents
def contents(lam):
    return [j - i for (i, j) in boxes(lam)]

def C(lam, k):
    return sum(c ** k for c in contents(lam))

def r_of(lam):
    return len(removable_corners(lam))

# ---------------------------------------------------------------- classes
def z(rho):
    res = 1
    for k, m in Counter(rho).items():
        res *= k ** m * factorial(m)
    return res

def is_even_type(rho):
    return all(m % 2 == 0 for k, m in Counter(rho).items() if k % 2 == 0)

def dfact(m):
    """(m-1)!! = number of perfect matchings of m points (m even), 1 for m = 0."""
    r = 1
    for t in range(m - 1, 0, -2):
        r *= t
    return r

def sigma(rho):
    """Number of square roots of a permutation of cycle type rho.  A square root pairs up some cycles of equal
    length (a pair of k-cycles is the square of a 2k-cycle, k ways) and takes a square root inside each
    remaining odd cycle (unique); even cycles must all be paired."""
    res = 1
    for k, m in Counter(rho).items():
        if k % 2 == 0:
            if m % 2:
                return 0
            res *= dfact(m) * k ** (m // 2)
        else:
            res *= sum(comb(m, 2 * j) * dfact(2 * j) * k ** j for j in range(m // 2 + 1))
    return res

def strip_ones(rho):
    return tuple(p for p in rho if p > 1)

def class_size(rhobar, n):
    s = sum(rhobar)
    if n < s:
        return 0
    return factorial(n) // (z(rhobar) * factorial(n - s))

def omega(lam, rhobar):
    """Central character omega_lam(K_rhobar) = |K| chi^lam(rhobar u 1^{n-|rhobar|}) / f^lam  (0 if n < |rhobar|)."""
    n = sum(lam)
    if n < sum(rhobar):
        return Fraction(0)
    return Fraction(class_size(rhobar, n) * chi(lam, rhobar), f_hook(lam))

def even_types_no_ones(max_size):
    return [rho for s in range(max_size + 1) for rho in partitions(s) if all(p > 1 for p in rho) and is_even_type(rho)]

# ---------------------------------------------------------------- u-vectors by DP on Young's lattice
def u_table(N, I):
    """dict lam -> tuple (u_0(lam), ..., u_I(lam)) for all partitions lam of all n <= N, with
    u_i(lam) = d_{n-i}(lam) = sum_{nu |- i} f^{lam/nu} computed by u_i(lam) = sum_c u_i(lam - c) for n > i,
    u_i(lam) = 1 for n = i, u_i(lam) = 0 for n < i."""
    U = {}
    for n in range(N + 1):
        for lam in partitions(n):
            if n == 0:
                U[lam] = tuple([1] + [0] * I); continue
            subs = [U[remove_box(lam, c)] for c in removable_corners(lam)]
            U[lam] = tuple(1 if i == n else (0 if i > n else sum(s[i] for s in subs)) for i in range(I + 1))
    return U

if __name__ == "__main__":
    # tiny self-tests
    assert chi((2, 1), (3,)) == -1 and chi((2, 1), (2, 1)) == 0 and chi((2, 1), (1, 1, 1)) == 2
    assert chi((2, 2), (2, 2)) == 2 and chi((3, 1), (4,)) == -1 and chi((2, 2), (4,)) == 0 and chi((3, 1), (2, 2)) == -1
    assert chi((5,), (2, 2)) == 1 and chi((1, 1, 1, 1, 1), (2, 2)) == 1 and chi((3, 2), (2, 2, 1)) == 1
    assert chi((2, 2, 1), (5,)) == 0 and chi((3, 1, 1), (5,)) == 1 and chi((4, 1), (5,)) == -1
    assert sigma((2, 2)) == 2 and sigma((3,)) == 1 and sigma((1, 1, 1)) == 4 and sigma((2, 2, 2, 2)) == 12
    print("topend_mn self-tests OK")
