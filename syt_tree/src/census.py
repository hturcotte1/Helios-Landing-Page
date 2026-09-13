"""Up-census s_k(lambda), down-census d_j(lambda) and the generating-function
identity  E_lambda(t) = exp(t + t^2/2) P_lambda(t).

All arithmetic is exact.
"""
from __future__ import annotations

from functools import lru_cache
from math import comb, factorial

from young import (Partition, add_box, addable_corners, contains, down_set,
                   f_skew, involutions, partitions, size, up_set)


# ---------------------------------------------------------------------------
# Down-census  d_j(lambda) = # ways to remove j boxes one at a time
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def d_vector(lam: Partition) -> tuple[int, ...]:
    """``(d_0(lam), ..., d_n(lam))`` with ``d_j = sum_c d_{j-1}(lam - c)``."""
    n = size(lam)
    if n == 0:
        return (1,)
    d = [0] * (n + 1)
    d[0] = 1
    for mu in down_set(lam):
        dm = d_vector(mu)
        for j in range(1, n + 1):
            d[j] += dm[j - 1]
    return tuple(d)


def d_vector_by_skew(lam: Partition) -> tuple[int, ...]:
    """Independent computation: d_j = sum_{nu subset lam, |lam/nu| = j} f^{lam/nu}."""
    n = size(lam)
    d = [0] * (n + 1)
    for m in range(n + 1):
        for nu in partitions(m):
            if contains(lam, nu):
                d[n - m] += f_skew(lam, nu)
    return tuple(d)


# ---------------------------------------------------------------------------
# Up-census  s_k(lambda) = # saturated chains of length k going up from lambda
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def s_up(lam: Partition, k: int) -> int:
    """Brute force by the definition: s_k(lam) = sum_{c addable} s_{k-1}(lam + c)."""
    if k == 0:
        return 1
    return sum(s_up(mu, k - 1) for mu in up_set(lam))


def s_up_by_skew(lam: Partition, k: int) -> int:
    """Independent computation: s_k = sum_{mu contains lam, |mu| = |lam| + k} f^{mu/lam}."""
    n = size(lam)
    return sum(f_skew(mu, lam) for mu in partitions(n + k) if contains(mu, lam))


def s_vector(lam: Partition, K: int) -> tuple[int, ...]:
    return tuple(s_up(lam, k) for k in range(K + 1))


# ---------------------------------------------------------------------------
# The identity  E_lambda(t) = exp(t + t^2/2) P_lambda(t), coefficientwise:
#     s_k(lam) = sum_j C(k, j) I_{k-j} d_j(lam)
# ---------------------------------------------------------------------------

def s_from_d(lam: Partition, k: int) -> int:
    d = d_vector(lam)
    return sum(comb(k, j) * involutions(k - j) * d[j] for j in range(0, min(k, len(d) - 1) + 1))


def check_identity(max_n: int, max_k: int, verbose: bool = False) -> bool:
    """Check s_k(lam) (brute force) == sum_j C(k,j) I_{k-j} d_j(lam) for all
    |lam| <= max_n and k <= max_k."""
    ok = True
    for n in range(max_n + 1):
        for lam in partitions(n):
            for k in range(max_k + 1):
                lhs = s_up(lam, k)
                rhs = s_from_d(lam, k)
                if lhs != rhs:
                    ok = False
                    print("IDENTITY FAILS", lam, k, lhs, rhs)
        if verbose:
            print(f"n={n}: identity holds for all partitions and k<={max_k}")
    return ok


# ---------------------------------------------------------------------------
# Polynomials as coefficient lists (exact rationals avoided: work with
# "exponential" coefficient lists, i.e. lists of d_j, s_k)
# ---------------------------------------------------------------------------

def d_polynomial_sympy(lam: Partition):
    """P_lambda(t) as a sympy polynomial."""
    import sympy as sp
    t = sp.symbols('t')
    d = d_vector(lam)
    return sum(sp.Rational(dj, factorial(j)) * t ** j for j, dj in enumerate(d)), t
