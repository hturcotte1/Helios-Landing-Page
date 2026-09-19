"""Shifted Young lattice, shifted standard Young tableaux and their up/down
censuses (Worley Problem 2, shifted case).

Conventions
-----------
* A strict partition is a tuple of *strictly* decreasing positive integers;
  the empty partition is ``()``.  Row ``i`` (0-indexed) of the shifted diagram
  occupies the boxes ``(i, j)`` for ``j = i, i+1, ..., i + lam[i] - 1``
  (0-indexed row/column).  A box ``(i, j)`` is *diagonal* iff ``j == i``.
* Cover relations of the shifted Young lattice (1-indexed statement, as in
  the task; here translated to 0-indexed rows):
    - a box can be added to row ``i < l`` iff ``i == 0`` or
      ``lam[i-1] >= lam[i] + 2``;  a new row ``l`` (a single diagonal box)
      can be added iff ``lam`` is empty or ``lam[l-1] >= 2``;
    - the last box of row ``i`` can be removed iff ``i == l-1`` or
      ``lam[i] >= lam[i+1] + 2``.
  The added box of row ``i < l`` sits at column ``i + lam[i] >= i + 1`` and is
  never diagonal; the new-row box ``(l, l)`` is always diagonal.  A removed
  box ``(i, i + lam[i] - 1)`` is diagonal iff ``lam[i] == 1`` (only possible
  for the last row).
* A shifted SYT of shape ``lam`` is a saturated chain from ``()`` to ``lam``;
  ``g^lam`` is their number.
* Up-census ``s_k(lam)`` = number of saturated chains of length ``k`` going
  up from ``lam``; down-census ``d_j(lam)`` = number of chains of length ``j``
  going down; weighted down-census ``d'_j(lam)`` = same with each step that
  removes an off-diagonal box weighted 2 and a diagonal box weighted 1.
* All arithmetic is exact (Python integers / ``fractions.Fraction``).
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import factorial
from typing import Iterator, Tuple

StrictPartition = Tuple[int, ...]
Box = Tuple[int, int]


# ---------------------------------------------------------------------------
# Strict partitions
# ---------------------------------------------------------------------------

def strict_partitions(n: int, max_part: int | None = None) -> Iterator[StrictPartition]:
    """Yield all strict partitions of ``n`` (parts <= ``max_part``) in reverse
    lexicographic order, ``(n,)`` first."""
    if max_part is None or max_part > n:
        max_part = n
    if n == 0:
        yield ()
        return
    for first in range(max_part, 0, -1):
        for rest in strict_partitions(n - first, first - 1):
            yield (first,) + rest


def is_strict_partition(lam) -> bool:
    return all(isinstance(x, int) and x > 0 for x in lam) and all(
        lam[i] > lam[i + 1] for i in range(len(lam) - 1)
    )


def size(lam: StrictPartition) -> int:
    return sum(lam)


def boxes(lam: StrictPartition) -> list[Box]:
    """Boxes of the shifted diagram (0-indexed)."""
    return [(i, j) for i, row in enumerate(lam) for j in range(i, i + row)]


# ---------------------------------------------------------------------------
# Covers (with box coordinates and diagonal flag)
# ---------------------------------------------------------------------------

def addable_boxes(lam: StrictPartition) -> list[tuple[int, Box, bool]]:
    """List of ``(row, (r, c), diagonal)`` for every box that can be added to
    ``lam`` leaving a strict partition.  ``row == len(lam)`` means a new row."""
    ell = len(lam)
    res = []
    for i in range(ell):
        if i == 0 or lam[i - 1] >= lam[i] + 2:
            res.append((i, (i, i + lam[i]), False))
    if ell == 0 or lam[ell - 1] >= 2:
        res.append((ell, (ell, ell), True))
    return res


def removable_boxes(lam: StrictPartition) -> list[tuple[int, Box, bool]]:
    """List of ``(row, (r, c), diagonal)`` for every box that can be removed
    from ``lam`` leaving a strict partition."""
    ell = len(lam)
    res = []
    for i in range(ell):
        if i == ell - 1 or lam[i] >= lam[i + 1] + 2:
            c = i + lam[i] - 1
            res.append((i, (i, c), c == i))
    return res


def add_box(lam: StrictPartition, i: int) -> StrictPartition:
    """Add a box to row ``i`` (``i == len(lam)`` starts a new row)."""
    assert i in [r for r, _, _ in addable_boxes(lam)], (lam, i)
    new = list(lam)
    if i == len(new):
        new.append(1)
    else:
        new[i] += 1
    return tuple(new)


def remove_box(lam: StrictPartition, i: int) -> StrictPartition:
    """Remove the last box of row ``i``."""
    assert i in [r for r, _, _ in removable_boxes(lam)], (lam, i)
    new = list(lam)
    new[i] -= 1
    if new[i] == 0:
        new.pop()
    return tuple(new)


def up_set(lam: StrictPartition) -> list[StrictPartition]:
    """All strict partitions covering ``lam``."""
    return [add_box(lam, i) for i, _, _ in addable_boxes(lam)]


def down_set(lam: StrictPartition) -> list[StrictPartition]:
    """All strict partitions covered by ``lam``."""
    return [remove_box(lam, i) for i, _, _ in removable_boxes(lam)]


def down_set_weighted(lam: StrictPartition) -> list[tuple[StrictPartition, int]]:
    """``(mu, w)`` for every ``mu`` covered by ``lam``; ``w = 1`` if the removed
    box is diagonal, ``w = 2`` otherwise (Fomin's weighted D')."""
    return [(remove_box(lam, i), 1 if diag else 2) for i, _, diag in removable_boxes(lam)]


# ---------------------------------------------------------------------------
# g^lambda: recursion and shifted hook length formula
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def g_rec(lam: StrictPartition) -> int:
    """Number of shifted SYT (saturated chains from ``()``) by the recursion
    ``g^lam = sum_{mu covered by lam} g^mu``."""
    if not lam:
        return 1
    return sum(g_rec(mu) for mu in down_set(lam))


def g_hook(lam: StrictPartition) -> int:
    """Shifted hook length formula
    ``g^lam = n!/(prod_i lam_i!) * prod_{i<j} (lam_i - lam_j)/(lam_i + lam_j)``."""
    n = size(lam)
    num = factorial(n)
    den = 1
    for x in lam:
        den *= factorial(x)
    for a in range(len(lam)):
        for b in range(a + 1, len(lam)):
            num *= lam[a] - lam[b]
            den *= lam[a] + lam[b]
    assert num % den == 0, lam
    return num // den


def shifted_chains(lam: StrictPartition) -> list[tuple[StrictPartition, ...]]:
    """All saturated chains ``() = nu_0 < nu_1 < ... < nu_n = lam`` (explicit
    enumeration, for cross-checks on small shapes)."""
    if not lam:
        return [((),)]
    res = []
    for mu in down_set(lam):
        for ch in shifted_chains(mu):
            res.append(ch + (lam,))
    return res


def chain_to_tableau(chain) -> dict[Box, int]:
    """Filling box -> entry (1..n) of the shifted tableau of a chain."""
    T: dict[Box, int] = {}
    for k in range(1, len(chain)):
        new = set(boxes(chain[k])) - set(boxes(chain[k - 1]))
        assert len(new) == 1, (chain[k - 1], chain[k])
        T[new.pop()] = k
    return T


def is_shifted_syt(T: dict[Box, int], lam: StrictPartition) -> bool:
    """Entries ``1..n`` on the shifted diagram of ``lam``, increasing along
    rows (left to right) and columns (top to bottom)."""
    if set(T) != set(boxes(lam)) or sorted(T.values()) != list(range(1, len(T) + 1)):
        return False
    for (i, j), v in T.items():
        if (i, j + 1) in T and T[(i, j + 1)] < v:
            return False
        if (i + 1, j) in T and T[(i + 1, j)] < v:
            return False
    return True


@lru_cache(maxsize=None)
def total_shifted_syt(k: int) -> int:
    """Number of shifted SYT with ``k`` boxes = sum over strict ``lam |- k`` of ``g^lam``
    = ``s_k(())``."""
    return sum(g_hook(lam) for lam in strict_partitions(k))


# ---------------------------------------------------------------------------
# Up-census (brute force) and down-censuses
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def s_up(lam: StrictPartition, k: int) -> int:
    """``s_k(lam)`` = number of saturated chains of length ``k`` going up from
    ``lam``, by the recursion ``s_k(lam) = sum_{mu covers lam} s_{k-1}(mu)``."""
    if k == 0:
        return 1
    return sum(s_up(mu, k - 1) for mu in up_set(lam))


def s_vector(lam: StrictPartition, K: int) -> tuple[int, ...]:
    return tuple(s_up(lam, k) for k in range(K + 1))


def s_up_by_g(lam: StrictPartition, k: int) -> int:
    """Independent computation: ``s_k(lam) = sum_{mu strict, |mu| = |lam|+k} g^{mu/lam}``
    where ``g^{mu/lam}`` counts chains from ``lam`` to ``mu`` (skew shifted SYT)."""
    n = size(lam)
    return sum(g_skew(mu, lam) for mu in strict_partitions(n + k))


def contains(lam: StrictPartition, nu: StrictPartition) -> bool:
    """Shifted diagram of ``nu`` inside that of ``lam`` (iff ``nu_i <= lam_i``)."""
    return len(nu) <= len(lam) and all(nu[i] <= lam[i] for i in range(len(nu)))


@lru_cache(maxsize=None)
def g_skew(lam: StrictPartition, nu: StrictPartition) -> int:
    """Number of saturated chains from ``nu`` up to ``lam`` (0 if not contained)."""
    if not contains(lam, nu):
        return 0
    if lam == nu:
        return 1
    return sum(g_skew(mu, nu) for mu in down_set(lam))


@lru_cache(maxsize=None)
def d_vector(lam: StrictPartition) -> tuple[int, ...]:
    """``(d_0(lam), ..., d_n(lam))``, ``d_j = sum_{mu covered} d_{j-1}(mu)``."""
    n = size(lam)
    d = [0] * (n + 1)
    d[0] = 1
    for mu in down_set(lam):
        dm = d_vector(mu)
        for j in range(1, n + 1):
            d[j] += dm[j - 1]
    return tuple(d)


@lru_cache(maxsize=None)
def dw_vector(lam: StrictPartition) -> tuple[int, ...]:
    """Weighted down-census ``(d'_0, ..., d'_n)``: removing an off-diagonal
    box weighs 2, a diagonal box weighs 1."""
    n = size(lam)
    d = [0] * (n + 1)
    d[0] = 1
    for mu, w in down_set_weighted(lam):
        dm = dw_vector(mu)
        for j in range(1, n + 1):
            d[j] += w * dm[j - 1]
    return tuple(d)


# ---------------------------------------------------------------------------
# Exponential generating functions as lists of exact rationals
# ---------------------------------------------------------------------------

def egf(coeffs, K: int) -> list[Fraction]:
    """``[c_0/0!, c_1/1!, ..., c_K/K!]`` (missing coefficients are 0)."""
    return [Fraction(coeffs[k], factorial(k)) if k < len(coeffs) else Fraction(0)
            for k in range(K + 1)]


def E_series(lam: StrictPartition, K: int) -> list[Fraction]:
    """``E_lam(t) = sum_k s_k(lam) t^k/k!`` truncated at ``t^K``."""
    return egf(s_vector(lam, K), K)


def P_series(lam: StrictPartition, K: int) -> list[Fraction]:
    """``P_lam(t) = sum_j d_j(lam) t^j/j!`` truncated at ``t^K``."""
    return egf(d_vector(lam), K)


def Pw_series(lam: StrictPartition, K: int) -> list[Fraction]:
    """``P'_lam(t) = sum_j d'_j(lam) t^j/j!`` truncated at ``t^K``."""
    return egf(dw_vector(lam), K)


def series_mul(a, b, K: int) -> list[Fraction]:
    c = [Fraction(0)] * (K + 1)
    for i, ai in enumerate(a[: K + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: K + 1 - i]):
            c[i + j] += ai * bj
    return c


def series_div(a, b, K: int) -> list[Fraction]:
    """``a / b`` as a power series mod ``t^{K+1}``; requires ``b[0] != 0``."""
    assert b[0] != 0
    q = [Fraction(0)] * (K + 1)
    for k in range(K + 1):
        s = Fraction(a[k]) if k < len(a) else Fraction(0)
        for j in range(1, k + 1):
            if j < len(b):
                s -= b[j] * q[k - j]
        q[k] = s / b[0]
    return q


def series_scale(a, c, K: int) -> list[Fraction]:
    """``a(c t)``."""
    return [Fraction(a[k]) * Fraction(c) ** k if k < len(a) else Fraction(0) for k in range(K + 1)]


def series_exp_poly(a, b, K: int) -> list[Fraction]:
    """``exp(a t + b t^2)`` mod ``t^{K+1}`` via ``f' = (a + 2 b t) f``."""
    f = [Fraction(0)] * (K + 1)
    f[0] = Fraction(1)
    a, b = Fraction(a), Fraction(b)
    for k in range(K):
        # (k+1) f_{k+1} = a f_k + 2 b f_{k-1}
        f[k + 1] = (a * f[k] + (2 * b * f[k - 1] if k >= 1 else 0)) / (k + 1)
    return f


def series_log(a, K: int) -> list[Fraction]:
    """``log a`` mod ``t^{K+1}``, ``a[0] == 1``."""
    assert a[0] == 1
    # (log a)' = a'/a
    da = [(k + 1) * a[k + 1] for k in range(K)] + [Fraction(0)]
    q = series_div(da, a, K)
    return [Fraction(0)] + [q[k - 1] / k for k in range(1, K + 1)]


# ---------------------------------------------------------------------------
# Fomin's commutation relations on the shifted lattice
# ---------------------------------------------------------------------------

def commutator_DU_UD(lam: StrictPartition, weighted: bool) -> dict[StrictPartition, int]:
    """Coefficients of ``(D U - U D) lam`` (``D = D'`` weighted if ``weighted``)
    as a dict strict partition -> integer (zero entries dropped)."""
    res: dict[StrictPartition, int] = {}
    for mu in up_set(lam):
        if weighted:
            for nu, w in down_set_weighted(mu):
                res[nu] = res.get(nu, 0) + w
        else:
            for nu in down_set(mu):
                res[nu] = res.get(nu, 0) + 1
    if weighted:
        for mu, w in down_set_weighted(lam):
            for nu in up_set(mu):
                res[nu] = res.get(nu, 0) - w
    else:
        for mu in down_set(lam):
            for nu in up_set(mu):
                res[nu] = res.get(nu, 0) - 1
    return {k: v for k, v in res.items() if v != 0}
