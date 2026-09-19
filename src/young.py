"""Basic combinatorics of partitions, Young's lattice and standard Young tableaux.

Conventions
-----------
* A partition is a tuple of positive integers in weakly decreasing order; the
  empty partition is ``()``.  Row ``i`` (0-indexed) has ``lam[i]`` boxes.
* A *box* is a pair ``(i, j)`` of 0-indexed (row, column) coordinates.
* All arithmetic is exact (Python integers).
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations
from math import comb, factorial
from typing import Iterator, Tuple

Partition = Tuple[int, ...]


# ---------------------------------------------------------------------------
# Partitions
# ---------------------------------------------------------------------------

def partitions(n: int, max_part: int | None = None) -> Iterator[Partition]:
    """Yield all partitions of ``n`` (parts <= ``max_part``) in reverse
    lexicographic order, i.e. ``(n)`` first and ``(1,)*n`` last."""
    if max_part is None or max_part > n:
        max_part = n
    if n == 0:
        yield ()
        return
    for first in range(max_part, 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def is_partition(lam) -> bool:
    return all(isinstance(x, int) and x > 0 for x in lam) and all(
        lam[i] >= lam[i + 1] for i in range(len(lam) - 1)
    )


def size(lam: Partition) -> int:
    return sum(lam)


def conjugate(lam: Partition) -> Partition:
    """Transpose (conjugate) partition."""
    if not lam:
        return ()
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0]))


def is_symmetric(lam: Partition) -> bool:
    return lam == conjugate(lam)


def contains(lam: Partition, nu: Partition) -> bool:
    """True iff the diagram of ``nu`` is contained in that of ``lam``."""
    if len(nu) > len(lam):
        return False
    return all(nu[i] <= lam[i] for i in range(len(nu)))


def boxes(lam: Partition) -> list[tuple[int, int]]:
    return [(i, j) for i, row in enumerate(lam) for j in range(row)]


# ---------------------------------------------------------------------------
# Corners
# ---------------------------------------------------------------------------

def removable_corners(lam: Partition) -> list[int]:
    """Rows ``i`` whose last box ``(i, lam[i]-1)`` can be removed leaving a
    partition, i.e. ``lam[i] > lam[i+1]`` (with ``lam[len] = 0``)."""
    ell = len(lam)
    return [i for i in range(ell) if lam[i] > (lam[i + 1] if i + 1 < ell else 0)]


def addable_corners(lam: Partition) -> list[int]:
    """Rows ``i`` (``0 <= i <= len(lam)``) to which a box can be added at
    position ``(i, lam[i])`` leaving a partition."""
    ell = len(lam)
    res = []
    for i in range(ell + 1):
        cur = lam[i] if i < ell else 0
        prev = lam[i - 1] if i > 0 else None
        if i == 0 or prev > cur:
            res.append(i)
    return res


def remove_box(lam: Partition, i: int) -> Partition:
    """Remove the last box of row ``i`` (must be a removable corner)."""
    assert i in removable_corners(lam), (lam, i)
    new = list(lam)
    new[i] -= 1
    if new[i] == 0:
        new.pop()
    return tuple(new)


def add_box(lam: Partition, i: int) -> Partition:
    """Add a box at the end of row ``i`` (must be an addable corner)."""
    assert i in addable_corners(lam), (lam, i)
    new = list(lam)
    if i == len(new):
        new.append(1)
    else:
        new[i] += 1
    return tuple(new)


def corner_boxes_removable(lam: Partition) -> list[tuple[int, int]]:
    return [(i, lam[i] - 1) for i in removable_corners(lam)]


def corner_boxes_addable(lam: Partition) -> list[tuple[int, int]]:
    return [(i, lam[i] if i < len(lam) else 0) for i in addable_corners(lam)]


def down_set(lam: Partition) -> list[Partition]:
    """All partitions covered by ``lam`` (one box removed)."""
    return [remove_box(lam, i) for i in removable_corners(lam)]


def up_set(lam: Partition) -> list[Partition]:
    """All partitions covering ``lam`` (one box added)."""
    return [add_box(lam, i) for i in addable_corners(lam)]


# ---------------------------------------------------------------------------
# Corner-run parametrisation  lam <-> ((a_1,b_1),...,(a_r,b_r))
# ---------------------------------------------------------------------------

def corner_runs(lam: Partition) -> list[tuple[int, int]]:
    """Return ``[(a_1,b_1),...,(a_r,b_r)]`` where the distinct part sizes are
    ``alpha_1 > ... > alpha_r``, ``a_i = alpha_i - alpha_{i+1}`` (``alpha_{r+1}=0``)
    and ``b_i`` is the multiplicity of ``alpha_i``.  Transposition reverses
    the list and swaps the two entries of each pair."""
    if not lam:
        return []
    runs = []
    i = 0
    while i < len(lam):
        j = i
        while j < len(lam) and lam[j] == lam[i]:
            j += 1
        nxt = lam[j] if j < len(lam) else 0
        runs.append((lam[i] - nxt, j - i))
        i = j
    return runs


def from_corner_runs(runs: list[tuple[int, int]]) -> Partition:
    lam: list[int] = []
    for k, (a, b) in enumerate(runs):
        alpha = sum(aa for aa, _ in runs[k:])
        lam.extend([alpha] * b)
    return tuple(lam)


# ---------------------------------------------------------------------------
# Hook lengths and f^lambda
# ---------------------------------------------------------------------------

def hook_lengths(lam: Partition) -> dict[tuple[int, int], int]:
    lt = conjugate(lam)
    return {(i, j): lam[i] - j + lt[j] - i - 1 for (i, j) in boxes(lam)}


def f_hook(lam: Partition) -> int:
    """Number of SYT of shape ``lam`` by the hook length formula."""
    n = size(lam)
    prod = 1
    for h in hook_lengths(lam).values():
        prod *= h
    num = factorial(n)
    assert num % prod == 0
    return num // prod


@lru_cache(maxsize=None)
def f_rec(lam: Partition) -> int:
    """Number of SYT of shape ``lam`` by the recursion f^lam = sum_c f^(lam-c)."""
    if not lam:
        return 1
    return sum(f_rec(mu) for mu in down_set(lam))


@lru_cache(maxsize=None)
def f_skew(lam: Partition, nu: Partition) -> int:
    """Number of standard skew tableaux of shape lam/nu, i.e. the number of
    saturated chains from ``nu`` up to ``lam`` in Young's lattice.  Returns 0
    if ``nu`` is not contained in ``lam``."""
    if not contains(lam, nu):
        return 0
    if lam == nu:
        return 1
    return sum(f_skew(mu, nu) for mu in down_set(lam))


# ---------------------------------------------------------------------------
# Standard Young tableaux as explicit fillings
# ---------------------------------------------------------------------------

def syt_of_shape(lam: Partition) -> list[dict[tuple[int, int], int]]:
    """All SYT of shape ``lam`` as dicts box -> entry (entries 1..n)."""
    n = size(lam)
    if n == 0:
        return [{}]
    res = []
    for i in removable_corners(lam):
        mu = remove_box(lam, i)
        for T in syt_of_shape(mu):
            T2 = dict(T)
            T2[(i, lam[i] - 1)] = n
            res.append(T2)
    return res


def shape_of(T: dict[tuple[int, int], int]) -> Partition:
    rows: dict[int, int] = {}
    for (i, j) in T:
        rows[i] = rows.get(i, 0) + 1
    return tuple(rows[i] for i in range(len(rows)))


def transpose_tableau(T: dict[tuple[int, int], int]) -> dict[tuple[int, int], int]:
    return {(j, i): v for (i, j), v in T.items()}


def is_syt(T: dict[tuple[int, int], int]) -> bool:
    lam = shape_of(T)
    if not is_partition(lam) or set(T.values()) != set(range(1, len(T) + 1)):
        return False
    for (i, j), v in T.items():
        if (i, j + 1) in T and T[(i, j + 1)] < v:
            return False
        if (i + 1, j) in T and T[(i + 1, j)] < v:
            return False
    return True


# ---------------------------------------------------------------------------
# Involution numbers (number of SYT of size k) and binomials
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def involutions(k: int) -> int:
    """I_k = number of involutions of [k] = number of SYT with k boxes."""
    if k <= 1:
        return 1
    return involutions(k - 1) + (k - 1) * involutions(k - 2)


def binomial(n: int, k: int) -> int:
    return comb(n, k)
