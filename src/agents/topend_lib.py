"""topend angle: library.
Independent Murnaghan-Nakayama characters (explicit rim-hook removal on the diagram, not beta-numbers),
content power sums, square-root counts sigma(rho), even-type partitions, and the
class-eigenvalue omega_lam(K_rho) = |K_rho| chi^lam(rho u 1^{n-|rho|}) / f^lam.
Exact integer / Fraction arithmetic only."""
import os, sys
from fractions import Fraction
from functools import lru_cache
from math import factorial, comb
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from young import partitions, conjugate, boxes, f_hook, f_skew, involutions

# ---------------------------------------------------------------- rim hooks
def rim_hooks(lam, k):
    """All ways to remove a rim hook (border strip) of size k from lam.
    Yields (mu, leg) where leg = (number of rows of the strip) - 1.
    Implementation: a rim hook of size k is determined by its lowest-leftmost box; we walk the rim.
    We use the standard description: strips correspond to pairs (i, j) of rows with
    lam_i - i - (lam_j - j) ... but here we do it by explicit diagram surgery for independence."""
    ell = len(lam)
    lamc = conjugate(lam)
    res = []
    # the rim (boundary) boxes in order from top-right to bottom-left:
    rim = []
    for i in range(ell):
        # boxes (i, j) with j from lam[i]-1 down to (lam[i+1]-1 if i+1<ell else 0)
        lo = lam[i + 1] - 1 if i + 1 < ell else 0
        lo = max(lo, 0)
        for j in range(lam[i] - 1, lo - 1, -1):
            rim.append((i, j))
    # a rim hook is a contiguous segment of the rim (in this order) of length k whose removal leaves a partition
    for start in range(len(rim) - k + 1):
        seg = rim[start:start + k]
        S = set(seg)
        # removal leaves a partition iff the segment is "closed": for every box (i,j) removed, (i, j+1) is
        # not in lam or removed, and (i+1, j) not in lam or removed.
        ok = True
        for (i, j) in seg:
            if j + 1 < lam[i] and (i, j + 1) not in S: ok = False; break
            if i + 1 < ell and j < lam[i + 1] and (i + 1, j) not in S: ok = False; break
        if not ok: continue
        # contiguity of segment as a connected strip: consecutive rim boxes are adjacent iff same row or
        # (i,j)->(i+1, j) ... rim order guarantees adjacency only if the segment does not "jump".
        for a, b in zip(seg, seg[1:]):
            if not ((a[0] == b[0] and a[1] == b[1] + 1) or (a[1] == b[1] and b[0] == a[0] + 1)):
                ok = False; break
        if not ok: continue
        mu = list(lam)
        for (i, j) in seg:
            mu[i] -= 1
        while mu and mu[-1] == 0: mu.pop()
        assert all(mu[t] >= mu[t + 1] for t in range(len(mu) - 1)), (lam, seg, mu)
        rows = {i for (i, _) in seg}
        res.append((tuple(mu), len(rows) - 1))
    return res

@lru_cache(maxsize=None)
def chi(lam, rho):
    """Murnaghan-Nakayama: chi^lam(rho); rho a tuple (any order) of cycle lengths, sum = |lam|."""
    if sum(rho) != sum(lam):
        raise ValueError((lam, rho))
    if not rho:
        return 1
    rho = tuple(sorted(rho, reverse=True))
    if rho[0] == 1:
        return f_hook(lam)
    k, rest = rho[0], rho[1:]
    return sum((-1) ** leg * chi(mu, rest) for mu, leg in rim_hooks(lam, k))

def chi_padded(lam, rho):
    """chi^lam(rho u 1^{n-|rho|})."""
    n = sum(lam); s = sum(rho)
    return chi(lam, tuple(rho) + (1,) * (n - s))

# ------------------------------------------------------------ contents etc.
def contents(lam):
    return [j - i for (i, j) in boxes(lam)]

def C(lam, k):
    return sum(c ** k for c in contents(lam))

def z(rho):
    from collections import Counter
    res = 1
    for k, m in Counter(rho).items():
        res *= k ** m * factorial(m)
    return res

def is_even_type(rho):
    from collections import Counter
    return all(m % 2 == 0 for k, m in Counter(rho).items() if k % 2 == 0)

def dfact(m):  # (m-1)!! for even m (number of perfect matchings on m points), = 1 for m = 0
    r = 1
    for t in range(m - 1, 0, -2): r *= t
    return r

def sigma(rho):
    """Number of square roots of a permutation of cycle type rho (formula; brute-force checked in topend_verify)."""
    from collections import Counter
    res = 1
    for k, m in Counter(rho).items():
        if k % 2 == 0:
            if m % 2: return 0
            res *= dfact(m) * k ** (m // 2)
        else:
            res *= sum(comb(m, 2 * j) * dfact(2 * j) * k ** j for j in range(m // 2 + 1))
    return res

def strip_ones(rho):
    return tuple(p for p in rho if p > 1)

def class_size(rho, n):
    """|K_{rho u 1^{n-|rho|}}| = n!/(z_rho (n-|rho|)!)."""
    s = sum(rho)
    if n < s: return 0
    return factorial(n) // (z(rho) * factorial(n - s))

def omega(lam, rho):
    """omega_lam(K_rho) = |K_rho| chi^lam(rho)/f^lam as an exact Fraction (rho without 1's, padded)."""
    n = sum(lam)
    if n < sum(rho): return Fraction(0)
    return Fraction(class_size(rho, n) * chi_padded(lam, rho), f_hook(lam))

def even_type_no_ones(max_size):
    res = []
    for s in range(max_size + 1):
        for rho in partitions(s):
            if all(p > 1 for p in rho) and is_even_type(rho):
                res.append(rho)
    return res

def cycle_type(perm):
    n = len(perm); seen = [False] * n; ct = []
    for i in range(n):
        if not seen[i]:
            l = 0; j = i
            while not seen[j]:
                seen[j] = True; j = perm[j]; l += 1
            ct.append(l)
    return tuple(sorted(ct, reverse=True))
