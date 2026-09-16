"""topend (third pass) -- independent library.

Everything here is written from scratch (no import from the earlier topend_* scripts):
  * Murnaghan-Nakayama characters via explicit border-strip removal on the diagram,
  * content moments, hook lengths / f^lam,
  * square-root counts sigma(rho) by brute force in S_i (i <= 8) and by the closed formula,
  * Aitken's determinant for skew SYT counts f^{lam/nu} (exact Fractions),
  * u_i from characters and from the d-vector.
All arithmetic exact (int / Fraction).
"""
from __future__ import annotations
import sys, os
from fractions import Fraction
from functools import lru_cache
from math import factorial, comb
from itertools import permutations
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from young import partitions, conjugate, f_hook, f_skew, hook_lengths
from census import d_vector

# ---------------------------------------------------------------- contents
def contents(lam):
    return [j - i for i, row in enumerate(lam) for j in range(row)]

def C(lam, k):
    return sum(c ** k for c in contents(lam))

def ffact(n, i):
    r = 1
    for t in range(i):
        r *= n - t
    return r

# ---------------------------------------------------------------- MN characters
def _rim_hooks(lam, k):
    """All ways to remove a border strip of size k from lam.  Yields (mu, height).
    Implemented on the diagram: a border strip is a connected skew shape lam/mu with
    no 2x2 square, i.e. mu is obtained by choosing rows i..i' and, walking the rim.
    We simply enumerate all mu subset lam with |lam/mu| = k and test the strip condition."""
    n = sum(lam)
    res = []
    for mu in _subpartitions_of_size(lam, n - k):
        cells = [(i, j) for i in range(len(lam)) for j in range((mu[i] if i < len(mu) else 0), lam[i])]
        if not cells:
            continue
        # no 2x2 square
        S = set(cells)
        if any((i + 1, j) in S and (i, j + 1) in S and (i + 1, j + 1) in S for (i, j) in cells):
            continue
        # connected (rookwise)
        seen = {cells[0]}; stack = [cells[0]]
        while stack:
            (i, j) = stack.pop()
            for nb in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)):
                if nb in S and nb not in seen:
                    seen.add(nb); stack.append(nb)
        if len(seen) != len(cells):
            continue
        height = len({i for (i, _) in cells}) - 1
        res.append((mu, height))
    return res

def _subpartitions_of_size(lam, m):
    """All partitions mu subset lam with |mu| = m."""
    out = []
    def rec(i, prev, remaining, cur):
        if i == len(lam):
            if remaining == 0:
                mu = tuple(x for x in cur if x > 0)
                out.append(mu)
            return
        hi = min(lam[i], prev, remaining)
        for x in range(hi, -1, -1):
            # remaining rows can hold at most sum(min(lam[t], x)) ...
            cur.append(x); rec(i + 1, x, remaining - x, cur); cur.pop()
            if x == 0:
                break
    rec(0, 10 ** 9, m, [])
    return out

@lru_cache(maxsize=None)
def chi(nu, rho):
    """chi^nu(rho) by Murnaghan-Nakayama, rho a cycle type (tuple, any order)."""
    if not rho:
        return 1 if not nu else 0
    if sum(rho) != sum(nu):
        raise ValueError
    k = max(rho)
    rest = list(rho); rest.remove(k); rest = tuple(sorted(rest, reverse=True))
    tot = 0
    for mu, h in _rim_hooks(nu, k):
        tot += (-1) ** h * chi(mu, rest)
    return tot

def chi_padded(lam, rho):
    """chi^lam(rho cup 1^{n-|rho|}) for lam |- n, computed directly by MN."""
    n = sum(lam)
    full = tuple(sorted(list(rho) + [1] * (n - sum(rho)), reverse=True))
    return chi(lam, full)

def z(rho):
    from collections import Counter
    r = 1
    for k, m in Counter(rho).items():
        r *= k ** m * factorial(m)
    return r

def class_size(rho, n):
    """number of permutations in S_n of cycle type rho cup 1^{n-|rho|}."""
    m = n - sum(rho)
    return factorial(n) // (z(rho) * factorial(m))

def omega(lam, rho):
    """central character omega_lam(K_rho) = |K_rho| chi^lam(rho)/f^lam (Fraction)."""
    n = sum(lam)
    if sum(rho) > n:
        return Fraction(0)
    return Fraction(class_size(rho, n) * chi_padded(lam, rho), f_hook(lam))

# ---------------------------------------------------------------- square roots
def sigma_formula(rho):
    from collections import Counter
    r = 1
    for k, m in Counter(rho).items():
        if k % 2 == 1:
            s = 0
            for j in range(0, m // 2 + 1):
                s += comb(m, 2 * j) * _dfact(2 * j - 1) * k ** j
            r *= s
        else:
            if m % 2:
                return 0
            r *= _dfact(m - 1) * k ** (m // 2)
    return r

def _dfact(x):  # (x)!! for odd x, with (-1)!! = 1
    r = 1
    while x > 1:
        r *= x; x -= 2
    return r

def cycle_type(p):
    n = len(p); seen = [False] * n; ct = []
    for i in range(n):
        if not seen[i]:
            l = 0; j = i
            while not seen[j]:
                seen[j] = True; j = p[j]; l += 1
            ct.append(l)
    return tuple(sorted(ct, reverse=True))

def sigma_brute(rho):
    """number of square roots of a fixed permutation of cycle type rho (brute force)."""
    i = sum(rho)
    # build g of type rho
    g = list(range(i)); pos = 0
    for k in rho:
        for t in range(k):
            g[pos + t] = pos + (t + 1) % k
        pos += k
    cnt = 0
    for h in permutations(range(i)):
        if all(h[h[x]] == g[x] for x in range(i)):
            cnt += 1
    return cnt

# ---------------------------------------------------------------- Aitken
def f_skew_aitken(lam, nu):
    """f^{lam/nu} = |lam/nu|! det[ 1/(lam_a - nu_b - a + b)! ]  (0 if nu not subset lam)."""
    l = len(lam)
    nu = tuple(nu) + (0,) * (l - len(nu))
    if len(nu) > l or any(nu[a] > lam[a] for a in range(l)):
        return 0
    m = sum(lam) - sum(nu)
    M = [[Fraction(1, factorial(lam[a] - nu[b] - a + b)) if lam[a] - nu[b] - a + b >= 0 else Fraction(0)
          for b in range(l)] for a in range(l)]
    d = det(M)
    val = d * factorial(m)
    assert val.denominator == 1
    return int(val)

def det(M):
    M = [row[:] for row in M]; n = len(M); d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            if M[r][c] != 0:
                fac = M[r][c] / M[c][c]
                for k in range(c, n):
                    M[r][k] -= fac * M[c][k]
    return d

# ---------------------------------------------------------------- u_i
def u_from_d(lam, i):
    n = sum(lam)
    return d_vector(lam)[n - i]

def u_from_chars(lam, i):
    """u_i = sum_{rho |- i} sigma(rho)/z_rho * chi^lam(rho cup 1^{n-i})  (Fraction)."""
    tot = Fraction(0)
    for rho in partitions(i):
        s = sigma_formula(rho)
        if s:
            tot += Fraction(s, z(rho)) * chi_padded(lam, rho)
    return tot
