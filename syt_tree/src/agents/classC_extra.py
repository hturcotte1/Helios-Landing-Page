"""classC: independent rechecks.
(1) d-vectors of two-row shapes and hooks by direct enumeration of removal sequences (no recursion), n<=12,
    compared with Lemma 4 formulas.
(2) The decision procedure of Step 3 of agent_notes/classC.md, run as an algorithm on the d-vector alone,
    must recover every r=2 member of C (up to transpose) for n<=60.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from math import comb
from young import conjugate, removable_corners, remove_box, corner_runs
from census import d_vector

def binom(n, k): return comb(n, k) if 0 <= k <= n else 0

def d_enum(lam):
    """count removal sequences by explicit depth-first enumeration."""
    n = sum(lam); d = [0] * (n + 1)
    def rec(mu, j):
        d[j] += 1
        for i in removable_corners(mu):
            rec(remove_box(mu, i), j + 1)
    rec(lam, 0)
    return tuple(d)

def d_tworow(l1, l2, j):
    delta = l1 - l2; k0 = max(0, -((delta - j) // 2))
    return sum(binom(j, k) - binom(j, j - delta - 1 - k) for k in range(k0, min(j, l2) + 1))

def d_hook(p, q, j):
    return sum(comb(j, i) for i in range(max(0, j - q), min(j, p) + 1))

def check_enum(N):
    for n in range(1, N + 1):
        for l2 in range(0, n // 2 + 1):
            l1 = n - l2
            lam = (l1, l2) if l2 else (l1,)
            de = d_enum(lam)
            assert de == d_vector(lam)
            assert all(de[j] == d_tworow(l1, l2, j) for j in range(n + 1)), lam
        for q in range(0, n):
            p = n - 1 - q; lam = (p + 1,) + (1,) * q
            de = d_enum(lam)
            assert de == d_vector(lam)
            assert all(de[j] == d_hook(p, q, j) for j in range(n)) and de[n] == comb(n - 1, q), lam
    print(f"(1) direct enumeration agrees with recursion and with Lemma 4 formulas for two-row shapes and hooks, n<={N}")

def decide(dv):
    """Given the d-vector of an r=2 member of C, return the shape (normalised: two-row (l1,l2) with l1>l2>=1,
    or hook (p+1,1^q) with p>=q>=1) following Step 3 of the proof.  Returns None if the vector is not
    consistent with any C-shape (must never happen)."""
    n = len(dv) - 1
    assert dv[1] == 2
    M = 0
    while dv[M + 1] == 2 ** (M + 1): M += 1
    cands = []
    # two-row candidates (delta, l2) with n = 2 l2 + delta, min = M
    for (delta, l2) in {(M, (n - M) // 2) if (n - M) % 2 == 0 else None, (n - 2 * M, M)} - {None}:
        if delta >= 1 and l2 >= 1 and min(delta, l2) == M: cands.append((l2 + delta, l2))
    # hook candidates p>=q with q = M
    q = M; p = n - 1 - q
    if p >= q >= 1: cands.append((p + 1,) + (1,) * q)
    good = [c for c in set(cands) if d_vector(c) == dv]
    return good

def check_decision(N):
    for n in range(3, N + 1):
        shapes = [(n - l2, l2) for l2 in range(1, (n - 1) // 2 + 1)] + [(p + 1,) + (1,) * (n - 1 - p) for p in range((n - 1 + 1) // 2, n - 1)]
        for lam in shapes:
            good = decide(d_vector(lam))
            # normalise hooks/two-row overlaps: the answer set must be {lam} up to transpose
            assert all(g == lam or g == conjugate(lam) for g in good) and good, (lam, good)
    print(f"(2) decision procedure recovers every r=2 C-shape (up to transpose) from its d-vector, n<={N}")

if __name__ == '__main__':
    check_enum(12)
    check_decision(60)
