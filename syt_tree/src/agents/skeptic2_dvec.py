"""skeptic2_dvec.py -- independent d-vector code (no import from census.py / young.py / skeptic_*.py).

d_j(lam) = number of ways to remove j boxes one at a time from lam (order matters),
computed level by level from the recursion d_j(lam) = sum_{c removable} d_{j-1}(lam - c), d_0 = 1.
A second, definitionally different computation (sum over nu subset lam, |lam/nu| = j, of f^{lam/nu},
the latter by Aitken's determinant  f^{lam/nu} = |lam/nu|! det[ 1/(lam_i - nu_j - i + j)! ]) is used as a
cross-check on small shapes.

Usage: python3 skeptic2_dvec.py NMAX  -> per-n collision census + digest comparison with data/collisions_python.txt
"""
import sys, os
from fractions import Fraction
from math import factorial
from itertools import combinations

def partitions(n, maxpart=None):
    """All partitions of n as decreasing tuples."""
    if maxpart is None or maxpart > n:
        maxpart = n
    if n == 0:
        yield ()
        return
    for first in range(maxpart, 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest

def conjugate(lam):
    if not lam:
        return ()
    return tuple(sum(1 for p in lam if p > i) for i in range(lam[0]))

def removable(lam):
    """Indices i such that lam minus one box in row i is a partition."""
    out = []
    L = len(lam)
    for i in range(L):
        if i == L - 1 or lam[i] > lam[i + 1]:
            out.append(i)
    return out

def remove_box(lam, i):
    l = list(lam)
    l[i] -= 1
    if l[i] == 0:
        l.pop()
    return tuple(l)

def all_dvectors(nmax, callback=None):
    """Compute d-vectors for all partitions of n = 0..nmax level by level.
    Returns dict for level nmax; callback(n, dict) is called at each level if given."""
    prev = {(): [1]}
    if callback:
        callback(0, prev)
    for n in range(1, nmax + 1):
        cur = {}
        for lam in partitions(n):
            d = [0] * (n + 1)
            d[0] = 1
            for i in removable(lam):
                pd = prev[remove_box(lam, i)]
                for j in range(1, n + 1):
                    d[j] += pd[j - 1]
            cur[lam] = d
        prev = cur
        if callback:
            callback(n, cur)
    return prev

# ---------- independent cross-check via Aitken's determinant ----------
def det_frac(M):
    n = len(M)
    M = [row[:] for row in M]
    det = Fraction(1)
    for c in range(n):
        p = None
        for r in range(c, n):
            if M[r][c] != 0:
                p = r
                break
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            det = -det
        det *= M[c][c]
        for r in range(c + 1, n):
            if M[r][c] != 0:
                f = M[r][c] / M[c][c]
                for k in range(c, n):
                    M[r][k] -= f * M[c][k]
    return det

def f_skew_aitken(lam, nu):
    L = len(lam)
    nu = tuple(nu) + (0,) * (L - len(nu))
    m = sum(lam) - sum(nu)
    M = []
    for i in range(L):
        row = []
        for j in range(L):
            e = lam[i] - nu[j] - i + j
            row.append(Fraction(1, factorial(e)) if e >= 0 else Fraction(0))
        M.append(row)
    v = det_frac(M) * factorial(m)
    assert v.denominator == 1
    return int(v)

def subpartitions(lam):
    """All partitions nu contained in lam."""
    def rec(i, cap):
        if i == len(lam):
            yield ()
            return
        for p in range(min(lam[i], cap), -1, -1):
            for rest in rec(i + 1, p):
                yield (p,) + rest
    for nu in rec(0, lam[0] if lam else 0):
        yield tuple(p for p in nu if p > 0)

def dvector_by_definition(lam):
    n = sum(lam)
    d = [0] * (n + 1)
    for nu in subpartitions(lam):
        j = n - sum(nu)
        d[j] += f_skew_aitken(lam, nu)
    return d

def digest(vecs):
    P = (1 << 61) - 1
    s = 0
    for d in vecs.values():
        for j, x in enumerate(d):
            s = (s + x * pow(7, j, P)) % P
    return s

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    here = os.path.dirname(os.path.abspath(__file__))
    ref = {}
    try:
        with open(os.path.join(here, "..", "..", "data", "collisions_python.txt")) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                t = line.split()
                ref[int(t[0])] = (int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]), int(t[7]))
    except FileNotFoundError:
        pass

    def report(n, vecs):
        # transpose invariance
        for lam, d in vecs.items():
            assert vecs[conjugate(lam)] == d, ("transpose failure", lam)
        groups = {}
        for lam, d in vecs.items():
            groups.setdefault(tuple(d), []).append(lam)
        nsym = sum(1 for lam in vecs if conjugate(lam) == lam)
        npairs = (len(vecs) - nsym) // 2
        bad = [g for g in groups.values() if not (len(g) == 1 or (len(g) == 2 and conjugate(g[0]) == g[1]))]
        dg = digest(vecs)
        mine = (len(vecs), len(groups), nsym, npairs, len(bad), dg)
        status = "MATCH" if ref.get(n) == mine else ("MISMATCH " + str(ref.get(n)) if n in ref else "noref")
        print(n, *mine, status, flush=True)
        for g in bad:
            print("  COLLISION GROUP:", g, flush=True)

    # cross-check with Aitken on all partitions of n <= 12 and a few larger ones
    import random
    random.seed(1)
    checked = 0
    lv = all_dvectors(12)
    lvls = {}
    all_dvectors(12, callback=lambda n, v: lvls.__setitem__(n, v))
    for n in range(0, 13):
        for lam, d in lvls[n].items():
            assert dvector_by_definition(lam) == d, ("Aitken mismatch", lam)
            checked += 1
    print("Aitken cross-check: all partitions of n<=12 agree (%d shapes)" % checked, flush=True)
    L20 = all_dvectors(20)
    sample = random.sample(sorted(L20), 20)
    for lam in sample:
        assert dvector_by_definition(lam) == L20[lam], ("Aitken mismatch", lam)
    print("Aitken cross-check: 20 random partitions of 20 agree", flush=True)

    print("# n  #partitions  #distinct  #symmetric  #transpose_pairs  #collision_groups  digest  status")
    all_dvectors(nmax, callback=report)
