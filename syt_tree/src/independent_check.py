"""Independent reproduction of the d-vector scan (no imports from the library).
Route A: partitions in multiplicity representation (dict part-size -> multiplicity), corners = distinct
         part sizes; level-by-level recursion; collision report + the same mod-(2^61-1) digest as the scans.
Route B (n <= 24): d-vectors recovered from UP-chain counts s_k(lam) (computed with addable-corner
         recursion, a different primitive) via the proved identity P = exp(-t-t^2/2) E, i.e.
         d_j = sum_i C(j,i) J_{j-i} s_i with J_m = m! [t^m] exp(-t - t^2/2)."""
import sys
from math import comb, factorial
from fractions import Fraction
P = (1 << 61) - 1

def parts_of(n, maxp=None):
    if maxp is None or maxp > n: maxp = n
    if n == 0: yield (); return
    for f in range(maxp, 0, -1):
        for r in parts_of(n - f, f): yield (f,) + r

def to_mult(lam):
    m = {}
    for x in lam: m[x] = m.get(x, 0) + 1
    return m

def from_mult(m):
    out = []
    for x in sorted(m, reverse=True): out += [x] * m[x]
    return tuple(out)

def children_down(m):
    """multiplicity dicts of lam - c, one per distinct part size."""
    res = []
    for x in m:
        mm = dict(m); mm[x] -= 1
        if mm[x] == 0: del mm[x]
        if x - 1 > 0: mm[x - 1] = mm.get(x - 1, 0) + 1
        res.append(from_mult(mm))
    return res

def conj(lam):
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()

def route_A(N):
    prev = {(): (1,)}
    out = {}
    for n in range(1, N + 1):
        cur = {}
        for lam in parts_of(n):
            d = [0] * (n + 1); d[0] = 1
            for mu in children_down(to_mult(lam)):
                pm = prev[mu]
                for j in range(1, n + 1): d[j] += pm[j - 1]
            cur[lam] = tuple(d)
        groups = {}
        for lam, dv in cur.items(): groups.setdefault(dv, []).append(lam)
        bad = [g for g in groups.values() if not (len(g) == 1 and g[0] == conj(g[0]) or len(g) == 2 and g[0] == conj(g[1]))]
        w = [pow(7, j, P) for j in range(n + 1)]
        dig = sum(sum(a * b for a, b in zip(dv, w)) for dv in cur.values()) % P
        out[n] = (len(cur), len(groups), len(bad), dig)
        prev = cur
        if n <= 24: out[('vec', n)] = cur
    return out

def route_B(N):
    """s_k(lam) for |lam| + k <= 2N via addable-corner recursion, then d_j via the identity."""
    K = N
    s = {}  # (lam, k) -> s_k
    def addable(lam):
        res = []
        for i in range(len(lam) + 1):
            cur = lam[i] if i < len(lam) else 0
            if i == 0 or lam[i - 1] > cur:
                new = list(lam)
                if i == len(new): new.append(1)
                else: new[i] += 1
                res.append(tuple(new))
        return res
    for m in range(2 * N, -1, -1):
        for lam in parts_of(m):
            s[(lam, 0)] = 1
            for k in range(1, 2 * N - m + 1):
                s[(lam, k)] = sum(s[(mu, k - 1)] for mu in addable(lam))
    # J_m = m! [t^m] exp(-t - t^2/2)
    J = [0] * (2 * N + 1); J[0] = 1
    if 2 * N >= 1: J[1] = -1
    for m in range(2, 2 * N + 1): J[m] = -J[m - 1] - (m - 1) * J[m - 2]   # from (e^{-t-t^2/2})' = -(1+t) e^{-t-t^2/2}
    vecs = {}
    for n in range(0, N + 1):
        for lam in parts_of(n):
            vecs[lam] = tuple(sum(comb(j, i) * J[j - i] * s[(lam, i)] for i in range(j + 1)) for j in range(n + 1))
    return vecs

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    A = route_A(N)
    print("route A (multiplicity representation): n p(n) #distinct #collisions digest")
    for n in range(1, N + 1):
        print(n, *A[n])
    NB = min(N, 24)
    B = route_B(NB)
    mism = 0
    for n in range(0, NB + 1):
        for lam in parts_of(n):
            if n >= 1 and A[('vec', n)][lam] != B[lam]: mism += 1
            if n == 0 and B[lam] != (1,): mism += 1
    print(f"route B (up-chains + identity) agrees with route A for all partitions of n <= {NB}: {'YES' if mism == 0 else 'NO (%d)' % mism}")
