"""skeptic2_targeted.py -- targeted search for candidate counterexamples to Conjecture B, n <= NMAX.

One representative per transpose class.  Transpose-invariant prefilters (own code, exact integers):
  K0 = (n, multiset of unordered corner-box pairs {a_i, b_i})     [fixes d_0..d_{2m}, local product theorem]
  Kf = K0 + f^lam (hook length formula)                           [fixes d_n, d_{n-1}, d_{n-2}]
  KC = K0 + (C_2, C_1^2, C_4)                                     [with f: fixes d_{n-3}..d_{n-6}]
  K2 = K0 + f + (C_2, C_1^2, C_4)
  K3 = K2 + (C_6 - 16 C_1 C_3)                                    [level-7 invariant, checked in skeptic2_level7.py]
For every K2-group of size >= 2 the exact d-vectors are computed and compared.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from skeptic2_dvec import partitions, conjugate, removable, remove_box

def corner_runs(lam):
    """[(a_i, b_i)]: distinct part sizes alpha_1 > ... > alpha_r, a_i = alpha_i - alpha_{i+1}, b_i = multiplicity."""
    runs = []
    i = 0
    L = len(lam)
    while i < L:
        j = i
        while j < L and lam[j] == lam[i]:
            j += 1
        nxt = lam[j] if j < L else 0
        runs.append((lam[i] - nxt, j - i))
        i = j
    return runs

def hook_product(lam):
    c = conjugate(lam)
    p = 1
    for i, row in enumerate(lam):
        for j in range(row):
            p *= (row - j) + (c[j] - i) - 1
    return p

def content_sums(lam):
    C1 = C2 = C3 = C4 = C6 = 0
    for i, row in enumerate(lam):
        for j in range(row):
            c = j - i
            c2 = c * c
            C1 += c; C2 += c2; C3 += c2 * c; C4 += c2 * c2; C6 += c2 * c2 * c2
    return C1, C2, C3, C4, C6

_memo = {(): [1]}
def dvec(lam):
    lam = tuple(lam)
    if lam in _memo:
        return _memo[lam]
    n = sum(lam)
    d = [0] * (n + 1); d[0] = 1
    for i in removable(lam):
        pd = dvec(remove_box(lam, i))
        for j in range(1, n + 1):
            d[j] += pd[j - 1]
    _memo[lam] = d
    return d

def run(n):
    from math import factorial
    nf = factorial(n)
    K0 = {}; Kf = {}; KC = {}; K2 = {}; K3 = {}
    nclasses = 0
    for lam in partitions(n):
        c = conjugate(lam)
        if c < lam:
            continue
        nclasses += 1
        cr = corner_runs(lam)
        k0 = tuple(sorted(tuple(sorted(p)) for p in cr))
        C1, C2, C3, C4, C6 = content_sums(lam)
        f = nf // hook_product(lam)
        kc = (C2, C1 * C1, C4)
        K0.setdefault(k0, []).append(lam)
        Kf.setdefault((k0, f), []).append(lam)
        KC.setdefault((k0, kc), []).append(lam)
        K2.setdefault((k0, f, kc), []).append(lam)
        K3.setdefault((k0, f, kc, C6 - 16 * C1 * C3), []).append(lam)
    def npairs(G):
        return sum(len(g) * (len(g) - 1) // 2 for g in G.values())
    surv = [g for g in K2.values() if len(g) >= 2]
    coll = []
    for g in surv:
        vs = [dvec(l) for l in g]
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                if vs[i] == vs[j]:
                    coll.append((g[i], g[j]))
    print("n=%2d classes=%7d sameK0=%9d sameK0+f=%6d sameK0+C=%6d sameK2=%4d sameK3=%4d collisions=%d" % (
        n, nclasses, npairs(K0), npairs(Kf), npairs(KC), npairs(K2), npairs(K3), len(coll)), flush=True)
    for g in surv:
        vs = [dvec(l) for l in g]
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                D = [t for t in range(n + 1) if vs[i][t] != vs[j][t]]
                print("   K2 survivor:", g[i], g[j], "first diff", D[0], "last diff", D[-1], "ham", len(D), flush=True)
    for cpair in coll:
        print("   !!! COLLISION:", cpair, flush=True)
    _memo.clear(); _memo[()] = [1]

if __name__ == "__main__":
    n0 = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    n1 = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    for n in range(n0, n1 + 1):
        t = time.time()
        run(n)
