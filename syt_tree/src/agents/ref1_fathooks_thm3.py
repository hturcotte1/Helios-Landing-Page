"""Referee #1 check of fathooks.md Theorem 3 (local multiset reading) and Corollary 3.1.
Independent implementation: own d-vector recursion, own hook-length f, own series arithmetic (Fractions).
Usage: python3 ref1_fathooks_thm3.py NMAX
"""
import sys
from fractions import Fraction
from functools import lru_cache
from math import factorial, comb
sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src')

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40

# ---------- own primitives ----------
def parts(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for p in range(min(n, maxp), 0, -1):
        for rest in parts(n - p, p): yield (p,) + rest

def f_hook_own(lam):
    n = sum(lam)
    if n == 0: return 1
    conj = [sum(1 for r in lam if r > j) for j in range(lam[0])]
    H = 1
    for i, r in enumerate(lam):
        for j in range(r):
            H *= (r - j - 1) + (conj[j] - i - 1) + 1
    assert factorial(n) % H == 0
    return factorial(n) // H

@lru_cache(maxsize=None)
def dvec_own(lam):
    """(d_0,...,d_n) via d_j(lam) = sum_c d_{j-1}(lam - c)."""
    n = sum(lam)
    if n == 0: return (1,)
    res = [0] * (n + 1); res[0] = 1
    for i in range(len(lam)):
        if i == len(lam) - 1 or lam[i] > lam[i + 1]:
            mu = list(lam); mu[i] -= 1
            mu = tuple(x for x in mu if x > 0)
            sub = dvec_own(mu)
            for j in range(1, n + 1): res[j] += sub[j - 1]
    return tuple(res)

def shape(a, b, c, d): return tuple([a + c] * b + [c] * d)

# series as list of Fractions, index = degree, truncated at DEG
DEG = 60
def ser_from_counts(counts):  # counts[j] -> sum counts[j] t^j/j!
    return [Fraction(counts[j], factorial(j)) if j < len(counts) else Fraction(0) for j in range(DEG + 1)]
def mul(A, B):
    C = [Fraction(0)] * (DEG + 1)
    for i, x in enumerate(A):
        if x == 0: continue
        for j in range(DEG + 1 - i): C[i + j] += x * B[j]
    return C
def inv(A):
    assert A[0] == 1
    B = [Fraction(0)] * (DEG + 1); B[0] = Fraction(1)
    for k in range(1, DEG + 1):
        B[k] = -sum(A[i] * B[k - i] for i in range(1, k + 1))
    return B

# involution numbers, D_x^{(j)} = sum_{rho |- j, rho_1 > x} f^rho
@lru_cache(maxsize=None)
def fsum_by_first_part(j):
    """dict rho_1 -> sum of f^rho over rho |- j with that first part"""
    out = {}
    for rho in parts(j):
        p1 = rho[0] if rho else 0
        out[p1] = out.get(p1, 0) + f_hook_own(rho)
    return out
JMAX = 42
Dtab = {}   # Dtab[x][j]
for x in range(0, JMAX + 1):
    Dtab[x] = [0] * (JMAX + 1)
for j in range(0, JMAX + 1):
    fb = fsum_by_first_part(j)
    for x in range(0, JMAX + 1):
        Dtab[x][j] = sum(v for p, v in fb.items() if p > x)
I_counts = [sum(fsum_by_first_part(j).values()) for j in range(JMAX + 1)]
I_ser = ser_from_counts(I_counts)
Iinv = inv(I_ser)
# sanity: involution recursion
for k in range(2, JMAX + 1): assert I_counts[k] == I_counts[k - 1] + (k - 1) * I_counts[k - 2]
assert Dtab[3][4] == 1 and Dtab[3][5] == 1 + 4  # (5), (4,1)

def read_algorithm(dv):
    """Theorem 3 algorithm exactly as stated in the text (uses only the d-vector and universal D_x^{(j)})."""
    assert dv[1] == 2
    P = ser_from_counts(dv)
    Q = [I_ser[j] - v for j, v in enumerate(mul(P, Iinv))]
    q = [Q[j] * factorial(j) for j in range(DEG + 1)]
    assert all(x.denominator == 1 for x in q)
    q = [int(x) for x in q]
    j = 1
    while q[j] == 0: j += 1
    x1 = j - 1
    N = {x1: q[x1 + 1]}
    assert N[x1] >= 1
    def peel(k):  # N_k = q_{k+1} - sum_{x<k} N_x D_x^{(k+1)}
        return q[k + 1] - sum(N[x] * Dtab[x][k + 1] for x in N if x < k)
    if N[x1] >= 2:
        x2 = x1
    else:
        k = x1 + 1
        while True:
            N[k] = peel(k)
            if N[k] != 0: x2 = k; break
            k += 1
    for k in range(x1 + 1, x1 + x2):
        if k not in N: N[k] = peel(k)
    return x1, x2, N, q

bad = 0; cnt = 0; tight = 0; balanced_cnt = 0; ncase_x1eq1 = 0
maxdeg_seen = 0
for n in range(3, NMAX + 1):
    for a in range(1, n):
        for b in range(1, n):
            for c in range(1, n):
                if b * (a + c) >= n: break
                rem = n - b * (a + c)
                if rem % c: continue
                d = rem // c
                if d < 1: continue
                lam = shape(a, b, c, d); assert sum(lam) == n
                cnt += 1
                dv = dvec_own(lam)
                if n <= 22:
                    from census import d_vector
                    assert dv == d_vector(lam), (lam)
                xs = sorted([a, b, c, d]); x1, x2 = xs[0], xs[1]
                Nk_true = {k: xs.count(k) for k in range(x1, x1 + x2)}
                # (3.1): q_j == sum_k N_k D_k^{(j)} for all j <= x1+x2
                P = ser_from_counts(dv)
                Q = [I_ser[j] - v for j, v in enumerate(mul(P, Iinv))]
                for j in range(0, x1 + x2 + 1):
                    lhs = Q[j] * factorial(j)
                    rhs = sum(xs.count(k) * Dtab[k][j] for k in set(xs))
                    if lhs != rhs: bad += 1; print("(3.1) FAIL", (a, b, c, d), j, lhs, rhs)
                # is it tight at j = x1+x2+1?
                j = x1 + x2 + 1
                if j <= n:
                    lhs = Q[j] * factorial(j); rhs = sum(xs.count(k) * Dtab[k][j] for k in set(xs))
                    if lhs != rhs: tight += 1
                # reading algorithm
                r1, r2, N, q = read_algorithm(dv)
                if (r1, r2) != (x1, x2) or N != Nk_true:
                    bad += 1; print("READ FAIL", (a, b, c, d), (r1, r2, N), (x1, x2, Nk_true))
                # corollary: balanced recognition
                bal = xs[3] <= x1 + x2 - 1
                read_bal = sum(N.values()) == 4
                if bal != read_bal: bad += 1; print("BAL FAIL", (a, b, c, d))
                if bal:
                    balanced_cnt += 1
                    ms = sorted(k for k, v in N.items() for _ in range(v))
                    if ms != xs: bad += 1; print("MULTISET FAIL", (a, b, c, d), ms)
                if x1 == 1: ncase_x1eq1 += 1
print(f"n<={NMAX}: {cnt} two-corner shapes, failures = {bad}; balanced = {balanced_cnt}; x1=1 cases = {ncase_x1eq1}; "
      f"shapes where congruence fails at degree x1+x2+1 = {tight}")
