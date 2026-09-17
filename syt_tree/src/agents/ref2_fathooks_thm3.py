"""Referee #2 independent check of fathooks.md section 2 (Theorem 3 + Corollary 3.1).
Own primitives throughout (no import from young/census except for a final cross-check of the d-vector).
Usage: python3 ref2_fathooks_thm3.py NMAX
Checks:
 (A) own d-vector recursion == brute-force filter/linear-extension count for all partitions n <= 9
 (B) G_{x,y} = I - D_x - D_y + E_{x,y} as series to degree DEG (inclusion-exclusion identity used in the proof)
 (C) (3.1): q_j = j![t^j](I - P/I) == sum_k N_k D_k^{(j)} for all j <= x1+x2, every 2-corner shape n <= NMAX;
     also tightness: how often it FAILS at j = x1+x2+1 (proof only claims mod t^{x1+x2+1})
 (D) own reading algorithm (exactly as stated in Thm 3) recovers x1, x2, N_k (k <= x1+x2-1), never reads q_j with
     j > x1+x2, and Corollary 3.1(b): balanced <=> sum N_k == 4, multiset recovered for balanced shapes
 (E) sanity: Lemma 1(i) window a+d >= x1+x2 and n >= x1+x2+1 always
"""
import sys
from fractions import Fraction
from functools import lru_cache
from math import factorial, comb
from itertools import permutations

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
DEG = 46

def parts(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for p in range(min(n, maxp), 0, -1):
        for rest in parts(n - p, p): yield (p,) + rest

def conj(lam):
    return tuple(sum(1 for r in lam if r > j) for j in range(lam[0])) if lam else ()

def fhook(lam):
    n = sum(lam)
    if n == 0: return 1
    ct = conj(lam); H = 1
    for i, r in enumerate(lam):
        for j in range(r):
            H *= (r - j) + (ct[j] - i) - 1
    assert factorial(n) % H == 0
    return factorial(n) // H

def corners(lam):
    return [i for i in range(len(lam)) if lam[i] > (lam[i + 1] if i + 1 < len(lam) else 0)]

def rm(lam, i):
    l = list(lam); l[i] -= 1
    if l[i] == 0: l.pop()
    return tuple(l)

@lru_cache(maxsize=None)
def dvec(lam):
    n = sum(lam)
    if n == 0: return (1,)
    d = [0] * (n + 1); d[0] = 1
    for i in corners(lam):
        dm = dvec(rm(lam, i))
        for j in range(1, n + 1): d[j] += dm[j - 1]
    return tuple(d)

# (A) brute force: filters S of lam (closed under right/down), e(S) = # linear extensions (box removed after all
# boxes weakly right/below it), d_j = sum_{|S|=j} e(S).
def dvec_brute(lam):
    bx = [(i, j) for i, r in enumerate(lam) for j in range(r)]
    n = len(bx)
    d = [0] * (n + 1)
    # enumerate filters as subsets: a subset is a filter iff closed under (i,j)->(i,j+1),(i+1,j) within lam
    bset = set(bx)
    def is_filter(S):
        for (i, j) in S:
            if (i, j + 1) in bset and (i, j + 1) not in S: return False
            if (i + 1, j) in bset and (i + 1, j) not in S: return False
        return True
    @lru_cache(maxsize=None)
    def e(S):  # number of removal orders of filter S (frozenset)
        if not S: return 1
        tot = 0
        for (i, j) in S:
            if (i, j + 1) not in S and (i + 1, j) not in S:  # removable now
                tot += e(S - {(i, j)})
        return tot
    from itertools import combinations
    for k in range(n + 1):
        for S in combinations(bx, k):
            S = frozenset(S)
            if is_filter(S): d[k] += e(S)
    return tuple(d)

# ---- series ----
def inv_num(k):
    return 1 if k <= 1 else inv_num(k - 1) + (k - 1) * inv_num(k - 2)
PART = [[(fhook(r), r[0] if r else 0, len(r)) for r in parts(k)] for k in range(DEG + 1)]
def egf(seq): return [Fraction(x, factorial(j)) for j, x in enumerate(seq)]
def I_ser(): return egf([sum(f for f, _, _ in PART[k]) for k in range(DEG + 1)])
def D_ser(x): return egf([sum(f for f, r1, l in PART[k] if r1 > x) for k in range(DEG + 1)])
def Dt_ser(y): return egf([sum(f for f, r1, l in PART[k] if l > y) for k in range(DEG + 1)])
def E_ser(x, y): return egf([sum(f for f, r1, l in PART[k] if r1 > x and l > y) for k in range(DEG + 1)])
def G_ser(x, y): return egf([sum(f for f, r1, l in PART[k] if r1 <= x and l <= y) for k in range(DEG + 1)])
def mul(p, q):
    r = [Fraction(0)] * (DEG + 1)
    for i, x in enumerate(p):
        if x == 0 or i > DEG: continue
        for j, y in enumerate(q):
            if i + j > DEG: break
            r[i + j] += x * y
    return r
def inv(p):
    r = [Fraction(0)] * (DEG + 1); r[0] = 1 / p[0]
    for k in range(1, DEG + 1):
        r[k] = -sum(p[i] * r[k - i] for i in range(1, k + 1)) / p[0]
    return r
I = I_ser(); Iinv = inv(I)
Dk = {}
def Dcoef(x, j):  # D_x^{(j)} = j! [t^j] D_x = sum_{rho |- j, rho_1 > x} f^rho
    if x not in Dk: Dk[x] = D_ser(x)
    return Dk[x][j] * factorial(j)

def fathook(a, b, c, d): return tuple([a + c] * b + [c] * d)

def shapes(N):
    for b in range(1, N + 1):
        for d in range(1, N + 1):
            for c in range(1, N + 1):
                if b * c + c * d > N: break
                for a in range(1, N + 1):
                    if a * b + b * c + c * d > N: break
                    yield (a, b, c, d)

# ---- (D) reading algorithm exactly as in Theorem 3, instrumented ----
def read_local(dv):
    P = egf(list(dv)) + [Fraction(0)] * (DEG + 1 - len(dv))
    PI = mul(P, Iinv)
    q = [(I[j] - PI[j]) * factorial(j) for j in range(DEG + 1)]
    maxread = [0]
    def Q(j):
        maxread[0] = max(maxread[0], j); return q[j]
    j = 1
    while Q(j) == 0: j += 1
    x1 = j - 1
    def toint(v):
        assert v.denominator == 1; return int(v)
    N = {x1: toint(Q(x1 + 1))}
    assert N[x1] >= 1
    if N[x1] >= 2:
        x2 = x1
    else:
        k = x1 + 1
        while True:
            Nk = toint(Q(k + 1) - sum(N[x] * Dcoef(x, k + 1) for x in N if x < k))
            N[k] = Nk
            if Nk != 0: x2 = k; break
            k += 1
    for k in range(x2 + 1, x1 + x2):
        N[k] = toint(Q(k + 1) - sum(N[x] * Dcoef(x, k + 1) for x in N if x < k))
    return x1, x2, N, maxread[0]

def main():
    fails = 0
    # (A)
    for n in range(0, 10):
        for lam in parts(n):
            if dvec(lam) != dvec_brute(lam): fails += 1; print("A FAIL", lam)
    print("(A) own recursion vs brute-force filters/linear extensions, n<=9: done, fails so far", fails)
    # cross-check with project census
    sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src')
    from census import d_vector
    # (B)
    for x in range(1, 8):
        for y in range(1, 8):
            G = G_ser(x, y); rhs = [I[k] - Dk.setdefault(x, D_ser(x))[k] - D_ser(y)[k] + E_ser(x, y)[k] for k in range(DEG + 1)]
            if G != rhs: fails += 1; print("B FAIL", x, y)
            if D_ser(y) != Dt_ser(y): fails += 1; print("B' FAIL (D_y != transpose sum)", y)
            # order of vanishing claims
            E = E_ser(x, y)
            if any(E[k] != 0 for k in range(x + y + 1)) or E[x + y + 1] == 0: fails += 1; print("B'' E order FAIL", x, y)
            D = D_ser(x)
            if any(D[k] != 0 for k in range(x + 1)) or D[x + 1] != Fraction(1, factorial(x + 1)): fails += 1; print("B''' D order FAIL", x)
    print("(B) inclusion-exclusion + orders of vanishing: fails so far", fails)
    # (C),(D),(E)
    cnt = 0; tight_fail = 0; tight_hold = 0; balanced = 0; maxread_viol = 0
    for (a, b, c, d) in shapes(NMAX):
        lam = fathook(a, b, c, d); n = sum(lam); cnt += 1
        dv = dvec(lam)
        if dv != d_vector(lam): fails += 1; print("dvec mismatch with census", lam)
        xs = sorted([a, b, c, d]); x1, x2 = xs[0], xs[1]
        # (E)
        if a + d < x1 + x2 or n < x1 + x2 + 1: fails += 1; print("E FAIL", (a, b, c, d))
        # (C)
        P = egf(list(dv)) + [Fraction(0)] * (DEG + 1 - len(dv))
        PI = mul(P, Iinv)
        q = [(I[j] - PI[j]) * factorial(j) for j in range(DEG + 1)]
        Nk = {k: xs.count(k) for k in set(xs)}
        for j in range(0, x1 + x2 + 1):
            rhs = sum(Nk[k] * Dcoef(k, j) for k in Nk)
            if q[j] != rhs: fails += 1; print("C FAIL (3.1)", (a, b, c, d), j, q[j], rhs)
        j = x1 + x2 + 1
        if j <= n:
            rhs = sum(Nk[k] * Dcoef(k, j) for k in Nk)
            if q[j] != rhs: tight_fail += 1
            else: tight_hold += 1
        # (D)
        r1, r2, N, mr = read_local(dv)
        truth = {k: xs.count(k) for k in range(x1, x1 + x2)}
        if (r1, r2) != (x1, x2) or N != truth: fails += 1; print("D FAIL read", (a, b, c, d), (r1, r2, N), truth)
        if mr > x1 + x2: maxread_viol += 1; print("D FAIL maxread", (a, b, c, d), mr, x1 + x2)
        isbal = xs[3] <= x1 + x2 - 1
        rec = sum(N.values()) == 4
        if isbal != rec: fails += 1; print("D FAIL balanced recognition", (a, b, c, d))
        if isbal:
            balanced += 1
            ms = sorted(k for k in N for _ in range(N[k]))
            if ms != xs: fails += 1; print("D FAIL multiset", (a, b, c, d), ms)
    print(f"(C),(D),(E) n<={NMAX}: {cnt} shapes, balanced = {balanced}, maxread violations = {maxread_viol}; "
          f"(3.1) at degree x1+x2+1: holds {tight_hold}, fails {tight_fail}")
    print("TOTAL FAILS =", fails)

if __name__ == '__main__':
    main()
