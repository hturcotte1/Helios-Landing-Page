"""Referee #1, independent check of fathooks.md Lemma 2 (box / up-census window).
No use of young.py / census.py; exact integers and Fractions only.

For lam = lam(a,b,c,d) = ((a+c)^b, c^d), B = (b+d) x (a+c) box, m = min(b,c):
 (R)  rotation bijection, per nu: for every nu <= lam, rho(nu) := rot(B \ nu) is a partition,
      a^d <= rho(nu) <= B, |rho| = ad + |lam/nu|, f^{lam/nu} = f^{rho/a^d}; the map is injective and hits
      every rho with a^d <= rho <= B.                                                  (n <= NR)
 (BOX) d_j(lam) = # j-step up-chains from a^d with top inside B, for ALL j (full vector) (n <= NBOX)
 (ii) d_j(lam) = s_j(a^d) (brute-force up-census, unbounded) = j![t^j] I(t) G_{a,d}(t) for j <= m (n <= NMAX)
 (iii) d_{m+1}(lam) = (m+1)![t^{m+1}] I G_{a,d} - [b=m] - [c=m]; and directly: # (m+1)-chains from a^d
      leaving B equals [b=m]+[c=m], and the escaping tops are exactly (a^d,1^{m+1}) / (a+m+1,a^{d-1})  (n <= NMAX)
 (S)  sharpness: (ii) fails at j = m+1 for every shape; (iii)'s correction is not enough at j = m+2 (report count)
 (V)  the hypothesis 'm+1 <= n' is vacuous (m+1 <= n always).
d_j(lam) is computed here by the removal recursion AND (for n <= NR) by brute-force sum over nu of f^{lam/nu}.
"""
import sys
from functools import lru_cache
from math import comb, factorial
from fractions import Fraction

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
NBOX = int(sys.argv[2]) if len(sys.argv) > 2 else 22
NR   = int(sys.argv[3]) if len(sys.argv) > 3 else 16

# ---------------- own primitives ----------------
def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for p in range(min(n, maxp), 0, -1):
        for rest in partitions(n - p, p): yield (p,) + rest

def conj(lam):
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()

def f_hook(lam):
    n = sum(lam)
    if n == 0: return 1
    lt = conj(lam); h = 1
    for i, r in enumerate(lam):
        for j in range(r): h *= (r - j) + (lt[j] - i) - 1
    assert factorial(n) % h == 0
    return factorial(n) // h

def corners(lam):
    return [i for i in range(len(lam)) if i == len(lam) - 1 or lam[i] > lam[i + 1]]

def remove(lam, i):
    l = list(lam); l[i] -= 1
    if l[i] == 0: l.pop()
    return tuple(l)

def addables(lam):
    return [i for i in range(len(lam) + 1) if i == 0 or (i == len(lam) and lam[-1] >= 1) or (i < len(lam) and lam[i - 1] > lam[i])]

def add(lam, i):
    l = list(lam)
    if i == len(l): l.append(1)
    else: l[i] += 1
    return tuple(l)

def contains(lam, nu):
    return len(nu) <= len(lam) and all(nu[i] <= lam[i] for i in range(len(nu)))

@lru_cache(maxsize=None)
def dvec(lam):
    n = sum(lam)
    if n == 0: return (1,)
    d = [0] * (n + 1); d[0] = 1
    for i in corners(lam):
        sub = dvec(remove(lam, i))
        for j, v in enumerate(sub): d[j + 1] += v
    return tuple(d)

@lru_cache(maxsize=None)
def f_skew(lam, nu):
    """number of saturated chains nu -> lam (down-recursion on lam)."""
    if not contains(lam, nu): return 0
    if lam == nu: return 1
    return sum(f_skew(remove(lam, i), nu) for i in corners(lam))

@lru_cache(maxsize=None)
def s_up(lam, k):
    if k == 0: return 1
    return sum(s_up(add(lam, i), k - 1) for i in addables(lam))

@lru_cache(maxsize=None)
def chains_in_box(rho, k, R, C):
    """# k-step up-chains from rho staying inside the R x C box."""
    if k == 0: return 1
    tot = 0
    for i in addables(rho):
        mu = add(rho, i)
        if len(mu) <= R and mu[0] <= C:
            tot += chains_in_box(mu, k - 1, R, C)
    return tot

@lru_cache(maxsize=None)
def tops_escaping(rho, k, R, C):
    """multiset (as sorted tuple) of tops of k-step chains from rho whose top lies outside the box."""
    if k == 0: return ()
    out = []
    for i in addables(rho):
        mu = add(rho, i)
        if k == 1:
            if not (len(mu) <= R and mu[0] <= C): out.append(mu)
        else:
            out.extend(tops_escaping(mu, k - 1, R, C))
    return tuple(sorted(out))

# series (EGF coefficients as Fractions)
def I_ser(N):   # exp(t + t^2/2)
    I = [0] * (N + 1); I[0] = 1
    if N >= 1: I[1] = 1
    for k in range(2, N + 1): I[k] = I[k - 1] + (k - 1) * I[k - 2]
    return [Fraction(I[k], factorial(k)) for k in range(N + 1)]

def G_ser(x, y, N):  # P_{x^y}(t) = sum_{rho in y x x box} f^rho t^|rho|/|rho|!
    G = [Fraction(0)] * (N + 1)
    for k in range(min(N, x * y) + 1):
        G[k] = Fraction(sum(f_hook(r) for r in partitions(k) if len(r) <= y and (not r or r[0] <= x)), factorial(k))
    return G

def mul(p, q):
    N = len(p) - 1
    r = [Fraction(0)] * (N + 1)
    for i in range(N + 1):
        if p[i]:
            for j in range(N + 1 - i): r[i + j] += p[i] * q[j]
    return r

def coeff(ser, j):
    v = ser[j] * factorial(j); assert v.denominator == 1; return int(v)

def fathook(a, b, c, d): return tuple([a + c] * b + [c] * d)

def fathooks_upto(N):
    for a in range(1, N + 1):
        for b in range(1, N + 1):
            if a * b >= N: break
            for c in range(1, N + 1):
                if a * b + b * c >= N: break
                for d in range(1, N + 1):
                    n = a * b + b * c + c * d
                    if n > N: break
                    yield a, b, c, d, n

def rot_complement(nu, R, C):
    """rotation by 180 degrees of B \ nu, returned as a tuple of row lengths (raise if not a partition)."""
    rows = [C - (nu[i] if i < len(nu) else 0) for i in range(R)][::-1]
    while rows and rows[-1] == 0: rows.pop()
    assert all(rows[i] >= rows[i + 1] for i in range(len(rows) - 1)), (nu, rows)
    return tuple(rows)

def main():
    bad = 0; cnt = 0; cntbox = 0; cntR = 0; sharp_m2_fail = 0; sharp_m2_tot = 0
    for a, b, c, d, n in fathooks_upto(NMAX):
        cnt += 1
        lam = fathook(a, b, c, d); m = min(b, c); R, C = b + d, a + c
        rect = tuple([a] * d)
        dv = dvec(lam)
        assert len(dv) == n + 1
        # (V) hypothesis vacuous
        if not (m + 1 <= n): bad += 1; print("V FAIL", (a, b, c, d))
        # (ii)
        IG = mul(I_ser(m + 2), G_ser(a, d, m + 2))
        for j in range(m + 1):
            sj = s_up(rect, j)
            if dv[j] != sj or dv[j] != coeff(IG, j):
                bad += 1; print("ii FAIL", (a, b, c, d), j, dv[j], sj, coeff(IG, j))
        # (iii)
        pred = coeff(IG, m + 1) - (b == m) - (c == m)
        if dv[m + 1] != pred: bad += 1; print("iii FAIL", (a, b, c, d), dv[m + 1], pred)
        if s_up(rect, m + 1) != coeff(IG, m + 1): bad += 1; print("s_{m+1} series FAIL", (a, b, c, d))
        esc = tops_escaping(rect, m + 1, R, C)
        expected = []
        if b == m: expected.append(tuple([a] * d + [1] * (m + 1)))
        if c == m: expected.append(tuple([a + m + 1] + [a] * (d - 1)))
        if esc != tuple(sorted(expected)): bad += 1; print("escape FAIL", (a, b, c, d), esc, expected)
        # (S) sharpness
        if dv[m + 1] == s_up(rect, m + 1): bad += 1; print("S: (ii) does not fail at m+1 ?!", (a, b, c, d))
        if m + 2 <= n:
            sharp_m2_tot += 1
            if dv[m + 2] != s_up(rect, m + 2) - (b == m) - (c == m): sharp_m2_fail += 1
        # (BOX)
        if n <= NBOX:
            cntbox += 1
            bx = tuple(chains_in_box(rect, j, R, C) for j in range(n + 1))
            if bx != dv: bad += 1; print("BOX FAIL", (a, b, c, d), bx, dv)
        # (R) rotation bijection per nu, and brute-force d-vector
        if n <= NR:
            cntR += 1
            seen = set(); dbrute = [0] * (n + 1)
            for k in range(n + 1):
                for nu in partitions(k):
                    if not contains(lam, nu): continue
                    rho = rot_complement(nu, R, C)
                    if not (contains(rho, rect) and len(rho) <= R and (not rho or rho[0] <= C)):
                        bad += 1; print("R range FAIL", (a, b, c, d), nu, rho)
                    if sum(rho) != a * d + n - k: bad += 1; print("R size FAIL", (a, b, c, d), nu, rho)
                    fs = f_skew(lam, nu)
                    if fs != f_skew(rho, rect): bad += 1; print("R fskew FAIL", (a, b, c, d), nu, rho)
                    dbrute[n - k] += fs
                    seen.add(rho)
            target = set(rho for kk in range(a * d, R * C + 1) for rho in partitions(kk)
                         if contains(rho, rect) and len(rho) <= R and rho[0] <= C)
            if seen != target: bad += 1; print("R bijection FAIL", (a, b, c, d), len(seen), len(target))
            if tuple(dbrute) != dv: bad += 1; print("R dvec FAIL", (a, b, c, d))
    print(f"[ii, iii, escape, V, sharp: n<={NMAX}] {cnt} shapes; [BOX: n<={NBOX}] {cntbox} shapes; "
          f"[R: n<={NR}] {cntR} shapes; failures = {bad}")
    print(f"sharpness at m+2: correction -[b=m]-[c=m] wrong for {sharp_m2_fail} of {sharp_m2_tot} shapes (informational)")
    # all-parameters-1 and other boundary shapes explicitly
    for p in [(1,1,1,1),(1,1,1,2),(2,1,1,1),(1,2,1,1),(1,1,2,1),(3,1,1,1),(1,1,1,3),(1,3,1,1),(1,1,3,1),(2,2,2,2)]:
        a,b,c,d = p; lam = fathook(*p); m = min(b,c); rect = tuple([a]*d)
        print(p, lam, "m=",m, "d=",dvec(lam)[:m+3], "s(a^d)=",[s_up(rect,j) for j in range(m+3)], "corr=",(b==m)+(c==m))

if __name__ == "__main__":
    main()
