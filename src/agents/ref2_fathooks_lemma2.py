"""Referee #2 independent check of fathooks.md Lemma 2 (box / up-census window).
Everything below is written from scratch (no import of young/census), exact arithmetic."""
import sys
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from math import factorial, comb

# ---------- basic partition helpers ----------
def down(lam):
    out = []
    for i in range(len(lam)):
        nxt = lam[i+1] if i+1 < len(lam) else 0
        if lam[i] > nxt:
            new = list(lam); new[i] -= 1
            if new[i] == 0: new.pop()
            out.append(tuple(new))
    return out

def up(lam):
    out = []
    for i in range(len(lam)+1):
        cur = lam[i] if i < len(lam) else 0
        prev = lam[i-1] if i > 0 else None
        if i == 0 or prev > cur:
            new = list(lam)
            if i == len(new): new.append(1)
            else: new[i] += 1
            out.append(tuple(new))
    return out

@lru_cache(maxsize=None)
def dvec(lam):
    n = sum(lam)
    if n == 0: return (1,)
    d = [0]*(n+1); d[0] = 1
    for mu in down(lam):
        dm = dvec(mu)
        for j in range(1, n+1): d[j] += dm[j-1]
    return tuple(d)

def conj(lam):
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()

def f_hook(lam):
    n = sum(lam); lt = conj(lam); prod = 1
    for i, r in enumerate(lam):
        for j in range(r):
            prod *= lam[i]-j + lt[j]-i - 1
    assert factorial(n) % prod == 0
    return factorial(n)//prod

def parts_in_box(rows, cols):
    """all partitions with <= rows parts each <= cols"""
    def rec(rem_rows, maxp):
        if rem_rows == 0:
            yield (); return
        for p in range(0, maxp+1):
            if p == 0:
                yield ()
            else:
                for rest in rec(rem_rows-1, p):
                    yield (p,)+rest
    seen = set()
    for p in rec(rows, cols):
        if p not in seen:
            seen.add(p); yield p

# ---------- definition-level cross check of dvec: filters + linear extensions (tiny n) ----------
def dvec_by_filters(lam):
    bx = [(i, j) for i, r in enumerate(lam) for j in range(r)]
    n = len(bx)
    idx = {b: k for k, b in enumerate(bx)}
    d = [0]*(n+1)
    # enumerate filters: subsets closed under moving right/down inside lam
    from itertools import product
    for mask in range(1 << n):
        S = [bx[k] for k in range(n) if mask >> k & 1]
        Sset = set(S)
        ok = True
        for (i, j) in S:
            if (i, j+1) in idx and (i, j+1) not in Sset: ok = False; break
            if (i+1, j) in idx and (i+1, j) not in Sset: ok = False; break
        if not ok: continue
        # linear extensions: count removal orders; box (i,j) removable when (i,j+1),(i+1,j) already removed
        @lru_cache(maxsize=None)
        def cnt(fs):
            if not fs: return 1
            tot = 0
            for b in fs:
                i, j = b
                if (i, j+1) not in fs and (i+1, j) not in fs:
                    tot += cnt(fs - {b})
            return tot
        d[len(S)] += cnt(frozenset(S))
    return tuple(d)

# ---------- Lemma 2 (i): chains from a^d inside the (b+d) x (a+c) box ----------
def box_chain_vector(a, b, c, d):
    R, C = b+d, a+c
    start = tuple([a]*d)
    n = a*b + b*c + c*d
    cur = {start: 1}
    out = [1]
    for j in range(1, n+1):
        nxt = {}
        for rho, w in cur.items():
            for mu in up(rho):
                if len(mu) <= R and mu[0] <= C:
                    nxt[mu] = nxt.get(mu, 0) + w
        cur = nxt
        out.append(sum(cur.values()))
    return tuple(out)

# ---------- up-census by direct chain counting ----------
@lru_cache(maxsize=None)
def s_up(lam, k):
    if k == 0: return 1
    return sum(s_up(mu, k-1) for mu in up(lam))

# ---------- series: I(t) G_{x,y}(t) ----------
def G_coeffs(x, y, N):
    """coefficient list of G_{x,y}(t) = sum_{rho in y x x box} f^rho t^|rho|/|rho|! up to t^N (Fractions)"""
    g = [Fraction(0)]*(N+1)
    for rho in parts_in_box(y, x):
        k = sum(rho)
        if k <= N: g[k] += Fraction(f_hook(rho), factorial(k))
    return g

def I_coeffs(N):
    # e^{t + t^2/2}: use involution numbers I_k
    I = [1, 1]
    for k in range(2, N+1): I.append(I[-1] + (k-1)*I[-2])
    return [Fraction(I[k], factorial(k)) for k in range(N+1)]

def mul(p, q):
    N = len(p)-1
    r = [Fraction(0)]*(N+1)
    for i in range(N+1):
        for j in range(N+1-i):
            r[i+j] += p[i]*q[j]
    return r

def lam_of(a, b, c, d):
    return tuple([a+c]*b + [c]*d)

def main(N=16, NBOX=20, NFILT=9):
    bad = 0; cnt = 0
    # 0. sanity of dvec against filter/linear-extension definition for all partitions n <= NFILT
    def all_parts(n, mx=None):
        if mx is None or mx > n: mx = n
        if n == 0: yield (); return
        for f in range(mx, 0, -1):
            for r in all_parts(n-f, f): yield (f,)+r
    for n in range(NFILT+1):
        for lam in all_parts(n):
            if dvec(lam) != dvec_by_filters(lam):
                bad += 1; print("DVEC DEF FAIL", lam)
    print(f"dvec vs filter definition: all partitions n<={NFILT} checked, failures so far {bad}")

    # 1. Lemma 2(i) full vector, all two-corner shapes n <= NBOX
    c1 = 0
    for a in range(1, NBOX):
        for b in range(1, NBOX):
            for c in range(1, NBOX):
                for d in range(1, NBOX):
                    n = a*b+b*c+c*d
                    if n > NBOX: continue
                    c1 += 1
                    if dvec(lam_of(a, b, c, d)) != box_chain_vector(a, b, c, d):
                        bad += 1; print("BOX FAIL", (a, b, c, d))
    print(f"Lemma 2(i) box formulation: {c1} shapes n<={NBOX}, failures so far {bad}")

    # 2. Lemma 2(ii),(iii) with direct s_up and with series, all two-corner shapes n <= N
    c2 = 0; window_tight = 0; m2_fail_ok = 0
    Iser = I_coeffs(N+3)
    for a in range(1, N):
        for b in range(1, N):
            for c in range(1, N):
                for d in range(1, N):
                    n = a*b+b*c+c*d
                    if n > N: continue
                    c2 += 1
                    lam = lam_of(a, b, c, d); dv = dvec(lam)
                    m = min(b, c); rect = tuple([a]*d)
                    ser = mul(Iser, G_coeffs(a, d, N+3))
                    for j in range(0, m+1):
                        su = s_up(rect, j)
                        sf = ser[j]*factorial(j)
                        if dv[j] != su: bad += 1; print("W2 s_up FAIL", (a,b,c,d), j, dv[j], su)
                        if sf != su: bad += 1; print("SERIES FAIL", (a,b,c,d), j, sf, su)
                    assert m+1 <= n
                    su = s_up(rect, m+1); sf = ser[m+1]*factorial(m+1)
                    if sf != su: bad += 1; print("SERIES FAIL", (a,b,c,d), m+1)
                    corr = (b == m) + (c == m)
                    if dv[m+1] != su - corr: bad += 1; print("W2' FAIL", (a,b,c,d), dv[m+1], su, corr)
                    if dv[m+1] == su: window_tight += 1  # would mean window extends without correction (impossible since corr>=1)
    print(f"Lemma 2(ii)/(iii): {c2} shapes n<={N}, failures so far {bad}; uncorrected m+1 agreements = {window_tight}")
    print("TOTAL FAILURES", bad)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    NBOX = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    main(N, NBOX)
