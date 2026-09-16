"""fathooks (v2): verify the three local 'window' formulas for the 2-corner shape
   lam(a,b,c,d) = ((a+c)^b, c^d),   a,b,c,d >= 1,   n = ab+bc+cd.
Notation (exponential generating functions, exact rationals):
   I(t)   = sum_k I_k t^k/k!            (involution numbers, = e^{t+t^2/2})
   G_{x,y}(t) = P_{x^y}(t) = sum_{rho in y x x box} f^rho t^|rho|/|rho|!   (G_{x,y} = G_{y,x})
   D_x(t) = sum_{rho : rho_1 > x} f^rho t^|rho|/|rho|!      so  sum_{l(rho)<=x} = I - D_x
Claims tested (all j within range 0..n):
 (W1)  d_j(lam) = j! [t^j] G_{a,b} G_{c,d}                          for j <= a+d,
       d_{a+d+1}(lam) = (a+d+1)! [t^{a+d+1}] G_{a,b} G_{c,d} + C(a+d, a).
 (W2)  d_j(lam) = j! [t^j] I G_{a,d}   ( = s_j(a^d), the up-census of the rectangle a^d )   for j <= min(b,c),
       d_{m+1}(lam) = (m+1)! [t^{m+1}] I G_{a,d} - 1 - [b = c],   m = min(b,c).
 (BOX) d_j(lam) = #{ j-step saturated chains a^d = rho^0 < rho^1 < ... < rho^j with every rho^i inside the (b+d) x (a+c) box }.
 (W3)  s_k(lam) = k! [t^k] (I - D_b) G_{a,d} (I - D_c)                for k <= min(a+b, c+d),
       s_{K+1}(lam) = (K+1)! [t^{K+1}] (I - D_b) G_{a,d} (I - D_c) + [a+b = K] C(a+b, a) + [c+d = K] C(c+d, c),  K = min(a+b,c+d).
Usage: python3 fathooks_v2_windows.py N   (all 2-corner shapes with n <= N; default 40).
"""
import os, sys
from fractions import Fraction
from math import factorial, comb
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from young import partitions, f_hook, involutions, up_set, size
from census import d_vector, s_up

DEG = 60


def fathook(a, b, c, d):
    return tuple([a + c] * b + [c] * d)


def fathooks_upto(N):
    for b in range(1, N + 1):
        for d in range(1, N + 1):
            for c in range(1, N + 1):
                if b * c + c * d > N:
                    break
                for a in range(1, N + 1):
                    if a * b + b * c + c * d > N:
                        break
                    yield (a, b, c, d)


def egf(seq):
    return [Fraction(x, factorial(j)) for j, x in enumerate(seq)]


def mul(p, q, deg=DEG):
    r = [Fraction(0)] * (deg + 1)
    for i, x in enumerate(p[:deg + 1]):
        if x == 0:
            continue
        for j, y in enumerate(q[:deg + 1 - i]):
            r[i + j] += x * y
    return r


I_ser = egf([involutions(k) for k in range(DEG + 1)])
_D_cache = {}


def D_ser(x):
    if x not in _D_cache:
        _D_cache[x] = egf([sum(f_hook(r) for r in partitions(k) if r and r[0] > x) for k in range(DEG + 1)])
    return _D_cache[x]


def IminusD(x):
    return [I_ser[k] - D_ser(x)[k] for k in range(DEG + 1)]


_G_cache = {}


def G_ser(x, y):
    key = (min(x, y), max(x, y))
    if key not in _G_cache:
        _G_cache[key] = egf([sum(f_hook(r) for r in partitions(k) if len(r) <= y and (not r or r[0] <= x)) for k in range(DEG + 1)])
    return _G_cache[key]


def coeff(ser, j):
    return ser[j] * factorial(j)


def box_chains(a, b, c, d, jmax):
    """number of j-step up-chains from a^d inside the (b+d) x (a+c) box, j = 0..jmax (exact, by recursion)."""
    from functools import lru_cache
    rows, cols = b + d, a + c

    @lru_cache(maxsize=None)
    def cnt(rho, j):
        if j == 0:
            return 1
        tot = 0
        for mu in up_set(rho):
            if len(mu) <= rows and mu[0] <= cols:
                tot += cnt(mu, j - 1)
        return tot
    return [cnt(tuple([a] * d), j) for j in range(jmax + 1)]


def main(N, do_box=True, do_s=True):
    bad = 0
    cnt = 0
    for (a, b, c, d) in fathooks_upto(N):
        lam = fathook(a, b, c, d)
        n = sum(lam)
        cnt += 1
        dv = d_vector(lam)
        # (W1)
        P1 = mul(G_ser(a, b), G_ser(c, d))
        for j in range(0, min(a + d, n) + 1):
            if dv[j] != coeff(P1, j):
                bad += 1
                print("W1 FAIL", (a, b, c, d), j)
        if a + d + 1 <= n and dv[a + d + 1] != coeff(P1, a + d + 1) + comb(a + d, a):
            bad += 1
            print("W1' FAIL", (a, b, c, d))
        # (W2)
        P2 = mul(I_ser, G_ser(a, d))
        m = min(b, c)
        for j in range(0, min(m, n) + 1):
            if dv[j] != coeff(P2, j):
                bad += 1
                print("W2 FAIL", (a, b, c, d), j)
        if m + 1 <= n and dv[m + 1] != coeff(P2, m + 1) - 1 - (b == c):
            bad += 1
            print("W2' FAIL", (a, b, c, d), dv[m + 1], coeff(P2, m + 1) - 1 - (b == c))
        # (BOX) full vector
        if do_box and n <= 24:
            bc = box_chains(a, b, c, d, n)
            if tuple(bc) != tuple(dv):
                bad += 1
                print("BOX FAIL", (a, b, c, d))
        # (W3)
        if do_s and n <= 22:
            K = min(a + b, c + d)
            P3 = mul(mul(IminusD(b), G_ser(a, d)), IminusD(c))
            for k in range(0, K + 2):
                sk = s_up(lam, k)
                pred = coeff(P3, k)
                if k == K + 1:
                    pred += (a + b == K) * comb(a + b, a) + (c + d == K) * comb(c + d, c)
                if sk != pred:
                    bad += 1
                    print("W3 FAIL", (a, b, c, d), k, sk, pred)
    print(f"2-corner shapes with n <= {N}: {cnt} shapes; failures: {bad}  (W1,W1',W2,W2' all n<={N}; BOX n<=24; W3,W3' n<=22)")


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 40)
