"""fathooks angle: verify the locality theorem for 2-corner shapes
   lam = ((a+c)^b, c^d):
   (T1) d_j(lam) = j! [t^j] P_{a^b}(t) P_{c^d}(t)  for all j <= a+d,
   (T2) d_{a+d+1}(lam) = (a+d+1)! [t^{a+d+1}] P_{a^b} P_{c^d} + C(a+d, a),
   and check Conjecture B among all 2-corner shapes with n <= N.
"""
import os, sys
from math import factorial, comb
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from young import conjugate, corner_runs
from census import d_vector

def fathook(a, b, c, d):
    return tuple([a + c] * b + [c] * d)

def fathooks_upto(N):
    for b in range(1, N + 1):
        for d in range(1, N + 1):
            for c in range(1, N + 1):
                if b * c + c * d > N: break
                for a in range(1, N + 1):
                    n = a * b + b * c + c * d
                    if n > N: break
                    yield (a, b, c, d)

def egf(dv, deg):
    return [Fraction(dv[j], factorial(j)) if j < len(dv) else Fraction(0) for j in range(deg + 1)]

def poly_mul(p, q, deg):
    r = [Fraction(0)] * (deg + 1)
    for i, x in enumerate(p[:deg + 1]):
        if x == 0: continue
        for j, y in enumerate(q[:deg + 1 - i]):
            r[i + j] += x * y
    return r

def main(N):
    bad = 0; cnt = 0
    dvs = {}
    for (a, b, c, d) in fathooks_upto(N):
        lam = fathook(a, b, c, d)
        assert corner_runs(lam) == [(a, b), (c, d)]
        n = sum(lam); cnt += 1
        dv = d_vector(lam)
        deg = min(n, a + d + 1)
        prod = poly_mul(egf(d_vector(tuple([a] * b)), deg), egf(d_vector(tuple([c] * d)), deg), deg)
        for j in range(0, min(a + d, n) + 1):
            if Fraction(dv[j], factorial(j)) != prod[j]:
                bad += 1; print("T1 FAILS", (a, b, c, d), j)
        if a + d + 1 <= n:
            if Fraction(dv[a + d + 1], factorial(a + d + 1)) != prod[a + d + 1] + Fraction(comb(a + d, a), factorial(a + d + 1)):
                bad += 1; print("T2 FAILS", (a, b, c, d))
        dvs.setdefault(dv, []).append((a, b, c, d))
    # Conjecture B among fathooks
    coll = 0
    for dv, L in dvs.items():
        if len(L) > 1:
            shapes = set(fathook(*x) for x in L)
            if len(shapes) > 2 or (len(shapes) == 2 and conjugate(list(shapes)[0]) not in shapes):
                coll += 1; print("COLLISION", L)
    print(f"checked {cnt} fathooks with n <= {N}: T1/T2 failures = {bad}, non-transpose collisions = {coll}")

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 40)
