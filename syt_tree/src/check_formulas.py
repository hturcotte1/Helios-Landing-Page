"""Numerical checks of two derived formulas (before they are used in results.md):
 (1) d_{n-3}(lam) = f^lam * ( 2/3 + (C2 - n(n-1)/2) / (n(n-1)(n-2)) ),  C2 = sum of squared contents,
     equivalently 3 d_{n-3} - 2 f = chi^lam(3-cycle) = f^{lam/3} - f^{lam/21} + f^{lam/111}.
 (2) Pfaffian formula: P_lam(t) = Pf(B) with B_ij = [x^{a_i} y^{a_j}] (x-y) e^{t(x+y)} / ((1-xy)(1-x)(1-y)),
     a = lam + delta, ell even.
"""
import os, sys
from fractions import Fraction
from math import factorial
sys.path.insert(0, os.path.dirname(__file__))
from young import partitions, boxes, f_hook, f_skew
from census import d_vector
import sympy as sp

def C2(lam):
    return sum((j - i) ** 2 for (i, j) in boxes(lam))

def check_C2(max_n):
    bad = 0
    for n in range(3, max_n + 1):
        for lam in partitions(n):
            dv = d_vector(lam); f = dv[n]; u3 = dv[n - 3]
            rhs = Fraction(f) * (Fraction(2, 3) + Fraction(C2(lam) - n * (n - 1) // 2, n * (n - 1) * (n - 2)))
            chi3 = f_skew(lam, (3,)) - f_skew(lam, (2, 1)) + f_skew(lam, (1, 1, 1))
            if Fraction(u3) != rhs or 3 * u3 - 2 * f != chi3:
                bad += 1; print("C2 formula FAILS", lam)
    print(f"C2 formula checked for all 3 <= |lam| <= {max_n}:", "OK" if bad == 0 else f"{bad} failures")

def pf_entries(alphas, t):
    x, y = sp.symbols('x y')
    deg = max(alphas)
    # series of (x-y) e^{t(x+y)} / ((1-xy)(1-x)(1-y)) up to x^deg y^deg
    ex = sum((t * x) ** k / sp.factorial(k) for k in range(deg + 1))
    ey = sum((t * y) ** k / sp.factorial(k) for k in range(deg + 1))
    gx = sum(x ** k for k in range(deg + 1)); gy = sum(y ** k for k in range(deg + 1))
    gxy = sum((x * y) ** k for k in range(deg + 1))
    G = sp.expand((x - y) * ex * ey * gx * gy * gxy)
    P = sp.Poly(G, x, y)
    coeff = {}
    for (i, j), c in P.terms():
        if i <= deg and j <= deg:
            coeff[(i, j)] = sp.expand(c)
    return coeff

def pfaffian(M, idx):
    if not idx: return sp.Integer(1)
    i = idx[0]; s = 0
    for k, j in enumerate(idx[1:]):
        rest = idx[1:k+1] + idx[k+2:]
        s += (-1) ** k * M[(i, j)] * pfaffian(M, rest)
    return sp.expand(s)

def check_pfaffian(max_n, ell=4):
    t = sp.symbols('t')
    bad = 0; cnt = 0
    for n in range(0, max_n + 1):
        for lam in partitions(n):
            if len(lam) > ell: continue
            lamp = list(lam) + [0] * (ell - len(lam))
            alphas = [lamp[i] + ell - 1 - i for i in range(ell)]
            coeff = pf_entries(alphas, t)
            M = {(i, j): coeff.get((alphas[i], alphas[j]), 0) for i in range(ell) for j in range(ell)}
            pf = pfaffian(M, list(range(ell)))
            dv = d_vector(lam)
            Plam = sum(sp.Rational(dj, factorial(j)) * t ** j for j, dj in enumerate(dv))
            if sp.expand(pf - Plam) != 0:
                bad += 1; print("Pfaffian FAILS", lam, sp.expand(pf), Plam)
            cnt += 1
    print(f"Pfaffian formula (ell={ell}) checked for {cnt} partitions with |lam| <= {max_n}, <= {ell} rows:", "OK" if bad == 0 else f"{bad} failures")

if __name__ == "__main__":
    check_C2(16)
    check_pfaffian(7, ell=2)
    check_pfaffian(7, ell=4)
