"""Exploration for Conjecture B:
 (a) verify the 'local product' theorem  d_j(lam) = j! [t^j] prod_i P_{a_i^{b_i}}(t)  for j <= 2m,
     m = min_i min(a_i,b_i), and the deficit formula at j = m+1;
 (b) for each n, the largest 'first differing index' between non-transpose partitions of n,
     i.e. how deep into the d-vector one must look to separate shapes.
"""
import os, sys
from math import factorial
sys.path.insert(0, os.path.dirname(__file__))
from young import partitions, conjugate, corner_runs, involutions
from census import d_vector
from fractions import Fraction

def poly_mul(p, q, deg):
    r = [Fraction(0)] * (deg + 1)
    for i, a in enumerate(p[:deg + 1]):
        if a == 0: continue
        for j, b in enumerate(q[:deg + 1 - i]):
            r[i + j] += a * b
    return r

def egf(dv, deg):
    return [Fraction(dv[j], factorial(j)) if j < len(dv) else Fraction(0) for j in range(deg + 1)]

def check_local_product(max_n):
    bad = 0
    for n in range(1, max_n + 1):
        for lam in partitions(n):
            runs = corner_runs(lam)
            m = min(min(a, b) for a, b in runs)
            r = len(runs)
            deg = 2 * m
            prod = [Fraction(1)] + [Fraction(0)] * deg
            for (a, b) in runs:
                rect = tuple([a] * b)
                prod = poly_mul(prod, egf(d_vector(rect), deg), deg)
            dv = d_vector(lam)
            for j in range(0, min(deg, n) + 1):
                if Fraction(dv[j], factorial(j)) != prod[j]:
                    bad += 1; print("local product FAILS", lam, j)
            # deficit at m+1 relative to exp(r(t+t^2/2)):  A_m + B_m
            if m + 1 <= n:
                # coefficient of t^{m+1}/(m+1)! in exp(r(t+t^2/2)) = sum over compositions... compute via r-fold convolution
                conv = [Fraction(1)] + [Fraction(0)] * (m + 1)
                inv = [Fraction(involutions(k), factorial(k)) for k in range(m + 2)]
                for _ in range(r):
                    conv = poly_mul(conv, inv, m + 1)
                full = conv[m + 1] * factorial(m + 1)
                deficit = full - dv[m + 1]
                AB = sum(1 for a, b in runs if a == m) + sum(1 for a, b in runs if b == m)
                if deficit != AB:
                    bad += 1; print("deficit formula FAILS", lam, deficit, AB)
                for j in range(0, m + 1):
                    # d_j equals the r-fold convolution for j <= m
                    convj = [Fraction(1)] + [Fraction(0)] * j
                    for _ in range(r):
                        convj = poly_mul(convj, inv, j)
                    if convj[j] * factorial(j) != dv[j]:
                        bad += 1; print("d_j = conv fails", lam, j)
    print(f"local product theorem checked for all |lam| <= {max_n}: {'OK' if bad == 0 else 'FAILURES: %d' % bad}")

def first_diff(dv1, dv2):
    for j in range(len(dv1)):
        if dv1[j] != dv2[j]: return j
    return None

def separation(max_n):
    print("n  max_first_differing_index  witnesses (up to 3)")
    for n in range(1, max_n + 1):
        reps = {}
        for lam in partitions(n):
            key = min(lam, conjugate(lam))
            reps[key] = d_vector(lam)
        keys = list(reps)
        worst = -1; wit = []
        for i in range(len(keys)):
            for k in range(i + 1, len(keys)):
                j = first_diff(reps[keys[i]], reps[keys[k]])
                assert j is not None, (keys[i], keys[k])
                if j > worst: worst, wit = j, [(keys[i], keys[k])]
                elif j == worst and len(wit) < 3: wit.append((keys[i], keys[k]))
        print(n, worst, wit)

if __name__ == "__main__":
    check_local_product(int(sys.argv[1]) if len(sys.argv) > 1 else 18)
    separation(int(sys.argv[2]) if len(sys.argv) > 2 else 24)
