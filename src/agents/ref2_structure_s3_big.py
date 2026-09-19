"""Referee #2: (a) cross-check own MN chi vs project's check_characters.chi (n<=10);
(b) Theorem S3 + derangement transform on selected large shapes (n = 24, 30) via class sums,
    comparing against the project's d_vector as well as own dvec."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))
from math import factorial, comb
from fractions import Fraction
from ref2_structure_s3 import chi, dvec, parts, z, sq_type, rhs_class, delta_class, f_dim, D
from check_characters import chi as chi_proj
from census import d_vector

bad = 0
for n in range(0, 11):
    for lam in parts(n):
        for rho in parts(n):
            if chi(lam, rho) != chi_proj(lam, rho): bad += 1; print("chi mismatch", lam, rho)
        if dvec(lam) != d_vector(lam): bad += 1; print("dvec mismatch", lam)
print("cross-check own chi/dvec vs project, n<=10: mismatches =", bad)

shapes = [(24,), (23,1), (12,12), (8,8,8), (6,6,6,6), (20,3,1), (13,1,1,1,1,1,1,1,1,1,1,1), (9,7,5,2,1),
          (30,), (29,1), (15,15), (10,10,10), (6,6,6,6,6), (25,3,2), (12,8,5,3,2), (16,7,4,2,1)]
bad = 0
for lam in shapes:
    n = sum(lam)
    d = dvec(lam); assert d == d_vector(lam)
    lhs = [factorial(n) * d[j] // factorial(j) for j in range(n + 1)]
    r = rhs_class(lam)
    ok1 = (lhs == r)
    delta = [delta_class(lam, m) for m in range(n + 1)]
    ok2 = all(factorial(i) * d[n - i] == sum(comb(i, m) * delta[m] for m in range(i + 1)) for i in range(n + 1))
    f = f_dim(lam)
    ok3 = (delta[0] == f and delta[1] == 0 and delta[2] == f and delta[3] == 2 * chi(lam, (3,) + (1,) * (n - 3))
           and delta[4] == 6 * chi(lam, (2, 2) + (1,) * (n - 4)) + 3 * f)
    val = sum((-1) ** j * lhs[j] for j in range(n + 1))
    ok4 = (val == delta[n]) and (lam != (n,) or val == D(n)) and (lam != (n-1,1) or val == (-1)**n*(n-1))
    print(lam, "S3:", ok1, "transform:", ok2, "small deltas:", ok3, "derangement value:", ok4, "n!P(-1) =", val, flush=True)
    if not (ok1 and ok2 and ok3 and ok4): bad += 1
print("large-shape mismatches =", bad)
