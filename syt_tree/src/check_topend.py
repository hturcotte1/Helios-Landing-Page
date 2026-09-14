"""Independent re-verification (own code only) of the key claims of agent_notes/topend.md:
 (A) omega_lam(K_3) = C2 - C(n,2);  omega(K_22) = C1^2/2 - 3C2/2 + C(n,2);
     omega(K_5) = C4 - (3n-10)C2 - 2C1^2 + n(n-1)(5n-19)/6,   omega(K_rho) = |K_rho| chi(rho)/f,  via Murnaghan-Nakayama;
 (C) level-6 redundancy: (n)_6 d_{n-6}/f = G6(n,C1^2,C2,C4) with the explicit G6, and functional dependence
     (n, d_n, d_{n-3}, d_{n-4}, d_{n-5}) -> d_{n-6} on all partitions of n <= 30;
 (S) (n, f, C2, C1^2, C4) separates transpose classes for n <= 48 and fails at n = 49 (pair given in the report);
 (M) the n=50 mirror-box-move pair has equal f, n, C1^2, C2, C4 (all even moments) and different d-vectors."""
import os, sys
from fractions import Fraction
from math import comb
sys.path.insert(0, os.path.dirname(__file__))
from young import partitions, conjugate, f_hook, boxes, f_skew
from census import d_vector
from check_characters import chi_lam

def C(lam, k): return sum((j - i) ** k for (i, j) in boxes(lam))

# (A)
bad = 0
for n in range(3, 13):
    for lam in partitions(n):
        f = f_hook(lam); c1, c2, c4 = C(lam, 1), C(lam, 2), C(lam, 4)
        w3 = Fraction(2 * comb(n, 3) * chi_lam(lam, (3,)), f)
        if w3 != c2 - comb(n, 2): bad += 1; print("A1 fails", lam)
        if n >= 4:
            w22 = Fraction(3 * comb(n, 4) * chi_lam(lam, (2, 2)), f)
            if w22 != Fraction(c1 * c1, 2) - Fraction(3 * c2, 2) + comb(n, 2): bad += 1; print("A2 fails", lam)
        if n >= 5:
            w5 = Fraction(24 * comb(n, 5) * chi_lam(lam, (5,)), f)
            if w5 != c4 - (3 * n - 10) * c2 - 2 * c1 * c1 + Fraction(n * (n - 1) * (5 * n - 19), 6): bad += 1; print("A3 fails", lam)
print("(A) Theorem A formulas via own MN code, 3<=n<=12:", "OK" if bad == 0 else f"{bad} failures")

# (C)
bad = 0
for n in range(6, 23):
    for lam in partitions(n):
        dv = d_vector(lam); f = dv[n]; c1, c2, c4 = C(lam, 1), C(lam, 2), C(lam, 4)
        G6 = (n * n - 11 * n + 42) * c1 * c1 + 2 * c2 * c2 + Fraction(2 * n ** 3, 3) * c2 - 16 * n * n * c2 + Fraction(328 * n, 3) * c2 - 210 * c2 + (n - 15) * c4
        lhs = Fraction(n * (n - 1) * (n - 2) * (n - 3) * (n - 4) * (n - 5) * dv[n - 6], f)
        if lhs != G6: bad += 1; print("G6 fails", lam, lhs, G6)
print("(C) explicit G6 formula, 6<=n<=22:", "OK" if bad == 0 else f"{bad} failures")
dep_bad = 0
for n in range(6, 31):
    seen = {}
    for lam in partitions(n):
        dv = d_vector(lam)
        key = (dv[n], dv[n - 3], dv[n - 4], dv[n - 5])
        if key in seen and seen[key] != dv[n - 6]: dep_bad += 1
        seen[key] = dv[n - 6]
print("(C) functional dependence d_{n-6} <- (n,d_n,d_{n-3},d_{n-4},d_{n-5}), n<=30:", "OK" if dep_bad == 0 else f"{dep_bad} violations")

# (S)
first_fail = None
for n in range(1, 51):
    groups = {}
    for lam in partitions(n):
        key = min(lam, conjugate(lam))
        k = (f_hook(lam), C(lam, 2), C(lam, 1) ** 2, C(lam, 4))
        groups.setdefault(k, set()).add(key)
    coll = [g for g in groups.values() if len(g) > 1]
    if coll and first_fail is None:
        first_fail = (n, coll); break
print("(S) (n,f,C2,C1^2,C4) first collision at n =", first_fail[0] if first_fail else None, ":", sorted(first_fail[1][0]) if first_fail else "")

# (M)
lam = (12, 10, 8, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1); mu = (12, 11, 8, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1)
assert sum(lam) == sum(mu) == 50
print("(M) f equal:", f_hook(lam) == f_hook(mu), "| even moments equal:", all(C(lam, k) == C(mu, k) for k in (2, 4, 6, 8)),
      "| C1:", C(lam, 1), C(mu, 1), "| C1*C3:", C(lam, 1) * C(lam, 3), C(mu, 1) * C(mu, 3))
dl, dm = d_vector(lam), d_vector(mu)
diff = [j for j in range(51) if dl[j] != dm[j]]
print("(M) d-vectors differ at indices (first/last):", diff[:3], diff[-3:], "| agree at n-3..n:", all(dl[50 - i] == dm[50 - i] for i in (0, 3, 4, 5, 6)), "| differ at n-7:", dl[43] != dm[43])
