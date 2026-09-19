"""Independent re-verification (third pass) of the character-theoretic core.
(a) MN sanity: orthogonality, transpose rule, chi(1^n)=f, n<=8/10.
(b) sigma(rho) closed formula vs brute force, |rho|<=8.
(c) Theorem B: u_i = sum_rho sigma/z chi^lam(rho cup 1^{n-i}) vs d-vector, i<=10, n<=13.
(d) Theorem A (A1)-(A4) via omega from MN, all lam |- n <= 16.
(e) Corollary B': C_2, C_1^2, C_4 recovered from (n,f,u_3,u_4,u_5) using only d-vectors, n<=20.
(f) Aitken determinant vs f_skew, all lam |- n<=9, all nu.
"""
import sys, time
from fractions import Fraction
from math import comb, factorial
from topend3_lib import *

t0 = time.time()
# (a)
for n in range(1, 9):
    P = list(partitions(n))
    for lam in P:
        assert chi(lam, (1,) * n) == f_hook(lam), (lam,)
        lt = conjugate(lam)
        for rho in P:
            sgn = (-1) ** (n - len(rho))
            assert chi(lt, rho) == sgn * chi(lam, rho)
    # row orthogonality  sum_rho chi^lam(rho) chi^mu(rho) / z_rho = delta
    for lam in P:
        for mu in P:
            s = sum(Fraction(chi(lam, rho) * chi(mu, rho), z(rho)) for rho in P)
            assert s == (1 if lam == mu else 0)
for n in range(9, 11):
    for lam in partitions(n):
        assert chi(lam, (1,) * n) == f_hook(lam)
print(f"(a) MN: orthogonality+transpose n<=8, chi(1^n)=f n<=10: OK [{time.time()-t0:.0f}s]")

# (b)
for i in range(1, 9):
    for rho in partitions(i):
        assert sigma_formula(rho) == sigma_brute(rho), rho
    # and sigma = sum_nu chi^nu(rho)  (Frobenius-Schur)
    for rho in partitions(i):
        assert sigma_formula(rho) == sum(chi(nu, rho) for nu in partitions(i)), rho
print(f"(b) sigma formula == brute force == sum_nu chi^nu, |rho|<=8: OK [{time.time()-t0:.0f}s]")

# (c)
cnt = 0
for n in range(1, 14):
    for lam in partitions(n):
        dv = d_vector(lam)
        for i in range(0, min(n, 10) + 1):
            assert u_from_chars(lam, i) == dv[n - i], (lam, i)
            cnt += 1
print(f"(c) Theorem B u_i (chars) == d_{{n-i}}, i<=10, n<=13: {cnt} checks OK [{time.time()-t0:.0f}s]")

# (d)
cnt = 0
for n in range(1, 17):
    for lam in partitions(n):
        c1, c2, c4 = C(lam, 1), C(lam, 2), C(lam, 4)
        A1 = c2 - comb(n, 2)
        A2 = Fraction(c1 * c1, 2) - Fraction(3 * c2, 2) + comb(n, 2)
        A3 = c4 - (3 * n - 10) * c2 - 2 * c1 * c1 + Fraction(n * (n - 1) * (5 * n - 19), 6)
        A4 = (Fraction(c2 * c2, 2) - Fraction(5 * c4, 2) + 3 * c1 * c1 - Fraction(n * n * c2, 2) + Fraction(13 * n * c2, 2)
              - 15 * c2 + Fraction(n ** 4, 8) - Fraction(7 * n ** 3, 4) + Fraction(47 * n * n, 8) - Fraction(17 * n, 4))
        assert omega(lam, (3,)) == A1, lam
        assert omega(lam, (2, 2)) == A2, lam
        assert omega(lam, (5,)) == A3, lam
        assert omega(lam, (3, 3)) == A4, lam
        cnt += 1
print(f"(d) Theorem A (A1)-(A4) via MN, all lam |- n<=16: {cnt} shapes OK [{time.time()-t0:.0f}s]")

# (e)
cnt = 0
for n in range(6, 21):
    for lam in partitions(n):
        dv = d_vector(lam); f = dv[n]; u3, u4, u5 = dv[n - 3], dv[n - 4], dv[n - 5]
        c2 = comb(n, 2) + Fraction(ffact(n, 3) * (3 * u3 - 2 * f), 3 * f)
        c1sq = 3 * c2 - n * (n - 1) + Fraction(6 * comb(n, 4) * (4 * u4 - 4 * u3 + f), f)
        c4 = Fraction(24 * comb(n, 5) * (5 * u5 - 5 * u4 + f), f) + (3 * n - 10) * c2 + 2 * c1sq - Fraction(n * (n - 1) * (5 * n - 19), 6)
        assert c2 == C(lam, 2) and c1sq == C(lam, 1) ** 2 and c4 == C(lam, 4), lam
        cnt += 1
print(f"(e) Corollary B' (C_2, C_1^2, C_4 from d-vector), all lam |- 6<=n<=20: {cnt} shapes OK [{time.time()-t0:.0f}s]")

# (f)
cnt = 0
for n in range(1, 10):
    for lam in partitions(n):
        for m in range(0, n + 1):
            for nu in partitions(m):
                assert f_skew_aitken(lam, nu) == f_skew(lam, nu), (lam, nu)
                cnt += 1
print(f"(f) Aitken determinant == f_skew, all lam |- n<=9, all nu: {cnt} checks OK [{time.time()-t0:.0f}s]")
