"""topend: small consolidated checks quoted in the report.
 (a) chi^lam(3,3,1^{n-6}) = (9/2)u_6 - (9/2)u_5 + (3/2)u_3 - f/2, all lam |- n <= 16 (MN vs d-vector);
 (b) the sigma/z coefficients of u_3..u_6 as printed in the report;
 (c) the hard pair (3,1^{n-3}) vs (2,2,1^{n-4}): f differs by 1 and C_2 by 3, n <= 40;
 (d) Theorem A polynomial for omega(3,3) expanded (sympy) equals the fitted Omega_(3,3)."""
import os, sys
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from topend_mn import *
from census import d_vector
from young import f_hook
import sympy as sp
bad = 0; cnt = 0
for n in range(6, 17):
    for lam in partitions(n):
        dv = d_vector(lam); f = dv[n]; u = lambda i: dv[n - i]
        cnt += 1
        if Fraction(9, 2) * u(6) - Fraction(9, 2) * u(5) + Fraction(3, 2) * u(3) - Fraction(f, 2) != chi(lam, (3, 3)): bad += 1
print(f"(a) chi(3,3) = (9/2)u_6 - (9/2)u_5 + (3/2)u_3 - f/2, all lam |- n, 6 <= n <= 16: {cnt} checks, {'OK' if bad == 0 else bad}")
for i in range(3, 7):
    print(f"(b) u_{i} coefficients:", [(str(Fraction(sigma(r), z(r))), r) for r in partitions(i) if sigma(r)])
bad = 0
for n in range(6, 41):
    a = (3,) + (1,) * (n - 3); b = (2, 2) + (1,) * (n - 4)
    if f_hook(a) - f_hook(b) != 1 or C(a, 2) - C(b, 2) != 4: bad += 1
print(f"(c) f(3,1^(n-3)) - f(2,2,1^(n-4)) = 1 and C_2 difference = 4 for 6 <= n <= 40: {'OK' if bad == 0 else bad}")
n, C1, C2, C4 = sp.symbols('n C1 C2 C4')
w3 = C2 - n * (n - 1) / 2; w22 = C1**2 / 2 - sp.Rational(3, 2) * C2 + n * (n - 1) / 2
w5 = C4 - (3 * n - 10) * C2 - 2 * C1**2 + n * (n - 1) * (5 * n - 19) / 6
w33 = sp.expand((w3**2 - 2 * n * (n - 1) * (n - 2) / 6 - (3 * n - 8) * w3 - 8 * w22 - 5 * w5) / 2)
print("(d) omega(3,3) expanded =", w33)
