"""Referee #2: Theorem 3 / Corollary 3.1 on 2-corner shapes with a+b+c+d <= S (large n, small perimeter), plus
the n = 5 remark (lambda(1,1,1,3) vs lambda(1,1,2,1) share (n, x1, x2, N_k)).  Reuses ref2_fathooks_thm3 primitives."""
import sys
sys.argv = [sys.argv[0]] + ['0']
import ref2_fathooks_thm3 as R
from fractions import Fraction
from math import factorial

S = 16
cnt = 0; fails = 0; maxn = 0; bal = 0
for a in range(1, S):
    for b in range(1, S):
        for c in range(1, S):
            for d in range(1, S):
                if a + b + c + d > S: continue
                lam = R.fathook(a, b, c, d); n = sum(lam)
                if n > 70: continue
                cnt += 1; maxn = max(maxn, n)
                dv = R.dvec(lam)
                xs = sorted([a, b, c, d]); x1, x2 = xs[0], xs[1]
                P = R.egf(list(dv)) + [Fraction(0)] * (R.DEG + 1 - len(dv))
                PI = R.mul(P, R.Iinv)
                q = [(R.I[j] - PI[j]) * factorial(j) for j in range(R.DEG + 1)]
                Nk = {k: xs.count(k) for k in set(xs)}
                for j in range(0, x1 + x2 + 1):
                    if q[j] != sum(Nk[k] * R.Dcoef(k, j) for k in Nk): fails += 1; print("3.1 FAIL", (a, b, c, d), j)
                r1, r2, N, mr = R.read_local(dv)
                truth = {k: xs.count(k) for k in range(x1, x1 + x2)}
                if (r1, r2) != (x1, x2) or N != truth or mr > x1 + x2: fails += 1; print("READ FAIL", (a, b, c, d))
                isbal = xs[3] <= x1 + x2 - 1
                if isbal != (sum(N.values()) == 4): fails += 1; print("BAL FAIL", (a, b, c, d))
                if isbal: bal += 1
print(f"perimeter <= {S}, n <= 70: {cnt} shapes (max n {maxn}), balanced {bal}, fails = {fails}")
# n = 5 remark
for p in [(1, 1, 1, 3), (1, 1, 2, 1)]:
    lam = R.fathook(*p); print(p, lam, R.read_local(R.dvec(lam))[:3])
