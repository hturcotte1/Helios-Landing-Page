"""Support for the class-C theorem (two-row, two-column, hooks):
 (1) among partitions with r=2, those with d_2 <= 4 are exactly the six families with A_1+B_1 >= 2;
 (2) among those, d_3 <= 8 holds exactly for C-members and the exceptions (3,2,2),(3,3,1),(3,3,2);
 (3) explicit formulas for two-row and hook d-vectors agree with d_vector."""
import os, sys
from math import comb
sys.path.insert(0, os.path.dirname(__file__))
from young import partitions, conjugate, corner_runs
from census import d_vector

def sgn(x): return (x > 0) - (x < 0)

def d_tworow(l1, l2, j):
    return sum(comb(j, b) * sgn(l1 - l2 + 1 - j + 2 * b) for b in range(0, min(j, l2) + 1) if j - b <= l1 + 1)

def d_hook(a, b, j):  # hook (a, 1^b), n = a+b
    n = a + b
    if j == n: return comb(n - 1, b)
    return sum(comb(j, i) for i in range(0, j + 1) if i <= a - 1 and j - i <= b)

def in_C(lam):
    return len(lam) <= 2 or lam[0] <= 2 or (len(lam) >= 2 and lam[1] <= 1)

bad = 0
for n in range(2, 41):
    for lam in partitions(n):
        dv = d_vector(lam); runs = corner_runs(lam)
        if len(lam) <= 2:
            l1, l2 = (lam + (0,))[:2]
            if any(d_tworow(l1, l2, j) != dv[j] for j in range(n + 1)): bad += 1; print("two-row formula fails", lam)
        if len(lam) >= 1 and all(x == 1 for x in lam[1:]) and lam[0] >= 2 and len(lam) >= 2:
            a, b = lam[0], len(lam) - 1
            if any(d_hook(a, b, j) != dv[j] for j in range(n + 1)): bad += 1; print("hook formula fails", lam)
        if dv[1] == 2:
            A1 = sum(1 for a, b in runs if a == 1); B1 = sum(1 for a, b in runs if b == 1)
            assert dv[2] == 6 - A1 - B1
            if dv[2] <= 4:
                assert A1 + B1 >= 2
                small = dv[3] <= 8
                exc = lam in [(3,2,2),(3,3,1),(3,3,2)]
                if small != (in_C(lam) or exc):
                    bad += 1; print("d_3 criterion fails", lam, dv[:4], in_C(lam))
            if in_C(lam):
                assert dv[2] <= 4 and dv[3] <= 8, (lam, dv[:4])
print("class-C support checks n<=40:", "OK" if bad == 0 else f"{bad} failures")
