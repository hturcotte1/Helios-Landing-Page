"""Referee #1: run the author's rec3 on all <=3-row shapes n<=300 (as claimed), and cross-check rec3 with
invariants derived from d-vectors (my own recursion) for n<=16."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from structure_recover_fast import rec3, shapes, Cpow
from ref1_structure_check import dvec, invariants_from_d, shapes3, conjugate
bad = []; cnt = 0; neg = []
for n in range(1, 301):
    for lam in shapes(n, 3):
        cnt += 1; C1 = Cpow(lam, 1)
        if C1 <= 0: neg.append(lam)
        if rec3(n, C1*C1, Cpow(lam, 2)) != [lam]: bad.append(lam)
print("author rec3, n<=300:", cnt, "shapes; not unique:", bad)
print("C1<=0:", neg)
# cross-check with d-vector-derived invariants
for n in range(4, 17):
    for lam in shapes3(n):
        C1sq, C2 = invariants_from_d(dvec(lam))
        assert C1sq.denominator == 1 and C2.denominator == 1
        r = rec3(n, int(C1sq), int(C2))
        if n >= 8: assert r == [lam], (lam, r)
        else: assert sorted(r) == sorted({lam, conjugate(lam)}), (lam, r)
print("rec3 fed with d-vector-derived (C1^2,C2): unique for 8<=n<=16, transpose pairs for n<=7: OK")
