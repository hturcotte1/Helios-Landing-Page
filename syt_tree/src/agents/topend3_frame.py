"""Single-box mirror moves inside the frame  lam = (12, 10, x_3, ..., x_10, 2, 1, 1)  (rows 3..10 free in [2,10]),
mu = lam - box(11,2) + box(2,11)  (content -9 -> +9).  Which inner shapes give f^lam = f^mu ?"""
from fractions import Fraction
from itertools import combinations_with_replacement
from topend3_lib import *
from young import hook_lengths
def hook_prod(lam):
    p = 1
    for h in hook_lengths(lam).values(): p *= h
    return p
good = []; total = 0
for inner in combinations_with_replacement(range(2, 11), 8):
    x = tuple(sorted(inner, reverse=True))
    lam = (12, 10) + x + (2, 1, 1)
    mu = (12, 11) + x + (2, 1, 1, 1)
    mu = tuple(list(mu[:10]) + [2, 1, 1, 1])  # rows: 12, 11, x..., 2 (row 11)?  careful
    # lam rows: 1:12, 2:10, 3..10: x, 11:2, 12:1, 13:1.  mu: row2 -> 11, row 11 -> 1: (12,11,x...,1,1,1)
    mu = (12, 11) + x + (1, 1, 1)
    total += 1
    if hook_prod(lam) == hook_prod(mu):
        good.append(x)
print(f"inner shapes tested: {total}; with f^lam == f^mu: {len(good)}")
for x in good: print(x, "n =", 12 + 10 + sum(x) + 4, " C_1 =", C((12, 10) + x + (2, 1, 1), 1))
