"""Independent checks of claims in agent_notes/fathooks.md and agent_notes/structure.md."""
import os, sys, itertools
from math import comb, factorial
from fractions import Fraction
sys.path.insert(0, os.path.dirname(__file__))
from young import partitions, involutions, boxes, f_hook, conjugate
from census import d_vector
from check_characters import chi
I = involutions
# --- fathooks closed forms (S1)-(S3) and Theorem F ingredients
bad = 0
for c in range(2, 20):
    for d in range(2, 20):
        lam = (c + 1,) + (c,) * d; n = sum(lam)
        if n > 44: continue
        dv = d_vector(lam); x = min(c, d)
        for j in range(0, x + 1):
            if dv[j] != I(j + 1): bad += 1; print("F1 window fails", lam, j)
        if dv[x + 1] != I(x + 2) - ((c == x) + (d == x)): bad += 1; print("F1 eps fails", lam)
        if c == x < d and dv[x + 2] != I(x + 3) - 2 * (x + 2) - (d == x + 1): bad += 1; print("F1 (c=x<d) d_{x+2} fails", lam)
        if d == x < c and dv[x + 2] != I(x + 3) - (x + 3) - (c == x + 1): bad += 1; print("F1 (d=x<c) d_{x+2} fails", lam)
for b in range(2, 20):
    for c in range(b, 20):
        lam = (c + 1,) * b + (c,); n = sum(lam)
        if n > 44: continue
        dv = d_vector(lam); x = b
        rect = d_vector((c + 1,) * (b + 1))
        if any(dv[j] != rect[j + 1] for j in range(n + 1)): bad += 1; print("(S1) fails", lam)
        for j in range(0, x + 1):
            if dv[j] != I(j + 1): bad += 1; print("F2 window fails", lam, j)
        if dv[x + 1] != I(x + 2) - 1 - (b == c): bad += 1; print("F2 eps fails", lam)
        if b < c and dv[x + 2] != I(x + 3) - (x + 3) - (c == x + 1): bad += 1; print("F2 d_{x+2} fails", lam)
print("fathooks closed forms (S1)-(S3), all F1/F2 shapes with parameters>=2 and n<=44:", "OK" if bad == 0 else f"{bad} failures")
# Theorem F itself against all partitions, n <= 26 (cheap replication)
bad = 0
for n in range(6, 27):
    groups = {}
    for lam in partitions(n): groups.setdefault(d_vector(lam), []).append(lam)
    for lam in partitions(n):
        isF = (len(lam) >= 3 and lam[0] == lam[1] + 1 and all(x == lam[1] for x in lam[1:])) or (len(lam) >= 2 and all(x == lam[0] for x in lam[:-1]) and lam[-1] == lam[0] - 1)
        if isF:
            g = groups[d_vector(lam)]
            if any(mu not in (lam, conjugate(lam)) for mu in g): bad += 1; print("Theorem F fails", lam, g)
print("Theorem F (F1/F2 shapes vs all partitions) n<=26:", "OK" if bad == 0 else f"{bad} failures")
# --- structure: Theorem S3  n! P_lam(t) = sum_tau chi^lam(tau^2) (1+t)^{fix tau}
def cycle_type(p):
    n = len(p); seen = [False] * n; ct = []
    for i in range(n):
        if not seen[i]:
            l = 0; j = i
            while not seen[j]: seen[j] = True; j = p[j]; l += 1
            ct.append(l)
    return tuple(sorted(ct, reverse=True))
bad = 0
for n in range(1, 8):
    perms = list(itertools.permutations(range(n)))
    for lam in partitions(n):
        dv = d_vector(lam)
        rhs = [0] * (n + 1)
        for p in perms:
            sq = tuple(p[p[i]] for i in range(n))
            fix = sum(1 for i in range(n) if sq[i] == i) if False else sum(1 for i in range(n) if p[i] == i)
            ch = chi(lam, cycle_type(sq))
            for j in range(n + 1): rhs[j] += ch * comb(fix, j)
        for j in range(n + 1):
            if Fraction(rhs[j] * factorial(j), factorial(n)) != dv[j]: bad += 1; print("S3 fails", lam, j); break
print("Theorem S3 (character sum over square roots) n<=7:", "OK" if bad == 0 else f"{bad} failures")
# --- sign lemma: C_1 > 0 for <=3 rows, n>=10; <=4 rows, n>=21
def C1(l): return sum(j - i for (i, j) in boxes(l))
worst3 = max(n for n in range(1, 41) for l in partitions(n) if len(l) <= 3 and C1(l) <= 0)
worst4 = max(n for n in range(1, 41) for l in partitions(n) if len(l) <= 4 and C1(l) <= 0)
print("largest n with a <=3-row shape having C_1<=0:", worst3, "| <=4 rows:", worst4)
