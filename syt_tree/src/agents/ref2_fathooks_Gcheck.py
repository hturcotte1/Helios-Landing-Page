"""Referee #2: the external fact used in Lemma 1's last line, G_{x,y} := sum_{rho in y x x box} f^rho t^|rho|/|rho|! = P_{x^y}
(Briefing fact 3 / Lemma 4.1).  Independent code, exact integers.  Rectangles with xy <= 40."""
from functools import lru_cache
from math import factorial
def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for p in range(min(n, maxp), 0, -1):
        for rest in partitions(n - p, p): yield (p,) + rest
def conj(lam): return tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()
def f_hook(lam):
    n = sum(lam)
    if n == 0: return 1
    lt = conj(lam); h = 1
    for i, r in enumerate(lam):
        for j in range(r): h *= (r - j) + (lt[j] - i) - 1
    return factorial(n) // h
def corners(lam): return [i for i in range(len(lam)) if i == len(lam) - 1 or lam[i] > lam[i + 1]]
def remove(lam, i):
    l = list(lam); l[i] -= 1
    if l[i] == 0: l.pop()
    return tuple(l)
@lru_cache(None)
def dvec(lam):
    n = sum(lam)
    if n == 0: return (1,)
    d = [0] * (n + 1); d[0] = 1
    for i in corners(lam):
        for j, v in enumerate(dvec(remove(lam, i))): d[j + 1] += v
    return tuple(d)
bad = 0; tot = 0
for x in range(1, 41):
    for y in range(1, 41):
        if x * y > 40: continue
        dv = dvec(tuple([x] * y)); tot += 1
        for j in range(x * y + 1):
            s = sum(f_hook(r) for r in partitions(j) if (not r or r[0] <= x) and len(r) <= y)
            if s != dv[j]: bad += 1; print("G FAIL", x, y, j)
print("G_{x,y} = P_{x^y} check, rectangles xy<=40:", tot, "rectangles, fails =", bad)
