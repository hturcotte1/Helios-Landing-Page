"""Fast, sympy-free version of structure_verify.py (3)-(4): recovery of <=3-row shapes from (n,C1^2,C2), n<=300,
and of <=4-row shapes from (n,C1^2,C2,C4), n<=150, by exact integer arithmetic; integer roots are found by
trial over the finite range [-L, n] (all x_i lie there)."""
import os, sys
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from young import boxes
def Cpow(lam, k): return sum((j-i)**k for (i, j) in boxes(lam))
def isq(v):
    if v < 0: return None
    r = int(v**0.5)
    while r*r > v: r -= 1
    while (r+1)*(r+1) <= v: r += 1
    return r if r*r == v else None
def int_roots(coefs, lo, hi):
    """monic integer poly (high->low). Return sorted-desc list of integer roots with multiplicity if it splits over
    integers in [lo,hi], else None."""
    deg = len(coefs)-1
    roots = []
    p = list(coefs)
    for r in range(lo, hi+1):
        while True:
            # evaluate
            v = 0
            for c in p: v = v*r + c
            if v != 0: break
            # divide by (T - r)
            q = [p[0]]
            for c in p[1:-1]: q.append(c + q[-1]*r)
            p = q; roots.append(r)
            if len(p) == 1: break
        if len(p) == 1: break
    return sorted(roots, reverse=True) if len(roots) == deg else None
def rec3(n, C1sq, C2):
    out = []; p1 = n-6; s = isq(C1sq)
    if s is None: return out
    for C1 in sorted({s, -s}):
        p2 = 2*(C1+4) - p1
        p3 = Fraction(6*(C2-6) - 3*p2 - p1, 2)
        if p3.denominator != 1: continue
        p3 = int(p3); e1 = p1; e2 = Fraction(e1*p1-p2, 2); e3 = Fraction(e2*p1 - e1*p2 + p3, 3)
        if e2.denominator != 1 or e3.denominator != 1: continue
        x = int_roots([1, -e1, int(e2), -int(e3)], -3, n)
        if x and x[0] > x[1] > x[2] >= -3:
            lam = tuple(v for v in (x[0]+1, x[1]+2, x[2]+3) if v > 0)
            if sum(lam) == n: out.append(lam)
    return out
def rec4(n, C1sq, C2, C4):
    out = []; p1 = n-10; s = isq(C1sq)
    if s is None: return out
    for C1 in sorted({s, -s}):
        p2 = 2*(C1+10) - p1
        p3 = Fraction(6*(C2-20) - 3*p2 - p1, 2)
        if p3.denominator != 1: continue
        p3 = int(p3); e1 = p1; e2 = Fraction(e1*p1-p2, 2); e3 = Fraction(e2*p1 - e1*p2 + p3, 3)
        if e2.denominator != 1 or e3.denominator != 1: continue
        e2 = int(e2); e3 = int(e3)
        A = e1*p3 - e2*p2 + e3*p1; B = e1*A - e2*p3 + e3*p2
        if e1 + 2 == 0: continue
        e4 = Fraction(30*(C4-116) - (6*B + 15*A + 10*p3 - p1), -30*(e1+2))
        if e4.denominator != 1: continue
        x = int_roots([1, -e1, e2, -e3, int(e4)], -4, n)
        if x and x[0] > x[1] > x[2] > x[3] >= -4:
            lam = tuple(v for v in (x[0]+1, x[1]+2, x[2]+3, x[3]+4) if v > 0)
            if sum(lam) == n: out.append(lam)
    return out
def shapes(n, L):
    def go(rem, maxp, k):
        if k == 0:
            if rem == 0: yield ()
            return
        for a in range(min(rem, maxp), -1, -1):
            for rest in go(rem-a, a, k-1): yield (a,)+rest
    for s in go(n, n, L):
        yield tuple(v for v in s if v > 0)
if __name__ == "__main__":
    bad = []; cnt = 0; neg = []
    for n in range(1, 301):
        for lam in shapes(n, 3):
            cnt += 1; C1 = Cpow(lam, 1)
            if C1 <= 0: neg.append(lam)
            if rec3(n, C1*C1, Cpow(lam, 2)) != [lam]: bad.append(lam)
    print(f"3 rows, n<=300: {cnt} shapes, not uniquely recovered: {bad}")
    print(f"3 rows: shapes with C1<=0: {neg}")
    bad = []; cnt = 0; neg = []
    for n in range(1, 151):
        for lam in shapes(n, 4):
            cnt += 1; C1 = Cpow(lam, 1)
            if C1 <= 0: neg.append(lam)
            r = rec4(n, C1*C1, Cpow(lam, 2), Cpow(lam, 4))
            if r != [lam]: bad.append((lam, r))
    print(f"4 rows, n<=150: {cnt} shapes, not uniquely recovered ({len(bad)}): {bad}")
    print(f"4 rows: {len(neg)} shapes with C1<=0, max n = {max(sum(l) for l in neg)}: {neg}")
