"""structure angle: exploratory checks.
1. Hermite/character expansion  P_lam(t) = sum_m w_m(lam) H_m(t+1),
   H_m(x) = sum_i C(m,2i)(2i-1)!! x^{m-2i}  (EGF e^{xu+u^2/2}),
   w_m(lam) = (1/m!) sum_{rho |- n-m, all parts >= 3} z_rho^{-1} chi^lam(sq(rho) cup 1^m).
2. n! P_lam(-1) table.
3. Three-row Pfaffian identity  P_(a,b,c) = P_(a+2,b+2)E_c - P_(a+2,c+1)E_{b+1} + P_(b+1,c+1)E_{a+2}.
4. Separation of 3-row shapes by polynomial invariants.
"""
import os, sys
from fractions import Fraction
from math import factorial, comb
from collections import Counter, defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from young import partitions, conjugate, f_hook, f_skew, boxes, involutions
from census import d_vector
from check_characters import chi

def Ppoly(lam):
    """coefficients of P_lam(t) as Fractions, index = power of t"""
    d = d_vector(tuple(lam))
    return [Fraction(dj, factorial(j)) for j, dj in enumerate(d)]

def H(m):
    """H_m(x) coefficients (integers), index = power of x"""
    c = [0]*(m+1)
    for i in range(0, m//2+1):
        c[m-2*i] = comb(m, 2*i) * (factorial(2*i)//(2**i*factorial(i)))
    return c

def shift_poly(c):
    """given coefficients in x, return coefficients in t where x = t+1"""
    out = [Fraction(0)]*len(c)
    for k, ck in enumerate(c):
        for i in range(k+1):
            out[i] += ck*comb(k, i)
    return out

def zee(rho):
    z = 1
    for k, m in Counter(rho).items():
        z *= k**m * factorial(m)
    return z

def sq(rho):
    out = []
    for k in rho:
        if k % 2: out.append(k)
        else: out += [k//2, k//2]
    return tuple(sorted(out, reverse=True))

def w(lam, m):
    n = sum(lam)
    tot = Fraction(0)
    for rho in partitions(n-m):
        if rho and rho[-1] < 3: continue
        rho2 = tuple(sorted(list(sq(rho)) + [1]*m, reverse=True))
        # chi^lam at cycle type rho2 (full type of S_n) via MN
        tot += Fraction(chi(tuple(lam), rho2), zee(rho))
    return tot/factorial(m)

def check_hermite(N):
    bad = 0
    for n in range(0, N+1):
        for lam in partitions(n):
            lam = tuple(lam)
            P = Ppoly(lam)
            Q = [Fraction(0)]*(n+1)
            for m in range(n+1):
                wm = w(lam, m)
                if wm == 0: continue
                hc = shift_poly([Fraction(c) for c in H(m)])
                for i, c in enumerate(hc):
                    Q[i] += wm*c
            if P != Q:
                bad += 1
                print("MISMATCH", lam, P, Q)
    print(f"Hermite/character expansion checked for all |lam| <= {N}: mismatches = {bad}")

def check_Pminus1(N):
    print("n! P_lam(-1):")
    for n in range(0, N+1):
        row = []
        for lam in partitions(n):
            d = d_vector(tuple(lam))
            v = sum((-1)**j * d[j] * factorial(n)//factorial(j) for j in range(n+1))
            row.append((lam, v))
        print(n, row)

def Epoly(k):
    return [Fraction(1, factorial(i)) for i in range(k+1)]

def pmul(a, b):
    out = [Fraction(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out

def padd(*ps):
    L = max(len(p) for p in ps)
    out = [Fraction(0)]*L
    for p in ps:
        for i, x in enumerate(p): out[i] += x
    return out

def strip(p):
    while len(p) > 1 and p[-1] == 0: p = p[:-1]
    return p

def two_row(p, q):
    """P_{(p,q)} for p >= q >= 0 (q=0 -> one row)"""
    if q == 0: return Ppoly((p,)) if p > 0 else [Fraction(1)]
    return Ppoly((p, q))

def check_three_row(N):
    bad = 0; cnt = 0
    for n in range(1, N+1):
        for lam in partitions(n):
            if len(lam) != 3: continue
            a, b, c = lam
            lhs = strip(Ppoly(lam))
            t1 = pmul(two_row(a+2, b+2), Epoly(c))
            t2 = [-x for x in pmul(two_row(a+2, c+1), Epoly(b+1))]
            t3 = pmul(two_row(b+1, c+1), Epoly(a+2))
            rhs = strip(padd(t1, t2, t3))
            cnt += 1
            if lhs != rhs:
                bad += 1; print("3-row MISMATCH", lam)
    print(f"three-row Pfaffian identity: {cnt} shapes with |lam| <= {N}, mismatches = {bad}")

def C(lam, k):
    return sum((j-i)**k for (i, j) in boxes(lam))

def check_three_row_separation(N):
    for name, key in [("(n,C2)", lambda l: (C(l,2),)),
                      ("(n,C2,C1^2)", lambda l: (C(l,2), C(l,1)**2)),
                      ("(n,C2,C1^2,C4)", lambda l: (C(l,2), C(l,1)**2, C(l,4))),
                      ("(n,C1^2)", lambda l: (C(l,1)**2,)),
                      ("(n,C2,r)", lambda l: (C(l,2), len(set(l))))]:
        first = None; ncoll = 0
        for n in range(1, N+1):
            groups = defaultdict(list)
            for a in range(1, n+1):
                for b in range(0, min(a, n-a)+1):
                    c = n-a-b
                    if c < 0 or c > b: continue
                    lam = tuple(x for x in (a, b, c) if x > 0)
                    groups[key(lam)].append(lam)
            for k, v in groups.items():
                if len(v) > 1:
                    # exclude transpose pairs
                    cls = set()
                    for l in v:
                        cls.add(min(l, conjugate(l)))
                    if len(cls) > 1:
                        ncoll += 1
                        if first is None: first = (n, v)
        print(f"3-row shapes (<=3 rows), n <= {N}, invariant {name}: collision groups = {ncoll}, first = {first}")

if __name__ == "__main__":
    check_hermite(9)
    check_Pminus1(9)
    check_three_row(16)
    check_three_row_separation(150)
