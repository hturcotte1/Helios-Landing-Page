"""Referee #2 independent check of Theorem S2 (structure.md section 1.3).
Own code: partitions, d-vector recursion, MN characters, content moments,
Faulhaber constants, Newton identities, recovery algorithm."""
import sys
from fractions import Fraction as Fr
from functools import lru_cache
from itertools import product

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k): yield (k,)+rest

def conj(lam):
    return tuple(sum(1 for p in lam if p > j) for j in range(lam[0])) if lam else ()

@lru_cache(None)
def dvec(lam):
    n = sum(lam)
    if n == 0: return (1,)
    # removable corners
    kids = []
    for i in range(len(lam)):
        if i == len(lam)-1 or lam[i] > lam[i+1]:
            mu = list(lam); mu[i] -= 1
            if mu[i] == 0: mu.pop()
            kids.append(tuple(mu))
    d = [0]*(n+1); d[0] = 1
    for mu in kids:
        dm = dvec(mu)
        for j in range(1, n+1): d[j] += dm[j-1]
    return tuple(d)

# Murnaghan-Nakayama via beta-sets / rim hooks
@lru_cache(None)
def chi(lam, rho):
    """character of irrep lam at cycle type rho (tuple, any order)."""
    if not rho:
        return 1 if sum(lam) == 0 else 0
    r = rho[0]; rest = rho[1:]
    L = len(lam)
    beta = [lam[i] + (L-1-i) for i in range(L)]  # distinct decreasing
    bs = set(beta)
    total = 0
    for b in beta:
        if b - r >= 0 and (b - r) not in bs:
            # height = number of beta elements strictly between b-r and b
            h = sum(1 for c in beta if b-r < c < b)
            nb = sorted([c for c in beta if c != b] + [b-r], reverse=True)
            newlam = tuple(nb[i] - (L-1-i) for i in range(L))
            newlam = tuple(p for p in newlam if p > 0)
            total += (-1)**h * chi(newlam, rest)
    return total

def contents(lam):
    return [j - i for i, row in enumerate(lam) for j in range(row)]
def Ck(lam, k):
    return sum(c**k for c in contents(lam))
def binom(n, k):
    from math import comb; return comb(n, k)

def S(k, m):  # Faulhaber as polynomial, integer m any sign
    m = Fr(m)
    if k == 1: return m*(m+1)/2
    if k == 2: return m*(m+1)*(2*m+1)/6
    if k == 3: return (m*(m+1)/2)**2
    if k == 4: return (6*m**5+15*m**4+10*m**3-m)/30
    raise ValueError

# ---------- (a) Lemma 1.1 and (1.2) on all shapes with <=L rows, n<=N ----------
def check_lemma11(N=22, Ls=(3,4,5)):
    bad = 0
    for L in Ls:
        for n in range(1, N+1):
            for lam in partitions(n):
                if len(lam) > L: continue
                lp = list(lam) + [0]*(L-len(lam))
                x = [lp[i]-(i+1) for i in range(L)]
                p = lambda k: sum(Fr(xi)**k for xi in x)
                assert all(x[i] > x[i+1] for i in range(L-1))
                for k in (1,2,3,4):
                    lhs = Ck(lam, k)
                    rhs = sum(S(k, xi) - S(k, -(i+1)) for i, xi in enumerate(x))
                    rhs2 = sum(S(k, xi) for xi in x) + (-1)**k * sum(S(k, i-1) for i in range(1, L+1))
                    if lhs != rhs or lhs != rhs2: bad += 1; print("L1.1 fail", lam, k, lhs, rhs, rhs2)
                const1 = sum(binom(i, 2) for i in range(1, L+1))
                const2 = sum(S(2, i-1) for i in range(1, L+1))
                const4 = sum(S(4, i-1) for i in range(1, L+1))
                if n != p(1) + L*(L+1)//2: bad += 1; print("n fail", lam)
                if Ck(lam,1) != (p(2)+p(1))/2 - const1: bad += 1; print("C1 fail", lam)
                if Ck(lam,2) != (2*p(3)+3*p(2)+p(1))/6 + const2: bad += 1; print("C2 fail", lam)
                if Ck(lam,4) != (6*p(5)+15*p(4)+10*p(3)-p(1))/30 + const4: bad += 1; print("C4 fail", lam)
        print(f"L={L}: constants n:{L*(L+1)//2} C1:{sum(binom(i,2) for i in range(1,L+1))} C2:{sum(S(2,i-1) for i in range(1,L+1))} C4:{sum(S(4,i-1) for i in range(1,L+1))}")
    # identity (ii)
    for k in (1,2,3,4):
        for m in range(-30, 31):
            assert S(k, -m) == (-1)**(k+1) * S(k, m-1), (k, m)
    print("Lemma 1.1 / (1.2) check done, failures:", bad)

# ---------- (b) sign lemma ----------
def check_sign(N=40):
    neg = []
    for n in range(1, N+1):
        for lam in partitions(n):
            if len(lam) > 4: continue
            a = lam[0]
            c1 = Ck(lam, 1)
            assert 2*c1 >= a*(a-1) - 20, lam
            if len(lam) <= 3: assert 2*c1 >= a*(a-1) - 8, lam
            if c1 <= 0: neg.append(lam)
    print("<=4-row shapes with C1<=0:", len(neg), "max n:", max(sum(l) for l in neg), "largest:", max(neg, key=sum))
    print(sorted(neg, key=sum))
    print("<=3-row with C1<=0:", [l for l in neg if len(l) <= 3])

# ---------- (c) d-vector determines C2, C1^2, C4 : Cor 4.10.2 formulas ----------
def invariants_from_d(d):
    n = len(d)-1
    f = d[n]; u = lambda i: d[n-i]
    # chi(3,1^{n-3}) = 3u3 - 2f ; chi(2,2,1^{n-4}) = 4u4-4u3+f ; chi(5,1^{n-5}) = 5u5-5u4+f
    chi3 = 3*u(3) - 2*f
    chi22 = 4*u(4) - 4*u(3) + f
    chi5 = 5*u(5) - 5*u(4) + f
    K3 = Fr(binom(n,3)*2 * chi3, f); K22 = Fr(3*binom(n,4) * chi22, f); K5 = Fr(24*binom(n,5)*chi5, f)
    C2 = K3 + binom(n,2)
    C1sq = 2*(K22 + Fr(3,2)*C2 - binom(n,2))
    C4 = K5 + (3*n-10)*C2 + 2*C1sq - Fr(n*(n-1)*(5*n-19), 6)
    return n, C1sq, C2, C4

def check_cor4102(N=16):
    bad = 0; cnt = 0
    for n in range(5, N+1):
        for lam in partitions(n):
            d = dvec(lam)
            # first check the character formulas from Lemma 4.5 with my own MN
            f = d[n]
            assert f == chi(lam, tuple([1]*n))
            u = lambda i: d[n-i]
            assert chi(lam, (3,)+ (1,)*(n-3)) == 3*u(3)-2*f, lam
            assert chi(lam, (2,2)+(1,)*(n-4)) == 4*u(4)-4*u(3)+f, lam
            assert chi(lam, (5,)+(1,)*(n-5)) == 5*u(5)-5*u(4)+f, lam
            nn, C1sq, C2, C4 = invariants_from_d(d)
            if (C1sq, C2, C4) != (Ck(lam,1)**2, Ck(lam,2), Ck(lam,4)):
                bad += 1; print("Cor4.10.2 fail", lam, C1sq, C2, C4, Ck(lam,1)**2, Ck(lam,2), Ck(lam,4))
            cnt += 1
    print("Cor 4.10.2 check on all partitions 5<=n<=%d: %d shapes, failures %d" % (N, cnt, bad))

# ---------- (d) recovery algorithm of Theorem S2 ----------
def int_roots_quartic(e1, e2, e3, e4, lo, hi):
    roots = []
    for T in range(lo, hi+1):
        if T**4 - e1*T**3 + e2*T**2 - e3*T + e4 == 0: roots.append(T)
    return roots

def recover4(n, C1sq, C2, C4):
    """All strictly-decreasing integer 4-tuples x (x4>=-4) consistent with invariants."""
    out = set()
    signs = {0} if C1sq == 0 else set()
    import math
    s = math.isqrt(int(C1sq))
    if s*s != C1sq: return out
    signs = {s, -s}
    for C1 in signs:
        p1 = Fr(n-10); p2 = 2*C1 + 20 - p1; p3 = 3*C2 - 60 - (3*p2+p1)/2
        e1 = p1; e2 = (e1*p1 - p2)/2; e3 = (e2*p1 - e1*p2 + p3)/3
        # 30(C4-116) = 6 p5 + 15 p4 + 10 p3 - p1, p4 = A - 4e4, p5 = B - 5 e1 e4 where
        A = e1*p3 - e2*p2 + e3*p1
        # p5 = e1 p4 - e2 p3 + e3 p2 - e4 p1 = e1(A-4e4) - e2 p3 + e3 p2 - e4 p1
        B = e1*A - e2*p3 + e3*p2
        # 30(C4-116) = 6(B - 5 e1 e4) + 15(A - 4 e4) + 10 p3 - p1
        rhs = 30*(C4-116) - 6*B - 15*A - 10*p3 + p1
        slope = -30*e1 - 60
        assert slope == -30*(n-8)
        if slope == 0: continue
        e4 = rhs/slope
        if any(v.denominator != 1 for v in (e1,e2,e3,e4)): continue
        roots = int_roots_quartic(int(e1), int(e2), int(e3), int(e4), -4, n)
        if len(roots) != 4: continue
        x = sorted(roots, reverse=True)
        lam = tuple(x[i] + i + 1 for i in range(4))
        if all(lam[i] >= lam[i+1] for i in range(3)) and lam[3] >= 0:
            out.add(tuple(p for p in lam if p > 0))
    return out

def check_recovery(N=45, from_d=True):
    amb = {}
    cnt = 0
    for n in range(1, N+1):
        for lam in partitions(n):
            if len(lam) > 4: continue
            cnt += 1
            if from_d and n >= 5:
                nn, C1sq, C2, C4 = invariants_from_d(dvec(lam))
                assert (C1sq, C2, C4) == (Ck(lam,1)**2, Ck(lam,2), Ck(lam,4))
            else:
                C1sq, C2, C4 = Ck(lam,1)**2, Ck(lam,2), Ck(lam,4)
            res = recover4(n, C1sq, C2, C4)
            if n == 8:
                assert res == set(), (lam, res); continue
            if lam not in res: print("RECOVERY MISSES", lam, res); amb.setdefault(n, []).append((lam, res))
            if res != {lam}:
                amb.setdefault(n, []).append((lam, res))
    print("recovery checked on", cnt, "shapes with <=4 rows, n<=", N)
    for n in sorted(amb):
        print(" n=%d ambiguous:" % n, amb[n])
        for lam, res in amb[n]:
            assert res == {lam, conj(lam)}, (lam, res)
    print("all ambiguities are exactly transpose pairs:", all(res == {lam, conj(lam)} for n in amb for lam, res in amb[n]))
    print("max ambiguous n:", max(amb) if amb else None)

if __name__ == "__main__":
    check_lemma11()
    check_sign()
    check_cor4102(16)
    check_recovery(int(sys.argv[1]) if len(sys.argv) > 1 else 30, from_d=False)
