"""Referee #1 independent checks of Theorem S1 (structure.md 1.1-1.2). Exact arithmetic, no imports from project."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import sys

# ---------- independent basics ----------
def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for a in range(min(n, maxp), 0, -1):
        for rest in partitions(n-a, a): yield (a,)+rest

def shapes3(n):
    for a in range(n, 0, -1):
        for b in range(min(a, n-a), -1, -1):
            c = n-a-b
            if 0 <= c <= b: yield tuple(v for v in (a,b,c) if v > 0)

def boxes(lam): return [(i, j) for i, l in enumerate(lam, 1) for j in range(1, l+1)]
def Ck(lam, k): return sum((j-i)**k for i, j in boxes(lam))
def conjugate(lam):
    if not lam: return ()
    return tuple(sum(1 for l in lam if l >= j) for j in range(1, lam[0]+1))

@lru_cache(None)
def dvec(lam):
    """(d_0..d_n): d_j = number of ways to remove j boxes one at a time."""
    n = sum(lam)
    if n == 0: return (1,)
    d = [0]*(n+1); d[0] = 1
    for i in range(len(lam)):
        if i == len(lam)-1 or lam[i] > lam[i+1]:
            nu = list(lam); nu[i] -= 1
            nu = tuple(v for v in nu if v > 0)
            sub = dvec(nu)
            for j in range(1, n+1): d[j] += sub[j-1]
    return tuple(d)

# ---------- Cor 4.10.2 route: C2 and C1^2 from the d-vector ----------
def invariants_from_d(d):
    n = len(d)-1
    f = d[n]; u3 = d[n-3]; u4 = d[n-4]
    chi3 = 3*u3 - 2*f                       # chi(3,1^{n-3})
    chi22 = 4*u4 - 4*u3 + f                 # chi(2,2,1^{n-4})
    C2 = (Fraction(chi3*n*(n-1)*(n-2), f) + Fraction(3*n*(n-1), 2))/3
    K22size = Fraction(n*(n-1)*(n-2)*(n-3), 8)
    om22 = K22size*Fraction(chi22, f)
    C1sq = 2*(om22 + Fraction(3, 2)*C2 - Fraction(n*(n-1), 2))
    return C1sq, C2

# ---------- Faulhaber / Lemma 1.1 ----------
def S(k, m):
    """Faulhaber polynomial S_k(m) evaluated at integer m (may be negative) using the polynomial identity."""
    # compute via Bernoulli-free approach: S_k(m) for m>=0 by summation; for m<0 use S_k(-m) = (-1)^{k+1} S_k(m-1)
    if m >= 0: return sum(u**k for u in range(1, m+1))
    return (-1)**(k+1) * S(k, -m-1)

def check_faulhaber_identity():
    # (ii) S_k(-m) = (-1)^{k+1} S_k(m-1): verify the polynomial identity independently by Lagrange fit degree k+1
    # We instead verify: the function S(k,.) defined above satisfies S(k,m)-S(k,m-1)=m^k for all integers m in [-30,30]
    for k in range(1, 7):
        for m in range(-30, 31):
            assert S(k, m) - S(k, m-1) == m**k, (k, m)
    return True

def check_lemma11(L, nmax):
    for n in range(1, nmax+1):
        for lam in partitions(n, n):
            if len(lam) > L: continue
            x = [ (lam[i] if i < len(lam) else 0) - (i+1) for i in range(L)]
            p = lambda k: sum(v**k for v in x)
            assert n == p(1) + L*(L+1)//2
            for k in (1, 2, 3, 4, 5, 6):
                lhs = Ck(lam, k)
                rhs = sum(S(k, xi) for xi in x) + (-1)**k * sum(S(k, i-1) for i in range(1, L+1))
                assert lhs == rhs, (lam, k, lhs, rhs)
            c1 = Fraction(p(2)+p(1), 2) - sum(i*(i-1)//2 for i in range(1, L+1))
            c2 = Fraction(2*p(3)+3*p(2)+p(1), 6) + sum(S(2, i-1) for i in range(1, L+1))
            c4 = Fraction(6*p(5)+15*p(4)+10*p(3)-p(1), 30) + sum(S(4, i-1) for i in range(1, L+1))
            assert (c1, c2, c4) == (Ck(lam,1), Ck(lam,2), Ck(lam,4)), lam
    return True

# ---------- independent recovery from (n, C1^2, C2), L=3 ----------
def recover3(n, C1sq, C2):
    """Return the set of all <=3-row partitions mu of n with C1(mu)^2==C1sq and C2(mu)==C2, via Newton."""
    out = set()
    s = None
    r = int(round(C1sq**0.5))
    for c in (r-1, r, r+1):
        if c >= 0 and c*c == C1sq: s = c
    if s is None: return out
    for C1 in {s, -s}:
        p1 = n - 6
        p2 = 2*C1 + 8 - p1
        p3 = Fraction(6*C2 - 36 - 3*p2 - p1, 2)
        e1 = Fraction(p1); e2 = (e1*p1 - p2)/2; e3 = (e2*p1 - e1*p2 + p3)/3
        # roots: brute force over triples? cheaper: test all integers in [-3, n] as roots
        roots = [t for t in range(-3, n+1) if t**3 - e1*t**2 + e2*t - e3 == 0]
        # need exactly 3 distinct roots (strictly decreasing) -> product check
        if len(roots) == 3:
            x = sorted(roots, reverse=True)
            # verify polynomial equals prod (T - x_i): compare coefficients
            if (x[0]+x[1]+x[2] == e1 and x[0]*x[1]+x[0]*x[2]+x[1]*x[2] == e2 and x[0]*x[1]*x[2] == e3):
                lam = tuple(v for v in (x[0]+1, x[1]+2, x[2]+3) if v > 0)
                if sum(lam) == n and all(lam[i] >= lam[i+1] for i in range(len(lam)-1)):
                    out.add(lam)
    return out

if __name__ == "__main__":
    print("Faulhaber (i)/(ii) consistency:", check_faulhaber_identity())
    for L in (3, 4, 5):
        print(f"Lemma 1.1 / (1.2) for L={L}, n<=24:", check_lemma11(L, 24))
    # Cor 4.10.2 route: from d-vector, all partitions n<=16 (all rows), compare with direct contents
    bad = 0; cnt = 0
    for n in range(4, 17):
        for lam in partitions(n):
            cnt += 1
            C1sq, C2 = invariants_from_d(dvec(lam))
            if C1sq != Ck(lam,1)**2 or C2 != Ck(lam,2):
                bad += 1; print("MISMATCH Cor4.10.2", lam, C1sq, C2, Ck(lam,1)**2, Ck(lam,2))
    print(f"Cor 4.10.2 (C1^2, C2 from d) checked on {cnt} partitions, 4<=n<=16: mismatches={bad}")
    # Lemma 1.3: <=3-row shapes with C1<=0, n<=80; and formula 2C1 = a(a-1)+b(b-3)+c(c-5)
    neg = []
    for n in range(1, 81):
        for lam in shapes3(n):
            a, b, c = (lam + (0,0,0))[:3]
            assert 2*Ck(lam,1) == a*(a-1)+b*(b-3)+c*(c-5), lam
            assert 2*Ck(lam,1) >= a*(a-1) - 8
            if Ck(lam,1) <= 0: neg.append(lam)
    print("<=3-row shapes with C1<=0 (n<=80):", neg, "max n =", max(map(sum, neg)))
    # Direct collision search on (n, C1^2, C2) among <=3-row shapes, n<=120, and recovery check
    from collections import defaultdict
    amb = []
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    for n in range(1, NMAX+1):
        groups = defaultdict(list)
        for lam in shapes3(n):
            groups[(Ck(lam,1)**2, Ck(lam,2))].append(lam)
            rec = recover3(n, Ck(lam,1)**2, Ck(lam,2))
            if rec != {lam}:
                amb.append((lam, sorted(rec)))
        for key, g in groups.items():
            if len(g) > 1: print("COLLISION (n,C1^2,C2):", n, key, g)
    print("Recovery ambiguities (n<=%d):" % NMAX)
    for lam, rec in amb:
        ok = set(rec) == {lam, conjugate(lam)}
        print("  ", lam, rec, "== {lam,lam^t}" if ok else "  !!! NOT transpose pair")
    print("max n with ambiguity:", max((sum(l) for l, _ in amb), default=None))
    # Also: with equal d-vectors among <=3-row shapes n<=16 (sanity: full conjecture within class)
    for n in range(1, 17):
        seen = {}
        for lam in partitions(n):
            if len(lam) > 3: continue
            key = dvec(lam)
            if key in seen and seen[key] != conjugate(lam): print("d-COLLISION", n, lam, seen[key])
            seen[key] = lam
    print("d-vector collision check among <=3-row shapes n<=16 done")
