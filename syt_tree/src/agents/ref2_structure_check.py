"""Referee #2, structure angle, Theorem S1 (<=3 rows). Fully independent code: own MN characters, own d-vector
recursion, own hook-length f, own content power sums. Exact integer arithmetic throughout."""
import sys
from fractions import Fraction
from functools import lru_cache
from math import comb
from collections import defaultdict

# ---------- partitions ----------
def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for a in range(min(n, maxp), 0, -1):
        for rest in partitions(n-a, a): yield (a,)+rest
def shapes3(n, L=3):
    """all partitions of n with <= L rows"""
    if L == 1:
        if n >= 1: yield (n,)
        return
    for a in range((n + L - 1)//L, n+1):
        for rest in (shapes3(n-a, L-1) if n-a > 0 else [()]):
            if not rest or rest[0] <= a: yield (a,)+rest
def conjugate(lam):
    return tuple(sum(1 for p in lam if p > j) for j in range(lam[0])) if lam else ()
def contents(lam):
    return [j - i for i, p in enumerate(lam) for j in range(p)]
def C(lam, k): return sum(c**k for c in contents(lam))
def hook_f(lam):
    lt = conjugate(lam); n = sum(lam)
    from math import factorial
    h = 1
    for i, p in enumerate(lam):
        for j in range(p): h *= (p - j) + (lt[j] - i) - 1
    assert factorial(n) % h == 0
    return factorial(n)//h

# ---------- d-vector: d_j(lam) = number of ways to remove j boxes one at a time ----------
@lru_cache(maxsize=None)
def dvec(lam):
    n = sum(lam)
    if n == 0: return (1,)
    d = [0]*(n+1); d[0] = 1
    for i, p in enumerate(lam):
        if p > 0 and (i+1 == len(lam) or lam[i+1] < p):
            child = list(lam); child[i] -= 1
            child = tuple(v for v in child if v > 0)
            dc = dvec(child)
            for j in range(1, n+1): d[j] += dc[j-1]
    return tuple(d)

# ---------- Murnaghan-Nakayama ----------
def mn_chi(lam, rho):
    """chi^lam(rho) via rim hooks; lam, rho tuples."""
    @lru_cache(maxsize=None)
    def rec(lam, idx):
        if idx == len(rho): return 1 if sum(lam) == 0 else 0
        r = rho[idx]; total = 0
        # remove a rim hook of size r from lam: use beta-numbers
        L = len(lam)
        beta = [lam[i] + (L-1-i) for i in range(L)]
        bset = set(beta)
        for b in beta:
            nb = b - r
            if nb >= 0 and nb not in bset:
                # sign = (-1)^{number of beta numbers strictly between nb and b}
                s = sum(1 for c in beta if nb < c < b)
                newb = sorted([c for c in beta if c != b] + [nb], reverse=True)
                newlam = tuple(newb[i] - (L-1-i) for i in range(L))
                newlam = tuple(v for v in newlam if v > 0)
                total += (-1)**s * rec(newlam, idx+1)
        return total
    return rec(tuple(lam), 0)

def check_char_and_content(N):
    """For all lam |- n <= N: chi(3,1^{n-3}) = 3u_3 - 2f, chi(2,2,1^{n-4}) = 4u_4 - 4u_3 + f, and the
    content formulas C_2 = C(n,2) + omega(K_3), C_1^2 = 2 omega(K_22) + 3 C_2 - 2 C(n,2)."""
    bad = 0; cnt = 0
    for n in range(4, N+1):
        for lam in partitions(n):
            cnt += 1
            d = dvec(lam); f = hook_f(lam)
            assert d[n] == f == d[n-1] == d[n-2]
            u3, u4 = d[n-3], d[n-4]
            chi3 = mn_chi(lam, (3,)+(1,)*(n-3)); chi22 = mn_chi(lam, (2,2)+(1,)*(n-4))
            if chi3 != 3*u3 - 2*f: bad += 1; print("chi3 mismatch", lam)
            if chi22 != 4*u4 - 4*u3 + f: bad += 1; print("chi22 mismatch", lam)
            omK3 = Fraction(2*comb(n,3)*chi3, f); omK22 = Fraction(3*comb(n,4)*chi22, f)   # |K_22| = n!/(8 (n-4)!) = 3 C(n,4)
            C2 = omK3 + comb(n,2)
            C1sq = 2*omK22 + 3*C2 - 2*comb(n,2)
            if C2 != C(lam,2) or C1sq != C(lam,1)**2:
                bad += 1; print("content mismatch", lam, C2, C(lam,2), C1sq, C(lam,1)**2)
    print(f"[A] chi(3), chi(2,2) from u3,u4 and C_2, C_1^2 from them: {cnt} shapes n<={N}, mismatches={bad}")

def S(k, m):
    # Faulhaber S_k(m)=sum_{u=1}^m u^k as polynomial; evaluate for any integer m via closed forms
    m = Fraction(m)
    if k == 1: return m*(m+1)/2
    if k == 2: return m*(m+1)*(2*m+1)/6
    if k == 3: return (m*(m+1)/2)**2
    if k == 4: return (6*m**5+15*m**4+10*m**3-m)/30
def check_faulhaber_facts():
    bad = 0
    for k in (1,2,3,4):
        for m in range(-15, 16):
            if S(k,m) - S(k,m-1) != Fraction(m)**k: bad += 1
            if S(k,-m) != (-1)**(k+1)*S(k,m-1): bad += 1
        for m in range(0, 15):
            if S(k,m) != sum(Fraction(u)**k for u in range(1,m+1)): bad += 1
    print(f"[B] Faulhaber facts (i),(ii) and closed forms, k<=4, |m|<=15: bad={bad}")

def check_lemma11_and_12(N, L):
    bad = 0; cnt = 0
    for n in range(1, N+1):
        for lam in shapes3(n, L):
            cnt += 1
            lp = list(lam) + [0]*(L-len(lam)); x = [lp[i]-(i+1) for i in range(L)]
            assert all(x[i] > x[i+1] for i in range(L-1)) and x[-1] >= -L
            p = lambda k: sum(v**k for v in x)
            for k in (1,2,3,4):
                if C(lam,k) != sum(S(k,xi) - S(k,-(i+1)) for i, xi in enumerate(x)): bad += 1
                if C(lam,k) != sum(S(k,xi) for xi in x) + (-1)**k*sum(S(k,i-1) for i in range(1,L+1)): bad += 1
            c1 = sum(comb(i,2) for i in range(1,L+1)); c2 = sum(S(2,i-1) for i in range(1,L+1)); c4 = sum(S(4,i-1) for i in range(1,L+1))
            if n != p(1) + L*(L+1)//2: bad += 1
            if C(lam,1) != Fraction(p(2)+p(1),2) - c1: bad += 1
            if C(lam,2) != Fraction(2*p(3)+3*p(2)+p(1),6) + c2: bad += 1
            if C(lam,4) != Fraction(6*p(5)+15*p(4)+10*p(3)-p(1),30) + c4: bad += 1
            if L == 3 and (c1,c2,c4) != (4,6,18): bad += 1
            if L == 4 and (c1,c2,c4) != (10,20,116): bad += 1
    print(f"[C] Lemma 1.1 and (1.2) with L={L}: {cnt} shapes n<={N}, bad={bad}")

def check_lemma13(N):
    bad = 0; neg = []
    for n in range(1, N+1):
        for lam in shapes3(n):
            a, b, c = (list(lam)+[0,0,0])[:3]
            c1 = C(lam,1)
            if 2*c1 != a*(a-1)+b*(b-3)+c*(c-5): bad += 1
            if 2*c1 < a*(a-1) - 8: bad += 1
            # row contribution formula
            if c1 != sum(Fraction(l*(l+1-2*(i+1)),2) for i,l in enumerate(lam)): bad += 1
            if c1 <= 0:
                neg.append(lam)
                if a > 3 or n > 9: bad += 1; print("Lemma 1.3 violated", lam)
    print(f"[D] Lemma 1.3, <=3 rows, n<={N}: bad={bad}; shapes with C_1<=0: {neg}")

def check_S1_invariants(N):
    """Group <=3-row shapes by (n, C_1^2, C_2); n>=10 -> singletons; n<=9 -> subsets of {lam, lam^t}."""
    bad = []; amb = []
    for n in range(1, N+1):
        g = defaultdict(list)
        for lam in shapes3(n):
            g[(C(lam,1)**2, C(lam,2))].append(lam)
        for key, cl in g.items():
            if len(cl) > 1:
                amb.append((n, cl))
                if n >= 10 or any(set(cl) - {l, conjugate(l)} for l in cl): bad.append((n, cl))
    print(f"[E] (n,C_1^2,C_2) classes among <=3-row shapes, n<={N}: non-singleton classes at n={sorted(set(a[0] for a in amb))}; violations={bad}")
    print("    non-singleton classes:", amb)

def check_recovery_newton(N):
    """Implement the proof's recovery independently: from (n, C1sq, C2), both signs, Newton -> cubic; brute-force
    which strictly decreasing integer triples x1>x2>x3>=-3 satisfy it (search all triples with lam>=0)."""
    bad = []; amb = []
    for n in range(1, N+1):
        for lam in shapes3(n):
            C1sq, C2 = C(lam,1)**2, C(lam,2)
            found = set()
            s = int(round(C1sq**0.5)); assert s*s == C1sq
            for C1 in {s, -s}:
                p1 = n - 6; p2 = 2*C1 + 8 - p1; p3 = Fraction(3*C2 - 18) - Fraction(3*p2 + p1, 2)
                e1 = p1; e2 = Fraction(e1*p1 - p2, 2); e3 = Fraction(e2*p1 - e1*p2 + p3, 3)
                # brute force triples
                for mu in shapes3(n):
                    mp = list(mu)+[0,0,0]; x = [mp[i]-(i+1) for i in range(3)]
                    if all(x[i] > x[i+1] for i in range(2)) and x[0]-e1 == 0 and x[0]*x[1]+x[0]*x[2]+x[1]*x[2] == e2 and x[0]*x[1]*x[2] == e3 and False: pass
                    if sum(x) == e1 and x[0]*x[1]+x[0]*x[2]+x[1]*x[2] == e2 and x[0]*x[1]*x[2] == e3:
                        found.add(mu)
            if found != {lam}:
                amb.append((lam, found))
                if n >= 8 or found != {lam, conjugate(lam)}: bad.append((lam, found))
    print(f"[F] recovery via Newton, brute-force root check, <=3 rows n<={N}: ambiguous={[a[0] for a in amb]}; violations={bad}")

def check_dvector_collisions(N, rows=None):
    tot_bad = []
    for n in range(1, N+1):
        g = defaultdict(list)
        for lam in partitions(n):
            if rows and len(lam) > rows: continue
            g[dvec(lam)].append(lam)
        for cl in g.values():
            if len(cl) > 1 and any(set(cl) - {l, conjugate(l)} for l in cl): tot_bad.append((n, cl))
    print(f"[G] d-vector collisions ({'all' if not rows else '<=%d rows'%rows}, n<={N}) other than transpose pairs: {tot_bad}")

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "1"):
        check_faulhaber_facts()
        check_char_and_content(16)
        check_lemma11_and_12(40, 3); check_lemma11_and_12(30, 4); check_lemma11_and_12(24, 5)
        check_lemma13(120)
    if which in ("all", "2"):
        check_S1_invariants(200)
        check_recovery_newton(40)
    if which in ("all", "3"):
        check_dvector_collisions(14)
        check_dvector_collisions(30, rows=3)
