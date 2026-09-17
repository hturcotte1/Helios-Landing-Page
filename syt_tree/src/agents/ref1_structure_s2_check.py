"""Referee #1, Theorem S2 (structure.md 1.3): fully independent re-derivation. No imports from the project.
Exact integer / Fraction arithmetic throughout.
Parts:
 (A) group-algebra identities of Thm 4.10.1 in Z[S_n], n<=7 (brute force over S_n)
 (B) own Murnaghan-Nakayama characters + own d-vector: Lemma 4.5 (u_3,u_4,u_5), sigma(rho) by brute force,
     and the full chain d -> (C_2, C_1^2, C_4) of Cor. 4.10.2, all lambda |- n <= NB
 (C) Lemma 1.1 / (1.2) for L=3,4,5, k<=6, n<=NC ; Lemma 1.3 (sign of C_1) by enumeration
 (D) Newton identities + slope claim on random integer 4-tuples and symbolically (Fractions)
 (E) own recovery algorithm from (n,C1^2,C2,C4) on all <=4-row shapes n<=NE ; brute-force collisions of
     (n,C1^2,C2,C4) among <=4-row shapes (n<=NE) and among ALL shapes (n<=NF)
"""
import sys, itertools, random
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
from collections import defaultdict

NB = int(sys.argv[1]) if len(sys.argv) > 1 else 12
NE = int(sys.argv[2]) if len(sys.argv) > 2 else 60
NF = int(sys.argv[3]) if len(sys.argv) > 3 else 18

# ---------------- basics ----------------
def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for a in range(min(n, maxp), 0, -1):
        for rest in partitions(n-a, a): yield (a,)+rest

def shapes_rows(n, L):
    """partitions of n with at most L parts (direct enumeration)"""
    def go(rem, maxp, k):
        if k == 0:
            if rem == 0: yield ()
            return
        for a in range(min(rem, maxp), -1, -1):
            for rest in go(rem-a, a, k-1): yield (a,)+rest
    for s in go(n, n, L):
        yield tuple(v for v in s if v > 0)

def conjugate(lam):
    if not lam: return ()
    return tuple(sum(1 for l in lam if l >= j) for j in range(1, lam[0]+1))

def contents(lam): return [j-i for i, l in enumerate(lam, 1) for j in range(1, l+1)]
def Ck(lam, k): return sum(c**k for c in contents(lam))

@lru_cache(None)
def dvec(lam):
    n = sum(lam)
    if n == 0: return (1,)
    # children: remove a removable corner
    kids = []
    for i in range(len(lam)):
        if i == len(lam)-1 or lam[i] > lam[i+1]:
            mu = list(lam); mu[i] -= 1
            kids.append(tuple(v for v in mu if v > 0))
    d = [0]*(n+1); d[0] = 1
    for mu in kids:
        dm = dvec(mu)
        for j in range(1, n+1): d[j] += dm[j-1]
    return tuple(d)

def f_of(lam): return dvec(lam)[-1]

# ---------------- own Murnaghan-Nakayama via beta-sets ----------------
@lru_cache(None)
def chi(lam, rho):
    """chi^lam at cycle type rho (tuple, any order), beta-set rim-hook removal."""
    if not rho:
        return 1 if sum(lam) == 0 else 0
    r = rho[0]; rest = rho[1:]
    L = len(lam)
    beta = [lam[i] + (L - 1 - i) for i in range(L)]  # distinct, decreasing
    bset = set(beta)
    tot = 0
    for b in beta:
        nb = b - r
        if nb < 0 or nb in bset: continue
        sign = (-1) ** sum(1 for x in beta if nb < x < b)
        newbeta = sorted([x for x in beta if x != b] + [nb], reverse=True)
        # convert back to partition
        L2 = len(newbeta)
        mu = tuple(v for v in (newbeta[i] - (L2-1-i) for i in range(L2)) if v > 0)
        tot += sign * chi(mu, rest)
    return tot

def chi_cls(lam, rho):
    n = sum(lam); k = sum(rho)
    return chi(lam, tuple(rho) + (1,)*(n-k))

# ---------------- (A) group algebra identities ----------------
def perm_mul(p, q):  # (p*q)(x) = p(q(x))
    return tuple(p[q[x]] for x in range(len(p)))
def cycle_type(p):
    n = len(p); seen = [False]*n; ct = []
    for i in range(n):
        if not seen[i]:
            j = i; c = 0
            while not seen[j]: seen[j] = True; j = p[j]; c += 1
            ct.append(c)
    return tuple(sorted(ct, reverse=True))
def transp(n, i, k):
    p = list(range(n)); p[i], p[k] = p[k], p[i]; return tuple(p)
def ga_mul(A, B):
    C = defaultdict(int)
    for p, a in A.items():
        for q, b in B.items():
            C[perm_mul(p, q)] += a*b
    return C
def ga_add(A, B, cb=1):
    C = defaultdict(int, A)
    for q, b in B.items(): C[q] += cb*b
    return C
def class_sum(n, rho):
    rho = tuple(sorted(rho, reverse=True)); full = rho + (1,)*(n-sum(rho))
    return {p: 1 for p in itertools.permutations(range(n)) if cycle_type(p) == full}
def check_A(nmax):
    ok = True
    for n in range(2, nmax+1):
        e = {tuple(range(n)): 1}
        J = [ {transp(n, i, k): 1 for i in range(k)} for k in range(n) ]  # J_0 = 0 (empty)
        K2 = class_sum(n, (2,)); K3 = class_sum(n, (3,)) if n >= 3 else {}
        K22 = class_sum(n, (2, 2)) if n >= 4 else {}; K5 = class_sum(n, (5,)) if n >= 5 else {}
        S1 = defaultdict(int); S2 = defaultdict(int); S4 = defaultdict(int)
        for k in range(1, n):
            Jk = J[k]; S1 = ga_add(S1, Jk)
            Jk2 = ga_mul(Jk, Jk); S2 = ga_add(S2, Jk2)
            S4 = ga_add(S4, ga_mul(Jk2, Jk2))
        def eq(X, Y):
            keys = set(X) | set(Y)
            return all(X.get(p, 0) == Y.get(p, 0) for p in keys)
        t1 = eq(S1, K2)
        t2 = eq(S2, ga_add({tuple(range(n)): comb(n, 2)}, K3))
        K2sq = ga_mul(K2, K2)
        rhs = ga_add(ga_add({tuple(range(n)): comb(n, 2)}, K3, 3), K22, 2)
        t3 = eq(K2sq, rhs)
        rhs4 = ga_add(ga_add(ga_add({tuple(range(n)): Fraction(n*(n-1)*(4*n-5), 6)}, K5), K3, 3*n-4), K22, 4)
        t4 = eq(S4, rhs4)
        print(f"  (A) n={n}: sumJ=K2 {t1}; sumJ^2=C(n,2)+K3 {t2}; K2^2 {t3}; sumJ^4 {t4}")
        ok &= t1 and t2 and t3 and t4
    return ok

# ---------------- (B) Lemma 4.5 and Cor 4.10.2 chain ----------------
def sigma_brute(rho):
    i = sum(rho); rho = tuple(sorted(rho, reverse=True))
    # number of square roots of ONE fixed permutation pi of type rho (Frobenius-Schur: sum_nu chi^nu(pi))
    pi = []; start = 0
    for part in rho: pi += [start + (j+1) % part for j in range(part)]; start += part
    pi = tuple(pi)
    return sum(1 for p in itertools.permutations(range(i)) if perm_mul(p, p) == pi)
def z_of(rho):
    z = 1
    for part, m in itertools.groupby(sorted(rho)):
        m = len(list(m)); z *= part**m * factorial(m)
    return z
def check_B(nmax):
    # sigma values used
    sig = {rho: sigma_brute(rho) for i in range(1, 6) for rho in partitions(i)}
    print("  (B) sigma(rho) brute force:", {r: s for r, s in sig.items() if s})
    ok = True
    for n in range(5, nmax+1):
        for lam in partitions(n):
            d = dvec(lam); f = d[n]
            u = lambda i: d[n-i]
            # Lemma 4.5 general formula for u_i, i<=5
            for i in range(1, 6):
                val = sum(Fraction(sig[rho], z_of(rho)) * chi_cls(lam, rho) for rho in partitions(i))
                if val != u(i): ok = False; print("   Lemma 4.5 FAIL", lam, i, val, u(i))
            c3 = 3*u(3) - 2*f; c22 = 4*u(4) - 4*u(3) + f; c5 = 5*u(5) - 5*u(4) + f
            if (c3, c22, c5) != (chi_cls(lam, (3,)), chi_cls(lam, (2, 2)), chi_cls(lam, (5,))):
                ok = False; print("   chi-from-u FAIL", lam)
            # central characters
            w3 = Fraction(n*(n-1)*(n-2)//3 * c3, f)
            w22 = Fraction(n*(n-1)*(n-2)*(n-3)//8 * c22, f)
            w5 = Fraction(n*(n-1)*(n-2)*(n-3)*(n-4)//5 * c5, f)
            C2 = w3 + comb(n, 2)
            C1sq = comb(n, 2) + 3*w3 + 2*w22
            C4 = w5 + (3*n-4)*w3 + 4*w22 + Fraction(n*(n-1)*(4*n-5), 6)
            # Thm 4.10.1 printed inversions
            w22b = Fraction(1, 2)*Ck(lam, 1)**2 - Fraction(3, 2)*Ck(lam, 2) + comb(n, 2)
            w5b = Ck(lam, 4) - (3*n-10)*Ck(lam, 2) - 2*Ck(lam, 1)**2 + Fraction(n*(n-1)*(5*n-19), 6)
            if (C2, C1sq, C4) != (Ck(lam, 2), Ck(lam, 1)**2, Ck(lam, 4)) or w22 != w22b or w5 != w5b:
                ok = False; print("   Cor 4.10.2 chain FAIL", lam, (C2, C1sq, C4), (Ck(lam, 2), Ck(lam, 1)**2, Ck(lam, 4)))
    return ok

# ---------------- (C) Lemma 1.1 / (1.2) / Lemma 1.3 ----------------
def S(k, m):
    """Faulhaber polynomial S_k(m)=sum_{u=1}^m u^k extended to all integers via S_k(m)-S_k(m-1)=m^k."""
    if m >= 0: return sum(u**k for u in range(1, m+1))
    # S_k(m) = S_k(m+1) - (m+1)^k, going down from 0
    return -sum(u**k for u in range(m+1, 1))  # sum over u = m+1..0 of u^k, negated
def check_C(nmax):
    ok = True
    for L in (3, 4, 5):
        for n in range(1, nmax+1):
            for lam in shapes_rows(n, L):
                pad = list(lam) + [0]*(L-len(lam)); x = [pad[i-1]-i for i in range(1, L+1)]
                p = lambda k: sum(v**k for v in x)
                for k in range(1, 7):
                    lhs = Ck(lam, k)
                    rhs = sum(S(k, xi) for xi in x) + (-1)**k * sum(S(k, i-1) for i in range(1, L+1))
                    if lhs != rhs: ok = False; print("   Lemma 1.1 FAIL", L, lam, k)
                consts = {3: (6, 4, 6, 18), 4: (10, 10, 20, 116), 5: (15, 20, 50, 470)}[L]
                f12 = (p(1) + consts[0],
                       Fraction(p(2)+p(1), 2) - consts[1],
                       Fraction(2*p(3)+3*p(2)+p(1), 6) + consts[2],
                       Fraction(6*p(5)+15*p(4)+10*p(3)-p(1), 30) + consts[3])
                if f12 != (n, Ck(lam, 1), Ck(lam, 2), Ck(lam, 4)):
                    ok = False; print("   (1.2) FAIL", L, lam, f12)
    # constants recomputed
    for L in (3, 4):
        print(f"  (C) L={L} constants: L(L+1)/2={L*(L+1)//2}, sum C(i,2)={sum(comb(i,2) for i in range(1,L+1))}, "
              f"sum S2(i-1)={sum(S(2,i-1) for i in range(1,L+1))}, sum S4(i-1)={sum(S(4,i-1) for i in range(1,L+1))}")
    # Lemma 1.3: row contribution and the bound
    for i in range(1, 5):
        for l in range(0, 30):
            if sum(j-i for j in range(1, l+1)) * 2 != l*(l+1-2*i): ok = False; print("   row contrib FAIL")
    mins = [min(l*(l-3) for l in range(0, 100)), min(l*(l-5) for l in range(0, 100)), min(l*(l-7) for l in range(0, 100))]
    print("  (C) minima of l(l-3), l(l-5), l(l-7):", mins, "sum", sum(mins))
    neg4 = [lam for n in range(1, 61) for lam in shapes_rows(n, 4) if Ck(lam, 1) <= 0]
    neg3 = [lam for n in range(1, 61) for lam in shapes_rows(n, 3) if Ck(lam, 1) <= 0]
    print(f"  (C) <=4-row shapes with C1<=0 (n<=60): count {len(neg4)}, max n {max(map(sum, neg4))}, max first part {max(l[0] for l in neg4)}")
    print("      list:", neg4)
    print(f"  (C) <=3-row shapes with C1<=0 (n<=60): count {len(neg3)}, max n {max(map(sum, neg3))}:", neg3)
    return ok

# ---------------- (D) Newton ----------------
def check_D():
    ok = True
    rng = random.Random(1)
    for _ in range(2000):
        x = [rng.randint(-30, 30) for _ in range(4)]
        p = lambda k: sum(v**k for v in x)
        e1 = sum(x); e2 = sum(a*b for a, b in itertools.combinations(x, 2))
        e3 = sum(a*b*c for a, b, c in itertools.combinations(x, 3)); e4 = x[0]*x[1]*x[2]*x[3]
        if p(4) != e1*p(3) - e2*p(2) + e3*p(1) - 4*e4: ok = False
        if p(5) != e1*p(4) - e2*p(3) + e3*p(2) - e4*p(1): ok = False
        if 2*e2 != e1*p(1) - p(2) or 3*e3 != e2*p(1) - e1*p(2) + p(3): ok = False
    # slope: vary e4 alone with p1,p2,p3 fixed (formal)
    for _ in range(200):
        p1, p2, p3, e4a, e4b = (rng.randint(-50, 50) for _ in range(5))
        e1 = p1; e2 = Fraction(e1*p1-p2, 2); e3 = Fraction(e2*p1-e1*p2+p3, 3)
        def F(e4):
            p4 = e1*p3 - e2*p2 + e3*p1 - 4*e4; p5 = e1*p4 - e2*p3 + e3*p2 - e4*p1
            return 6*p5 + 15*p4 + 10*p3 - p1
        if e4a != e4b and Fraction(F(e4a)-F(e4b), e4a-e4b) != -30*(e1+2): ok = False; print("   slope FAIL")
    print("  (D) Newton identities and slope -30(e1+2):", ok)
    return ok

# ---------------- (E) recovery + collisions ----------------
def isqrt_exact(v):
    if v < 0: return None
    r = int(round(v**0.5))
    while r*r > v: r -= 1
    while (r+1)*(r+1) <= v: r += 1
    return r if r*r == v else None
def int_roots_monic(coefs, lo, hi):
    """all integer roots in [lo,hi] with multiplicity; returns None if polynomial does not split there."""
    deg = len(coefs)-1; roots = []
    for r in range(lo, hi+1):
        while True:
            v = 0
            for c in coefs: v = v*r + c
            if v != 0: break
            q = [coefs[0]]
            for c in coefs[1:-1]: q.append(c + q[-1]*r)
            coefs = q; roots.append(r)
            if len(coefs) == 1: break
        if len(coefs) == 1: break
    return roots if len(roots) == deg else None
def recover4(n, C1sq, C2, C4, use_sign_lemma):
    out = []
    s = isqrt_exact(C1sq)
    if s is None: return out
    signs = [s] if (use_sign_lemma and n >= 21) else sorted({s, -s})
    p1 = n - 10
    for C1 in signs:
        p2 = 2*(C1 + 10) - p1
        p3 = Fraction(6*(C2 - 20) - 3*p2 - p1, 2)
        e1 = p1; e2 = Fraction(e1*p1 - p2, 2); e3 = Fraction(e2*p1 - e1*p2 + p3, 3)
        if n == 8: continue
        A = e1*p3 - e2*p2 + e3*p1; B = e1*A - e2*p3 + e3*p2   # p4, p5 at e4 = 0
        e4 = Fraction(30*(C4 - 116) - (6*B + 15*A + 10*p3 - p1), -30*(e1 + 2))
        cs = [1, -e1, e2, -e3, e4]
        if any(c.denominator != 1 for c in map(Fraction, cs)): continue
        cs = [int(c) for c in cs]
        x = int_roots_monic(cs, -4, n)
        if x is None: continue
        x = sorted(x, reverse=True)
        if not (x[0] > x[1] > x[2] > x[3] >= -4): continue
        lam = tuple(v for v in (x[0]+1, x[1]+2, x[2]+3, x[3]+4) if v > 0)
        if sum(lam) == n and all(lam[i] >= lam[i+1] for i in range(len(lam)-1)): out.append(lam)
    return out
def check_E(nE, nF):
    ok = True
    # (E1) recovery with the proof's argument (sign lemma for n>=21), all <=4-row shapes n<=nE
    bad = []; cnt = 0
    for n in range(1, nE+1):
        for lam in shapes_rows(n, 4):
            cnt += 1
            r = recover4(n, Ck(lam, 1)**2, Ck(lam, 2), Ck(lam, 4), True)
            if r != [lam]: bad.append((lam, r))
    ns = sorted({sum(l) for l, _ in bad})
    print(f"  (E1) proof's recovery (sign lemma at n>=21) on {cnt} <=4-row shapes, n<={nE}: failures at n in {ns}")
    badbig = [(l, r) for l, r in bad if sum(l) >= 15 and sum(l) != 8]
    print(f"       failures with n>=15: {badbig}")
    nontr = [(l, r) for l, r in bad if sum(l) != 8 and sorted(r) != sorted({l, conjugate(l)})]
    print(f"       failures (n!=8) not equal to {{lam, lam^t}}: {nontr}")
    ok &= not badbig and not nontr
    # (E2) brute-force collisions of (n,C1^2,C2,C4) within <=4-row shapes, n<=nE
    for n in range(1, nE+1):
        groups = defaultdict(list)
        for lam in shapes_rows(n, 4): groups[(Ck(lam, 1)**2, Ck(lam, 2), Ck(lam, 4))].append(lam)
        coll = [v for v in groups.values() if len(v) > 1]
        if coll and (n >= 15 or any(sorted(v) != sorted({v[0], conjugate(v[0])}) for v in coll)):
            print(f"  (E2) n={n}: collisions of (C1^2,C2,C4) among <=4-row shapes: {coll}")
            if n >= 21: ok = False
        elif coll: print(f"  (E2) n={n}: only transpose-pair collisions: {coll}")
    # (E3) collisions of (n,C1^2,C2,C4) among ALL shapes (not claimed; boundary info), n<=nF
    for n in range(1, nF+1):
        groups = defaultdict(list)
        for lam in partitions(n): groups[(Ck(lam, 1)**2, Ck(lam, 2), Ck(lam, 4))].append(lam)
        coll = [v for v in groups.values() if len(v) > 1 and sorted(v) != sorted({v[0], conjugate(v[0])})]
        if coll: print(f"  (E3) n={n}: (C1^2,C2,C4) collisions among ALL shapes beyond transposes (count {len(coll)}), e.g. {coll[:3]}")
    # (E4) d-vector collisions among all shapes n<=nF, own d-vector (sanity for the hypothesis d(lam)=d(mu))
    for n in range(1, nF+1):
        groups = defaultdict(list)
        for lam in partitions(n): groups[dvec(lam)].append(lam)
        coll = [v for v in groups.values() if len(v) > 1 and sorted(v) != sorted({v[0], conjugate(v[0])})]
        if coll: ok = False; print(f"  (E4) n={n}: d-vector collisions beyond transposes: {coll}")
    print(f"  (E4) no non-transpose d-vector collisions for n<={nF}: checked")
    return ok

if __name__ == "__main__":
    print("(A) group algebra identities, n<=7"); a = check_A(7); print("  (A) OK:", a)
    print(f"(B) Lemma 4.5 + Cor 4.10.2 chain, n<={NB}"); b = check_B(NB); print("  (B) OK:", b)
    print("(C) Lemma 1.1, (1.2), Lemma 1.3"); c = check_C(24); print("  (C) OK:", c)
    print("(D) Newton"); d = check_D()
    print(f"(E) recovery n<={NE}, collisions"); e = check_E(NE, NF); print("  (E) OK:", e)
    print("ALL OK:", a and b and c and d and e)
