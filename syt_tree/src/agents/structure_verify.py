"""structure angle: verification of the claims in agent_notes/structure.md.

(1) Theorem S3: n! P_lam(t) = sum_{tau in S_n} chi^lam(tau^2) (1+t)^{fix(tau)}   (brute force over S_n, n<=7)
    and the derangement value n! P_lam(-1) = sum_{tau derangement} chi^lam(tau^2).
(2) Content power-sum formulas for shapes padded to L rows, x_i = lam_i - i:
      n   = p_1 + L(L+1)/2
      C_1 = (p_2+p_1)/2 - sum_{i<=L} i(i-1)/2
      C_2 = (2p_3+3p_2+p_1)/6 + sum_{i<=L} S_2(i-1)
      C_4 = (6p_5+15p_4+10p_3-p_1)/30 + sum_{i<=L} S_4(i-1)
    checked for all lam with <= L rows, n <= 30, L = 3,4,5.
(3) Recovery algorithm: from (n, C_1^2, C_2, C_4) recover a <=4-row shape (Theorem S2) - tested on all
    <=4-row shapes n <= 120; and <=3-row shapes from (n, C_1^2, C_2) (Theorem S1), n <= 200.
    Lower bound C_1 >= 1 for n >= 10 (3 rows) / n >= 21 (4 rows) checked exhaustively.
(4) Sign-of-C_1 lemma exhaustively: list all <=3-row shapes with C_1 <= 0 and all <=4-row shapes with C_1 <= 0.
(5) Pfaffian expansion along the last index for 5-row shapes (ell = 6), n <= 12.
(6) Question (c): no constant-coefficient bilinear differential identity
      sum_c P_{lam-c}(t) P_{lam-c}(u) = sum_{a,b} beta_{ab} P^{(a)}(t) P^{(b)}(u)   (n = 7..10).
(7) Psi_lam(s) = sum_m s^m I_m w_m(lam) takes non-negative integer values at s = 0..5, n <= 10.
"""
import os, sys, itertools
from fractions import Fraction
from math import factorial, comb
from collections import Counter, defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from young import partitions, conjugate, f_hook, boxes, involutions
from census import d_vector
from check_characters import chi

def cycle_type(perm):
    n = len(perm); seen = [False]*n; ct = []
    for i in range(n):
        if seen[i]: continue
        L = 0; j = i
        while not seen[j]:
            seen[j] = True; j = perm[j]; L += 1
        ct.append(L)
    return tuple(sorted(ct, reverse=True))

def check_S3(N):
    bad = 0
    for n in range(1, N+1):
        perms = list(itertools.permutations(range(n)))
        # chi^lam on classes
        classes = {}
        for lam in partitions(n):
            lam = tuple(lam)
            tab = {}
            for rho in partitions(n):
                tab[tuple(rho)] = chi(lam, tuple(rho))
            classes[lam] = tab
        for lam in partitions(n):
            lam = tuple(lam)
            tab = classes[lam]
            # LHS: n! * d_j/j! coefficients
            d = d_vector(lam)
            lhs = [Fraction(factorial(n)*d[j], factorial(j)) for j in range(n+1)]
            rhs = [Fraction(0)]*(n+1)
            der = 0
            for p in perms:
                sq = tuple(p[p[i]] for i in range(n))
                fx = sum(1 for i in range(n) if p[i] == i)
                c = tab[cycle_type(sq)]
                for j in range(fx+1):
                    rhs[j] += c*comb(fx, j)
                if fx == 0: der += c
            val_m1 = sum((-1)**j * lhs[j] for j in range(n+1))
            if lhs != rhs or val_m1 != der:
                bad += 1; print("S3 MISMATCH", lam, lhs, rhs, val_m1, der)
    print(f"(1) Theorem S3 (n! P_lam(t) = sum_tau chi(tau^2)(1+t)^fix) brute force n <= {N}: mismatches = {bad}")

def Cpow(lam, k): return sum((j-i)**k for (i, j) in boxes(lam))
def S2(m): return m*(m+1)*(2*m+1)//6
def S4(m): return (6*m**5+15*m**4+10*m**3-m)//30

def check_content_formulas(N, L):
    bad = 0; cnt = 0
    for n in range(1, N+1):
        for lam in partitions(n):
            if len(lam) > L: continue
            lp = list(lam) + [0]*(L-len(lam))
            x = [lp[i]-(i+1) for i in range(L)]
            p = lambda k: sum(xi**k for xi in x)
            cnt += 1
            ok = (n == p(1) + L*(L+1)//2)
            ok &= Fraction(p(2)+p(1), 2) - sum(i*(i-1)//2 for i in range(1, L+1)) == Cpow(lam, 1)
            ok &= Fraction(2*p(3)+3*p(2)+p(1), 6) + sum(S2(i-1) for i in range(1, L+1)) == Cpow(lam, 2)
            ok &= Fraction(6*p(5)+15*p(4)+10*p(3)-p(1), 30) + sum(S4(i-1) for i in range(1, L+1)) == Cpow(lam, 4)
            if not ok: bad += 1; print("content formula MISMATCH", lam)
    print(f"(2) content power-sum formulas, L={L}, all shapes with <= {L} rows and n <= {N}: {cnt} shapes, mismatches = {bad}")

def recover3(n, C1sq, C2):
    """Theorem S1 algorithm. Returns list of candidate 3-row shapes (a,b,c) (padded) consistent with the data,
    trying both signs of C_1."""
    out = []
    L = 3
    p1 = n - 6
    for C1 in sorted({r for r in (isqrt_exact(C1sq), -isqrt_exact(C1sq)) if r is not None}):
        p2 = 2*(C1 + 4) - p1          # C_1 = (p2+p1)/2 - 4
        p3 = Fraction(6*(C2 - 6) - 3*p2 - p1, 2)   # C_2 = (2p3+3p2+p1)/6 + 6
        if p3.denominator != 1: continue
        p3 = int(p3)
        # elementary symmetric functions
        e1 = p1; e2 = Fraction(e1*p1 - p2, 2); e3 = Fraction(e2*p1 - e1*p2 + p3, 3)
        if e2.denominator != 1 or e3.denominator != 1: continue
        roots = integer_roots([1, -e1, int(e2), -int(e3)])
        if roots is None: continue
        x = sorted(roots, reverse=True)
        if x[0] > x[1] > x[2] >= -3:
            lam = tuple(x[i] + (i+1) for i in range(3))
            lam = tuple(v for v in lam if v > 0)
            if sum(lam) == n: out.append(lam)
    return out

def isqrt_exact(v):
    if v < 0: return None
    r = int(round(v**0.5))
    while r*r > v: r -= 1
    while (r+1)*(r+1) <= v: r += 1
    return r if r*r == v else None

def integer_roots(coeffs):
    """monic integer polynomial (list high->low); return all integer roots with multiplicity if it splits, else None"""
    from sympy import Poly, symbols, roots as sroots
    T = symbols('T')
    P = Poly(sum(c*T**(len(coeffs)-1-i) for i, c in enumerate(coeffs)), T)
    rs = sroots(P)
    tot = sum(rs.values())
    if tot != P.degree(): return None
    out = []
    for r, m in rs.items():
        if not (r.is_integer): return None
        out += [int(r)]*m
    return out

def recover4(n, C1sq, C2, C4):
    """Theorem S2 algorithm: 4 rows (padded)."""
    out = []
    p1 = n - 10
    for C1 in sorted({r for r in (isqrt_exact(C1sq), -isqrt_exact(C1sq)) if r is not None}):
        p2 = 2*(C1 + 10) - p1
        p3 = Fraction(6*(C2 - 20) - 3*p2 - p1, 2)
        if p3.denominator != 1: continue
        p3 = int(p3)
        e1 = p1; e2 = Fraction(e1*p1 - p2, 2); e3 = Fraction(e2*p1 - e1*p2 + p3, 3)
        if e2.denominator != 1 or e3.denominator != 1: continue
        e2 = int(e2); e3 = int(e3)
        # p4 = e1 p3 - e2 p2 + e3 p1 - 4 e4 ; p5 = e1 p4 - e2 p3 + e3 p2 - e4 p1
        # C4 = (6 p5 + 15 p4 + 10 p3 - p1)/30 + 116
        # => 30(C4-116) = 6 p5 + 15 p4 + 10 p3 - p1, linear in e4 with coefficient 6*(-5e1) + 15*(-4) = -30(e1+2)
        A = e1*p3 - e2*p2 + e3*p1           # p4 = A - 4 e4
        B = e1*A - e2*p3 + e3*p2            # p5 = B - 4 e1 e4 - e4 p1 = B - 5 e1 e4
        rhs = 30*(C4 - 116) - (6*B + 15*A + 10*p3 - p1)
        if e1 + 2 == 0: continue
        e4 = Fraction(rhs, -30*(e1+2))
        if e4.denominator != 1: continue
        roots = integer_roots([1, -e1, e2, -e3, int(e4)])
        if roots is None: continue
        x = sorted(roots, reverse=True)
        if x[0] > x[1] > x[2] > x[3] >= -4:
            lam = tuple(x[i] + (i+1) for i in range(4))
            lam = tuple(v for v in lam if v > 0)
            if sum(lam) == n: out.append(lam)
    return out

def check_recovery():
    # 3 rows
    bad3 = []; minC1 = {}
    for n in range(1, 201):
        for a in range(1, n+1):
            for b in range(0, min(a, n-a)+1):
                c = n-a-b
                if c < 0 or c > b: continue
                lam = tuple(v for v in (a, b, c) if v > 0)
                C1 = Cpow(lam, 1)
                if C1 <= 0: minC1.setdefault(3, []).append((lam, C1))
                cands = recover3(n, C1*C1, Cpow(lam, 2))
                if cands != [lam]:
                    bad3.append((lam, cands))
    print(f"(3) Theorem S1 recovery, all <=3-row shapes n <= 200: shapes not uniquely recovered = {len(bad3)}: {bad3}")
    print(f"(4) <=3-row shapes with C_1 <= 0: {minC1.get(3)}")
    bad4 = []; low4 = []
    for n in range(1, 121):
        for lam in partitions(n):
            pass
        # enumerate <=4 rows directly
        for a in range(1, n+1):
            for b in range(0, min(a, n-a)+1):
                for c in range(0, min(b, n-a-b)+1):
                    d = n-a-b-c
                    if d < 0 or d > c: continue
                    lam = tuple(v for v in (a, b, c, d) if v > 0)
                    C1 = Cpow(lam, 1)
                    if C1 <= 0: low4.append((lam, C1))
                    cands = recover4(n, C1*C1, Cpow(lam, 2), Cpow(lam, 4))
                    if cands != [lam]:
                        bad4.append((lam, cands))
    print(f"(3) Theorem S2 recovery, all <=4-row shapes n <= 120: shapes not uniquely recovered = {len(bad4)}")
    for b in bad4: print("     ", b)
    print(f"(4) <=4-row shapes with C_1 <= 0 (max n = {max(sum(l) for l, _ in low4)}): {len(low4)} shapes:", low4)

# ---------- (5) Pfaffian expansion, ell = 6 ----------
def Bij(ai, aj):
    """B_ij(t) coefficient list (Fractions)"""
    out = [Fraction(0)]*(ai+aj+1)
    for a in range(ai+1):
        for b in range(aj+1):
            s = ai - a - aj + b
            sg = (s > 0) - (s < 0)
            if sg: out[a+b] += Fraction(sg, factorial(a)*factorial(b))
    return out

def pmul(a, b):
    out = [Fraction(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b): out[i+j] += x*y
    return out
def padd(*ps):
    L = max(len(p) for p in ps); out = [Fraction(0)]*L
    for p in ps:
        for i, x in enumerate(p): out[i] += x
    return out
def pscale(c, p): return [c*x for x in p]
def strip(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p

def pfaff(M):
    n = len(M)
    if n == 0: return [Fraction(1)]
    tot = [Fraction(0)]
    for j in range(1, n):
        idx = [k for k in range(n) if k not in (0, j)]
        sub = [[M[a][b] for b in idx] for a in idx]
        term = pmul(M[0][j], pfaff(sub))
        tot = padd(tot, pscale((-1)**(j+1), term))
    return tot

def check_pf_expansion(N):
    bad = 0; cnt = 0
    for n in range(1, N+1):
        for lam in partitions(n):
            if len(lam) != 5: continue
            ell = 6
            lp = list(lam)+[0]*(ell-len(lam))
            alpha = [lp[i] + ell-1-i for i in range(ell)]
            M = [[Bij(alpha[i], alpha[j]) if i != j else [Fraction(0)] for j in range(ell)] for i in range(ell)]
            P = strip(pfaff(M))
            # expansion along last index: P = sum_{j<ell} (-1)^{j+1} (1-indexed) B_{j,ell} Pf(M without j,ell)
            tot = [Fraction(0)]
            for j in range(ell-1):
                idx = [k for k in range(ell) if k not in (j, ell-1)]
                sub = [[M[a][b] for b in idx] for a in idx]
                E = [Fraction(1, factorial(i)) for i in range(alpha[j])]   # E_{alpha_j - 1}
                assert strip(M[j][ell-1]) == strip(E)
                tot = padd(tot, pscale((-1)**j, pmul(E, pfaff(sub))))
            cnt += 1
            if strip(tot) != P or P != strip([Fraction(dj, factorial(j)) for j, dj in enumerate(d_vector(tuple(lam)))]):
                bad += 1; print("PF EXPANSION MISMATCH", lam)
    print(f"(5) Pfaffian (ell=6) equals P_lam and its expansion along index 6 for all 5-row shapes n <= {N}: {cnt} shapes, mismatches = {bad}")

# ---------- (6) bilinear test ----------
def Ppoly(lam):
    d = d_vector(tuple(lam)); return [Fraction(dj, factorial(j)) for j, dj in enumerate(d)]
def deriv(p, a):
    for _ in range(a):
        p = [i*p[i] for i in range(1, len(p))] or [Fraction(0)]
    return p

def check_bilinear(n):
    from sympy import Matrix, Rational
    from young import down_set
    lams = [tuple(l) for l in partitions(n)]
    K = n
    unknowns = [(a, b) for a in range(K+1) for b in range(K+1)]
    rows = []; rhs = []
    for lam in lams:
        P = Ppoly(lam)
        Pd = [deriv(P, a) for a in range(K+1)]
        target = [[Fraction(0)]*(n+1) for _ in range(n+1)]
        for mu in down_set(lam):
            Q = Ppoly(mu)
            for i, x in enumerate(Q):
                for j, y in enumerate(Q): target[i][j] += x*y
        for i in range(n+1):
            for j in range(n+1):
                row = []
                for (a, b) in unknowns:
                    va = Pd[a][i] if i < len(Pd[a]) else 0
                    vb = Pd[b][j] if j < len(Pd[b]) else 0
                    row.append(Rational(va*vb))
                rows.append(row); rhs.append(Rational(target[i][j]))
    A = Matrix(rows); bvec = Matrix(rhs)
    rA = A.rank(); rAb = A.row_join(bvec).rank()
    print(f"(6) n={n}: bilinear system sum_c P_(lam-c)(t)P_(lam-c)(u) = sum beta_ab P^(a)(t)P^(b)(u): rank A = {rA}, rank [A|b] = {rAb} -> {'CONSISTENT' if rA == rAb else 'NO SOLUTION'}")

# ---------- (7) Psi ----------
def Hcoef(m):
    c = [0]*(m+1)
    for i in range(0, m//2+1):
        c[m-2*i] = comb(m, 2*i) * (factorial(2*i)//(2**i*factorial(i)))
    return c

def w_from_d(lam):
    """solve P_lam(s-1) = sum_m w_m H_m(s) for w (triangular, H_m monic of degree m)"""
    n = sum(lam)
    P = Ppoly(lam)
    # P(s-1) coefficients in s
    Q = [Fraction(0)]*(n+1)
    for j, c in enumerate(P):
        for i in range(j+1):
            Q[i] += c*comb(j, i)*(-1)**(j-i)
    w = [Fraction(0)]*(n+1)
    for m in range(n, -1, -1):
        w[m] = Q[m]
        for i, c in enumerate(Hcoef(m)):
            Q[i] -= w[m]*c
    return w

def check_psi(N):
    bad = 0
    for n in range(1, N+1):
        for lam in partitions(n):
            w = w_from_d(tuple(lam))
            for s in range(0, 6):
                val = sum(Fraction(s**m * involutions(m)) * w[m] for m in range(n+1))
                if val.denominator != 1 or val < 0:
                    bad += 1; print("PSI not nonneg integer", lam, s, val)
    print(f"(7) Psi_lam(s) = sum_m s^m I_m w_m nonneg integer for s=0..5, all n <= {N}: violations = {bad}")

if __name__ == "__main__":
    check_S3(7)
    for L in (3, 4, 5): check_content_formulas(30, L)
    check_recovery()
    check_pf_expansion(12)
    for n in (7, 8, 9, 10): check_bilinear(n)
    check_psi(10)
