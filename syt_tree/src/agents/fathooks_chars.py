"""fathooks: verify the character identities used for the polynomial invariants.
 (1) chi^lam((123)) = 3 u_3 - 2 f          (briefing fact 5)
 (2) chi^lam((12)(34)) = 4 u_4 - 4 u_3 + f
 (3) chi^lam((12345)) = 5 u_5 - (13/12) f - (5/3) chi(3) - (5/4) chi(2,2)   [from u_5 = 13/60 f + 1/3 chi3 + 1/4 chi22 + 1/5 chi5]
 (4) group algebra: sum_k J_k^2 = C(n,2) e + C_3 ;  sum_k J_k^3 = (2n-3) C_2 + C_4 ;
     sum_k J_k^4 = C_5 + 4 C_22 + (3n-4) C_3 + n(n-1)(4n-5)/6 e    (C_rho = sum of all elements of cycle type rho)
 (5) consequences on V^lam (eigenvalues): p1^2 = 3 p2 - 2C(n,2) + 2 w22,  p4 = w5 + 4 w22 + (3n-4) w3 + n(n-1)(4n-5)/6,
     with w_rho = |C_rho| chi(rho)/f, p_k = content power sums.
 (6) u_6 is a function of (u_3,u_4,u_5,f,n): checked by deriving w33 from w3^2 via C_3*C_3 expansion (numerically).
Characters by Murnaghan-Nakayama (own implementation, exact).
"""
import sys, os
from functools import lru_cache
from math import comb, factorial
from fractions import Fraction
from itertools import permutations
from collections import Counter
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from young import partitions, f_hook, conjugate, boxes
from census import d_vector

# ---- Murnaghan-Nakayama via beta-numbers ----
@lru_cache(maxsize=None)
def mn_char(lam, rho):
    """chi^lam(rho), lam, rho partitions (tuples) of the same n."""
    if not rho:
        return 1 if not lam else 0
    r = rho[0]; rest = rho[1:]
    L = len(lam)
    beta = [lam[i] + (L - 1 - i) for i in range(L)]   # distinct
    bset = set(beta); total = 0
    for idx, bnum in enumerate(beta):
        nb = bnum - r
        if nb < 0 or nb in bset: continue
        # sign = (-1)^(number of beta numbers strictly between nb and bnum)
        sgn = (-1) ** sum(1 for x in beta if nb < x < bnum)
        newbeta = sorted([x for x in beta if x != bnum] + [nb], reverse=True)
        newlam = tuple(x for x in (newbeta[i] - (L - 1 - i) for i in range(L)) if x > 0)
        total += sgn * mn_char(newlam, rest)
    return total

def cycle_type(p):
    n = len(p); seen = [False]*n; ct = []
    for i in range(n):
        if not seen[i]:
            l = 0; j = i
            while not seen[j]:
                seen[j] = True; j = p[j]; l += 1
            ct.append(l)
    return tuple(sorted(ct, reverse=True))

def compose(p, q):  # (p*q)(x) = p(q(x))
    return tuple(p[q[x]] for x in range(len(p)))

def transp(n, i, k):
    p = list(range(n)); p[i], p[k] = p[k], p[i]; return tuple(p)

def check_group_algebra(n):
    """Compute sum_k J_k^m (m=2,3,4) as class functions (element -> coefficient) and compare with claimed class-sum expansions."""
    e = tuple(range(n))
    J = []
    for k in range(n):
        J.append([transp(n, i, k) for i in range(k)])   # J_k = sum_{i<k} (i k), J_0 = 0
    def power_sum(m):
        tot = Counter()
        for k in range(n):
            # J_k^m: all products of m transpositions from J[k]
            cur = Counter({e: 1})
            for _ in range(m):
                nxt = Counter()
                for g, cg in cur.items():
                    for t in J[k]:
                        nxt[compose(g, t)] += cg
                cur = nxt
            tot.update(cur)
        return tot
    def class_sum_coeffs(vec):
        """check vec is a class function; return dict cycle_type -> coefficient"""
        out = {}
        for g in permutations(range(n)):
            ct = cycle_type(g); v = vec.get(g, 0)
            if ct in out:
                assert out[ct] == v, ("not central", ct, out[ct], v)
            else:
                out[ct] = v
        return {ct: v for ct, v in out.items() if v != 0}
    def pad(rho): return tuple(sorted(list(rho) + [1]*(n - sum(rho)), reverse=True))
    ok = True
    e_ct = (1,)*n
    got2 = class_sum_coeffs(power_sum(2)); exp2 = {e_ct: comb(n,2), pad((3,)): 1}
    got3 = class_sum_coeffs(power_sum(3)); exp3 = {pad((2,)): 2*n-3, pad((4,)): 1}
    got4 = class_sum_coeffs(power_sum(4)); exp4 = {pad((5,)): 1, pad((2,2)): 4, pad((3,)): 3*n-4, e_ct: n*(n-1)*(4*n-5)//6}
    for got, exp, name in ((got2,exp2,'J^2'),(got3,exp3,'J^3'),(got4,exp4,'J^4')):
        exp = {k: v for k, v in exp.items() if v != 0 and sum(k) == n}
        if got != exp:
            ok = False; print("MISMATCH", name, "n=",n, got, exp)
    print(f"group algebra identities (J^2,J^3,J^4) n={n}: {'OK' if ok else 'FAIL'}")
    return ok

def content_power(lam, k):
    return sum((j - i) ** k for (i, j) in boxes(lam))

def check_characters(N):
    bad = 0
    for n in range(1, N + 1):
        pad = lambda rho: tuple(sorted(list(rho) + [1]*(n - sum(rho)), reverse=True))
        for lam in partitions(n):
            dv = d_vector(lam); f = f_hook(lam)
            u = lambda i: dv[n - i] if i <= n else 0
            chi3 = mn_char(lam, pad((3,))) if n >= 3 else None
            chi22 = mn_char(lam, pad((2,2))) if n >= 4 else None
            chi5 = mn_char(lam, pad((5,))) if n >= 5 else None
            chi33 = mn_char(lam, pad((3,3))) if n >= 6 else None
            p1 = content_power(lam, 1); p2 = content_power(lam, 2); p4 = content_power(lam, 4)
            if n >= 3 and chi3 != 3*u(3) - 2*f: bad += 1; print("FAIL(1)", lam)
            if n >= 4 and chi22 != 4*u(4) - 4*u(3) + f: bad += 1; print("FAIL(2)", lam)
            if n >= 5:
                rhs = Fraction(5)*u(5) - Fraction(13,12)*f - Fraction(5,3)*chi3 - Fraction(5,4)*chi22
                if rhs != chi5: bad += 1; print("FAIL(3)", lam, rhs, chi5)
            # eigenvalue identities
            w3 = Fraction(2*comb(n,3)*chi3, f) if n >= 3 else Fraction(0)
            w22 = Fraction(3*comb(n,4)*chi22, f) if n >= 4 else Fraction(0)
            w5 = Fraction(24*comb(n,5)*chi5, f) if n >= 5 else Fraction(0)
            if p2 != comb(n,2) + w3: bad += 1; print("FAIL(p2)", lam)
            if p1*p1 != 3*p2 - 2*comb(n,2) + 2*w22: bad += 1; print("FAIL(p1sq)", lam)
            if p4 != w5 + 4*w22 + (3*n-4)*w3 + Fraction(n*(n-1)*(4*n-5), 6): bad += 1; print("FAIL(p4)", lam, p4, w5 + 4*w22 + (3*n-4)*w3 + Fraction(n*(n-1)*(4*n-5), 6))
    print(f"character identities (1)-(3),(5) for all partitions n<={N}: {'OK' if bad == 0 else 'FAIL %d' % bad}")
    return bad == 0

def check_u6_redundant(N):
    """Is u_6 determined by (u_3,u_4,u_5,f,n)? Test: fit u_6 as a linear combination of f, chi3, chi22, chi5, chi3^2/f, n-dependent... 
    Simpler: check that (u3,u4,u5,f,n) equal => u6 equal, over all pairs of partitions of the same n."""
    from collections import defaultdict
    coll = 0
    for n in range(6, N+1):
        g = defaultdict(set)
        for lam in partitions(n):
            dv = d_vector(lam)
            g[(dv[n-3], dv[n-4], dv[n-5], dv[n])].add(dv[n-6])
        coll += sum(1 for v in g.values() if len(v) > 1)
    print(f"u_6 determined by (u3,u4,u5,f) within each n<={N}: {'yes' if coll == 0 else 'NO (%d)' % coll}")

if __name__ == '__main__':
    for n in (4, 5, 6, 7): check_group_algebra(n)
    check_characters(int(sys.argv[1]) if len(sys.argv) > 1 else 16)
    check_u6_redundant(22)
