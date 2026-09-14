"""topend task (2): the level polynomials.
  g_i(lam) := (n)_i u_i(lam)/f^lam = sum_{rhobar even-type, no 1's, |rhobar|<=i} sigma(rhobar) I_{i-|rhobar|} C(n-|rhobar|, i-|rhobar|) Omega_rhobar(n;C)
with Omega from topend_omega_results.pkl (PROVED there for all n).
 (1) print the expansion of u_i as a rational combination of even-type character values (coefficients sigma/z), i <= 10;
 (2) build G_i(n; C_1..C_8) exactly (sympy Rationals), i <= 10;
 (3) verify G_i(lam) == (n)_i d_{n-i}(lam)/f^lam for all lam |- n <= NV using the d-vector only (no characters);
 (4) reduction: A_5 = Q[n, C_2, C_1^2, C_4]; show G_6 in A_5; write G_7..G_10 = (new part) + (element of A_{i-1}),
     using the successive substitutions C_6 -> h7 + 16 C_1 C_3 etc.; print the new generators h_i;
 (5) witness pairs: lam, mu |- n (n <= NW) with equal (C_2, C_1^2, C_4) but different h_7 (so u_7/f is not a function of
     (n, C_2, C_1^2, C_4)); similarly for h_8 given h_7, etc.;
 (6) sanity: sigma(rho) == sum_{mu |- i} chi^mu(rho) for all rho |- i <= 10 (Frobenius-Schur count, checked by MN)."""
import os, sys, pickle, itertools
from fractions import Fraction
from math import comb, factorial
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from topend_mn import *
from young import partitions, f_hook, involutions, conjugate
from census import d_vector

HERE = os.path.dirname(os.path.abspath(__file__))
n = sp.Symbol('n')
Cs = {k: sp.Symbol(f'C{k}') for k in range(1, 11)}
def R(fr): return sp.Rational(fr.numerator, fr.denominator)

def load_omegas():
    with open(os.path.join(HERE, 'topend_omega_results.pkl'), 'rb') as fh:
        terms = pickle.load(fh)
    Om = {}
    for rho, ts in terms.items():
        e = sp.Integer(0)
        for a, (j, mu) in ts:
            m = R(a) * n ** j
            for p in mu: m *= Cs[p]
            e += m
        Om[rho] = sp.expand(e)
    return Om

def falling(x, k):
    r = sp.Integer(1)
    for t in range(k): r *= (x - t)
    return r
def binom_poly(x, k): return falling(x, k) / sp.factorial(k)

def build_G(imax, Om):
    G = {}
    for i in range(imax + 1):
        e = sp.Integer(0)
        for rho, om in Om.items():
            s = sum(rho)
            if s > i: continue
            e += sp.Integer(sigma(rho)) * sp.Integer(involutions(i - s)) * binom_poly(n - s, i - s) * om
        G[i] = sp.expand(e)
    return G

def eval_poly(expr, lam):
    subs = {n: sum(lam)}
    for k, sym in Cs.items(): subs[sym] = C(lam, k)
    v = sp.Rational(expr.subs(subs))
    return Fraction(int(v.p), int(v.q))

def verify_G(G, NV):
    bad = 0; cnt = 0
    for N in range(NV + 1):
        for lam in partitions(N):
            dv = d_vector(lam); f = dv[N]
            for i in G:
                if i > N: continue
                cnt += 1
                if eval_poly(G[i], lam) != Fraction(factorial(N) // factorial(N - i) * dv[N - i], f):
                    bad += 1; print("G mismatch", lam, i)
    print(f"(3) G_i(lam) == (n)_i d_(n-i)(lam)/f^lam for all lam |- n <= {NV}, i <= {max(G)} [d-vector only]: {cnt} checks, {'OK' if bad == 0 else str(bad)+' FAIL'}")
    return bad == 0

def u_expansion(i):
    """u_i = sum_{rho |- i even type} (sigma(rho)/z(rho)) chi^lam(rho u 1^{n-i}); return list of (coef, rho)."""
    return [(Fraction(sigma(rho), z(rho)), rho) for rho in partitions(i) if sigma(rho)]

if __name__ == "__main__":
    NV = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    NW = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    # (6)
    bad = 0
    for i in range(11):
        for rho in partitions(i):
            if sigma(rho) != sum(chi(mu, rho) for mu in partitions(i)): bad += 1
    print(f"(6) sigma(rho) == sum_mu chi^mu(rho) for all rho |- i <= 10: {'OK' if bad == 0 else bad}")
    # (1)
    print("(1) u_i as combinations of character values chi^lam(rho u 1^(n-i)), rho even type:")
    for i in range(3, 11):
        print(f"   u_{i} = " + " + ".join(f"({c})chi{list(rho)}" for c, rho in u_expansion(i)))
    Om = load_omegas()
    print("(2) Omega polynomials (proved for all n by Theorem P + rank check):")
    for rho, om in Om.items(): print(f"   Omega_{list(rho)} =", om)
    G = build_G(10, Om)
    for i in range(3, 11): print(f"   G_{i} =", G[i])
    with open(os.path.join(HERE, 'topend_G_results.pkl'), 'wb') as fh:
        pickle.dump({i: sp.srepr(G[i]) for i in G}, fh)
    verify_G({i: G[i] for i in range(0, 11)}, NV)
    # (4) reduction
    C1, C2, C3, C4, C5, C6, C8 = [Cs[k] for k in (1, 2, 3, 4, 5, 6, 8)]
    print("(4) reduction modulo A_5 = Q[n, C_2, C_1^2, C_4]:")
    free5 = {C3, C5, C6, C8}
    def uses(expr, syms): return bool(expr.free_symbols & syms)
    print("   G_3 - C_2 =", sp.factor(G[3] - C2))
    print("   G_4 - C_1^2 =", sp.expand(G[4] - C1**2), "   [element of Q[n, C_2]]")
    print("   G_5 - C_4 =", sp.expand(G[5] - C4), "   [element of Q[n, C_2, C_1^2]]")
    print("   G_6 uses only", G[6].free_symbols, "=> G_6 in A_5:", not uses(G[6], free5))
    # all monomials of G_6 have even degree in C_1?
    P6 = sp.Poly(G[6], C1); print("   G_6 as polynomial in C_1 has only even powers:", all(m[0] % 2 == 0 for m in P6.monoms()))
    # level 7: new part
    h7 = C6 - 16 * C1 * C3
    rest7 = sp.expand(G[7] - h7)
    print("   h_7 := C_6 - 16 C_1 C_3 ;  G_7 - h_7 uses", rest7.free_symbols, "=> in A_5:", not uses(rest7, free5))
    # substitution C6 -> h7s + 16 C1 C3 for the later levels
    H7 = sp.Symbol('h7'); H8 = sp.Symbol('h8'); H9 = sp.Symbol('h9')
    G8s = sp.expand(G[8].subs(C6, H7 + 16 * C1 * C3))
    print("   G_8 with C_6 = h_7 + 16 C_1 C_3:", G8s)
    # collect the part not in A_7 = Q[n, C_2, C_1^2, C_4, h_7]: monomials involving C3, C5, C8 or odd power of C1
    def split(expr, allowed):
        """return (new, old): old = sum of monomials whose symbols are within `allowed` and have even C1-degree."""
        new = sp.Integer(0); old = sp.Integer(0)
        for term in sp.Add.make_args(sp.expand(expr)):
            d = term.as_powers_dict()
            syms = set(term.free_symbols)
            if syms <= allowed and d.get(C1, 0) % 2 == 0: old += term
            else: new += term
        return sp.expand(new), sp.expand(old)
    A7 = {n, C1, C2, C4, H7}
    new8, old8 = split(G8s, A7)
    print("   new part of G_8 (mod A_7):", sp.factor(new8))
    h8 = new8
    G9s = sp.expand(G[9].subs(C6, H7 + 16 * C1 * C3))
    # eliminate C3^2 using h8: h8 = 2 C3^2 + (a n + b) C1 C3  => C3^2 = (h8 - (a n + b) C1 C3)/2
    P = sp.Poly(h8, C3)
    coeffs = {m[0]: c for m, c in zip(P.monoms(), P.coeffs())}
    print("   h_8 as polynomial in C_3:", coeffs)
    c2, c1 = coeffs.get(2, 0), coeffs.get(1, 0)
    G9s = sp.expand(G9s.subs(C3**2, (H8 - c1 * C3) / c2))
    A8 = A7 | {H8}
    new9, old9 = split(G9s, A8)
    print("   new part of G_9 (mod A_8, after C_3^2 -> (h_8 - ...)/2):", sp.factor(new9))
    h9 = new9
    G10s = sp.expand(G[10].subs(C6, H7 + 16 * C1 * C3))
    # C3^2 -> ..., C3^3 etc: reduce polynomial in C3 modulo the quadratic relation
    def reduce_C3(expr):
        expr = sp.expand(expr)
        while sp.Poly(expr, C3).degree() >= 2:
            expr = sp.expand(expr.subs(C3**2, (H8 - c1 * C3) / c2))  # may leave higher powers; loop
            # sympy subs of C3**2 in C3**k works for k>=2 (C3**3 = C3*C3**2) — ensure by rewriting
            Pq = sp.Poly(expr, C3)
            if Pq.degree() >= 2:
                q, r = sp.div(Pq, sp.Poly(C3**2 + (c1 / c2) * C3 - H8 / c2, C3))
                expr = sp.expand(r.as_expr())
        return expr
    G10s = reduce_C3(G10s)
    # eliminate C8 using h9: h9 = C8 + (..) C1 C5 + ... => C8 = h9 - rest
    P9 = sp.Poly(h9, C8); assert P9.degree() == 1
    c8 = P9.coeffs()[0] if P9.monoms()[0][0] == 1 else None
    C8_expr = sp.expand((H9 - (h9 - sp.Poly(h9, C8).coeff_monomial(C8) * C8)) / sp.Poly(h9, C8).coeff_monomial(C8))
    G10s = sp.expand(G10s.subs(C8, C8_expr))
    G10s = reduce_C3(G10s)
    A9 = A8 | {H9}
    new10, old10 = split(G10s, A9)
    print("   new part of G_10 (mod A_9):", sp.factor(new10))
    with open(os.path.join(HERE, 'topend_levels_newparts.pkl'), 'wb') as fh:
        pickle.dump({7: sp.srepr(h7), 8: sp.srepr(h8), 9: sp.srepr(h9), 10: sp.srepr(new10)}, fh)
    # (5) witness pairs
    print(f"(5) witness pairs (n <= {NW}):")
    def h7f(lam): return C(lam, 6) - 16 * C(lam, 1) * C(lam, 3)
    for N in range(4, NW + 1):
        reps = {}
        for lam in partitions(N):
            key = min(lam, conjugate(lam))
            if key not in reps: reps[key] = lam
        groups = {}
        for lam in reps.values():
            groups.setdefault((C(lam, 2), C(lam, 1) ** 2, C(lam, 4)), []).append(lam)
        found = [(g, [h7f(l) for l in g]) for g in groups.values() if len(g) > 1 and len({h7f(l) for l in g}) > 1]
        if found:
            g, hs = found[0]
            print(f"   n={N}: equal (C_2, C_1^2, C_4) but different h_7: {g}  h_7 = {hs}   (u_7 values: {[d_vector(l)[N-7] for l in g]})")
            break
