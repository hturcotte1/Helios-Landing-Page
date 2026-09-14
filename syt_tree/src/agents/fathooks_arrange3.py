"""fathooks: arrangement lemma in the balanced region  x1<=x2<=x3<=x4<=x1+x2-1.
For every pair (A,B) of distinct arrangements (12 up to reversal) of the sorted values:
  Delta n := n(A)-n(B), Delta p2 := p2(A)-p2(B).
  If Delta n has only linear factors x_i-x_j: on each such tie the tuples must coincide/reverse (checked).
  Else Delta n is linear in x4 (checked); solve x4 = rational(x1,x2,x3); Num := numerator of Delta p2 there.
  Region parametrization: x1 = v+w+1, x2 = x1+u, x3 = x2+v  (u,v,w >= 0 integers)  <=>  1<=x1<=x2<=x3<=x1+x2-1.
  For each irreducible factor F of Num: (i) if F = x_i - x_j: check that the tie together with Delta n = 0 forces trivial pair;
  (ii) else substitute and check all coefficients of F(u,v,w) have the same sign and F is not identically zero;
  also record the x4 formula's requirement x3 <= x4 <= x1+x2-1 as a fallback.
Everything exact (sympy, integer polynomials)."""
import sympy as sp
from itertools import permutations, combinations
x1,x2,x3,x4 = X = sp.symbols('x1 x2 x3 x4')
u,v,w = sp.symbols('u v w')
i,j = sp.symbols('i j')
def rect(p, q, shift, k):
    return sp.expand(sp.summation(sp.summation((shift + j - i)**k, (j, 0, q-1)), (i, 0, p-1)))
def n_of(a,b,c,d): return sp.expand(a*b + b*c + c*d)
def p2_of(a,b,c,d): return sp.expand(rect(b, a+c, 0, 2) + rect(d, c, -b, 2))
sub = {x1: v+w+1, x2: v+w+1+u, x3: 2*v+w+1+u}
arrs = []
for perm in permutations(X):
    if perm[::-1] not in arrs: arrs.append(perm)
def trivial_on_tie(A, B, dn, f):
    """tie f = x_i - x_j = 0 : is (A==B or A==rev B) forced, possibly using Delta n = 0 further?"""
    fs = [s for s in X if f.has(s)]
    if len(fs) != 2: return None
    v_, s_ = fs[1], fs[0]
    A2 = tuple(sp.sympify(t).subs(v_, s_) for t in A); B2 = tuple(sp.sympify(t).subs(v_, s_) for t in B)
    if A2 == B2 or A2 == B2[::-1]: return True
    # use Delta n = 0 on the tie
    dn2 = sp.factor(dn.subs(v_, s_))
    return ('needs', dn2, A2, B2)
def sign_definite(F):
    G = sp.Poly(sp.expand(F.subs(sub)), u, v, w)
    cs = G.coeffs()
    if all(c > 0 for c in cs): return '+'
    if all(c < 0 for c in cs): return '-'
    return None
report = []
for A, B in combinations(arrs, 2):
    dn = sp.expand(n_of(*A) - n_of(*B)); dp = sp.expand(p2_of(*A) - p2_of(*B))
    fl = sp.factor_list(dn)[1]
    if all(sp.Poly(f, *X).total_degree() == 1 for f, m in fl):
        for f, m in fl:
            r = trivial_on_tie(A, B, dn, f)
            if r is not True: report.append(("LINEAR-NONTRIVIAL", A, B, f, r))
        continue
    assert sp.Poly(dn, x4).degree() == 1, (A, B, dn)
    sol = sp.solve(dn, x4)[0]
    num = sp.numer(sp.together(dp.subs(x4, sol)))
    den = sp.denom(sp.together(sol))
    verdicts = []
    for F, m in sp.factor_list(sp.expand(num))[1]:
        if sp.Poly(F, *X).total_degree() == 1 and len([s for s in X if F.has(s)]) == 2 and all(abs(c) == 1 for c in sp.Poly(F, *X).coeffs()):
            r = trivial_on_tie(A, B, dn, F)
            verdicts.append(("tie", F, "trivial" if r is True else r))
        else:
            sd = sign_definite(F)
            verdicts.append(("factor", F, sd if sd else "INDEFINITE"))
    bad = [vd for vd in verdicts if vd[2] not in ('+', '-', 'trivial')]
    report.append(("QUAD", A, B, sol, verdicts, bad))
nbad = 0
for r in report:
    if r[0] == "LINEAR-NONTRIVIAL":
        print("LINEAR case not settled:", r); nbad += 1
    else:
        _, A, B, sol, verdicts, bad = r
        print(f"pair {A} vs {B}: x4 = {sol}")
        for vd in verdicts: print("     ", vd[0], vd[1], "->", vd[2])
        if bad: nbad += 1; print("   *** UNSETTLED")
print("unsettled pairs:", nbad)
