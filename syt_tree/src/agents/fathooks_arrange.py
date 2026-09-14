"""fathooks: arrangement lemma. Given the multiset {x1,x2,x3,x4} of parameters, do n and p2 determine the
arrangement (a,b,c,d) up to reversal?  Symbolic analysis of Delta n, Delta p2 for pi=(x1,x2,x3,x4) vs each other ordering."""
import sympy as sp
from itertools import permutations
x1,x2,x3,x4 = X = sp.symbols('x1 x2 x3 x4')
i,j = sp.symbols('i j')
def rect(p, q, shift, k):
    return sp.expand(sp.summation(sp.summation((shift + j - i)**k, (j, 0, q-1)), (i, 0, p-1)))
def n_of(a,b,c,d): return sp.expand(a*b + b*c + c*d)
def p2_of(a,b,c,d): return sp.expand(rect(b, a+c, 0, 2) + rect(d, c, -b, 2))
def p1_of(a,b,c,d): return sp.expand(rect(b, a+c, 0, 1) + rect(d, c, -b, 1))
pi = (x1,x2,x3,x4)
n0, p20, p10 = n_of(*pi), p2_of(*pi), p1_of(*pi)
seen = set()
for perm in permutations(X):
    if perm == pi or perm == pi[::-1]: continue
    if perm[::-1] in seen: continue
    seen.add(perm)
    dn = sp.factor(n_of(*perm) - n0)
    dp2 = sp.factor(p2_of(*perm) - p20)
    print("pi' =", perm)
    print("   Delta n  =", dn)
    print("   Delta p2 =", dp2)
    # on each linear factor of Delta n, restrict Delta p2
    num = sp.factor_list(sp.expand(n_of(*perm) - n0))[1]
    for fac, mult in num:
        if sp.degree(fac, gen=None) if False else sp.Poly(fac, *X).total_degree() == 1:
            # solve fac = 0 for the last variable appearing
            v = [s for s in X if fac.has(s)][-1]
            sol = sp.solve(fac, v)[0]
            r = sp.factor(sp.expand((p2_of(*perm) - p20).subs(v, sol)))
            print("     on", fac, "= 0  (", v, "=", sol, "):  Delta p2 =", r)
