"""fathooks: arrangement lemma, all 66 pairs of the 12 arrangements (up to reversal) of sorted values x1<=x2<=x3<=x4.
For each pair: factor Delta n; if Delta n is an irreducible quadric, eliminate one variable and factor Delta p2 on it."""
import sympy as sp
from itertools import permutations, combinations
x1,x2,x3,x4 = X = sp.symbols('x1 x2 x3 x4', positive=True)
i,j = sp.symbols('i j')
def rect(p, q, shift, k):
    return sp.expand(sp.summation(sp.summation((shift + j - i)**k, (j, 0, q-1)), (i, 0, p-1)))
def n_of(a,b,c,d): return sp.expand(a*b + b*c + c*d)
def p2_of(a,b,c,d): return sp.expand(rect(b, a+c, 0, 2) + rect(d, c, -b, 2))
arrs = []
for perm in permutations(X):
    if perm[::-1] not in arrs: arrs.append(perm)
assert len(arrs) == 12
hard = []
for A, B in combinations(arrs, 2):
    dn = sp.expand(n_of(*A) - n_of(*B)); dp = sp.expand(p2_of(*A) - p2_of(*B))
    fl = sp.factor_list(dn)[1]
    lin = [f for f, m in fl if sp.Poly(f, *X).total_degree() == 1]
    quad = [f for f, m in fl if sp.Poly(f, *X).total_degree() >= 2]
    if not quad:
        # all factors linear: each factor x_i - x_j = 0 or x_i = 0; check tuples coincide
        ok = True
        for f in lin:
            fs = [s for s in X if f.has(s)]
            if len(fs) == 1: continue   # x_i = 0 impossible
            v = fs[-1]; sol = sp.solve(f, v)[0]
            A2 = tuple(sp.sympify(t).subs(v, sol) for t in A); B2 = tuple(sp.sympify(t).subs(v, sol) for t in B)
            if not (A2 == B2 or A2 == B2[::-1]):
                ok = False; print("LINEAR CASE NOT TRIVIAL", A, B, f)
        print("pair", A, B, ": Delta n =", sp.factor(dn), " -> settled by n" if ok else "")
    else:
        hard.append((A, B, dn, dp))
print(len(hard), "pairs need p2")
for A, B, dn, dp in hard:
    print("="*80); print("pair", A, "vs", B); print("  Delta n =", dn)
    # eliminate a variable in which dn is linear
    for v in X[::-1]:
        if sp.Poly(dn, v).degree() == 1:
            sol = sp.solve(dn, v)[0]
            num = sp.factor(sp.numer(sp.together(dp.subs(v, sol))))
            print(f"  solve for {v} = {sp.factor(sol)}")
            print(f"  Delta p2 numerator on Delta n = 0: {num}")
            break
