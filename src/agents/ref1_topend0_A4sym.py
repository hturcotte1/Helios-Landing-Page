import sympy as sp
n, C1, C2, C4 = sp.symbols('n C1 C2 C4')
binom = lambda a, k: sp.binomial(a, k)
A1 = C2 - sp.expand_func(binom(n, 2))
A2 = C1**2/2 - sp.Rational(3,2)*C2 + sp.expand_func(binom(n,2))
A3 = C4 - (3*n-10)*C2 - 2*C1**2 + n*(n-1)*(5*n-19)/6
A4 = sp.Rational(1,2)*(A1**2 - 2*sp.expand_func(binom(n,3)) - (3*n-8)*A1 - 8*A2 - 5*A3)
A4b = (C2**2/2 - sp.Rational(5,2)*C4 + 3*C1**2 - n**2*C2/2 + sp.Rational(13,2)*n*C2 - 15*C2
       + n**4/8 - sp.Rational(7,4)*n**3 + sp.Rational(47,8)*n**2 - sp.Rational(17,4)*n)
print("A4 expansion difference:", sp.simplify(sp.expand(A4 - A4b)))
print("A4 expanded:", sp.expand(A4))
