"""Symbolic content power sums p_k(a,b,c,d) of the fat hook ((a+c)^b, c^d) and their transpose behaviour."""
import sympy as sp
a,b,c,d,i,j = sp.symbols('a b c d i j')
def rect(p, q, shift, k):   # rows i<p, cols j<q, content shift + j - i
    return sp.expand(sp.summation(sp.summation((shift + j - i)**k, (j, 0, q-1)), (i, 0, p-1)))
def pk(k):
    return sp.expand(rect(b, a + c, 0, k) + rect(d, c, -b, k))
T = {a: d, b: c, c: b, d: a}
for k in range(5):
    e = pk(k)
    print(f"p{k} = {sp.factor(e)}")
    print("   transpose:", "invariant" if sp.expand(e.subs(T, simultaneous=True) - e) == 0 else ("anti" if sp.expand(e.subs(T, simultaneous=True) + e) == 0 else "neither"))
# express p2 in the invariant coordinates s=a+d, u=b+c, P=ab+cd, Q=ad, R=bc
s,u,P,Q,R = sp.symbols('s u P Q R')
p2 = pk(2); p1 = pk(1)
# try: p2 as polynomial in n=P+R, s, u, Q, R ... use Groebner-free approach: substitute d = s - a, c = u - b and check symmetric structure
print("p1 =", sp.factor(p1))
print("2*p1 =", sp.expand(2*p1))
