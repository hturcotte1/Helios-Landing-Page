import sympy as sp
a,b,c,d,u,v,w,s = sp.symbols('a b c d u v w s')
def rect_pow(p, q, shift, k):
    # sum_{i<p rows, j<q cols} (shift + j - i)^k
    i, j = sp.symbols('i j')
    return sp.expand(sp.summation(sp.summation((shift + j - i)**k, (j, 0, q - 1)), (i, 0, p - 1)))
def pk(k):
    return sp.expand(rect_pow(b, c, 0, k) + rect_pow(b, a, c, k) + rect_pow(d, c, -b, k))
n = pk(0); p1 = pk(1); p2 = pk(2); p4 = pk(4)
print("n =", sp.factor(n))
print("p1 =", sp.factor(p1))
print("p2 =", sp.factor(p2))
print("p2 in terms of S=a+b+c+d etc:")
S = sp.symbols('S')
print(sp.factor(p2.subs(d, S - a - b - c)))
print("p4 =", sp.factor(p4))
# check transpose antisymmetry / symmetry
T = {a: d, b: c, c: b, d: a}
print("p1 anti:", sp.expand(p1.subs(T, simultaneous=True) + p1) == 0, " p2 sym:", sp.expand(p2.subs(T, simultaneous=True) - p2) == 0)
# content generating function check: sum z^content
z = sp.symbols('z')
def rect_gf(p, q, shift):
    i, j = sp.symbols('i j')
    return sp.summation(sp.summation(z**(shift + j - i), (j, 0, q - 1)), (i, 0, p - 1))
for vals in [(1,1,1,1),(2,1,3,2),(3,2,1,4),(1,3,2,2)]:
    aa,bb,cc,dd = vals
    SS = aa+bb+cc+dd; vv = bb+dd
    gf = sum(z**(j-i) for i in range(bb) for j in range(aa+cc)) + sum(z**(j-i) for i in range(bb, bb+dd) for j in range(cc))
    formula = (z**SS - z**(SS-bb) + z**(cc+dd) - z**cc + 1 - z**vv) / (z**(vv-1) * (z-1)**2)
    print(vals, sp.simplify(gf - formula) == 0)
