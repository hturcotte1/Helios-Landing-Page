"""fathooks: closed forms of n, C_1, C_2 for lam(a,b,c,d) = ((a+c)^b, c^d) (sympy, exact), checked against direct
summation for all 2-corner shapes with n <= 60."""
import sympy as sp
a,b,c,d,i,j = sp.symbols('a b c d i j')
def rect(p, q, shift, k):   # rows i<p, cols j<q, content shift + j - i
    return sp.expand(sp.summation(sp.summation((shift + j - i)**k, (j, 0, q-1)), (i, 0, p-1)))
def pk(k): return sp.expand(rect(b, a + c, 0, k) + rect(d, c, -b, k))
n_, C1, C2 = pk(0), pk(1), pk(2)
print("n  =", sp.factor(n_)); print("C1 =", sp.factor(C1)); print("6*C2 =", sp.expand(6*C2))
T = {a: d, b: c, c: b, d: a}
print("C1 transpose-antisymmetric:", sp.expand(C1.subs(T, simultaneous=True) + C1) == 0)
print("C2 transpose-symmetric:", sp.expand(C2.subs(T, simultaneous=True) - C2) == 0)
def direct(A,B,C,D):
    return sum((jj-ii)**2 for ii in range(B) for jj in range(A+C)) + sum((jj-ii)**2 for ii in range(B,B+D) for jj in range(C))
bad = 0; cnt = 0
f2 = sp.lambdify((a,b,c,d), C2)
for B in range(1, 61):
    for D in range(1, 61):
        for C in range(1, 61):
            if B*C + C*D > 60: break
            for A in range(1, 61):
                if A*B+B*C+C*D > 60: break
                cnt += 1
                if int(C2.subs({a:A,b:B,c:C,d:D})) != direct(A,B,C,D): bad += 1
print(f"C2 formula checked on {cnt} shapes with n <= 60: failures = {bad}")
