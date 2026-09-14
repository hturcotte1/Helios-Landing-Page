"""fathooks: arrangement lemma in the balanced region, improved certificates.
Region B := {1 <= x1 <= x2 <= x3 <= x4 <= x1+x2-1} (integers).  For each pair (A,B) of the 12 arrangements:
 * Delta n factors into linear forms -> tie resolution (recursive): on each tie, the remaining Delta n=0 must force
   another tie ... until tuples coincide or reverse.  Single-variable factors are nonzero.
 * else Delta n linear in x4: x4 = N/D.  Numerator Num of Delta p2 on it.  Each factor F:
     - tie factor: recursive tie resolution using Delta n = 0
     - sign-definite on region via substitution x1=v+w+1, x2=x1+u, x3=x2+v (u,v,w>=0): done
     - else: try certificate using the constraint x3 <= x4  <=>  N - x3 D >= 0 (D>0 on region) and x4 <= x1+x2-1:
       numeric scan over the integer region up to bound to report the sign pattern, and attempt: F*D expressed as
       (sign-definite) + c*(N - x3*D) with c rational found by sympy (linear in unknown c)."""
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

def sign_definite(F):
    G = sp.Poly(sp.expand(F.subs(sub)), u, v, w); cs = G.coeffs()
    if all(c > 0 for c in cs): return '+'
    if all(c < 0 for c in cs): return '-'
    return None

def resolve_ties(A, B, dn, depth=0):
    """dn = 0 must hold. Return True if every solution branch forces A==B or A==rev(B) as tuples."""
    A = tuple(sp.sympify(t) for t in A); B = tuple(sp.sympify(t) for t in B)
    if A == B or A == B[::-1]: return True
    dn = sp.expand(dn)
    if dn == 0: return False
    fl = sp.factor_list(dn)[1]
    for f, m in fl:
        syms = [s for s in X if f.has(s)]
        P = sp.Poly(f, *X)
        if P.total_degree() != 1:
            return False   # cannot resolve
        if len(syms) == 1: continue  # x_i = 0 impossible; or x_i + x_j = 0 impossible
        if len(syms) == 2 and sorted(abs(c) for c in P.coeffs()) == [1,1] and sum(P.coeffs()) == 0:
            s_, v_ = syms   # f = ±(s - v): tie v = s
            A2 = tuple(t.subs(v_, s_) for t in A); B2 = tuple(t.subs(v_, s_) for t in B)
            if A2 == B2 or A2 == B2[::-1]: continue
            # further: on this tie, is there any freedom left? The tie alone satisfies dn=0; need more invariants.
            return ('tie-open', f, A2, B2)
        else:
            if all(c > 0 for c in P.coeffs()) or all(c < 0 for c in P.coeffs()): continue  # never zero
            return False
    return True

def region_scan(F, Nx4, Dx4, bound=40):
    """numeric: signs of F at integer points of region B (with x4 = N/D integral)."""
    signs = set(); pts = []
    for a in range(1, bound+1):
        for b in range(a, bound+1):
            for c in range(b, a+b):
                D = Dx4.subs({x1:a,x2:b,x3:c}); N = Nx4.subs({x1:a,x2:b,x3:c})
                if D == 0: continue
                if N % D != 0: continue
                d = N // D
                if not (c <= d <= a+b-1): continue
                val = F.subs({x1:a,x2:b,x3:c})
                signs.add(sp.sign(val)); 
                if val == 0: pts.append((a,b,c,d))
    return signs, pts

unsettled = []
for A, B in combinations(arrs, 2):
    dn = sp.expand(n_of(*A) - n_of(*B)); dp = sp.expand(p2_of(*A) - p2_of(*B))
    fl = sp.factor_list(dn)[1]
    if all(sp.Poly(f, *X).total_degree() == 1 for f, m in fl):
        r = resolve_ties(A, B, dn)
        if r is not True: unsettled.append((A, B, 'linear', r))
        continue
    sol = sp.solve(dn, x4)[0]; Nx4, Dx4 = sp.fraction(sp.together(sol))
    # normalise D > 0 on region
    if sign_definite(Dx4) == '-': Nx4, Dx4 = -Nx4, -Dx4
    assert sign_definite(Dx4) == '+', (A, B, Dx4)
    num = sp.numer(sp.together(dp.subs(x4, sol)))
    problems = []
    for F, m in sp.factor_list(sp.expand(num))[1]:
        P = sp.Poly(F, *X); syms = [s for s in X if F.has(s)]
        if len(syms) == 1: continue
        if P.total_degree() == 1 and len(syms) == 2 and sum(P.coeffs()) == 0:
            r = resolve_ties(A, B, dn * 0 + F) if False else None
            s_, v_ = syms
            A2 = tuple(sp.sympify(t).subs(v_, s_) for t in A); B2 = tuple(sp.sympify(t).subs(v_, s_) for t in B)
            r = resolve_ties(A2, B2, dn.subs(v_, s_))
            if r is not True: problems.append(('tie', F, r))
            continue
        sd = sign_definite(F)
        if sd: continue
        # certificate attempt: F*D^k + c*(N - x3*D)*M sign-definite?  Try F + c*(N - x3 D)/D ... use F*D + c*(N - x3*D)
        found = None
        for c in [sp.Rational(k, 2) for k in range(-40, 41)]:
            G = sp.expand(F * Dx4 + c * (Nx4 - x3 * Dx4))
            sd2 = sign_definite(G)
            if sd2: found = (c, sd2); break
        if found is None:
            # try with the other constraint  (x1+x2-1) D - N >= 0
            for c in [sp.Rational(k, 2) for k in range(-40, 41)]:
                G = sp.expand(F * Dx4 + c * ((x1 + x2 - 1) * Dx4 - Nx4))
                sd2 = sign_definite(G)
                if sd2: found = ('upper', c, sd2); break
        if found is None:
            signs, pts = region_scan(F, Nx4, Dx4, 40)
            problems.append(('indefinite', F, signs, pts[:5]))
        else:
            problems.append(('certified', F, found))
    hard = [p for p in problems if p[0] != 'certified']
    if hard: unsettled.append((A, B, sol, problems))
    else:
        cert = [p for p in problems if p[0] == 'certified']
        if cert: print("pair", A, B, "certified via constraint:", [(str(p[1]), p[2]) for p in cert])
print("="*60); print("UNSETTLED:", len(unsettled))
for item in unsettled:
    print(item)
