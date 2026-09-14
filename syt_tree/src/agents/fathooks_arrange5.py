"""fathooks: arrangement lemma in the balanced region with invariants (n, p2, h), h = number of hook pairs
{a,b},{c,d},{a,d} equal (as multisets of values) to {x1,x2}.  Split by tie pattern (composition of 4), strict variables y1<...<yk.
Region parametrization (all s_i >= 0 integers):  y_{i+1} = y_i + 1 + s_{i+1};  y_1 chosen so that x4 <= x1+x2-1 holds.
For each (pair, pattern): trivial (tuples equal/reversed) | h differs | Delta n never zero | eliminate a variable via Delta n=0,
then every factor of the Delta p2 numerator must be sign-definite (all coefficients same sign, nonzero) in the s-parametrization,
or the eliminated variable's formula violates the region (checked by sign-definiteness of N - y*D style expressions).
Prints every case with its verdict; the final line counts unsettled cases."""
import sympy as sp
from itertools import permutations, combinations
x1,x2,x3,x4 = X = sp.symbols('x1 x2 x3 x4')
i,j = sp.symbols('i j')
def rect(p, q, shift, k):
    return sp.expand(sp.summation(sp.summation((shift + j - i)**k, (j, 0, q-1)), (i, 0, p-1)))
def n_of(a,b,c,d): return sp.expand(a*b + b*c + c*d)
def p2_of(a,b,c,d): return sp.expand(rect(b, a+c, 0, 2) + rect(d, c, -b, 2))
N_expr = n_of(x1,x2,x3,x4); P2_expr = p2_of(x1,x2,x3,x4)
def n_arr(A): return N_expr.subs(dict(zip(X, A)), simultaneous=True)
def p2_arr(A): return P2_expr.subs(dict(zip(X, A)), simultaneous=True)
arrs = []
for perm in permutations(X):
    if perm[::-1] not in arrs: arrs.append(perm)
def compositions(n):
    if n == 0: yield (); return
    for f in range(1, n+1):
        for r in compositions(n-f): yield (f,) + r
def all_pos(P):
    cs = sp.Poly(P).coeffs() if P != 0 else [0]
    if P == 0: return None
    if all(c > 0 for c in cs): return '+'
    if all(c < 0 for c in cs): return '-'
    return None
def hcount(A, x1v, x2v):
    a,b,c,d = A; target = sorted([x1v, x2v], key=lambda t: sp.default_sort_key(t))
    cnt = 0
    for P in ((a,b),(c,d),(a,d)):
        if sorted(P, key=lambda t: sp.default_sort_key(t)) == target: cnt += 1
    return cnt
unsettled = []; total = 0
for pat in compositions(4):
    k = len(pat); Y = sp.symbols('y1:%d' % (k+1)); S = sp.symbols('s1:%d' % (k+1))
    # x -> y map
    xy = {}; idx = 0
    for gi, g in enumerate(pat):
        for _ in range(g): xy[X[idx]] = Y[gi]; idx += 1
    x1v, x2v, x4v = xy[x1], xy[x2], xy[x4]
    # parametrization
    ys = {}
    # y_i = y_1 + (i-1) + sum_{j=2..i} s_j
    def y_expr(iidx, y1): return y1 + (iidx) + sum(S[jj] for jj in range(1, iidx+1))
    if pat[0] >= 2:  # x2 = y1: need y_k <= 2 y1 - 1  -> y1 >= k + sum_{j>=2} s_j
        y1 = k + sum(S[1:]) + S[0]
    else:            # x2 = y2 = y1+1+s2: need y_k <= 2y1 + s2 -> y1 >= (k-1) + sum_{j>=3} s_j
        y1 = (k-1) + sum(S[2:]) + S[0]
    for ii in range(k): ys[Y[ii]] = sp.expand(y_expr(ii, y1))
    # sanity: region constraints hold identically with nonneg s: y1>=1, y_{i+1}-y_i>=1, x1+x2-1-x4>=0
    assert all_pos(ys[Y[0]]) == '+'
    assert all_pos(sp.expand(x1v.subs(ys)+x2v.subs(ys)-1-x4v.subs(ys))) in ('+', None) and sp.Poly(sp.expand(x1v.subs(ys)+x2v.subs(ys)-1-x4v.subs(ys))).coeffs() and all(c >= 0 for c in sp.Poly(sp.expand(x1v.subs(ys)+x2v.subs(ys)-1-x4v.subs(ys)), *S).coeffs())
    for A, B in combinations(arrs, 2):
        total += 1
        Ay = tuple(sp.sympify(t).subs(xy) for t in A); By = tuple(sp.sympify(t).subs(xy) for t in B)
        if Ay == By or Ay == By[::-1]: continue
        hA, hB = hcount(Ay, x1v, x2v), hcount(By, x1v, x2v)
        if hA != hB: continue
        dn = sp.expand(n_arr(Ay) - n_arr(By)); dp = sp.expand(p2_arr(Ay) - p2_arr(By))
        dn_s = sp.expand(dn.subs(ys)); dp_s = sp.expand(dp.subs(ys))
        if all_pos(dn_s): continue                     # Delta n never vanishes in region
        if dn_s == 0:
            if all_pos(dp_s): continue
            unsettled.append((pat, A, B, 'dn==0', sp.factor(dp))); continue
        # eliminate a variable appearing linearly in dn (prefer largest y)
        done = False
        for yv in Y[::-1]:
            if sp.Poly(dn, yv).degree() == 1:
                sol = sp.solve(dn, yv)[0]; Nn, Dd = sp.fraction(sp.together(sol))
                # region check of the eliminated variable: is yv = Nn/Dd compatible?  Use: substitute all OTHER y's by parametrization
                others = {yy: ys[yy] for yy in Y if yy != yv}
                Dd_s = sp.expand(Dd.subs(others)); Nn_s = sp.expand(Nn.subs(others))
                sD = all_pos(Dd_s); sN = all_pos(Nn_s)
                if sD and sN and sD != sN:   # yv would be negative
                    done = True; break
                num = sp.numer(sp.together(dp.subs(yv, sol)))
                fl = sp.factor_list(sp.expand(num))[1]
                allok = True; bad = []
                for F, m in fl:
                    Fs = sp.expand(F.subs(others))
                    if all_pos(Fs): continue
                    # try the constraint  yv >= y_{prev}+1  and  yv <= (upper) via  F*Dd + c*(Nn - (y_prev+1)*Dd)
                    ok2 = False
                    yi = list(Y).index(yv)
                    lower = (Y[yi-1] + 1) if yi > 0 else 1
                    # upper bound: balanced  x4 <= x1+x2-1  and ordering yv <= y_{next}-1
                    uppers = []
                    if yi < k-1: uppers.append(Y[yi+1] - 1)
                    # balanced bound in terms of yv if yv is x4 or contributes
                    bal = x1v + x2v - 1 - x4v   # >= 0
                    if sD == '+':
                        cons = [sp.expand(Nn - lower*Dd)]  # >= 0
                        for up in uppers: cons.append(sp.expand(up*Dd - Nn))  # >= 0
                        if bal.has(yv):
                            cons.append(sp.expand((bal.subs(yv, sol))*Dd))
                        for con in cons:
                            for c in [sp.Rational(kk, 4) for kk in range(-80, 81)]:
                                G = sp.expand((F*Dd + c*con).subs(others))
                                if all_pos(G): ok2 = True; break
                            if ok2: break
                    if not ok2: allok = False; bad.append(F)
                if allok: done = True
                else: unsettled.append((pat, A, B, yv, sol, bad))
                break
        if not done and not any(u[1] == A and u[2] == B and u[0] == pat for u in unsettled):
            unsettled.append((pat, A, B, 'no linear variable', dn))
print("cases:", total, " unsettled:", len(unsettled))
for u in unsettled: print(u)
