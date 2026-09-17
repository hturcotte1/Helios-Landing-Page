"""Referee #1: line-by-line checks of the proof of Theorem 3 (fathooks.md §2) as exact series, plus large boundary shapes."""
import sys, time
sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src/agents')
from ref1_fathooks_thm3 import *   # runs the n<=NMAX scan (NMAX from argv[1]); import with small NMAX
from fractions import Fraction
from math import factorial, comb

DEGS = 14
def series_over(pred):  # sum_{rho: pred} f^rho t^|rho|/|rho|!  up to degree DEGS
    S = [Fraction(0)] * (DEG + 1)
    for j in range(DEGS + 1):
        for rho in parts(j):
            if pred(rho): S[j] += Fraction(f_hook_own(rho), factorial(j))
    return S
def first(rho): return rho[0] if rho else 0
def length(rho): return len(rho)
def order(S):
    for j, v in enumerate(S):
        if v != 0: return j
    return None

bad = 0
# Step A: G_{x,y} = I - D_x - D_y + E_{x,y}; G_{x,y} = P_{x^y}; G symmetric
for x in range(1, 6):
    for y in range(1, 6):
        G = series_over(lambda r: first(r) <= x and length(r) <= y)
        Dx = series_over(lambda r: first(r) > x)
        Dy = series_over(lambda r: first(r) > y)
        Dy_via_length = series_over(lambda r: length(r) > y)
        E = series_over(lambda r: first(r) > x and length(r) > y)
        I = series_over(lambda r: True)
        rhs = [I[j] - Dx[j] - Dy[j] + E[j] for j in range(DEG + 1)]
        if any(G[j] != rhs[j] for j in range(DEGS + 1)): bad += 1; print("G identity FAIL", x, y)
        if any(Dy[j] != Dy_via_length[j] for j in range(DEGS + 1)): bad += 1; print("D_y via length FAIL", y)
        Prect = ser_from_counts(dvec_own(tuple([x] * y)))
        if any(G[j] != Prect[j] for j in range(min(DEGS, x * y) + 1)): bad += 1; print("G = P_rect FAIL", x, y)
        # vanishing orders
        if order(Dx) != x + 1 or Dx[x + 1] != Fraction(1, factorial(x + 1)): bad += 1; print("ord D FAIL", x)
        if x + y + 1 <= DEGS:
            if order(E) != x + y + 1 or E[x + y + 1] != Fraction(comb(x + y, y), factorial(x + y + 1)): bad += 1; print("ord E FAIL", x, y)
print("Step A (G = I - D_x - D_y + E, orders): failures =", bad)

# Step B: the product expansion.  For all (a,b,c,d) with values <= 4: check
#   G_{a,b} G_{c,d} == I^2 - I (D_a+D_b+D_c+D_d)  mod t^{x1+x2+1}   (pure series identity, no shapes involved)
#   and P_lam == G_{a,b} G_{c,d} mod t^{a+d+1} (Lemma 1(i)) and a+d >= x1+x2.
badB = 0
Dser = {x: series_over(lambda r, x=x: first(r) > x) for x in range(1, 9)}
Gser = {}
for x in range(1, 9):
    for y in range(1, 9):
        Gser[(x, y)] = series_over(lambda r, x=x, y=y: first(r) <= x and length(r) <= y)
Iser = series_over(lambda r: True)
I2 = mul(Iser, Iser)
for a in range(1, 6):
    for b in range(1, 6):
        for c in range(1, 6):
            for d in range(1, 6):
                xs = sorted([a, b, c, d]); x1, x2 = xs[0], xs[1]
                assert a + d >= x1 + x2 and min(a, b) + min(c, d) >= x1 + x2 and a + b >= x1 + x2 and c + d >= x1 + x2
                lhs = mul(Gser[(a, b)], Gser[(c, d)])
                Dsum = [Dser[a][j] + Dser[b][j] + Dser[c][j] + Dser[d][j] for j in range(DEG + 1)]
                rhs = [I2[j] - v for j, v in enumerate(mul(Iser, Dsum))]
                M = min(x1 + x2, DEGS)
                if any(lhs[j] != rhs[j] for j in range(M + 1)): badB += 1; print("expansion FAIL", (a, b, c, d))
                # dividing by I:  I - lhs/I == sum N_k D_k  mod t^{x1+x2+1}
                Qs = [Iser[j] - v for j, v in enumerate(mul(lhs, inv(Iser)))]
                if any(Qs[j] != Dsum[j] for j in range(M + 1)): badB += 1; print("division FAIL", (a, b, c, d))
print("Step B (product expansion & division by I, all params <= 5, deg <= 14): failures =", badB)

# Step C: larger boundary shapes beyond n=40
cases = [(4, 4, 4, 4), (3, 4, 5, 6), (3, 4, 5, 7), (4, 3, 6, 5), (2, 5, 6, 3), (5, 3, 4, 4), (3, 3, 3, 5), (3, 3, 5, 3),
         (2, 2, 3, 12), (2, 12, 3, 2), (1, 8, 8, 8), (1, 7, 8, 8), (8, 1, 8, 7), (6, 4, 2, 9), (5, 5, 5, 5)]
badC = 0
for (a, b, c, d) in cases:
    t0 = time.time()
    lam = shape(a, b, c, d); n = sum(lam)
    dv = dvec_own(lam)
    xs = sorted([a, b, c, d]); x1, x2 = xs[0], xs[1]
    Nk_true = {k: xs.count(k) for k in range(x1, x1 + x2)}
    r1, r2, N, q = read_algorithm(dv)
    P = ser_from_counts(dv); Q = [I_ser[j] - v for j, v in enumerate(mul(P, Iinv))]
    cong = all(Q[j] * factorial(j) == sum(xs.count(k) * Dtab[k][j] for k in set(xs)) for j in range(x1 + x2 + 1))
    ok = cong and (r1, r2) == (x1, x2) and N == Nk_true and ((sum(N.values()) == 4) == (xs[3] <= x1 + x2 - 1))
    if not ok: badC += 1
    print(f"  {(a,b,c,d)} n={n} x1,x2={x1},{x2} read=({r1},{r2},{N}) balanced={xs[3] <= x1+x2-1} ok={ok} [{time.time()-t0:.1f}s]")
print("Step C (large boundary shapes): failures =", badC)
