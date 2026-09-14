"""fathooks: the local multiset-reading algorithm applied to actual d-vectors of all 2-corner shapes, n <= N.
 Input: d-vector only.  Output: (x1, x2, {N_k : k <= x1+x2-1}).  Check against the true parameters.
 Universal series: I(t) = sum_k I_k t^k/k!, D_x(t) = sum_{rho: rho_1 > x} f^rho t^|rho|/|rho|!.
 Reading: (I - P/I) = sum_x D_x  mod t^{x1+x2+1}."""
import sys, os
from fractions import Fraction
from math import factorial
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from young import partitions, f_hook, involutions
from census import d_vector
from fathooks_local import fathook, fathooks_upto

K = 45
Dx = {}   # Dx[x][k] = (k!)[t^k] D_x = sum_{rho |- k, rho_1 > x} f^rho
for x in range(0, K):
    Dx[x] = [0]*(K+1)
    for k in range(K+1):
        Dx[x][k] = sum(f_hook(r) for r in partitions(k) if r and r[0] > x)
Iser = [Fraction(involutions(k), factorial(k)) for k in range(K+1)]
def ser_mul(p, q, K):
    r = [Fraction(0)]*(K+1)
    for a, pa in enumerate(p):
        if pa == 0: continue
        for b in range(K+1-a): r[a+b] += pa*q[b]
    return r
def ser_inv(p, K):
    r = [Fraction(0)]*(K+1); r[0] = 1/p[0]
    for k in range(1, K+1):
        r[k] = -sum(p[i]*r[k-i] for i in range(1, k+1))/p[0]
    return r
Iinv = ser_inv(Iser, K)

def decode(dv):
    """returns (x1, x2, Ndict) read from the d-vector alone."""
    n = len(dv)-1
    P = [Fraction(dv[k], factorial(k)) if k <= n else Fraction(0) for k in range(K+1)]
    S = ser_mul(P, Iinv, K)                       # P/I
    T = [Iser[k] - S[k] for k in range(K+1)]      # I - P/I = sum_x D_x  (mod t^{x1+x2+1})
    N = {}; found = []
    for k in range(1, K):                          # coefficient of t^k: (1/k!) [ sum_x N_x * Dx[x][k] ]
        known = sum(N[x]*Dx[x][k] for x in N)
        val = T[k]*factorial(k) - known           # = N_{k-1} * Dx[k-1][k] = N_{k-1}
        assert val.denominator == 1
        Nk = int(val)
        if Nk < 0: raise ValueError("negative")
        if Nk: N[k-1] = Nk; found += [k-1]*Nk
        if len(found) >= 2 and k-1 >= found[0] + found[1] - 1:
            break
    x1, x2 = found[0], found[1]
    return x1, x2, {k: v for k, v in N.items() if k <= x1+x2-1}

def main(N):
    bad = 0; cnt = 0
    for (a,b,c,d) in fathooks_upto(N):
        cnt += 1
        xs = sorted([a,b,c,d])
        x1, x2, Nd = decode(d_vector(fathook(a,b,c,d)))
        truth = {k: xs.count(k) for k in range(1, xs[0]+xs[1]) if xs.count(k)}
        if (x1, x2) != (xs[0], xs[1]) or Nd != truth:
            bad += 1; print("DECODE FAIL", (a,b,c,d), x1, x2, Nd, truth)
    print(f"local decoding of (x1,x2,N_k for k<=x1+x2-1) from the d-vector: {cnt} 2-corner shapes n<={N}, {'OK' if bad==0 else 'FAIL %d'%bad}")

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 40)
