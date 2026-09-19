"""fathooks: verify the explicit formulas used in the six-families theorem.
 F1(c,d) = (c+1, c^d) = runs (1,1,c,d);  F2(b,c) = ((c+1)^b, c) = runs (1,b,c,1);  I_k = involution numbers.
 Claims (c,d,b >= 1 unless stated):
  (S1) d_j(F2(b,c)) = d_{j+1}((c+1)^{b+1}) for all j.
  (S2) for b <= c: d_j(F2) = I_{j+1} (j<=b); d_{b+1}(F2) = I_{b+2} - 1 - [c=b]; d_{b+2}(F2) = I_{b+3} - (b+3) - [c=b+1] - (b+3)[c=b]  (needs b+2 <= n).
  (S3) x = min(c,d): d_j(F1) = I_{j+1} (j<=x); d_{x+1}(F1) = I_{x+2} - [c=x] - [d=x];
       if c < d: d_{c+2}(F1) = I_{c+3} - 2(c+2) - [d=c+1];   if d <= c: d_{d+2}(F1) = I_{d+3} - (d+3) - [c=d+1] - 2(d+2)[c=d].
  (S4) rectangle: d_j(c^d) = I_j for j <= x=min(c,d); d_{x+1} = I_{x+1} - [c=x] - [d=x]; d_{x+2}(c^d) = I_{x+2} - (x+2)([c=x]+[d=x]) - [c=x+1] - [d=x+1]  (x+2 <= cd).
 Also: brute-force Conjecture B for all lam with r=2 and A_1+B_1>=2 against ALL partitions of n, n <= 30.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from young import partitions, corner_runs, conjugate, involutions as I
from census import d_vector

def F1(c, d): return tuple([c+1] + [c]*d)
def F2(b, c): return tuple([c+1]*b + [c])
def rect(c, d): return tuple([c]*d)

def main(N):
    bad = 0
    for c in range(1, N):
        for d in range(1, N):
            if c*d + c + 1 > N: continue
            # S1 with (b,c) := (d,c)
            b = d
            if (c+1)*b + c <= N:
                dv = d_vector(F2(b, c)); dR = d_vector(rect(c+1, b+1))
                if any(dv[j] != dR[j+1] for j in range(len(dv))): bad += 1; print("S1 FAIL", b, c)
                if b <= c:
                    n = (c+1)*b + c
                    if any(dv[j] != I(j+1) for j in range(b+1)): bad += 1; print("S2a FAIL", b, c)
                    if dv[b+1] != I(b+2) - 1 - (c == b): bad += 1; print("S2b FAIL", b, c)
                    if b+2 <= n and dv[b+2] != I(b+3) - (b+3) - (c == b+1) - (b+3)*(c == b): bad += 1; print("S2c FAIL", b, c, dv[b+2])
            # S3
            dv = d_vector(F1(c, d)); x = min(c, d); n = c*d + c + 1
            if any(dv[j] != I(j+1) for j in range(x+1)): bad += 1; print("S3a FAIL", c, d)
            if dv[x+1] != I(x+2) - (c == x) - (d == x): bad += 1; print("S3b FAIL", c, d)
            if c < d and c+2 <= n and dv[c+2] != I(c+3) - 2*(c+2) - (d == c+1): bad += 1; print("S3c FAIL", c, d)
            if d <= c and d+2 <= n and dv[d+2] != I(d+3) - (d+3) - (c == d+1) - 2*(d+2)*(c == d): bad += 1; print("S3d FAIL", c, d, dv[d+2])
            # S4
            if c*d <= N:
                dr = d_vector(rect(c, d))
                if any(dr[j] != I(j) for j in range(x+1)): bad += 1; print("S4a FAIL", c, d)
                if x+1 <= c*d and dr[x+1] != I(x+1) - (c == x) - (d == x): bad += 1; print("S4b FAIL", c, d)
                if x+2 <= c*d and dr[x+2] != I(x+2) - (x+2)*((c == x)+(d == x)) - (c == x+1) - (d == x+1): bad += 1; print("S4c FAIL", c, d)
    print(f"formulas S1-S4 checked for all parameter values with n <= {N}: {'OK' if bad == 0 else 'FAIL %d' % bad}")

def brute(N):
    bad = 0; cnt = 0
    for n in range(1, N+1):
        groups = {}
        for lam in partitions(n):
            groups.setdefault(d_vector(lam), []).append(lam)
        for lam in partitions(n):
            runs = corner_runs(lam)
            if len(runs) != 2: continue
            (a, b), (c, d) = runs
            if [a, b, c, d].count(1) < 2: continue
            cnt += 1
            same = groups[d_vector(lam)]
            if any(mu != lam and mu != conjugate(lam) for mu in same): bad += 1; print("B FAILS", lam, same)
    print(f"Conjecture B for r=2 shapes with >=2 unit parameters vs all partitions, n <= {N}: {cnt} shapes, {'OK' if bad == 0 else 'FAIL'}")

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 60)
    brute(int(sys.argv[2]) if len(sys.argv) > 2 else 30)
