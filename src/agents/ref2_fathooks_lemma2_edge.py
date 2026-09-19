"""Referee #2: (a) push Lemma 2(ii)/(iii) to larger n, (b) probe j = m+2 to see the window is genuinely exact
(i.e. the lemma's claimed range is not silently larger/smaller), (c) spot-check tiny/boundary shapes by hand."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref2_fathooks_lemma2 import dvec, s_up, lam_of, G_coeffs, I_coeffs, mul
from math import factorial

N = int(sys.argv[1]) if len(sys.argv) > 1 else 34
bad = 0; cnt = 0; m2_diff = {}
Iser = I_coeffs(N+3)
for a in range(1, N):
    for b in range(1, N):
        for c in range(1, N):
            for d in range(1, N):
                n = a*b+b*c+c*d
                if n > N: continue
                cnt += 1
                dv = dvec(lam_of(a, b, c, d)); m = min(b, c); rect = tuple([a]*d)
                ser = mul(Iser, G_coeffs(a, d, N+3))
                for j in range(m+1):
                    if dv[j] != ser[j]*factorial(j): bad += 1; print("W2 FAIL", (a,b,c,d), j)
                s1 = ser[m+1]*factorial(m+1)
                if dv[m+1] != s1 - (b == m) - (c == m): bad += 1; print("W2' FAIL", (a,b,c,d))
                if m+2 <= n:
                    diff = ser[m+2]*factorial(m+2) - dv[m+2]
                    key = (b == m, c == m, b == m+1, c == m+1)
                    m2_diff.setdefault(key, set()).add((m, int(diff)))
print(f"n<={N}: {cnt} shapes, failures = {bad}")
print("j=m+2 deficits s_{m+2}(a^d)-d_{m+2}(lam) by (b=m,c=m,b=m+1,c=m+1), as set of (m,deficit):")
for k, v in sorted(m2_diff.items()): print("  ", k, sorted(v)[:12])
# hand boundary cases
for (a,b,c,d) in [(1,1,1,1),(2,1,1,1),(1,1,1,2),(1,2,2,1),(3,1,2,1),(1,3,1,3),(2,2,2,2)]:
    lam = lam_of(a,b,c,d); dv = dvec(lam); m = min(b,c); rect = tuple([a]*d)
    print((a,b,c,d), lam, "d=", dv[:m+3], " s_j(a^d)=", [s_up(rect,j) for j in range(m+3)], " corr=", (b==m)+(c==m))
