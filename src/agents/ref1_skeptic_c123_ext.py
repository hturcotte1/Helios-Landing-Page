"""ref1_skeptic_c123_ext.py -- extension: C2 for all 3-row nu with |nu| <= 16 (indeed nu_1 <= 16), j <= 24;
C1 direct on a larger box (lam in [-2,8]^3 dominant, mu in [-1,6]^3 dominant, only lam with lam_3 = -2..8 sampled every case)."""
import os
import sys, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref1_skeptic_c123 import Nch, char_GL3, decompose, pmul, dual_weight, sgn3, log
bad = cnt = 0
for n1 in range(0, 17):
    for n2 in range(0, n1+1):
        for n3 in range(0, n2+1):
            nu = (n1, n2, n3); dual = (n1-n3, n1-n2, 0)
            for j in range(0, 25):
                cnt += 1
                if Nch(nu, j) != Nch(dual, j): bad += 1; log("FAIL", nu, j)
log("EXT: C2 for all nu with nu_1 <= 16 (covers all 3-row partitions of n <= 16), j <= 24: %d cases, %d failures" % (cnt, bad))
bad = cnt = 0
doms = [l for l in itertools.product(range(-2, 9), repeat=3) if l[0] >= l[1] >= l[2] and l[0] - l[2] <= 8]
mus = [m for m in itertools.product(range(-1, 7), repeat=3) if m[0] >= m[1] >= m[2] and m[0] - m[2] <= 6]
for lam in doms:
    chl = char_GL3(lam)
    for mu in mus:
        chw = char_GL3(mu); chwd = char_GL3(dual_weight(mu))
        n1 = sum(decompose(pmul(chl, chw)).values()); n2 = sum(decompose(pmul(chl, chwd)).values())
        cnt += 1
        if n1 != n2: bad += 1; log("C1 FAIL", lam, mu, n1, n2)
log("EXT: C1 direct, lam in [-2,8]^3 dominant (width<=8), mu in [-1,6]^3 dominant (width<=6): %d cases, %d failures" % (cnt, bad))
