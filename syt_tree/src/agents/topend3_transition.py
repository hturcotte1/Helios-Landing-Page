"""For nu |- m with addable corners at contents c and -c (c > 0): compare
   D := f^{nu+(c)} - f^{nu+(-c)}   (equivalently the transition probabilities)  with  C_1(nu).
Records the joint sign distribution.  Also checks Kerov's formula f^{nu+x}/f^nu = (m+1) P_Y(x)/P_X'(x) exactly."""
import sys
from fractions import Fraction
from collections import Counter
from topend3_lib import *
from young import addable_corners, add_box, removable_corners
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 26
signs = Counter(); kerov_checks = 0
def corners(nu):
    X = [(nu[i] if i < len(nu) else 0) - i for i in addable_corners(nu)]
    Y = [nu[i] - 1 - i for i in removable_corners(nu)]
    return X, Y
for m in range(1, MMAX + 1):
    for nu in partitions(m):
        X, Y = corners(nu); f = f_hook(nu)
        addc = {(nu[i] if i < len(nu) else 0) - i: i for i in addable_corners(nu)}
        # Kerov formula check for every addable corner
        for x, i in addc.items():
            num = 1
            for y in Y: num *= (x - y)
            den = 1
            for x2 in X:
                if x2 != x: den *= (x - x2)
            assert Fraction(f_hook(add_box(nu, i)), f) == Fraction((m + 1) * num, den), (nu, x)
            kerov_checks += 1
        c1 = C(nu, 1)
        for c in range(1, m + 2):
            if c in addc and -c in addc:
                D = f_hook(add_box(nu, addc[c])) - f_hook(add_box(nu, addc[-c]))
                signs[(1 if D > 0 else (-1 if D < 0 else 0), 1 if c1 > 0 else (-1 if c1 < 0 else 0))] += 1
print(f"Kerov transition formula verified exactly for all addable corners of all nu |- m <= {MMAX}: {kerov_checks} checks")
print("joint distribution of (sign(f^{nu+c} - f^{nu-c}), sign(C_1(nu))) over all (nu, c) with +-c addable, m <=", MMAX)
for k, v in sorted(signs.items()): print("  ", k, v)
