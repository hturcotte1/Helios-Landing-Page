"""Hypothesis test: for lam |- n, a removable corner y of content -c and an addable corner x of content +c
(c != 0), mu := lam - y + x.  Is  f^mu == f^lam  <=>  C_1(lam) == -c ?   All lam |- n <= NMAX."""
import sys
from topend3_lib import *
from young import removable_corners, addable_corners, remove_box, add_box, hook_lengths
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
def hp(lam):
    p = 1
    for h in hook_lengths(lam).values(): p *= h
    return p
tested = 0; eq_f = 0; eq_c = 0; both = 0; bad = []
for n in range(2, NMAX + 1):
    for lam in partitions(n):
        c1 = C(lam, 1)
        rem = {(lam[i] - 1 - i): i for i in removable_corners(lam)}   # content -> row
        addc = {(lam[i] if i < len(lam) else 0) - i: i for i in addable_corners(lam)}
        for cy, i in rem.items():
            c = -cy
            if c == 0 or c not in addc: continue
            mu = add_box(remove_box(lam, i), addc[c] if addc[c] < len(remove_box(lam, i)) + 1 else addc[c])
            # careful: after removing y the addable corners may shift; x (content c) is still addable unless x was created/destroyed
            # by the removal; check directly:
            lam2 = remove_box(lam, i)
            j = addc[c]
            if j not in addable_corners(lam2): continue
            mu = add_box(lam2, j)
            tested += 1
            e1 = (hp(mu) == hp(lam)); e2 = (c1 == -c)
            eq_f += e1; eq_c += e2; both += (e1 and e2)
            if e1 != e2: bad.append((lam, c, mu, e1, e2))
print(f"n <= {NMAX}: single-box mirror moves tested: {tested}; f equal: {eq_f}; C_1 = -c: {eq_c}; both: {both}; mismatches: {len(bad)}")
for b in bad:
    if b[3]: print("  f equal but C_1 != -c:", b)
