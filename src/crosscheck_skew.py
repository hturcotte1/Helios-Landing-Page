"""Independent cross-check of the d-vector recursion against the definition
d_j = sum_{nu subset lam, |lam/nu| = j} f^{lam/nu} (skew chain counting), for all
partitions of n <= 14 and a random sample of partitions of n = 20..24."""
import os, sys, random
sys.path.insert(0, os.path.dirname(__file__))
from young import partitions
from census import d_vector, d_vector_by_skew
random.seed(1)
bad = 0; cnt = 0
for n in range(0, 15):
    for lam in partitions(n):
        cnt += 1
        if d_vector(lam) != d_vector_by_skew(lam): bad += 1; print("MISMATCH", lam)
print(f"all partitions n<=14: {cnt} checked, {bad} mismatches")
for n in range(20, 25):
    ps = list(partitions(n)); random.shuffle(ps)
    for lam in ps[:5]:
        cnt += 1
        if d_vector(lam) != d_vector_by_skew(lam): bad += 1; print("MISMATCH", lam)
        else: print("ok", lam)
print(f"total {cnt} checked, {bad} mismatches")
