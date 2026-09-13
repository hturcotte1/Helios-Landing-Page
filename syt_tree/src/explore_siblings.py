"""Do simple invariants distinguish non-mirror siblings lam+c, lam+c'?
For all lam |- n <= N and all pairs of distinct addable corners c != c' that are not mirror images
(i.e. not (lam symmetric and c' = c^t)), compare f, (f,r), (f,C2), (r,d2,f), and the full d-vector."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from young import partitions, addable_corners, add_box, conjugate, f_hook, boxes, corner_boxes_addable, is_symmetric
from census import d_vector

def C2(lam): return sum((j - i) ** 2 for (i, j) in boxes(lam))

N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
coll = {"f": 0, "f,r": 0, "f,C2": 0, "r,d2,f": 0, "dvec": 0}
examples = {k: [] for k in coll}
npairs = 0
for n in range(1, N + 1):
    for lam in partitions(n):
        sym = is_symmetric(lam)
        kids = []
        for i in addable_corners(lam):
            mu = add_box(lam, i)
            box = (i, lam[i] if i < len(lam) else 0)
            kids.append((box, mu))
        for a in range(len(kids)):
            for b in range(a + 1, len(kids)):
                (c, mu), (c2, mu2) = kids[a], kids[b]
                if sym and c2 == (c[1], c[0]):
                    continue  # mirror pair: d-vectors equal by symmetry
                npairs += 1
                f1, f2 = f_hook(mu), f_hook(mu2)
                inv = {"f": (f1, f2), "f,r": ((f1, len(mu) and d_vector(mu)[1]), (f2, d_vector(mu2)[1])),
                       "f,C2": ((f1, C2(mu)), (f2, C2(mu2))),
                       "r,d2,f": ((d_vector(mu)[1], d_vector(mu)[2], f1), (d_vector(mu2)[1], d_vector(mu2)[2], f2)),
                       "dvec": (d_vector(mu), d_vector(mu2))}
                for k, (x, y) in inv.items():
                    if x == y:
                        coll[k] += 1
                        if len(examples[k]) < 6: examples[k].append((lam, mu, mu2))
    print(n, npairs, coll, flush=True)
for k, ex in examples.items():
    print(k, ex)
