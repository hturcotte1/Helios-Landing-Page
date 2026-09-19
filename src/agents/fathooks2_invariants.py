"""fathooks: which PROVED-determined closed-form invariants separate 2-corner shapes lam(a,b,c,d) = ((a+c)^b, c^d)
up to transpose (a,b,c,d) -> (d,c,b,a)?  No d-vectors are computed here; everything is closed form (exact integers):
  n = ab+bc+cd;  H = product of hook lengths (H <-> f^lam given n);  C_k = sum of k-th powers of contents;
  LOC = (x1, x2, (N_k)_{k <= x1+x2-1}) the local multiset data of Theorem 3;  N1 = #unit parameters.
Prints, for each invariant set, the number of transpose classes of 2-corner shapes with n <= N that collide, and the
first few collisions.  Usage: python3 fathooks2_invariants.py N
"""
import sys
from collections import defaultdict

def fathooks_upto(N):
    for b in range(1, N + 1):
        for d in range(1, N + 1):
            for c in range(1, N + 1):
                if b * c + c * d > N: break
                for a in range(1, N + 1):
                    if a * b + b * c + c * d > N: break
                    yield (a, b, c, d)

def hook_product(a, b, c, d):
    H = 1
    for i in range(b):
        for j in range(a): H *= (a + b - 1 - i - j)
    for i in range(d):
        for j in range(c): H *= (c + d - 1 - i - j)
    S = a + b + c + d
    for i in range(b):
        for j in range(c): H *= (S - 1 - i - j)
    return H

def content_pows(a, b, c, d, K=4):
    ps = [0] * (K + 1)
    for i in range(b):
        for j in range(a + c):
            for k in range(K + 1): ps[k] += (j - i) ** k
    for i in range(b, b + d):
        for j in range(c):
            for k in range(K + 1): ps[k] += (j - i) ** k
    return ps

def loc(a, b, c, d):
    xs = sorted([a, b, c, d]); x1, x2 = xs[0], xs[1]
    return (x1, x2, tuple(xs.count(k) for k in range(x1, x1 + x2)))

def tclass(q): return min(q, q[::-1])

KEYS = {
 'n,H': lambda a,b,c,d,ps: (hook_product(a,b,c,d),),
 'n,C2': lambda a,b,c,d,ps: (ps[2],),
 'n,H,C2': lambda a,b,c,d,ps: (hook_product(a,b,c,d), ps[2]),
 'n,C2,C1^2,C4': lambda a,b,c,d,ps: (ps[2], ps[1]**2, ps[4]),
 'n,H,C2,C1^2,C4': lambda a,b,c,d,ps: (hook_product(a,b,c,d), ps[2], ps[1]**2, ps[4]),
 'n,LOC': lambda a,b,c,d,ps: (loc(a,b,c,d),),
 'n,LOC,H': lambda a,b,c,d,ps: (loc(a,b,c,d), hook_product(a,b,c,d)),
 'n,LOC,C2': lambda a,b,c,d,ps: (loc(a,b,c,d), ps[2]),
 'n,LOC,H,C2': lambda a,b,c,d,ps: (loc(a,b,c,d), hook_product(a,b,c,d), ps[2]),
 'n,LOC,H,C2,C1^2,C4': lambda a,b,c,d,ps: (loc(a,b,c,d), hook_product(a,b,c,d), ps[2], ps[1]**2, ps[4]),
 'n,multiset': lambda a,b,c,d,ps: (tuple(sorted([a,b,c,d])),),
 'n,multiset,H': lambda a,b,c,d,ps: (tuple(sorted([a,b,c,d])), hook_product(a,b,c,d)),
 'n,multiset,C2': lambda a,b,c,d,ps: (tuple(sorted([a,b,c,d])), ps[2]),
}

def main(N):
    shapes = list(fathooks_upto(N))
    data = [(q, content_pows(*q)) for q in shapes]
    print(f"{len(shapes)} 2-corner shapes with n <= {N}")
    for name, kf in KEYS.items():
        g = defaultdict(set)
        for q, ps in data:
            g[(ps[0],) + kf(*q, ps)].add(tclass(q))
        colls = [(k, sorted(v)) for k, v in g.items() if len(v) > 1]
        colls.sort(key=lambda kv: kv[0][0])
        nmin = colls[0][0][0] if colls else None
        print(f"{name:24s}: {len(colls):6d} colliding keys; first collision at n = {nmin}; e.g. {colls[:3]}")

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 120)
