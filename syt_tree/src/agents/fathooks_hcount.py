"""fathooks: (1) verify on exact d-vectors that for BALANCED 2-corner shapes (x4 <= x1+x2-1)
      o! [t^o] P_lam = o! [t^o] (I^2 - I*sum_x D_x) + h * C(x1+x2, x1),   o = x1+x2+1,
    where h = #{P in {{a,b},{c,d},{a,d}} : P = {x1,x2} as multisets}.
 (2) check that (sorted multiset, n, p2, h) separates balanced shapes up to transpose for n <= N2 (no d-vectors)."""
import sys, os
from fractions import Fraction
from math import factorial, comb
from collections import defaultdict
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from young import partitions, f_hook, involutions
from census import d_vector
from fathooks_local import fathook, fathooks_upto
from fathooks_invariants import content_sums, transpose_key
K = 45
Dx = {x: [sum(f_hook(r) for r in partitions(k) if r and r[0] > x) for k in range(K+1)] for x in range(K)}
Iser = [Fraction(involutions(k), factorial(k)) for k in range(K+1)]
def mul(p, q):
    r = [Fraction(0)]*(K+1)
    for a_, pa in enumerate(p):
        if pa == 0: continue
        for b_ in range(K+1-a_): r[a_+b_] += pa*q[b_]
    return r
def hcount(a,b,c,d):
    xs = sorted([a,b,c,d]); t = xs[:2]
    return sum(1 for P in ((a,b),(c,d),(a,d)) if sorted(P) == t)
def part1(N):
    bad = 0; cnt = 0
    for (a,b,c,d) in fathooks_upto(N):
        xs = sorted([a,b,c,d])
        if xs[3] > xs[0] + xs[1] - 1: continue
        cnt += 1
        o = xs[0] + xs[1] + 1
        dv = d_vector(fathook(a,b,c,d))
        sumD = [Fraction(sum(Dx[x][k] for x in (a,b,c,d)), factorial(k)) for k in range(K+1)]
        known = mul(Iser, Iser)[o] - mul(Iser, sumD)[o]
        lhs = Fraction(dv[o], factorial(o)) if o <= sum(fathook(a,b,c,d)) else Fraction(0)
        if lhs != known + Fraction(hcount(a,b,c,d)*comb(xs[0]+xs[1], xs[0]), factorial(o)):
            bad += 1; print("H-FORMULA FAIL", (a,b,c,d))
    print(f"(1) h-formula on balanced shapes n<={N}: {cnt} shapes, {'OK' if bad==0 else 'FAIL %d'%bad}")
def part2(N):
    g = defaultdict(set); cnt = 0
    for (a,b,c,d) in fathooks_upto(N):
        xs = sorted([a,b,c,d])
        if xs[3] > xs[0] + xs[1] - 1: continue
        cnt += 1
        p1, p2 = content_sums(a,b,c,d)
        g[(tuple(xs), a*b+b*c+c*d, p2, hcount(a,b,c,d))].add(transpose_key(a,b,c,d))
    coll = [(k, v) for k, v in g.items() if len(v) > 1]
    print(f"(2) balanced shapes n<={N}: {cnt} shapes; (multiset,n,p2,h) collisions: {len(coll)}")
    for k, v in coll[:10]: print("   ", k, v)
if __name__ == '__main__':
    part1(40); part2(int(sys.argv[1]) if len(sys.argv) > 1 else 400)
