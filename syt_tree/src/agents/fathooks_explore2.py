"""fathooks: (i) first differing index (bottom/top) of d-vectors for pairs colliding on cheap invariants;
(ii) does p4 help; (iii) same-multiset arrangement collisions."""
import sys, os
from collections import defaultdict
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from fathooks_local import fathook, fathooks_upto
from fathooks_invariants import hook_product, content_sums, local_key, transpose_key
from fathooks_inv2 import content_pows
from census import d_vector

def keyf(a,b,c,d, use):
    n = a*b+b*c+c*d; ps = content_pows(a,b,c,d,4); xs = sorted([a,b,c,d])
    key=[n]
    for u in use:
        if u=='p2': key.append(ps[2])
        elif u=='p1sq': key.append(ps[1]**2)
        elif u=='p4': key.append(ps[4])
        elif u=='m': key.append((xs[0], xs.count(xs[0])))
        elif u=='loc': key.append(local_key(a,b,c,d))
        elif u=='ms': key.append(tuple(xs))
    return tuple(key)

def collisions(N, use):
    g = defaultdict(set)
    for q in fathooks_upto(N): g[keyf(*q, use)].add(transpose_key(*q))
    return [(k, sorted(v)) for k, v in g.items() if len(v) > 1]

if __name__ == '__main__':
    N = int(sys.argv[1])
    for use in (['p2','p1sq','p4'], ['p2','p1sq','p4','loc'], ['p2','p1sq','p4','m'], ['ms'], ['ms','p2'], ['ms','p2','p1sq']):
        coll = collisions(N, use)
        print(f"N={N} {use}: {len(coll)} colliding keys")
        for k, v in sorted(coll)[:12]: print("   n=%d" % k[0], v)
    # first differing index for (n,p2,p1sq,loc) collisions, only where d-vector computable (n<=45)
    print("first-diff analysis for (n,p2,p1sq,m)-collisions with n<=44:")
    for k, v in sorted(collisions(44, ['p2','p1sq','m'])):
        n = k[0]
        dvs = [d_vector(fathook(*q)) for q in v]
        for i in range(1, len(v)):
            diff = [j for j in range(n+1) if dvs[0][j] != dvs[i][j]]
            print("   n=%d %s vs %s: first diff j=%d (x1+x2=%d, a+d=%d/%d), last diff j=%d (n-j=%d), #diff=%d" % (
                n, v[0], v[i], diff[0], sorted(v[0])[0]+sorted(v[0])[1], v[0][0]+v[0][3], v[i][0]+v[i][3], diff[-1], n-diff[-1], len(diff)))
