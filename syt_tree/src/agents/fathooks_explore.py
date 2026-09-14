"""fathooks angle (second pass): collision structure of cheap invariants on 2-corner shapes.
key options: n, p2, p1sq, m (min & multiplicity), loc (N_k for k <= x1+x2-1), H (hook product), s (a+d), u (b+c)
"""
import sys, os
from collections import defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from fathooks_local import fathook, fathooks_upto
from fathooks_invariants import hook_product, content_sums, local_key, transpose_key

def key_of(a, b, c, d, use):
    n = a*b + b*c + c*d
    p1, p2 = content_sums(a, b, c, d)
    xs = sorted([a, b, c, d])
    key = [n]
    for u in use:
        if u == 'p2': key.append(p2)
        elif u == 'p1sq': key.append(p1*p1)
        elif u == 'm': key.append((xs[0], xs.count(xs[0])))
        elif u == 'loc': key.append(local_key(a, b, c, d))
        elif u == 'H': key.append(hook_product(a, b, c, d))
        elif u == 's': key.append(a + d)
        elif u == 'u': key.append(b + c)
        elif u == 'ad': key.append(a * d)
        elif u == 'bc': key.append(b * c)
    return tuple(key)

def main(N, use, show=40):
    groups = defaultdict(set); cnt = 0
    for (a, b, c, d) in fathooks_upto(N):
        cnt += 1
        groups[key_of(a, b, c, d, use)].add(transpose_key(a, b, c, d))
    coll = [(k, sorted(v)) for k, v in groups.items() if len(v) > 1]
    print(f"N={N} invariants={use}: {cnt} shapes, {len(coll)} colliding keys, {sum(len(v) for _,v in coll)} shapes involved")
    for k, v in sorted(coll)[:show]:
        print("  n=%d" % k[0], v)
    return coll

if __name__ == '__main__':
    N = int(sys.argv[1]); use = sys.argv[2].split(',') if len(sys.argv) > 2 else []
    main(N, use)
