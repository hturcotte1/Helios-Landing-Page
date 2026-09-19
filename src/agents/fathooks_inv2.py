import sys
from collections import defaultdict
from fathooks_local import fathooks_upto
from fathooks_invariants import hook_product, local_key, transpose_key

def content_pows(a, b, c, d, K):
    ps = [0] * (K + 1)
    for i in range(b):
        for j in range(a + c):
            for k in range(K + 1): ps[k] += (j - i) ** k
    for i in range(b, b + d):
        for j in range(c):
            for k in range(K + 1): ps[k] += (j - i) ** k
    return ps

def main(N, use):
    groups = defaultdict(set); cnt = 0
    for (a, b, c, d) in fathooks_upto(N):
        cnt += 1
        ps = content_pows(a, b, c, d, 6)
        key = [ps[0]]
        for u in use:
            if u == 'H': key.append(hook_product(a, b, c, d))
            elif u == 'loc': key.append(local_key(a, b, c, d))
            elif u == 'p1sq': key.append(ps[1] ** 2)
            elif u == 'p3sq': key.append(ps[3] ** 2)
            elif u == 'p2': key.append(ps[2])
            elif u == 'p4': key.append(ps[4])
            elif u == 'p6': key.append(ps[6])
            elif u == 'p1p3': key.append(ps[1] * ps[3])
        groups[tuple(key)].add(transpose_key(a, b, c, d))
    coll = [(k, sorted(v)) for k, v in groups.items() if len(v) > 1]
    print(f"N={N} invariants={use}: {cnt} shapes, {len(coll)} colliding keys")
    for k, v in sorted(coll)[:30]:
        print("  n=%d" % k[0], v)

if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2].split(','))
