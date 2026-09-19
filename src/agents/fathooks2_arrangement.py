"""fathooks: the arrangement statement.  For every multiset 1 <= x1 <= x2 <= x3 <= x4 <= M, the 12 arrangements
(a,b,c,d) up to reversal of the multiset are separated by (n, C_2), and also by (n, H).  Exact integers.
Usage: python3 fathooks2_arrangement.py M"""
import sys
from itertools import permutations, combinations_with_replacement
from collections import defaultdict
from fathooks2_invariants import hook_product, content_pows, tclass
def main(M):
    badC = badH = 0; cnt = 0
    for xs in combinations_with_replacement(range(1, M + 1), 4):
        arrs = set(tclass(p) for p in permutations(xs))
        cnt += 1
        gC = defaultdict(set); gH = defaultdict(set)
        for q in arrs:
            ps = content_pows(*q, 2)
            gC[(ps[0], ps[2])].add(q); gH[(ps[0], hook_product(*q))].add(q)
        for k, v in gC.items():
            if len(v) > 1: badC += 1; print("(n,C2) collision", k, sorted(v))
        for k, v in gH.items():
            if len(v) > 1: badH += 1; print("(n,H) collision", k, sorted(v))
    print(f"multisets with values <= {M}: {cnt}; (n,C2)-collisions {badC}; (n,H)-collisions {badH}")
if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 30)
