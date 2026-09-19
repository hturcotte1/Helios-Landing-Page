"""fathooks angle: which closed-form invariants separate 2-corner shapes up to transpose?
Invariants (all PROVED to be determined by the d-vector, see report):
  n, H = product of hook lengths (= n!/f), p2 = sum of squared contents, p1sq = (sum of contents)^2,
  local data: x1 = min param, x2 = second smallest, N_k = #{params = k} for k <= x1+x2-1.
"""
import os, sys
from math import comb
from collections import defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from fathooks_local import fathook, fathooks_upto

def hook_product(a, b, c, d):
    # rows < b, cols >= c : b x a rectangle hooks;  rows >= b, cols < c: d x c rectangle;
    # rows < b, cols < c : hook = a+b+c+d-1-i-j
    H = 1
    for i in range(b):
        for j in range(a):
            H *= (a + b - 1 - i - j)
    for i in range(d):
        for j in range(c):
            H *= (c + d - 1 - i - j)
    S = a + b + c + d
    for i in range(b):
        for j in range(c):
            H *= (S - 1 - i - j)
    return H

def content_sums(a, b, c, d):
    p1 = p2 = 0
    for i in range(b):
        for j in range(a + c):
            p1 += j - i; p2 += (j - i) ** 2
    for i in range(b, b + d):
        for j in range(c):
            p1 += j - i; p2 += (j - i) ** 2
    return p1, p2

def p1_formula(a, b, c, d):
    return (b * (a + c) * (a + c - b) + c * d * (c - 2 * b - d)) // 2

def local_key(a, b, c, d):
    xs = sorted([a, b, c, d])
    x1, x2 = xs[0], xs[1]
    return (x1, x2, tuple(sum(1 for x in xs if x == k) for k in range(x1, x1 + x2)))

def transpose_key(a, b, c, d):
    return min((a, b, c, d), (d, c, b, a))

def main(N, use):
    groups = defaultdict(set)
    cnt = 0
    for (a, b, c, d) in fathooks_upto(N):
        cnt += 1
        n = a * b + b * c + c * d
        p1, p2 = content_sums(a, b, c, d)
        assert 2 * p1 == b * (a + c) * (a + c - b) + c * d * (c - 2 * b - d)
        key = [n]
        if 'H' in use: key.append(hook_product(a, b, c, d))
        if 'p2' in use: key.append(p2)
        if 'p1sq' in use: key.append(p1 * p1)
        if 'loc' in use: key.append(local_key(a, b, c, d))
        groups[tuple(key)].add(transpose_key(a, b, c, d))
    coll = [(k, sorted(v)) for k, v in groups.items() if len(v) > 1]
    print(f"N={N} invariants={use}: {cnt} shapes, {len(coll)} colliding keys")
    for k, v in sorted(coll)[:40]:
        print("  n=%d" % k[0], v)
    return coll

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    use = sys.argv[2].split(',') if len(sys.argv) > 2 else ['H', 'p2', 'p1sq', 'loc']
    main(N, use)
