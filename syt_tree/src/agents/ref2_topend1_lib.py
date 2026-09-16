"""ref2_topend1: independent library. Murnaghan-Nakayama by beta-numbers, contents, f^lam, partitions.
All exact integer arithmetic."""
from functools import lru_cache
from math import factorial, comb
from fractions import Fraction

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

def beta(lam):
    l = len(lam)
    return tuple(sorted((lam[i] + (l - 1 - i) for i in range(l)), reverse=True))

@lru_cache(maxsize=None)
def chi(lam, rho):
    """chi^lam(rho), lam, rho partitions (tuples) of the same n; rho may be given with 1's or padded."""
    if sum(lam) != sum(rho): raise ValueError
    if not lam: return 1
    k = rho[0]; rest = rho[1:]
    B = list(beta(lam)); S = set(B)
    tot = 0
    for b in B:
        nb = b - k
        if nb < 0 or nb in S: continue
        sign = (-1) ** sum(1 for x in B if nb < x < b)
        newB = sorted([x for x in B if x != b] + [nb], reverse=True)
        l = len(newB)
        newlam = tuple(newB[i] - (l - 1 - i) for i in range(l))
        newlam = tuple(x for x in newlam if x > 0)
        tot += sign * chi(newlam, rest)
    return tot

def f_hook(lam):
    n = sum(lam); conj = conjugate(lam)
    h = 1
    for i, li in enumerate(lam):
        for j in range(li):
            h *= (li - j) + (conj[j] - i) - 1
    return factorial(n) // h

def conjugate(lam):
    if not lam: return ()
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0]))

def C(lam, k):
    return sum((j - i) ** k for i, li in enumerate(lam) for j in range(li))

def z(rho):
    from collections import Counter
    r = 1
    for k, m in Counter(rho).items():
        r *= k ** m * factorial(m)
    return r

def class_size(rho, n):
    """|K_rho(n)| = n!/z(rho u 1^{n-|rho|})"""
    s = sum(rho)
    if n < s: return 0
    full = tuple(sorted(rho + (1,) * (n - s), reverse=True))
    return factorial(n) // z(full)

def omega(lam, rho):
    """omega_lam(K_rho) = |K_rho| chi^lam(rho u 1^{n-|rho|}) / f^lam, exact Fraction; 0 if n<|rho|"""
    n = sum(lam); s = sum(rho)
    if n < s: return Fraction(0)
    full = tuple(sorted(rho + (1,) * (n - s), reverse=True))
    return Fraction(class_size(rho, n) * chi(lam, full), f_hook(lam))

def even_type_no_ones(maxsize):
    from collections import Counter
    res = []
    for s in range(0, maxsize + 1):
        for rho in partitions(s):
            if 1 in rho: continue
            c = Counter(rho)
            if all(m % 2 == 0 for k, m in c.items() if k % 2 == 0):
                res.append(rho)
    return res

if __name__ == "__main__":
    # sanity: column orthogonality and f = chi(1^n) for n<=8; hook formula vs chi
    for n in range(0, 9):
        lams = list(partitions(n))
        for lam in lams:
            assert chi(lam, (1,) * n) == f_hook(lam), lam
        # row orthogonality sum_rho |K_rho| chi^lam chi^mu = n! delta
        for a in lams:
            for b in lams:
                s = sum(class_size(tuple(x for x in rho if x > 1), n) * chi(a, rho) * chi(b, rho) for rho in lams)
                assert s == (factorial(n) if a == b else 0)
    print("MN sanity OK n<=8; even-type rho |rho|<=10:", even_type_no_ones(10), len(even_type_no_ones(10)))
