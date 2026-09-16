"""Referee #2: independent check of Theorem A (A1)-(A4) for all lam |- n <= NMAX, incl. boundary n < |rho|.
Own Murnaghan-Nakayama via beta-numbers (memoised); cross-checked by column orthogonality and chi(1^n)=f (n<=9).
Also: symbolic verification of the algebra in the proof of Theorem A (sympy)."""
import sys
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src')
from young import partitions, f_hook

def ff(n, i):  # falling factorial
    r = 1
    for j in range(i): r *= (n - j)
    return r

@lru_cache(maxsize=None)
def chi(lam, rho):
    """MN: lam, rho tuples (partitions), |lam| = |rho|. rho processed largest part first."""
    if not rho: return 1
    r = rho[0]; rest = rho[1:]
    L = len(lam)
    beta = [lam[i] + (L - 1 - i) for i in range(L)]  # strictly decreasing
    bset = set(beta)
    total = 0
    for i, b in enumerate(beta):
        nb = b - r
        if nb < 0 or nb in bset: continue
        # sign: number of beta's strictly between nb and b
        h = sum(1 for c in beta if nb < c < b)
        newbeta = sorted([c for c in beta if c != b] + [nb], reverse=True)
        M = len(newbeta)
        newlam = tuple(x for x in (newbeta[j] - (M - 1 - j) for j in range(M)) if x > 0)
        total += (-1) ** h * chi(newlam, rest)
    return total

def z(rho):
    from collections import Counter
    r = 1
    for k, m in Counter(rho).items(): r *= k ** m * factorial(m)
    return r

def contents(lam):
    return [j - i for i, row in enumerate(lam) for j in range(row)]

def omega(lam, rhobar):
    """omega_lam(K_rho) = |K_rho| chi^lam(rho u 1^{n-|rho|}) / f^lam ; 0 if n < |rho|."""
    n = sum(lam)
    if n < sum(rhobar): return Fraction(0)
    rho = tuple(sorted(list(rhobar) + [1] * (n - sum(rhobar)), reverse=True))
    size = Fraction(ff(n, sum(rhobar)), z(rhobar))
    return size * chi(lam, rho) / f_hook(list(lam))

def sanity(nmax=9):
    bad = 0
    for n in range(nmax + 1):
        P = [tuple(l) for l in partitions(n)]
        for lam in P:
            if chi(lam, tuple([1] * n)) != f_hook(list(lam)): bad += 1
        for rho in P:
            for tau in P:
                s = sum(chi(l, rho) * chi(l, tau) for l in P)
                if s != (z(rho) if rho == tau else 0): bad += 1
    print(f'MN sanity (chi(1^n)=f, column orthogonality) n <= {nmax}:', 'OK' if bad == 0 else f'{bad} BAD')

def A_formulas(n, C1, C2, C4):
    A1 = C2 - comb(n, 2)
    A2 = Fraction(C1 * C1, 2) - Fraction(3 * C2, 2) + comb(n, 2)
    A3 = C4 - (3 * n - 10) * C2 - 2 * C1 * C1 + Fraction(n * (n - 1) * (5 * n - 19), 6)
    A4 = Fraction(C2 * C2, 2) - Fraction(5 * C4, 2) + 3 * C1 * C1 - Fraction(n * n * C2, 2) + Fraction(13 * n * C2, 2) \
         - 15 * C2 + Fraction(n ** 4, 8) - Fraction(7 * n ** 3, 4) + Fraction(47 * n * n, 8) - Fraction(17 * n, 4)
    A4b = Fraction(1, 2) * (A1 * A1 - 2 * comb(n, 3) - (3 * n - 8) * A1 - 8 * A2 - 5 * A3)
    return A1, A2, A3, A4, A4b

def checkA(nmax=16):
    bad = 0; cnt = 0
    for n in range(0, nmax + 1):
        for lam in partitions(n):
            lam = tuple(lam)
            c = contents(lam)
            C1, C2, C4 = sum(c), sum(x * x for x in c), sum(x ** 4 for x in c)
            A1, A2, A3, A4, A4b = A_formulas(n, C1, C2, C4)
            om = [omega(lam, r) for r in [(3,), (2, 2), (5,), (3, 3)]]
            cnt += 1
            if om != [A1, A2, A3, A4] or A4 != A4b:
                bad += 1; print('MISMATCH', n, lam, om, [A1, A2, A3, A4, A4b])
        print(f'n={n} done', flush=True)
    print(f'Theorem A (A1)-(A4) for all lam |- n <= {nmax}: {cnt} shapes,', 'OK' if bad == 0 else f'{bad} BAD')

def symbolic():
    import sympy as sp
    n, C1, C2, C4 = sp.symbols('n C1 C2 C4')
    Bn2 = n * (n - 1) / 2; Bn3 = n * (n - 1) * (n - 2) / 6
    w3 = C2 - Bn2                                    # from (I2)
    w22 = (C1 ** 2 - Bn2 - 3 * w3) / 2               # from (I3)
    w5 = C4 - (3 * n - 4) * w3 - 4 * w22 - n * (n - 1) * (4 * n - 5) / 6   # from (I5)
    w33 = (w3 ** 2 - 2 * Bn3 - (3 * n - 8) * w3 - 8 * w22 - 5 * w5) / 2    # from (I6)
    A2 = C1 ** 2 / 2 - sp.Rational(3, 2) * C2 + Bn2
    A3 = C4 - (3 * n - 10) * C2 - 2 * C1 ** 2 + n * (n - 1) * (5 * n - 19) / 6
    A4 = C2 ** 2 / 2 - sp.Rational(5, 2) * C4 + 3 * C1 ** 2 - n ** 2 * C2 / 2 + sp.Rational(13, 2) * n * C2 - 15 * C2 \
         + n ** 4 / 8 - sp.Rational(7, 4) * n ** 3 + sp.Rational(47, 8) * n ** 2 - sp.Rational(17, 4) * n
    print('symbolic (A2):', sp.expand(w22 - A2) == 0)
    print('symbolic (A3):', sp.expand(w5 - A3) == 0)
    print('symbolic (A4):', sp.expand(w33 - A4) == 0)
    print('  w33 expanded =', sp.expand(w33))

if __name__ == '__main__':
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    symbolic()
    sanity(9)
    checkA(NMAX)
