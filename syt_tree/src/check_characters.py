"""Verify Lemma 4.5's explicit formulas for u_3, u_4, u_5 in terms of characters,
computing chi^lam(rho cup 1^{n-i}) = sum_{nu |- i} chi^nu(rho) f^{lam/nu} with chi^nu(rho)
from an independent Murnaghan-Nakayama implementation (border strips)."""
import os, sys
from fractions import Fraction
sys.path.insert(0, os.path.dirname(__file__))
from young import partitions, f_skew, f_hook, conjugate

def border_strips(lam, k):
    """Yield (mu, height) for all border strips of size k removable from lam (mu = lam minus strip)."""
    # use beta-numbers: lam -> set {lam_i - i} ; removing a k-strip = moving a bead from x to x-k if x-k not occupied
    ell = len(lam)
    beta = [lam[i] - i for i in range(ell)]  # distinct decreasing
    S = set(beta)
    for idx, x in enumerate(beta):
        y = x - k
        if y in S: continue
        if y < -ell: continue  # would need more rows than allowed? handle via padding below
        newbeta = sorted([b for b in beta if b != x] + [y], reverse=True)
        # height = number of beads strictly between y and x
        height = sum(1 for b in beta if y < b < x)
        # convert back: mu_i = newbeta_i + i ; must be >= 0 and weakly decreasing; allow trailing zeros
        mu = [newbeta[i] + i for i in range(ell)]
        if any(m < 0 for m in mu): continue
        while mu and mu[-1] == 0: mu.pop()
        if any(mu[i] < mu[i+1] for i in range(len(mu)-1)): continue
        yield tuple(mu), height

def chi(nu, rho):
    """Murnaghan-Nakayama: chi^nu(rho), rho a cycle type (tuple)."""
    if not rho:
        return 1 if not nu else 0
    k = rho[0]; rest = rho[1:]
    # pad nu with zeros to allow strips extending rows beyond ell? rows cannot be created by removal, fine.
    total = 0
    for mu, h in border_strips(nu, k):
        total += (-1) ** h * chi(mu, rest)
    return total

def chi_lam(lam, rho):
    """chi^lam(rho cup 1^{n-i}) via sum_nu chi^nu(rho) f^{lam/nu}."""
    i = sum(rho)
    return sum(chi(nu, rho) * f_skew(lam, nu) for nu in partitions(i))

# sanity: character table of S_3, S_4
assert chi((2,1),(3,)) == -1 and chi((2,1),(2,1)) == 0 and chi((2,1),(1,1,1)) == 2
assert chi((2,2),(2,2)) == 2 and chi((3,1),(4,)) == -1 and chi((2,2),(4,)) == 0 and chi((3,1),(2,2)) == -1
from census import d_vector
bad = 0
for n in range(5, 15):
    for lam in partitions(n):
        dv = d_vector(lam); f = dv[n]
        u3, u4, u5 = dv[n-3], dv[n-4], dv[n-5]
        c3 = chi_lam(lam, (3,)); c22 = chi_lam(lam, (2,2)); c5 = chi_lam(lam, (5,))
        ok3 = Fraction(u3) == Fraction(2,3)*f + Fraction(1,3)*c3
        ok4 = Fraction(u4) == Fraction(5,12)*f + Fraction(1,4)*c22 + Fraction(1,3)*c3
        ok5 = Fraction(u5) == Fraction(13,60)*f + Fraction(1,4)*c22 + Fraction(1,3)*c3 + Fraction(1,5)*c5
        if not (ok3 and ok4 and ok5):
            bad += 1; print("FAIL", lam, ok3, ok4, ok5)
print("u_3,u_4,u_5 character formulas checked for all 5<=n<=14:", "OK" if bad == 0 else f"{bad} failures")
