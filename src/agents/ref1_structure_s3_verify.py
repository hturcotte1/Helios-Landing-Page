"""Referee #1: independent verification of Theorem S3 and Corollaries S3.1, S3.2 (structure.md section 2).
No imports from the project. Exact integer arithmetic only.
 (A) own Murnaghan-Nakayama characters (rim hooks by row-span formula), self-checked by column orthogonality,
     chi(1^n) = f (hook formula), and chi^{lam^t} = sgn * chi^lam.
 (B) own d-vector by the removal recursion d_j = sum_c d_{j-1}(lam - c), cross-checked (n<=9) with the
     order-filter definition sum over nu of f^{lam/nu}.
 (C) Theorem S3 by brute force over all permutations of S_n (n <= NB) : n! d_j / j! = sum_tau chi(tau^2) C(fix,j).
 (D) Theorem S3 by conjugacy classes for all lam |- n <= NC, both forms (coefficients of (1+t)^fix, and d_j form).
 (E) Cor S3.1: n! P(-1) = sum over derangements of chi(tau^2); = D_n for (n), = (-1)^n (n-1) for (n-1,1), n <= NC.
 (F) Cor S3.2: delta_m = sum_{Der_m} chi(tau^2 cup 1^{n-m}); i! u_i = sum_m C(i,m) delta_m, all lam |- n <= NC;
     delta_0=delta_2=f, delta_1=0, delta_3 = 2 chi(3,1^{n-3}), delta_4 = 6 chi(2,2,1^{n-4}) + 3f;
     also delta_m by brute force over derangement permutations for m <= 8 vs. class formula.
 (G) Hermite/involution expansion (iv): P(t) = sum_m w_m Hhat_m(1+t), n <= NG.
"""
import sys, itertools, time
from fractions import Fraction
from functools import lru_cache
from math import factorial, comb
from collections import Counter

def partitions(n, maxp=None):
    if maxp is None or maxp > n: maxp = n
    if n == 0: yield (); return
    for a in range(maxp, 0, -1):
        for rest in partitions(n-a, a): yield (a,)+rest

def conjugate(lam):
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()

def f_hook(lam):
    n = sum(lam); lt = conjugate(lam); prod = 1
    for i, l in enumerate(lam):
        for j in range(l):
            prod *= (l - j) + (lt[j] - i) - 1
    assert factorial(n) % prod == 0
    return factorial(n)//prod

def zee(rho):
    z = 1
    for k, m in Counter(rho).items(): z *= k**m * factorial(m)
    return z

# ---------- (A) characters ----------
def rim_hooks(lam, k):
    """yield (mu, height) for each border strip of size k removable from lam. Strip spans rows i..i2 (0-indexed)."""
    L = len(lam)
    lamp = list(lam) + [0]
    for i in range(L):
        for i2 in range(i, L):
            # mu_r = lam_{r+1} - 1 for i <= r < i2 ; mu_{i2} = lam_i + (i2 - i) - k
            mu = list(lam)
            ok = True
            for r in range(i, i2):
                mu[r] = lamp[r+1] - 1
                if mu[r] < 0: ok = False; break
            if not ok: continue
            last = lam[i] + (i2 - i) - k
            if last < lamp[i2+1] or last > lam[i2] - 1: continue
            mu[i2] = last
            # for r in [i, i2): need mu[r] >= mu[r+1]: lam_{r+1}-1 >= lam_{r+2}-1 ok; and mu[i2-1] = lam_{i2}-1 >= last ok
            while mu and mu[-1] == 0: mu.pop()
            yield tuple(mu), i2 - i

@lru_cache(maxsize=None)
def chi(lam, rho):
    """Murnaghan-Nakayama, rho a tuple (any order); remove parts from the largest first."""
    if not rho: return 1 if not lam else 0
    rho = tuple(sorted(rho, reverse=True))
    k, rest = rho[0], rho[1:]
    return sum((-1)**h * chi(mu, rest) for mu, h in rim_hooks(lam, k))

def check_characters(N):
    for n in range(1, N+1):
        parts = list(partitions(n))
        for lam in parts:
            assert chi(lam, (1,)*n) == f_hook(lam), ("f", lam)
            lt = conjugate(lam)
            for rho in parts:
                sgn = (-1)**(n - len(rho))
                assert chi(lt, rho) == sgn*chi(lam, rho), ("transpose", lam, rho)
        # column orthogonality: sum_rho (n!/z_rho) chi^lam(rho) chi^mu(rho) = n! delta
        for a in range(len(parts)):
            for b in range(a, len(parts)):
                s = sum(factorial(n)//zee(rho) * chi(parts[a], rho)*chi(parts[b], rho) for rho in parts)
                assert s == (factorial(n) if a == b else 0), ("orth", parts[a], parts[b])
    print(f"(A) own MN characters pass hook/transpose/orthogonality checks for n <= {N}")

# ---------- (B) d-vector ----------
def children(lam):
    out = []
    for i in range(len(lam)):
        if i == len(lam)-1 or lam[i] > lam[i+1]:
            mu = list(lam); mu[i] -= 1
            if mu[i] == 0: mu.pop()
            out.append(tuple(mu))
    return out

@lru_cache(maxsize=None)
def dvec(lam):
    n = sum(lam)
    if n == 0: return (1,)
    d = [0]*(n+1); d[0] = 1
    for mu in children(lam):
        dm = dvec(mu)
        for j in range(1, n+1): d[j] += dm[j-1]
    return tuple(d)

def contains(lam, nu):
    return len(nu) <= len(lam) and all(nu[i] <= lam[i] for i in range(len(nu)))

@lru_cache(maxsize=None)
def f_skew(lam, nu):
    if not contains(lam, nu): return 0
    if lam == nu: return 1
    return sum(f_skew(mu, nu) for mu in children(lam))

def check_dvec(N):
    for n in range(0, N+1):
        for lam in partitions(n):
            d = dvec(lam)
            alt = [0]*(n+1)
            for m in range(n+1):
                for nu in partitions(m):
                    alt[n-m] += f_skew(lam, nu)
            assert list(d) == alt, ("dvec", lam)
    print(f"(B) own d-vector recursion agrees with sum_nu f^(lam/nu) for n <= {N}")

# ---------- helpers ----------
def cycle_type(p):
    n = len(p); seen = [False]*n; ct = []
    for i in range(n):
        if seen[i]: continue
        L = 0; j = i
        while not seen[j]: seen[j] = True; j = p[j]; L += 1
        ct.append(L)
    return tuple(sorted(ct, reverse=True))

def sq_type(rho):
    out = []
    for k in rho:
        if k % 2: out.append(k)
        else: out += [k//2, k//2]
    return tuple(sorted(out, reverse=True))

def pad(rho, n):
    return tuple(sorted(list(rho) + [1]*(n - sum(rho)), reverse=True))

# ---------- (C) brute force over S_n ----------
def check_S3_brute(NB):
    for n in range(0, NB+1):
        # aggregate over permutations: count[(type of tau^2, fix)]
        agg = Counter()
        for p in itertools.permutations(range(n)):
            sq = tuple(p[p[i]] for i in range(n))
            fx = sum(1 for i in range(n) if p[i] == i)
            agg[(cycle_type(sq), fx)] += 1
        # sanity: types of squares from the class formula must agree with brute force
        agg2 = Counter()
        for rho in partitions(n):
            agg2[(sq_type(rho), rho.count(1))] += factorial(n)//zee(rho)
        assert agg == agg2, ("square types", n)
        for lam in partitions(n):
            d = dvec(lam)
            for j in range(n+1):
                rhs = sum(c*chi(lam, ty)*comb(fx, j) for (ty, fx), c in agg.items())
                assert rhs*factorial(j) == factorial(n)*d[j], ("S3 brute", lam, j)
    print(f"(C) Theorem S3 brute force over all permutations of S_n, n <= {NB}: OK")

# ---------- (D) class-based S3 ----------
def poly_1pt_pow(k, n):
    """coefficients of (1+t)^k as list length n+1"""
    return [comb(k, j) for j in range(n+1)]

def check_S3_classes(NC):
    for n in range(0, NC+1):
        parts = list(partitions(n))
        cls = [(factorial(n)//zee(rho), sq_type(rho), rho.count(1)) for rho in parts]
        for lam in parts:
            d = dvec(lam)
            lhs = [factorial(n)*d[j]//factorial(j) for j in range(n+1)]   # n! P(t) coefficients
            assert all(factorial(n)*d[j] % factorial(j) == 0 for j in range(n+1))
            rhs = [0]*(n+1)
            for size, ty, fx in cls:
                c = chi(lam, ty)
                if c == 0: continue
                for j in range(fx+1): rhs[j] += size*c*comb(fx, j)
            assert lhs == rhs, ("S3 classes", lam, lhs, rhs)
    print(f"(D) Theorem S3 by conjugacy classes, all lam |- n <= {NC}: OK")

# ---------- (E) derangement value ----------
def D(n):
    return sum((-1)**k * factorial(n)//factorial(k) for k in range(n+1))

def check_S31(NC):
    for n in range(0, NC+1):
        parts = list(partitions(n))
        for lam in parts:
            d = dvec(lam)
            val = sum((-1)**j*factorial(n)*d[j]//factorial(j) for j in range(n+1))
            der = sum(factorial(n)//zee(rho)*chi(lam, sq_type(rho)) for rho in parts if 1 not in rho)
            assert val == der, ("S3.1", lam, val, der)
            if lam == (n,): assert val == D(n), ("D_n", n, val)
            if n >= 2 and lam == tuple(x for x in (n-1, 1) if x): assert val == (-1)**n*(n-1), ("(n-1,1)", n, val)
    print(f"(E) Cor S3.1 derangement value, all lam |- n <= {NC}; (n) -> D_n and (n-1,1) -> (-1)^n (n-1): OK")
    # print table for hooks and two-row shapes n<=12 to compare with structure_extra (10)
    for n in range(2, 13):
        v = lambda lam: sum((-1)**j*factorial(n)*dvec(lam)[j]//factorial(j) for j in range(n+1))
        hooks = [v(tuple([n-k]+[1]*k)) for k in range(n)]
        rows = [v(tuple(x for x in (n-k, k) if x)) for k in range(n//2+1)]
        print("   n!P(-1)", n, "hooks:", hooks, "two-row:", rows)

# ---------- (F) derangement transform ----------
def delta_classes(lam, m):
    n = sum(lam)
    return sum(factorial(m)//zee(rho)*chi(lam, pad(sq_type(rho), n)) for rho in partitions(m) if 1 not in rho)

def delta_brute(lam, m):
    n = sum(lam); tot = 0
    for p in itertools.permutations(range(m)):
        if any(p[i] == i for i in range(m)): continue
        sq = tuple(p[p[i]] for i in range(m))
        tot += chi(lam, pad(cycle_type(sq), n))
    return tot

def check_S32(NC, MB):
    for n in range(0, NC+1):
        for lam in partitions(n):
            d = dvec(lam); f = d[n]
            delta = [delta_classes(lam, m) for m in range(n+1)]
            for i in range(n+1):
                u = d[n-i]
                assert factorial(i)*u == sum(comb(i, m)*delta[m] for m in range(i+1)), ("S3.2", lam, i)
            # polynomial form: n! P(s-1) = sum_m C(n,m) delta_m s^{n-m}; check at s = 0..n+1
            for s in range(0, n+2):
                lhsv = sum(factorial(n)*d[j]//factorial(j)*(s-1)**j for j in range(n+1))
                rhsv = sum(comb(n, m)*delta[m]*s**(n-m) for m in range(n+1))
                assert lhsv == rhsv, ("S3.2 poly", lam, s)
            assert delta[0] == f
            if n >= 1: assert delta[1] == 0
            if n >= 2: assert delta[2] == f
            if n >= 3: assert delta[3] == 2*chi(lam, pad((3,), n)), ("delta3", lam)
            if n >= 4: assert delta[4] == 6*chi(lam, pad((2, 2), n)) + 3*f, ("delta4", lam)
            if n <= MB:
                for m in range(n+1):
                    assert delta[m] == delta_brute(lam, m), ("delta brute", lam, m)
    print(f"(F) Cor S3.2 derangement transform + delta_0..delta_4 formulas, all lam |- n <= {NC} (brute-force deltas for n <= {MB}): OK")

# ---------- (G) involution expansion ----------
def invol_poly(m):
    # Hhat_m(x) = sum_i C(m,2i)(2i-1)!! x^{m-2i}, as coefficient list in x
    c = [0]*(m+1)
    for i in range(m//2+1):
        c[m-2*i] = comb(m, 2*i)*factorial(2*i)//(2**i*factorial(i))
    return c

def check_G(NG):
    for n in range(0, NG+1):
        for lam in partitions(n):
            d = dvec(lam)
            P = [Fraction(d[j], factorial(j)) for j in range(n+1)]
            Q = [Fraction(0)]*(n+1)
            for m in range(n+1):
                wm = Fraction(sum(factorial(n-m)//zee(rho)*chi(lam, pad(sq_type(rho), n))
                                  for rho in partitions(n-m) if not rho or rho[-1] >= 3), factorial(m)*factorial(n-m))
                if wm == 0: continue
                H = invol_poly(m)
                for k, hk in enumerate(H):
                    for j in range(k+1): Q[j] += wm*hk*comb(k, j)
            assert P == Q, ("G", lam)
    print(f"(G) involution/Hermite expansion P = sum_m w_m Hhat_m(1+t), all lam |- n <= {NG}: OK")

if __name__ == "__main__":
    NC = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    t0 = time.time()
    check_characters(9); print(" ", round(time.time()-t0, 1), "s")
    check_dvec(9); print(" ", round(time.time()-t0, 1), "s")
    check_S3_brute(8); print(" ", round(time.time()-t0, 1), "s")
    check_S3_classes(NC); print(" ", round(time.time()-t0, 1), "s")
    check_S31(NC); print(" ", round(time.time()-t0, 1), "s")
    check_S32(NC, 8); print(" ", round(time.time()-t0, 1), "s")
    check_G(min(NC, 12)); print(" ", round(time.time()-t0, 1), "s")
    print("ALL OK")
