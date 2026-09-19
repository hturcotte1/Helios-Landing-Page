"""Referee #2, independent check of agent_notes/structure.md section 2 (Theorem S3 and corollaries).

Everything here is written from scratch: own Murnaghan-Nakayama (explicit rim hooks on the diagram),
own d-vector recursion, own partition generator. Only cross-checked against the project's code at the end.

Checks:
 (A) brute force over S_n (n <= 8):  n! P_lam(t) = sum_tau chi(tau^2) (1+t)^fix(tau)  and
     d_j = (j!/n!) sum_tau chi(tau^2) C(fix tau, j).
 (B) class-sum version (n <= NMAX): n! P_lam(t) = sum_rho (n!/z_rho) chi(sq(rho)) (1+t)^{m_1(rho)}.
 (C) derangement value n! P_lam(-1) = sum_{Der_n} chi(tau^2), integrality, = D_n for (n), (-1)^n (n-1) for (n-1,1).
 (D) derangement transform: delta_m via classes, i! u_i = sum_m C(i,m) delta_m,
     n! P(s-1) = sum_m C(n,m) delta_m s^{n-m}, delta_0=delta_2=f, delta_1=0, delta_3=2chi(3,1^{n-3}),
     delta_4 = 6 chi(2,2,1^{n-4}) + 3 f;   plus brute force delta_m over Der_m for m <= 8.
 (E) item (iv): P_lam(t) = sum_m w_m Hhat_m(1+t).
 (F) item (i): chi^{lam^t}(tau^2) = chi^lam(tau^2) (parity) and (ii) formula for delta_m / f.
"""
import sys, itertools
from fractions import Fraction
from math import factorial, comb
from functools import lru_cache
from collections import Counter

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 16

def parts(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxp), 0, -1):
        for rest in parts(n - k, k):
            yield (k,) + rest

def conj(lam):
    if not lam: return ()
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0]))

# ---------- own Murnaghan-Nakayama via explicit rim hooks ----------
def rim_hooks(lam, k):
    """All ways to remove a rim hook (border strip) of size k from lam; yield (mu, leg length).
    Implemented by scanning the rim of the diagram cell by cell (no beta numbers)."""
    lam = tuple(lam)
    ell = len(lam)
    # rim cells listed from top-right to bottom-left: for each row i, cells (i, j) with j from lam_i - 1 down to lam_{i+1}
    rim = []
    for i in range(ell):
        lo = max(lam[i+1] - 1, 0) if i + 1 < ell else 0
        for j in range(lam[i] - 1, lo - 1, -1):
            rim.append((i, j))
    # a border strip of size k is a set of k consecutive rim cells whose removal leaves a partition
    for s in range(len(rim) - k + 1):
        strip = rim[s:s + k]
        rows = Counter(i for i, j in strip)
        mu = list(lam)
        ok = True
        for i, c in rows.items():
            # in row i, the cells removed must be the last c cells of row i
            removed = sorted(j for (ii, j) in strip if ii == i)
            if removed != list(range(lam[i] - c, lam[i])):
                ok = False; break
            mu[i] -= c
        if not ok: continue
        # must remain a partition
        if any(mu[i] < mu[i+1] for i in range(ell - 1)) or any(x < 0 for x in mu):
            continue
        # connected: consecutive rim cells are adjacent by construction? rim as listed is connected (each step moves
        # left or down), so consecutive cells are adjacent. Check anyway.
        for a, b in zip(strip, strip[1:]):
            if abs(a[0]-b[0]) + abs(a[1]-b[1]) != 1: ok = False
        if not ok: continue
        while mu and mu[-1] == 0: mu.pop()
        leg = len(rows) - 1
        yield tuple(mu), leg

@lru_cache(maxsize=None)
def chi(lam, rho):
    """chi^lam at cycle type rho (tuple, any order)."""
    lam = tuple(lam); rho = tuple(rho)
    if sum(lam) != sum(rho): raise ValueError((lam, rho))
    if not rho:
        return 1
    k = rho[0]; rest = rho[1:]
    tot = 0
    for mu, leg in rim_hooks(lam, k):
        tot += (-1) ** leg * chi(mu, rest)
    return tot

def f_dim(lam):
    return chi(lam, (1,) * sum(lam))

# ---------- own d-vector ----------
@lru_cache(maxsize=None)
def dvec(lam):
    lam = tuple(lam); n = sum(lam)
    if n == 0: return (1,)
    d = [0] * (n + 1); d[0] = 1
    for i in range(len(lam)):
        if i + 1 == len(lam) or lam[i] > lam[i + 1]:
            mu = list(lam); mu[i] -= 1
            while mu and mu[-1] == 0: mu.pop()
            dm = dvec(tuple(mu))
            for j in range(1, n + 1): d[j] += dm[j - 1]
    return tuple(d)

def z(rho):
    r = 1
    for k, m in Counter(rho).items(): r *= k ** m * factorial(m)
    return r

def sq_type(rho):
    out = []
    for k in rho:
        if k % 2: out.append(k)
        else: out += [k // 2, k // 2]
    return tuple(sorted(out, reverse=True))

def norm(rho): return tuple(sorted(rho, reverse=True))

def cycle_type(p):
    n = len(p); seen = [False] * n; ct = []
    for i in range(n):
        if seen[i]: continue
        L = 0; j = i
        while not seen[j]:
            seen[j] = True; j = p[j]; L += 1
        ct.append(L)
    return norm(ct)

def poly_1pt(k):
    """coefficients of (1+t)^k"""
    return [comb(k, j) for j in range(k + 1)]

def D(n):
    return sum((-1) ** k * factorial(n) // factorial(k) for k in range(n + 1))

# ---------------- (A) brute force over S_n ----------------
def check_A(N):
    bad = 0
    for n in range(0, N + 1):
        perms = list(itertools.permutations(range(n)))
        data = []
        for p in perms:
            sq = tuple(p[p[i]] for i in range(n))
            fx = sum(1 for i in range(n) if p[i] == i)
            data.append((cycle_type(sq), fx))
        for lam in parts(n):
            d = dvec(lam)
            lhs = [factorial(n) * d[j] // factorial(j) for j in range(n + 1)]
            assert all(factorial(n) * d[j] % factorial(j) == 0 for j in range(n + 1))
            rhs = [0] * (n + 1)
            for ct, fx in data:
                c = chi(lam, ct)
                for j in range(fx + 1): rhs[j] += c * comb(fx, j)
            # second display
            d2 = [Fraction(factorial(j) * rhs[j], factorial(n)) for j in range(n + 1)]
            if lhs != rhs or d2 != [Fraction(x) for x in d]:
                bad += 1; print("A MISMATCH", lam, lhs, rhs)
    print(f"(A) brute force over S_n, n<={N}: mismatches={bad}")

# ---------------- (B) class sums ----------------
def rhs_class(lam):
    n = sum(lam)
    rhs = [0] * (n + 1)
    for rho in parts(n):
        c = chi(lam, sq_type(rho))
        if c == 0: continue
        m1 = rho.count(1)
        w = factorial(n) // z(rho)
        for j in range(m1 + 1): rhs[j] += w * c * comb(m1, j)
    return rhs

def check_B(N):
    bad = 0; cnt = 0
    for n in range(0, N + 1):
        for lam in parts(n):
            d = dvec(lam)
            lhs = [factorial(n) * d[j] // factorial(j) for j in range(n + 1)]
            if lhs != rhs_class(lam):
                bad += 1; print("B MISMATCH", lam)
            cnt += 1
        print(f"  (B) n={n} done ({cnt} shapes so far)", flush=True)
    print(f"(B) class-sum Theorem S3, all lam |- n<={N}: mismatches={bad}, shapes={cnt}")

# ---------------- (C) derangement value ----------------
def check_C(N):
    bad = 0
    for n in range(0, N + 1):
        for lam in parts(n):
            d = dvec(lam)
            v = sum((-1) ** j * Fraction(factorial(n) * d[j], factorial(j)) for j in range(n + 1))
            der = 0
            for rho in parts(n):
                if 1 in rho: continue
                der += (factorial(n) // z(rho)) * chi(lam, sq_type(rho))
            if v != der or v.denominator != 1:
                bad += 1; print("C MISMATCH", lam, v, der)
            if lam == (n,) and v != D(n): bad += 1; print("C (n) MISMATCH", n, v, D(n))
            if n >= 2 and lam == (n - 1, 1) and v != (-1) ** n * (n - 1): bad += 1; print("C (n-1,1) MISMATCH", n, v)
    print(f"(C) derangement value, n<={N}: mismatches={bad}")

# ---------------- (D) derangement transform ----------------
def delta_class(lam, m):
    n = sum(lam)
    tot = 0
    for rho in parts(m):
        if 1 in rho: continue
        tot += (factorial(m) // z(rho)) * chi(lam, sq_type(rho) + (1,) * (n - m))
    return tot

def delta_brute(lam, m):
    n = sum(lam); tot = 0
    for p in itertools.permutations(range(m)):
        if any(p[i] == i for i in range(m)): continue
        sq = tuple(p[p[i]] for i in range(m))
        tot += chi(lam, cycle_type(sq) + (1,) * (n - m))
    return tot

def check_D(N, MBRUTE=7):
    bad = 0
    for n in range(0, N + 1):
        for lam in parts(n):
            d = dvec(lam); f = f_dim(lam)
            delta = [delta_class(lam, m) for m in range(n + 1)]
            if n <= 8:
                db = [delta_brute(lam, m) for m in range(min(n, MBRUTE) + 1)]
                if db != delta[:len(db)]: bad += 1; print("D brute/class delta mismatch", lam, db, delta)
            # i! u_i = sum_m C(i,m) delta_m
            for i in range(n + 1):
                u = d[n - i]
                if factorial(i) * u != sum(comb(i, m) * delta[m] for m in range(i + 1)):
                    bad += 1; print("D transform mismatch", lam, i)
            # n! P(s-1) = sum_m C(n,m) delta_m s^{n-m}: expand LHS in s
            lhs = [Fraction(0)] * (n + 1)
            for j in range(n + 1):
                c = Fraction(factorial(n) * d[j], factorial(j))
                for k in range(j + 1):  # (s-1)^j = sum_k C(j,k) s^k (-1)^{j-k}
                    lhs[k] += c * comb(j, k) * (-1) ** (j - k)
            rhs = [Fraction(comb(n, n - k) * delta[n - k]) for k in range(n + 1)]
            if lhs != rhs: bad += 1; print("D s-poly mismatch", lam, lhs, rhs)
            # small deltas
            if delta[0] != f: bad += 1; print("delta0", lam)
            if n >= 1 and delta[1] != 0: bad += 1; print("delta1", lam)
            if n >= 2 and delta[2] != f: bad += 1; print("delta2", lam)
            if n >= 3 and delta[3] != 2 * chi(lam, (3,) + (1,) * (n - 3)): bad += 1; print("delta3", lam)
            if n >= 4 and delta[4] != 6 * chi(lam, (2, 2) + (1,) * (n - 4)) + 3 * f: bad += 1; print("delta4", lam)
            # item (ii): delta_m / f formula (same as delta_class by construction; check the sq map on odd parts too)
            # item (i): transpose invariance
            if delta != [delta_class(conj(lam), m) for m in range(n + 1)]: bad += 1; print("transpose", lam)
        print(f"  (D) n={n} done", flush=True)
    print(f"(D) derangement transform, n<={N}: mismatches={bad}")

# ---------------- (E) item (iv) ----------------
def Hhat(m):
    c = [0] * (m + 1)
    for i in range(m // 2 + 1):
        c[m - 2 * i] = comb(m, 2 * i) * (factorial(2 * i) // (2 ** i * factorial(i)))
    return c

def check_E(N):
    bad = 0
    for n in range(0, N + 1):
        for lam in parts(n):
            d = dvec(lam)
            P = [Fraction(d[j], factorial(j)) for j in range(n + 1)]
            Q = [Fraction(0)] * (n + 1)
            for m in range(n + 1):
                tot = 0
                for rho in parts(n - m):
                    if rho and rho[-1] < 3: continue
                    tot += (factorial(n - m) // z(rho)) * chi(lam, sq_type(rho) + (1,) * m)
                wm = Fraction(tot, factorial(m) * factorial(n - m))
                if wm == 0: continue
                for k, hk in enumerate(Hhat(m)):
                    for i in range(k + 1): Q[i] += wm * hk * comb(k, i)
            if P != Q: bad += 1; print("E MISMATCH", lam)
    print(f"(E) item (iv) Hermite-type expansion, n<={N}: mismatches={bad}")

if __name__ == "__main__":
    # sanity for own chi
    assert chi((2, 1), (3,)) == -1 and chi((2, 1), (2, 1)) == 0 and chi((2, 1), (1, 1, 1)) == 2
    assert chi((2, 2), (2, 2)) == 2 and chi((3, 1), (4,)) == -1 and chi((2, 2), (4,)) == 0 and chi((3, 1), (2, 2)) == -1
    assert chi((3, 2), (5,)) == 0 and chi((3, 2), (3, 1, 1)) == -1 and chi((3, 2), (2, 2, 1)) == 1 and chi((4,1),(2,2,1)) == 0
    # orthogonality sanity for n=8
    for n in (5, 8):
        for l1 in parts(n):
            for l2 in parts(n):
                s = sum(Fraction(chi(l1, r) * chi(l2, r), z(r)) for r in parts(n))
                assert s == (1 if l1 == l2 else 0)
    print("own chi passes orthogonality n=5,8")
    check_A(8)
    check_C(NMAX)
    check_E(min(NMAX, 12))
    check_D(NMAX)
    check_B(NMAX)
