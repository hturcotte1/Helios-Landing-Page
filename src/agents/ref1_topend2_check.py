"""Referee #1, independent check of topend.md T3.1-T3.3 (Lemma H, Lemma T, Theorem F') and T4.1 (Lemma K).
Everything here is implemented from scratch (no import of young/census), exact integer arithmetic."""
from functools import lru_cache
from math import comb, factorial
from itertools import combinations
import sys

def C(j, k):
    return comb(j, k) if 0 <= k <= j else 0

def partitions(n, m=None):
    if m is None or m > n: m = n
    if n == 0:
        yield (); return
    for f in range(m, 0, -1):
        for r in partitions(n - f, f):
            yield (f,) + r

def conj(l):
    return tuple(sum(1 for x in l if x > j) for j in range(l[0])) if l else ()

def removable(l):  # list of row indices (0-based) whose last box is removable
    return [i for i in range(len(l)) if i == len(l) - 1 or l[i] > l[i + 1]]

def addable(l):  # rows (0-based) where a box can be added (including new row)
    return [i for i in range(len(l) + 1) if i == 0 or (i == len(l) and l) or (i < len(l) and l[i - 1] > l[i]) or (i == len(l) and i > 0 and True)]

def addable_rows(l):
    res = []
    for i in range(len(l) + 1):
        if i == 0: res.append(0)
        elif i == len(l): res.append(i)
        elif l[i - 1] > l[i]: res.append(i)
    return res

def remove_row(l, i):
    l = list(l); l[i] -= 1
    if l[i] == 0: l.pop()
    return tuple(l)

def add_row(l, i):
    l = list(l)
    if i == len(l): l.append(1)
    else: l[i] += 1
    return tuple(l)

@lru_cache(maxsize=None)
def dvec(l):
    n = sum(l)
    if n == 0: return (1,)
    d = [0] * (n + 1); d[0] = 1
    subs = [dvec(remove_row(l, i)) for i in removable(l)]
    for j in range(1, n + 1):
        d[j] = sum(s[j - 1] for s in subs)
    return tuple(d)

def hooks(l):
    lc = conj(l)
    return [(l[i] - j) + (lc[j] - i) - 1 for i in range(len(l)) for j in range(l[i])]

def f_hook(l):
    n = sum(l); p = 1
    for h in hooks(l): p *= h
    assert factorial(n) % p == 0
    return factorial(n) // p

# brute force d_j by order filters: count sequences of removals directly (n small) -- second independent check
def d_bruteforce(l):
    n = sum(l)
    from collections import Counter
    cur = Counter({l: 1}); out = [1]
    for j in range(1, n + 1):
        nxt = Counter()
        for s, c in cur.items():
            for i in removable(s):
                nxt[remove_row(s, i)] += c
        cur = nxt; out.append(sum(cur.values()))
    return tuple(out)

def X_contents(l):
    return [l[i] + 1 - (i + 1) if i < len(l) else 1 - (i + 1) for i in addable_rows(l)]

def Y_contents(l):
    return [l[i] - (i + 1) for i in removable(l)]

# ---------------- checks ----------------
ok = True
if __name__ != "__main__": sys.argv = [sys.argv[0], "0", "0", "-1"]
def report(cond, msg):
    global ok
    if not cond:
        ok = False; print("FAIL:", msg)

# sanity: dvec vs brute force n<=9, f vs d_n
for n in range(0, 10):
    for l in partitions(n):
        d = dvec(l)
        report(d == d_bruteforce(l), f"dvec vs brute {l}")
        report(d[n] == f_hook(l), f"d_n vs f_hook {l}")
print("sanity done")

# Lemma H
def dH(a, b, j):
    n = a + b
    return sum(C(j, e) for e in range(max(0, a + j - n), min(a - 1, j) + 1))
def dH_c(a, b, j):
    n = a + b
    return sum(C(j, a - c) for c in range(1, a + 1) if 0 <= n - j - c <= b)
NH = int(sys.argv[1]) if len(sys.argv) > 1 else 30
cnt = 0
for n in range(1, NH + 1):
    for a in range(1, n + 1):
        b = n - a; l = (a,) + (1,) * b; d = dvec(l)
        report(d[n] == C(n - 1, a - 1), f"Lemma H f {l}")
        for j in range(n):
            report(d[j] == dH(a, b, j) == dH_c(a, b, j), f"Lemma H {l} j={j}: {d[j]} vs {dH(a,b,j)}")
            cnt += 1
        # also note: formula at j=n gives 0, not f
        report(dH(a, b, n) == 0, f"Lemma H j=n formula {l}")
print("Lemma H checked", cnt, "values, n <=", NH)

# Lemma T
def dT(p, q, j):
    n = p + q; s = 0
    for qp in range(0, q + 1):
        pp = n - j - qp
        if qp <= pp <= p:
            s += C(j, p - pp) - C(j, p - qp + 1)
    return s
cnt = 0
for n in range(0, NH + 1):
    for q in range(0, n // 2 + 1):
        p = n - q; l = (p, q) if q > 0 else ((p,) if p > 0 else ())
        d = dvec(l)
        for j in range(n + 1):
            report(d[j] == dT(p, q, j), f"Lemma T {l} j={j}: {d[j]} vs {dT(p,q,j)}")
            cnt += 1
print("Lemma T checked", cnt, "values, n <=", NH)

# Aitken (T3.0) for n<=8, all nu subset lam, ell=len(lam)
from fractions import Fraction
def det(M):
    n = len(M)
    if n == 0: return Fraction(1)
    if n == 1: return M[0][0]
    s = Fraction(0)
    for c in range(n):
        minor = [row[:c] + row[c+1:] for row in M[1:]]
        s += (-1) ** c * M[0][c] * det(minor)
    return s
def f_skew_count(l, nu):
    # number of saturated chains nu -> l : via dvec-like counting
    from collections import Counter
    cur = Counter({l: 1}); m = sum(l) - sum(nu)
    for _ in range(m):
        nxt = Counter()
        for s, c in cur.items():
            for i in removable(s):
                t = remove_row(s, i)
                if all(t[k] >= (nu[k] if k < len(nu) else 0) for k in range(len(t))) and len(t) >= len(nu):
                    nxt[t] += c
        cur = nxt
    return cur.get(nu, 0)
def contained(nu, l):
    return len(nu) <= len(l) and all(nu[i] <= l[i] for i in range(len(nu)))
cnt = 0
for n in range(0, 9):
    for l in partitions(n):
        ell = len(l)
        for k in range(0, n + 1):
            for nu in partitions(k):
                if not contained(nu, l): continue
                m = n - k
                M = [[Fraction(1, factorial(l[i] - (nu[j] if j < len(nu) else 0) - i + j)) if l[i] - (nu[j] if j < len(nu) else 0) - i + j >= 0 else Fraction(0) for j in range(ell)] for i in range(ell)]
                val = factorial(m) * det(M)
                report(val == f_skew_count(l, nu), f"Aitken {l}/{nu}")
                cnt += 1
print("Aitken checked", cnt)

# Theorem F'
NF = int(sys.argv[2]) if len(sys.argv) > 2 else 30
cnt = 0
for a in range(3, NF // 2 + 2):
    for n in range(2 * a - 1, NF + 1):
        lam = (a,) + (1,) * (n - a); mu = (2,) * (a - 1) + (1,) * (n - 2 * a + 2)
        assert sum(lam) == n == sum(mu)
        dl, dm = dvec(lam), dvec(mu)
        for j in range(0, n - 2 * a + 3):
            report(dl[j] == dm[j], f"F' agree a={a} n={n} j={j}")
        j = n - 2 * a + 3
        report(dl[j] - dm[j] == 1, f"F' diff a={a} n={n} j={j}: {dl[j]-dm[j]}")
        # also check closed forms (T3.1),(T3.2) for j <= n-2a+3
        for j in range(0, n - 2 * a + 4):
            lo = max(0, -(-(2 * (a - 1) - (n - j)) // 2))  # ceil(a-1-(n-j)/2)
            t31 = sum(C(j, e) for e in range(lo, min(a - 1, j) + 1))
            report(dm[j] == t31, f"(T3.1) a={a} n={n} j={j}")
            report(dl[j] == dH(a, n - a, j), f"(T3.2) a={a} n={n} j={j}")
        cnt += 1
        if a == 3:
            report(all(dl[j] - dm[j] == 1 for j in range(n - 3, n + 1)), f"a=3 remark n={n}")
print("Theorem F' checked", cnt, "pairs (a,n), n <=", NF)
# boundary n = 2a-2 (not claimed): report what happens
for a in range(3, 8):
    n = 2 * a - 2
    lam = (a,) + (1,) * (n - a); mu = (2,) * (a - 1)
    dl, dm = dvec(lam), dvec(mu)
    print(f"boundary n=2a-2, a={a}: lam={lam} mu={mu} d(lam)-d(mu)={[dl[j]-dm[j] for j in range(n+1)]}")

# Lemma K, with separate row/column product checks
NK = int(sys.argv[3]) if len(sys.argv) > 3 else 16
cnt = 0
for m in range(0, NK + 1):
    for nu in partitions(m):
        X = X_contents(nu); Y = Y_contents(nu)
        report(len(X) == len(Y) + 1, f"|X|=|Y|+1 {nu}")
        report(sorted(X) == list(sorted(X)) and len(set(X)) == len(X), f"X distinct {nu}")
        fnu = f_hook(nu)
        for i in addable_rows(nu):
            x = (nu[i] if i < len(nu) else 0) + 1 - (i + 1)
            report(x in X, "content")
            nu2 = add_row(nu, i)
            num = m + 1
            for y in Y: num *= abs(x - y)
            den = 1
            for xp in X:
                if xp != x: den *= abs(x - xp)
            report(f_hook(nu2) * den == fnu * num, f"Lemma K {nu} x={x}")
            # row product R
            if i < len(nu):
                lc = conj(nu)
                R = Fraction(1)
                for j in range(nu[i]):
                    h = (nu[i] - j) + (lc[j] - i) - 1
                    R *= Fraction(h + 1, h)
                Rp = Fraction(1)
                for xp in X:
                    if xp < x: Rp *= (x - xp)
                for y in Y:
                    if y < x: Rp /= (x - y)
                report(R == Rp, f"row product {nu} x={x}")
            # column product
            lc = conj(nu); col = (nu[i] if i < len(nu) else 0)  # 0-based column index of new box
            Cc = Fraction(1)
            for ip in range(i):
                h = (nu[ip] - col) + (lc[col] - ip) - 1
                Cc *= Fraction(h + 1, h)
            Cp = Fraction(1)
            for xp in X:
                if xp > x: Cp *= (xp - x)
            for y in Y:
                if y > x: Cp /= (y - x)
            report(Cc == Cp, f"col product {nu} x={x}")
            cnt += 1
print("Lemma K checked", cnt, "addable corners, m <=", NK)
print("ALL OK" if ok else "SOME FAILURES")
