"""Closed-form d-vectors of hooks and two-row shapes (Aitken), Theorem F_a and its sharpness.

Lemma H (hooks).  lam = (a, 1^b), n = a + b, a >= 1, b >= 0.  For 0 <= j < n,
   d_j(lam) = sum_{c=1}^{a} [0 <= n-j-c <= b] * C(j, a-c),      d_n(lam) = C(n-1, a-1).
Lemma T (two rows).  lam = (p, q), p >= q >= 0, n = p+q.  For 0 <= j <= n,
   d_j(lam) = sum_{(p',q'): p'+q' = n-j, 0 <= q' <= q, q' <= p' <= p} [ C(j, p-p') - C(j, p-q'+1) ].
Theorem F_a + sharpness.  For a >= 3, n >= 2a-1, lam = (a,1^{n-a}), mu = (2^{a-1}, 1^{n-2a+2}) = (n-a+1, a-1)^t:
   d_j(lam) = d_j(mu) for j <= n-2a+2   and   d_{n-2a+3}(lam) - d_{n-2a+3}(mu) = 1.
"""
import sys, time
from math import comb
from topend3_lib import *
t0 = time.time()
def Cb(j, k): return comb(j, k) if 0 <= k <= j else 0
def d_hook(a, b):
    n = a + b; d = []
    for j in range(n):
        d.append(sum(Cb(j, a - c) for c in range(1, a + 1) if 0 <= n - j - c <= b))
    d.append(comb(n - 1, a - 1)); return tuple(d)
def d_tworow(p, q):
    n = p + q; d = []
    for j in range(n + 1):
        s = 0
        for qq in range(0, q + 1):
            pp = n - j - qq
            if qq <= pp <= p:
                s += Cb(j, p - pp) - Cb(j, p - qq + 1)
        d.append(s)
    return tuple(d)
# verify against the Young-lattice recursion
NV = int(sys.argv[1]) if len(sys.argv) > 1 else 30
cnt = 0
for n in range(1, NV + 1):
    for a in range(1, n + 1):
        assert d_hook(a, n - a) == d_vector((a,) + (1,) * (n - a)), (a, n); cnt += 1
    for q in range(0, n // 2 + 1):
        lam = (n - q, q) if q else (n,)
        assert d_tworow(n - q, q) == d_vector(lam), (n, q); cnt += 1
print(f"Lemmas H, T verified against d_vector for all hooks and two-row shapes with n <= {NV}: {cnt} shapes OK [{time.time()-t0:.0f}s]")
# Theorem F_a and sharpness via closed forms (large range), and via d_vector (moderate range)
AMAX, NMAX = 12, 150
cnt = 0
for a in range(3, AMAX + 1):
    for n in range(2 * a - 1, NMAX + 1):
        dl = d_hook(a, n - a); dm = d_tworow(n - a + 1, a - 1)
        s = n - 2 * a + 2
        assert dl[:s + 1] == dm[:s + 1], (a, n)
        assert dl[s + 1] - dm[s + 1] == 1, (a, n, dl[s + 1], dm[s + 1])
        cnt += 1
print(f"Theorem F_a (agreement j <= n-2a+2) and sharpness (difference exactly 1 at j = n-2a+3) via closed forms: a <= {AMAX}, n <= {NMAX}: {cnt} cases OK [{time.time()-t0:.0f}s]")
cnt = 0
for a in range(3, 8):
    for n in range(2 * a - 1, 41):
        lam = (a,) + (1,) * (n - a); mu = (2,) * (a - 1) + (1,) * (n - 2 * a + 2)
        dl, dm = d_vector(lam), d_vector(mu); s = n - 2 * a + 2
        assert dl[:s + 1] == dm[:s + 1] and dl[s + 1] - dm[s + 1] == 1, (a, n)
        cnt += 1
print(f"same via d_vector recursion directly: a <= 7, n <= 40: {cnt} cases OK [{time.time()-t0:.0f}s]")
# the closed-form values of the difference at all higher j (for the record) for a = 3
for a in (3, 4):
    for n in (10, 15, 20):
        dl = d_hook(a, n - a); dm = d_tworow(n - a + 1, a - 1)
        print(f"a={a} n={n}: d(lam)-d(mu) = {[x - y for x, y in zip(dl, dm)]}")
