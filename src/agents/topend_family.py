"""topend task (4): the families F_a: lam = (a, 1^{n-a}), mu = (2^{a-1}, 1^{n-2a+2}); claim d_j(lam) = d_j(mu) for
j <= n - 2a + 2 (PROVED by the poset isomorphism), and (observed) d_{n-2a+3} differ.  Check a = 3..7, n <= 44."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from census import d_vector
bad = 0
for a in range(3, 8):
    for n in range(2 * a - 1, 45):
        lam = (a,) + (1,) * (n - a); mu = (2,) * (a - 1) + (1,) * (n - 2 * a + 2)
        dl, dm = d_vector(lam), d_vector(mu)
        K = n - 2 * a + 2
        if any(dl[j] != dm[j] for j in range(K + 1)): bad += 1; print("FAIL agree", a, n)
        if K + 1 <= n and dl[K + 1] == dm[K + 1]: print("note: also agree at j =", K + 1, "for a,n =", a, n)
print("family F_a check (a=3..7, n<=44):", "OK" if bad == 0 else bad)
