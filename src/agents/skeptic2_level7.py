"""skeptic2_level7.py -- checks (own code) of the top-end structure used as prefilters:
 (a) d_{n-3} is a function of (n, f, C_2)                      [results.md Thm 4.6]
 (b) d_{n-4}, d_{n-5}, d_{n-6} are functions of (n, f, C_2, C_1^2, C_4)   [Cor 4.10.2, Thm 4.10.3]
 (c) d_{n-7} is a function of (n, f, C_2, C_1^2, C_4, C_6 - 16 C_1 C_3)   [Thm 4.10.4, level-7 invariant]
     and NOT of (n, f, C_2, C_1^2, C_4) alone (first failure n reported).
 (d) d_{n-8} is NOT a function of the level-7 data (first failure n reported), i.e. level 8 is new.
Checked for all partitions of n <= NMAX.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from skeptic2_dvec import all_dvectors
from skeptic2_targeted import content_sums, hook_product
from math import factorial

def check(n, vecs, state):
    if n < 8:
        return
    nf = factorial(n)
    tabs = {k: {} for k in ["a", "b4", "b5", "b6", "c", "c_weak", "d"]}
    for lam, d in vecs.items():
        C1, C2, C3, C4, C6 = content_sums(lam)
        f = nf // hook_product(lam)
        assert f == d[n]
        k3 = (f, C2); k5 = (f, C2, C1 * C1, C4); k7 = k5 + (C6 - 16 * C1 * C3,)
        for name, key, val in [("a", k3, d[n - 3]), ("b4", k5, d[n - 4]), ("b5", k5, d[n - 5]), ("b6", k5, d[n - 6]),
                               ("c", k7, d[n - 7]), ("c_weak", k5, d[n - 7]), ("d", k7, d[n - 8])]:
            tabs[name].setdefault(key, set()).add(val)
    for name in tabs:
        ok = all(len(s) == 1 for s in tabs[name].values())
        if not ok and name not in state:
            state[name] = n
        if name in ("a", "b4", "b5", "b6", "c"):
            assert ok, ("functional dependence fails", name, n)
    print("n=%d: (a),(b),(c) hold; d_{n-7} function of level-5 data: %s; d_{n-8} function of level-7 data: %s" % (
        n, "c_weak" not in state, "d" not in state), flush=True)

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    state = {}
    all_dvectors(nmax, callback=lambda n, v: check(n, v, state))
    print("first n where d_{n-7} is not a function of (n,f,C2,C1^2,C4):", state.get("c_weak"))
    print("first n where d_{n-8} is not a function of level-7 data:", state.get("d"))
