"""ref2_skeptic_top4.py -- referee #2: is d_{n-4} injective in C_1^2 given (n, f, C_2)?  (used for 'H_k disagrees at d_{n-4}')"""
import sys
sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src')
from young import partitions
from census import d_vector
def C(lam):
    c1=c2=0
    for i,l in enumerate(lam):
        for c in range(l): c1+=c-i; c2+=(c-i)**2
    return c1,c2
for n in range(4, 21):
    g = {}
    for lam in partitions(n):
        d = d_vector(lam); c1,c2 = C(lam)
        g.setdefault((n, d[n], c2), {}).setdefault(d[n-4], set()).add(c1*c1)
    bad = sum(1 for key, m in g.items() for dn4, s in m.items() if len(s) > 1)
    bad2 = sum(1 for key, m in g.items() if len({next(iter(s)) for s in m.values()}) < len(m))  # same C1^2 -> different d_{n-4}?
    print("n=%d: d_{n-4} <-> C_1^2 bijective given (n,f,C_2): %s" % (n, bad == 0 and bad2 == 0))
