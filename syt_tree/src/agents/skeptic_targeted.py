"""skeptic_targeted.py -- targeted search for n <= NMAX (default 50).
Prefilter keys (all transpose-invariant, all cheap exact integers):
  K0 = (n, multiset of unordered corner-box pairs {a_i,b_i})
  K1 = K0 + f^lam
  K2 = K1 + (C_2, C_1^2, C_4)
For each key report the number of surviving pairs of transpose classes; for K2 survivors
(and, up to a cap, K1 survivors) compute exact d-vectors with an independent memoised
recursion and report the Hamming distance and first/last differing index.
"""
import sys, time
from math import prod
sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src/agents')
from skeptic_dvec import gen_partitions, transpose, removable, remove_box
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 50
NMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 10
CAP = 60

def corner_multiset(lam):
    parts = sorted(set(lam), reverse=True) + [0]
    pairs = []
    for i in range(len(parts)-1):
        a, b = parts[i]-parts[i+1], lam.count(parts[i])
        pairs.append((min(a, b), max(a, b)))
    return tuple(sorted(pairs))

def contents(lam):
    c1 = c2 = c4 = 0
    for i, p in enumerate(lam):
        for j in range(p):
            c = j - i; c1 += c; c2 += c*c; c4 += c**4
    return (c2, c1*c1, c4)

def f_hook(lam):
    n = sum(lam); lt = transpose(lam)
    h = prod(lam[i] - j + lt[j] - i - 1 for i in range(len(lam)) for j in range(lam[i]))
    from math import factorial
    assert factorial(n) % h == 0
    return factorial(n) // h

memo = {(): (1,)}
def dvec(lam):
    v = memo.get(lam)
    if v is not None: return v
    n = sum(lam)
    vec = [1] + [0]*n
    for i in removable(lam):
        sub = dvec(remove_box(lam, i))
        for j in range(n): vec[j+1] += sub[j]
    memo[lam] = tuple(vec)
    return memo[lam]

for n in range(NMIN, NMAX+1):
    t0 = time.time()
    g0 = {}
    for lam in gen_partitions(n):
        if lam > transpose(lam): continue        # one rep per transpose class
        g0.setdefault((corner_multiset(lam), contents(lam)), []).append(lam)
    gC = {}   # corner multiset only
    for (cm, ct), L in g0.items(): gC.setdefault(cm, []).extend(L)
    n0 = sum(len(L)*(len(L)-1)//2 for L in gC.values())
    n0c = sum(len(L)*(len(L)-1)//2 for L in g0.values())
    # K1: corner multiset + f
    g1 = {}
    for cm, L in gC.items():
        if len(L) < 2: continue
        for lam in L: g1.setdefault((cm, f_hook(lam)), []).append(lam)
    n1 = sum(len(L)*(len(L)-1)//2 for L in g1.values())
    # K2: corner multiset + f + contents
    g2 = {}
    for (cm, f), L in g1.items():
        if len(L) < 2: continue
        for lam in L: g2.setdefault((cm, f, contents(lam)), []).append(lam)
    n2 = sum(len(L)*(len(L)-1)//2 for L in g2.values())
    print(f"n={n} classes={sum(len(L) for L in gC.values())} pairs: corner={n0} corner+C={n0c} corner+f={n1} corner+f+C={n2} ({time.time()-t0:.0f}s)", flush=True)
    survivors = [(L, 'K2') for L in g2.values() if len(L) >= 2]
    k1_only = [L for L in g1.values() if len(L) >= 2]
    shown = 0
    for L, tag in survivors + [(L, 'K1') for L in k1_only]:
        if tag == 'K1' and shown >= CAP: break
        for x in range(len(L)):
            for y in range(x+1, len(L)):
                lam, mu = L[x], L[y]
                if tag == 'K1' and contents(lam) == contents(mu): continue  # already in K2
                if tag == 'K1': shown += 1
                dl, dm = dvec(lam), dvec(mu)
                diff = [j for j in range(n+1) if dl[j] != dm[j]]
                print(f"   [{tag}] {lam} vs {mu}: ham={len(diff)} diff_idx={diff[:6]}..{diff[-4:]} "
                      f"m_bottom_agree={diff[0]} top_agree={n-diff[-1]} C={contents(lam)}|{contents(mu)}", flush=True)
    memo.clear()
