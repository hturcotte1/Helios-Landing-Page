"""skeptic_near.py -- structured near-collision scan, n <= NMAX (default 40), pure Python, exact integers.
For pairs {lam,mu} of transpose classes of partitions of n:
  ham   = #{j in [0,n] : d_j(lam) != d_j(mu)},  first/last = smallest/largest such j,
  a = first (prefix agreement length), b = n-last (suffix agreement length), window w = n+1-a-b.
(1) All pairs with ham <= K, found by pigeonhole: coordinates split into K+1 interleaved blocks
    (block t = {j : j = t mod K+1}); a pair with <= K mismatches agrees on a whole block; group by block, compare in groups.
(2) The pairs maximising a+b (= minimising the window): for each a, group by prefix (d_0..d_{a-1});
    inside a group the longest common suffix is attained by adjacent elements after sorting by reversed vector.
(3) Counts of pairs agreeing on the top s entries and on the bottom a entries (heuristic input).
"""
import sys, json
sys.path.insert(0, '/home/user/Helios-Landing-Page/syt_tree/src/agents')
from skeptic_dvec import levels, transpose
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
K = int(sys.argv[2]) if len(sys.argv) > 2 else 8
out = {}

def runs(lam):
    parts = sorted(set(lam), reverse=True) + [0]
    return [(parts[i]-parts[i+1], lam.count(parts[i])) for i in range(len(parts)-1)]

def boxes(lam):
    return {(i, j) for i, p in enumerate(lam) for j in range(p)}

def classify(lam, mu):
    tags = []
    for a, b in ((lam, mu), (lam, transpose(mu))):
        ba, bb = boxes(a), boxes(b)
        if len(ba - bb) == 1 and len(bb - ba) == 1:
            (i1, j1), = ba - bb; (i2, j2), = bb - ba
            tags.append(f"boxmove c{j1-i1}->{j2-i2}")
    for a, b in ((lam, mu), (lam, transpose(mu)), (transpose(lam), mu), (transpose(lam), transpose(mu))):
        if len(a) == 2 and a[1] >= 2 and b == (a[0],) + (1,) * a[1]:
            tags.append("two-row/hook (lam1,m)~(lam1,1^m)")
    return tags

def lcs_suffix(u, v):
    k = 0
    while k < len(u) and u[-1-k] == v[-1-k]: k += 1
    return k

for n, D in levels(NMAX):
    if n < 6: continue
    reps = sorted(lam for lam in D if lam <= transpose(lam))
    V = [tuple(D[lam]) for lam in reps]
    N = len(reps)
    # (1) Hamming <= K by pigeonhole
    found = {}
    for t in range(K + 1):
        idx = list(range(t, n + 1, K + 1))
        g = {}
        for i in range(N):
            g.setdefault(tuple(V[i][j] for j in idx), []).append(i)
        for L in g.values():
            if len(L) < 2: continue
            for x in range(len(L)):
                for y in range(x + 1, len(L)):
                    i, j = L[x], L[y]
                    if (i, j) in found: continue
                    diff = [q for q in range(n + 1) if V[i][q] != V[j][q]]
                    if len(diff) <= K: found[(i, j)] = diff
    near = sorted(((len(diff), diff[0], diff[-1], i, j, diff) for (i, j), diff in found.items()))
    # (2) max (a+b)
    best = (-1, None)
    for a in range(0, n + 1):
        g = {}
        for i in range(N):
            g.setdefault(V[i][:a], []).append(i)
        for L in g.values():
            if len(L) < 2: continue
            L.sort(key=lambda i: V[i][::-1])
            for x in range(len(L) - 1):
                i, j = L[x], L[x + 1]
                b = min(lcs_suffix(V[i], V[j]), n + 1 - a)
                if a + b > best[0]: best = (a + b, (i, j))
    ab, (i, j) = best
    diff = [q for q in range(n + 1) if V[i][q] != V[j][q]]
    # (3) counts
    top_counts, bot_counts = {}, {}
    for s in range(3, 13):
        g = {}
        for v in V: g[v[n - s + 1:]] = g.get(v[n - s + 1:], 0) + 1
        top_counts[s] = sum(c * (c - 1) // 2 for c in g.values())
    for a in range(2, 14):
        g = {}
        for v in V: g[v[:a]] = g.get(v[:a], 0) + 1
        bot_counts[a] = sum(c * (c - 1) // 2 for c in g.values())
    rec = {'N': N, 'npairs': N * (N - 1) // 2, 'top_counts': top_counts, 'bot_counts': bot_counts,
           'min_ham': near[0][0] if near else None, 'n_near': len(near),
           'near': [{'lam': reps[i], 'mu': reps[j], 'ham': h, 'diff_idx': d, 'runs': [runs(reps[i]), runs(reps[j])],
                     'tags': classify(reps[i], reps[j])} for h, f_, l_, i, j, d in near[:40]],
           'best_ends': {'a_plus_b': ab, 'lam': reps[i], 'mu': reps[j], 'ham': len(diff), 'diff_idx': diff,
                         'runs': [runs(reps[i]), runs(reps[j])], 'tags': classify(reps[i], reps[j])}}
    out[n] = rec
    print(f"n={n} classes={N} min_ham={rec['min_ham']} #pairs(ham<={K})={len(near)} "
          f"best a+b={ab} window={n+1-ab}: {reps[i]} {reps[j]} ham={len(diff)} diff={diff} tags={rec['best_ends']['tags']}", flush=True)
    print("   top_counts", top_counts, "\n   bot_counts", bot_counts, flush=True)
    for e in rec['near'][:8]:
        print("   ", e['ham'], e['diff_idx'], e['lam'], e['mu'], e['runs'], e['tags'], flush=True)
    json.dump(out, open('/home/user/Helios-Landing-Page/syt_tree/src/agents/skeptic_near.json', 'w'), default=str)
