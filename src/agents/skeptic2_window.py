"""skeptic2_window.py -- near-collision structure for n <= NMAX (one representative per transpose class).

For a pair (lam, mu) of transpose classes of size n let a = first index with d_a(lam) != d_a(mu),
b = n - (last differing index), w = n + 1 - a - b = length of the disagreement window.
Reports per n:
  * wmin_all : min w over all pairs, with the pairs attaining it,
  * wmin_f   : min w over pairs with equal f (b >= 3), the attaining pairs,
  * amax(b>=B): largest prefix agreement among pairs agreeing on the top B entries, B = 3,4,5,6,
  * the numbers of pairs agreeing on the bottom a0 entries / top b0 entries (for the heuristic),
  * V_j = number of distinct values of d_j over the classes.
Exact integer arithmetic; own d-vector code (skeptic2_dvec).
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from skeptic2_dvec import all_dvectors, conjugate

def reps(vecs):
    out = {}
    for lam, d in vecs.items():
        c = conjugate(lam)
        key = min(lam, c)
        out[key] = d
    return out

def max_common_prefix_pairs(keys, vecs, reverse=False, top=3):
    """Among vecs restricted to keys, find the maximal common prefix length (of the reversed vectors if reverse)
    between two distinct classes, and up to `top` pairs attaining it. Exact via sorting."""
    items = []
    for k in keys:
        v = vecs[k]
        items.append((tuple(reversed(v)) if reverse else tuple(v), k))
    items.sort()
    best, pairs = -1, []
    for i in range(len(items) - 1):
        u, v = items[i][0], items[i + 1][0]
        c = 0
        L = len(u)
        while c < L and u[c] == v[c]:
            c += 1
        if c > best:
            best, pairs = c, [(items[i][1], items[i + 1][1])]
        elif c == best and len(pairs) < top:
            pairs.append((items[i][1], items[i + 1][1]))
    return best, pairs

def diffs(u, v):
    return [j for j in range(len(u)) if u[j] != v[j]]

def analyse(n, vecs, out):
    R = reps(vecs)
    keys = list(R)
    N = len(keys)
    res = {"n": n, "classes": N}
    # V_j
    res["V"] = [len(set(R[k][j] for k in keys)) for j in range(n + 1)]
    # pairs agreeing on bottom a0 entries, top b0 entries
    bot, topc = [], []
    for a0 in range(1, n + 2):
        g = {}
        for k in keys:
            g[tuple(R[k][:a0])] = g.get(tuple(R[k][:a0]), 0) + 1
        bot.append(sum(s * (s - 1) // 2 for s in g.values()))
    for b0 in range(1, n + 2):
        g = {}
        for k in keys:
            g[tuple(R[k][n + 1 - b0:])] = g.get(tuple(R[k][n + 1 - b0:]), 0) + 1
        topc.append(sum(s * (s - 1) // 2 for s in g.values()))
    res["pairs_agree_bottom"] = bot   # index a0-1 : pairs agreeing on d_0..d_{a0-1}
    res["pairs_agree_top"] = topc     # index b0-1 : pairs agreeing on d_{n+1-b0}..d_n
    # min window over all pairs: for each prefix group compute max suffix agreement
    wmin, wpairs = None, []
    for a0 in range(1, n + 1):
        groups = {}
        for k in keys:
            groups.setdefault(tuple(R[k][:a0]), []).append(k)
        bestb = -1; bp = []
        for g in groups.values():
            if len(g) < 2:
                continue
            b, prs = max_common_prefix_pairs(g, R, reverse=True)
            if b > bestb:
                bestb, bp = b, prs
        if bestb < 0:
            break
        w = n + 1 - a0 - bestb
        if wmin is None or w < wmin:
            wmin, wpairs = w, bp
    res["wmin_all"] = wmin
    res["wmin_all_pairs"] = [(p, q, diffs(R[p], R[q])) for p, q in wpairs[:2]]
    # pairs with equal f (b >= 3): exact list of (a, b, w) for all such pairs
    byf = {}
    for k in keys:
        byf.setdefault(R[k][n], []).append(k)
    best = None; bestpairs = []
    amax_top = {3: (-1, None), 4: (-1, None), 5: (-1, None), 6: (-1, None), 7: (-1, None)}
    nf = 0
    for g in byf.values():
        if len(g) < 2:
            continue
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                u, v = R[g[i]], R[g[j]]
                nf += 1
                D = diffs(u, v)
                a, b = D[0], n - D[-1]
                w = n + 1 - a - b
                if best is None or w < best:
                    best, bestpairs = w, [(g[i], g[j], a, b)]
                elif w == best and len(bestpairs) < 3:
                    bestpairs.append((g[i], g[j], a, b))
                for B in amax_top:
                    if b >= B and a > amax_top[B][0]:
                        amax_top[B] = (a, (g[i], g[j], a, b, D))
    res["pairs_equal_f"] = nf
    res["wmin_f"] = best
    res["wmin_f_pairs"] = bestpairs
    res["amax_top"] = {B: amax_top[B] for B in amax_top}
    out.append(res)
    print("n=%2d classes=%6d  wmin_all=%s %s | equal-f pairs=%d wmin_f=%s %s | amax(b>=4)=%s amax(b>=5)=%s amax(b>=6)=%s" % (
        n, N, wmin, [(p, q, (D[0], n - D[-1])) for p, q, D in res["wmin_all_pairs"][:1]], nf, best,
        [(p, q, a, b) for p, q, a, b in bestpairs[:1]],
        amax_top[4][0], amax_top[5][0], amax_top[6][0]), flush=True)
    if amax_top[4][1]:
        print("    best (f,C2)-pair:", amax_top[4][1][:4], "diff idx", amax_top[4][1][4][:3], "...", amax_top[4][1][4][-3:], flush=True)
    if amax_top[5][1]:
        print("    best top-5 pair:", amax_top[5][1][:4], flush=True)

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    out = []
    all_dvectors(nmax, callback=lambda n, v: analyse(n, v, out) if n >= 6 else None)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "skeptic2_window.json"), "w") as fh:
        json.dump(out, fh, default=str)
