"""skeptic_dvec.py -- independent d-vector code (shares no code with census.py / young.py).

Partition = tuple of weakly decreasing positive ints.
d_j(lam) = number of ways to remove j boxes one at a time (d_0 = 1),
computed level by level:  d_j(lam) = sum_{c removable} d_{j-1}(lam - c).
Only levels k-1 and k are kept in memory.

Usage:  python skeptic_dvec.py NMAX  -> for each n <= NMAX checks that
  (i)  d(lam) == d(lam^t) for every lam,
  (ii) no two partitions of n with mu not in {lam, lam^t} share a d-vector,
  and prints per-n counts + a digest.
Also exposes  levels(NMAX)  generator yielding (n, dict lam -> d-vector) for reuse.
"""
import sys, hashlib

def gen_partitions(n, maxpart=None):
    if maxpart is None or maxpart > n:
        maxpart = n
    if n == 0:
        yield ()
        return
    for first in range(maxpart, 0, -1):
        for rest in gen_partitions(n - first, first):
            yield (first,) + rest

def transpose(lam):
    if not lam:
        return ()
    return tuple(sum(1 for p in lam if p > i) for i in range(lam[0]))

def removable(lam):
    """indices i such that removing a box from row i leaves a partition"""
    out = []
    L = len(lam)
    for i in range(L):
        if i == L - 1 or lam[i] > lam[i + 1]:
            out.append(i)
    return out

def remove_box(lam, i):
    if lam[i] == 1:
        return lam[:i] + lam[i + 1:]
    return lam[:i] + (lam[i] - 1,) + lam[i + 1:]

def levels(nmax):
    """yield (n, D) with D[lam] = [d_0,...,d_n] for all lam |- n, n = 0..nmax"""
    prev = {(): [1]}
    yield 0, prev
    for n in range(1, nmax + 1):
        cur = {}
        for lam in gen_partitions(n):
            vec = [1] + [0] * n
            for i in removable(lam):
                sub = prev[remove_box(lam, i)]
                for j in range(n):
                    vec[j + 1] += sub[j]
            cur[lam] = vec
        yield n, cur
        prev = cur

def main(nmax):
    P = (1 << 61) - 1
    for n, D in levels(nmax):
        # (i) transpose invariance
        for lam, v in D.items():
            if D[transpose(lam)] != v:
                print("TRANSPOSE FAILURE", lam); sys.exit(1)
        # (ii) collisions
        groups = {}
        for lam, v in D.items():
            groups.setdefault(tuple(v), []).append(lam)
        nsym = sum(1 for lam in D if transpose(lam) == lam)
        npairs = (len(D) - nsym) // 2
        bad = [g for g in groups.values() if len(g) > 1 and not (len(g) == 2 and transpose(g[0]) == g[1])]
        digest = 0
        for lam, v in D.items():
            for j, x in enumerate(v):
                digest = (digest + x * pow(7, j, P)) % P
        sha = hashlib.sha256(repr(sorted(groups)).encode()).hexdigest()[:16]
        print(f"n={n:2d} p(n)={len(D):6d} distinct={len(groups):6d} sym={nsym:3d} pairs={npairs:6d} "
              f"collision_groups={len(bad)} digest={digest} sha={sha}")
        for g in bad:
            print("   COLLISION:", g)
        sys.stdout.flush()

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 30)
