"""Phase 2: level-by-level exact computation of the down-census vectors
(d_0(lam),...,d_n(lam)) for ALL partitions lam of n, n = 0..N, with a
collision report per n (groups of partitions with equal d-vector that are
not of the form {lam} (lam symmetric) or {lam, lam^t}).

Usage: python3 dcensus_scan.py N [dump_upto]
Writes data/collisions_python.txt and data/dvectors_n{n}.txt.gz for n <= dump_upto.
"""
import gzip
import hashlib
import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))
from young import partitions, removable_corners, remove_box, conjugate

DATA = os.path.join(os.path.dirname(__file__), "..", "data")


def level_dvectors(n, prev):
    """prev: dict partition-of-(n-1) -> d-vector (tuple).  Returns dict for n."""
    cur = {}
    for lam in partitions(n):
        d = [0] * (n + 1)
        d[0] = 1
        for i in removable_corners(lam):
            dm = prev[remove_box(lam, i)]
            for j in range(1, n + 1):
                d[j] += dm[j - 1]
        cur[lam] = tuple(d)
    return cur


def collision_report(n, cur):
    groups = {}
    for lam, dv in cur.items():
        groups.setdefault(dv, []).append(lam)
    bad = []
    n_sym = 0
    for dv, lams in groups.items():
        if len(lams) == 1:
            assert lams[0] == conjugate(lams[0]) or True  # a lone non-symmetric lam would mean d(lam) != d(lam^t): impossible
            if lams[0] == conjugate(lams[0]):
                n_sym += 1
            else:
                raise AssertionError(("d(lam) != d(lam^t)?!", lams[0]))
        elif len(lams) == 2 and lams[0] == conjugate(lams[1]):
            pass
        else:
            bad.append((dv, lams))
    # digest of the multiset of d-vectors of this level (reproducibility check)
    h = hashlib.sha256()
    for dv in sorted(groups):
        h.update((" ".join(map(str, dv)) + "\n").encode())
    # order-independent digest mod P = 2^61-1 shared with the C++ scanner:
    #   sum over all lam, j of d_j(lam) * 7^j  (mod P)
    P = (1 << 61) - 1
    w = [pow(7, j, P) for j in range(n + 1)]
    dig = 0
    for dv in cur.values():
        dig += sum(dj * wj for dj, wj in zip(dv, w))
    return groups, bad, n_sym, h.hexdigest(), dig % P


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    dump_upto = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    os.makedirs(DATA, exist_ok=True)
    rep = open(os.path.join(DATA, "collisions_python.txt"), "w")
    rep.write("# n  p(n)  #distinct_dvectors  #symmetric  #transpose_pairs  #collision_groups  sha256(sorted dvectors)  digest_mod_2^61-1\n")
    prev = {(): (1,)}
    t0 = time.time()
    for n in range(0, N + 1):
        cur = level_dvectors(n, prev) if n > 0 else prev
        groups, bad, n_sym, digest, dig_p = collision_report(n, cur)
        npairs = (len(cur) - n_sym) // 2
        line = f"{n} {len(cur)} {len(groups)} {n_sym} {npairs} {len(bad)} {digest} {dig_p}"
        rep.write(line + "\n")
        rep.flush()
        print(line, f"[{time.time()-t0:.1f}s]", flush=True)
        for dv, lams in bad:
            msg = f"COLLISION n={n}: partitions {lams} share d-vector {dv}"
            rep.write(msg + "\n")
            print(msg, flush=True)
        if n <= dump_upto:
            with gzip.open(os.path.join(DATA, f"dvectors_n{n:02d}.txt.gz"), "wt") as f:
                for lam in partitions(n):
                    f.write(",".join(map(str, lam)) + ": " + " ".join(map(str, cur[lam])) + "\n")
        prev = cur
    rep.close()


if __name__ == "__main__":
    main()
