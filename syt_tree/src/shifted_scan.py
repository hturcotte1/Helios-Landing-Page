"""Scans for the shifted case (Worley Problem 2).

Usage (run from anywhere; all paths are absolute):
    python3 shifted_scan.py identity [NMAX=8] [K=10]
        Task 2: tests E_lam(t) against P_lam, P'_lam and Phi(t) P_lam(c t),
        c in {1, 1/2, 2}; verifies the corrected differential relation;
        writes data/shifted_identity.txt.
    python3 shifted_scan.py fomin [NMAX=10]
        Verifies D'U - UD' = I (weighted D') and DU - UD = I - Pi (unweighted D,
        Pi = projection onto strict partitions with last part 1) for |lam| <= NMAX;
        writes data/shifted_fomin.txt.
    python3 shifted_scan.py census [N=40] [K=14]
        Task 3: truncated up-census (s_0..s_K) of every strict partition of n <= N
        by level DP; per-n collision report, smallest separating K_n, cross-size
        collisions; writes data/shifted_census_N{N}_K{K}.txt.
    python3 shifted_scan.py siblings [N=39] [K=14]
        Task 4: for every strict lam with |lam| <= N and every pair of distinct
        covers mu, mu', check that the truncated censuses differ; writes
        data/shifted_siblings_N{N}_K{K}.txt.
Exact integer / rational arithmetic throughout.
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction
from functools import lru_cache
from math import factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shifted import (E_series, P_series, Pw_series, commutator_DU_UD,
                     down_set_weighted, s_up, series_div, series_exp_poly,
                     series_log, series_mul, series_scale, strict_partitions,
                     total_shifted_syt, up_set, addable_boxes)

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def first_diff(a, b):
    for k in range(min(len(a), len(b))):
        if a[k] != b[k]:
            return k
    return None


def fmt(series):
    return "[" + ", ".join(str(x) for x in series) + "]"


# ---------------------------------------------------------------------------
# Task 2: identity search
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def c_up(lam, k):
    """Number of chains of length k up from lam ending at a strict partition
    whose last part is 1 (k = 0: indicator of lam itself)."""
    if k == 0:
        return 1 if (lam and lam[-1] == 1) else 0
    return sum(c_up(mu, k - 1) for mu in up_set(lam))


def run_identity(NMAX=8, K=10):
    out = []
    lams = [lam for n in range(NMAX + 1) for lam in strict_partitions(n)]
    E = {lam: E_series(lam, K) for lam in lams}
    P = {lam: P_series(lam, K) for lam in lams}
    Pw = {lam: Pw_series(lam, K) for lam in lams}
    out.append(f"# identity search: all strict partitions with |lam| <= {NMAX}, series mod t^{K+1}")
    out.append(f"E_empty coefficients s_k(()) = total number of shifted SYT of size k: "
               f"{[total_shifted_syt(k) for k in range(K + 1)]}")
    out.append(f"log E_empty = {fmt(series_log(E[()], K))}")
    out.append(f"  (not of the form a t + b t^2: coefficient of t^3 is {series_log(E[()], K)[3]} != 0)")
    # (a) ratio E/Q independent of lam?
    for name, Q in [("P_lam", P), ("P'_lam", Pw)]:
        R = {lam: series_div(E[lam], Q[lam], K) for lam in lams}
        ref = R[()]
        bad = [(lam, first_diff(R[lam], ref)) for lam in lams if R[lam] != ref]
        out.append(f"E_lam / {name} independent of lam? {'YES' if not bad else 'NO'}; "
                   f"{len(bad)}/{len(lams)} shapes differ from lam=(); first failures (lam, first differing t-power): {bad[:6]}")
        if bad:
            lam0 = bad[0][0]
            out.append(f"    E_()/{name}(()) = {fmt(ref[:6])} ...;  E_{lam0}/{name}({lam0}) = {fmt(R[lam0][:6])} ...")
    # (b) E_lam == Phi(t) Q_lam(c t) with Phi forced = E_empty (Q_empty = 1)
    for name, Q in [("P_lam", P), ("P'_lam", Pw)]:
        for c in [Fraction(1), Fraction(1, 2), Fraction(2)]:
            rhs = {lam: series_mul(E[()], series_scale(Q[lam], c, K), K) for lam in lams}
            bad = [(lam, first_diff(E[lam], rhs[lam])) for lam in lams if E[lam] != rhs[lam]]
            out.append(f"E_lam == E_empty(t) * {name}({c} t)? {'YES' if not bad else 'NO'}; "
                       f"{len(bad)}/{len(lams)} failures; first: {bad[:6]}")
            if bad:
                lam0, k0 = bad[0]
                out.append(f"    lam={lam0}: coefficient of t^{k0}: LHS {E[lam0][k0]}, RHS {rhs[lam0][k0]}")
            # also allow an arbitrary universal factor: is E_lam/Q_lam(ct) lam-independent?
            R = {lam: series_div(E[lam], series_scale(Q[lam], c, K), K) for lam in lams}
            bad2 = [(lam, first_diff(R[lam], R[()])) for lam in lams if R[lam] != R[()]]
            out.append(f"    ... with an arbitrary universal factor Phi (ratio E_lam/{name}({c} t) lam-independent)? "
                       f"{'YES' if not bad2 else 'NO'}; first: {bad2[:4]}")
    # (c) structural obstruction: E_(1) = E_()' so E_(1)/E_() = (log E_())' is not a polynomial
    R1 = series_div(E[(1,)], E[()], K)
    dlog = [(k + 1) * series_log(E[()], K)[k + 1] for k in range(K)] + [None]
    assert R1[:K] == dlog[:K]
    out.append(f"E_(1)/E_() = (log E_())' = {fmt(R1[:K])} ... (E_(1) = E_()' since () has the single cover (1)).")
    out.append("  This is not a polynomial (nonzero through t^%d), so NO identity E_lam = Phi(t) Q_lam(t) with Q_lam a polynomial" % K)
    out.append("  of degree <= |lam| (whatever Q_lam is) can hold with a universal Phi: Phi would have to equal E_() (Q_() = 1)")
    out.append("  and then Q_(1) = E_(1)/E_() is not a polynomial.")
    # (d) exp(a t + b t^2) comparison for Phi := E_empty
    e = series_exp_poly(1, Fraction(1, 2), K)
    out.append(f"E_empty vs exp(t + t^2/2) (involutions): first differing t-power {first_diff(E[()], e)}; "
               f"k! * coefficients: {[E[()][k] * factorial(k) for k in range(K + 1)]} vs {[e[k] * factorial(k) for k in range(K + 1)]}")
    # (e) the corrected differential relation that DOES hold:
    #   2 s_{k+1}(lam) = sum_{mu<lam} w(mu) s_k(mu) + 2 s_k(lam) + k s_{k-1}(lam) - c_k(lam)
    ok = True
    first_bad = None
    for lam in lams:
        for k in range(K):
            lhs = 2 * s_up(lam, k + 1)
            rhs = (sum(w * s_up(mu, k) for mu, w in down_set_weighted(lam)) + 2 * s_up(lam, k)
                   + (k * s_up(lam, k - 1) if k else 0) - c_up(lam, k))
            if lhs != rhs:
                ok = False
                first_bad = first_bad or (lam, k, lhs, rhs)
    out.append("Corrected relation  A_lam' = (1/2) A_{D' lam} + (1 + t/2) A_lam - (1/2) C_lam  "
               "(C_lam = EGF of chains up from lam ending at a partition with last part 1), i.e.")
    out.append("  2 s_{k+1}(lam) = sum_w w s_k(mu) + 2 s_k(lam) + k s_{k-1}(lam) - c_k(lam):  "
               f"{'HOLDS' if ok else 'FAILS ' + str(first_bad)} for all |lam| <= {NMAX}, k < {K}")
    # (f) the unshifted-style relation without the correction term, for the record
    bad = []
    for lam in lams:
        for k in range(K):
            lhs = 2 * s_up(lam, k + 1)
            rhs = (sum(w * s_up(mu, k) for mu, w in down_set_weighted(lam)) + 2 * s_up(lam, k)
                   + (k * s_up(lam, k - 1) if k else 0))
            if lhs != rhs:
                bad.append((lam, k, lhs, rhs))
                break
    out.append(f"  without the -c_k(lam) term it fails for {len(bad)}/{len(lams)} shapes; first: {bad[:4]}")
    text = "\n".join(out)
    print(text)
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, "shifted_identity.txt"), "w") as f:
        f.write(text + "\n")


# ---------------------------------------------------------------------------
# Fomin relations
# ---------------------------------------------------------------------------

def run_fomin(NMAX=10):
    out = [f"# Fomin commutation relations on the shifted Young lattice, all strict partitions of n <= {NMAX}"]
    n_ok = 0
    bad_w = []
    bad_u = []
    for n in range(NMAX + 1):
        for lam in strict_partitions(n):
            cw = commutator_DU_UD(lam, True)
            if cw != {lam: 1}:
                bad_w.append((lam, cw))
            cu = commutator_DU_UD(lam, False)
            expected = {} if (lam and lam[-1] == 1) else {lam: 1}
            if cu != expected:
                bad_u.append((lam, cu))
            n_ok += 1
    out.append(f"D'U - UD' = I (D' weighted: off-diagonal box 2, diagonal box 1): "
               f"{'HOLDS' if not bad_w else 'FAILS'} on all {n_ok} strict partitions; failures: {bad_w[:5]}")
    out.append(f"DU - UD (D unweighted) = I - Pi, Pi = projection on strict partitions with last part 1: "
               f"{'HOLDS' if not bad_u else 'FAILS'}; failures: {bad_u[:5]}")
    ex = [(lam, commutator_DU_UD(lam, False)) for lam in [(1,), (2,), (2, 1), (3, 1), (3, 2), (3, 2, 1)]]
    out.append(f"  examples (DU-UD) lam: {ex}")
    out.append("  so DU - UD = I fails exactly on the strict partitions ending in 1 (where (DU-UD) lam = 0),")
    out.append("  because #addable(lam) - #removable(lam) = 1 - [lam_l = 1] while all cross terms cancel.")
    text = "\n".join(out)
    print(text)
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, "shifted_fomin.txt"), "w") as f:
        f.write(text + "\n")


# ---------------------------------------------------------------------------
# Truncated census by level DP
# ---------------------------------------------------------------------------

def truncated_census(N, K):
    """dict lam -> (s_0(lam), ..., s_K(lam)) for all strict partitions of n <= N.
    s_k(lam) = sum_{mu covers lam} s_{k-1}(mu); level k needs all |lam| <= N + K - k."""
    covers = {}
    for m in range(N + K + 1):
        for lam in strict_partitions(m):
            covers[lam] = up_set(lam) if m < N + K else []
    cur = {lam: 1 for lam in covers}                 # s_0
    vec = {lam: [1] for lam in covers if sum(lam) <= N}
    for k in range(1, K + 1):
        nxt = {}
        limit = N + K - k
        for lam, ups in covers.items():
            if sum(lam) <= limit:
                nxt[lam] = sum(cur[mu] for mu in ups)
        cur = nxt
        for lam in vec:
            vec[lam].append(cur[lam])
    return {lam: tuple(v) for lam, v in vec.items()}


def smallest_separating_K(vecs, K):
    """Smallest K' <= K such that the prefixes (s_0..s_K') are pairwise distinct; None if not separated."""
    for Kp in range(K + 1):
        seen = set()
        ok = True
        for v in vecs:
            p = v[: Kp + 1]
            if p in seen:
                ok = False
                break
            seen.add(p)
        if ok:
            return Kp
    return None


def run_census(N=40, K=14):
    t0 = time.time()
    cen = truncated_census(N, K)
    t1 = time.time()
    out = [f"# truncated up-census (s_0..s_{K}) of all strict partitions of n <= {N}  [{t1 - t0:.1f}s DP]",
           "# n  q(n)  #distinct_censuses  #collision_groups  K_n(smallest K separating all strict partitions of n)  max_s_K"]
    all_ok = True
    by_n = {}
    for lam, v in cen.items():
        by_n.setdefault(sum(lam), []).append((lam, v))
    Kn = {}
    for n in range(N + 1):
        items = by_n[n]
        groups = {}
        for lam, v in items:
            groups.setdefault(v, []).append(lam)
        coll = [g for g in groups.values() if len(g) > 1]
        kn = smallest_separating_K([v for _, v in items], K)
        Kn[n] = kn
        if coll:
            all_ok = False
        out.append(f"{n} {len(items)} {len(groups)} {len(coll)} {kn if kn is not None else '>' + str(K)} "
                   f"{max(v[K] for _, v in items)}")
        for g in coll[:5]:
            out.append(f"   COLLISION at n={n}: {g}")
    # cross-size collisions
    groups = {}
    for lam, v in cen.items():
        groups.setdefault(v, []).append(lam)
    cross = [g for g in groups.values() if len(g) > 1]
    out.append(f"within-size collisions: {'NONE' if all_ok else 'FOUND'} for all n <= {N}, K = {K}")
    out.append(f"cross-size collisions (equal truncated census, different |lam|, all |lam| <= {N}): {len(cross)}; examples: {cross[:5]}")
    # does s_1 (=#addable) alone, or (s_1, s_2) separate? summary of K_n
    out.append("K_n table: " + " ".join(f"{n}:{Kn[n]}" for n in range(N + 1)))
    out.append(f"max K_n over n <= {N}: {max(k for k in Kn.values() if k is not None)}")
    out.append(f"time {time.time() - t0:.1f}s")
    text = "\n".join(out)
    print(text)
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, f"shifted_census_N{N}_K{K}.txt"), "w") as f:
        f.write(text + "\n")
    return cen


def run_siblings(N=39, K=14):
    t0 = time.time()
    cen = truncated_census(N + 1, K)
    out = [f"# sibling rigidity: for every strict lam with |lam| <= {N} and every pair of distinct covers mu != mu',",
           f"# the truncated censuses (s_0..s_{K}) of mu and mu' differ.",
           "# n  q(n)  #lam_with_one_addable_box  #sibling_pairs  #equal_census_pairs  K_sib(n)(smallest K separating all sibling pairs)"]
    total_pairs = 0
    total_bad = 0
    for n in range(N + 1):
        pairs = 0
        bad = []
        one = 0
        ksib = 0
        for lam in strict_partitions(n):
            ups = up_set(lam)
            if len(ups) == 1:
                one += 1
                assert lam == tuple(range(len(lam), 0, -1))   # staircase or empty
            for a in range(len(ups)):
                for b in range(a + 1, len(ups)):
                    pairs += 1
                    va, vb = cen[ups[a]], cen[ups[b]]
                    fd = first_diff(va, vb)
                    if fd is None:
                        bad.append((lam, ups[a], ups[b]))
                    else:
                        ksib = max(ksib, fd)
        total_pairs += pairs
        total_bad += len(bad)
        out.append(f"{n} {len(list(strict_partitions(n)))} {one} {pairs} {len(bad)} {ksib if not bad else '>' + str(K)}")
        for b in bad[:5]:
            out.append(f"   EQUAL SIBLING CENSUS: {b}")
    out.append(f"total sibling pairs checked: {total_pairs}; pairs with equal truncated census: {total_bad}")
    out.append(f"sibling rigidity for all |lam| <= {N} with K = {K}: {'HOLDS' if total_bad == 0 else 'FAILS'}")
    out.append(f"time {time.time() - t0:.1f}s")
    text = "\n".join(out)
    print(text)
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, f"shifted_siblings_N{N}_K{K}.txt"), "w") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "identity"
    args = [int(x) for x in sys.argv[2:]]
    if mode == "identity":
        run_identity(*args)
    elif mode == "fomin":
        run_fomin(*args)
    elif mode == "census":
        run_census(*args)
    elif mode == "siblings":
        run_siblings(*args)
    else:
        raise SystemExit(__doc__)
