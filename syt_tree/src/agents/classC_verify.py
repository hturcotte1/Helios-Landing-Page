"""classC: computational verification of every lemma used in agent_notes/classC.md.

All arithmetic exact.  Run:  python3 src/agents/classC_verify.py [NMAX_FINAL]
Each section prints the range checked and asserts.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from math import comb
from young import partitions, conjugate, corner_runs, from_corner_runs, f_hook, boxes, removable_corners, remove_box, down_set
from census import d_vector

def binom(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def r_of(lam):
    return len(removable_corners(lam))

def A1B1(lam):
    runs = corner_runs(lam)
    return sum(1 for a, b in runs if a == 1), sum(1 for a, b in runs if b == 1)

def in_C(lam):
    # avoids (3,2,1): two rows, or two columns, or hook
    return len(lam) <= 2 or lam[0] <= 2 or (len(lam) >= 2 and lam[1] <= 1)

# ---------------------------------------------------------------- Lemma 1: r(lam - c_i) = r + 1 - [a_i=1] - [b_i=1]; d_2 = r(r+1) - A_1 - B_1
def check_lemma1(N):
    for n in range(1, N + 1):
        for lam in partitions(n):
            runs = corner_runs(lam); r = len(runs)
            corners = removable_corners(lam)
            assert len(corners) == r
            for i, row in enumerate(corners):
                a, b = runs[i]
                nu = remove_box(lam, row)
                assert r_of(nu) == r + 1 - (a == 1) - (b == 1), (lam, i)
            A1, B1 = A1B1(lam)
            dv = d_vector(lam)
            dj = lambda j: dv[j] if j < len(dv) else 0
            assert dj(1) == r
            assert dj(2) == r * (r + 1) - A1 - B1, lam
    print(f"Lemma 1 (corner count after removal; d_1=r; d_2=r(r+1)-A_1-B_1): OK for all partitions of n<={N}")

# ---------------------------------------------------------------- Lemma 2: six families
def family_members(lam):
    runs = corner_runs(lam)
    if len(runs) != 2: return set()
    (a1, b1), (a2, b2) = runs
    fam = set()
    if b1 == 1 and b2 == 1: fam.add('tworow')
    if a1 == 1 and a2 == 1: fam.add('twocol')
    if b1 == 1 and a2 == 1: fam.add('hook')
    if a1 == 1 and b2 == 1: fam.add('F2')
    if a1 == 1 and b1 == 1: fam.add('F1')
    if a2 == 1 and b2 == 1: fam.add('F4')
    return fam

def check_lemma2(N):
    for n in range(1, N + 1):
        for lam in partitions(n):
            if r_of(lam) != 2: continue
            A1, B1 = A1B1(lam)
            fam = family_members(lam)
            assert (A1 + B1 >= 2) == (len(fam) > 0), lam
            # explicit shapes
            (a1, b1), (a2, b2) = corner_runs(lam)
            if 'tworow' in fam: assert lam == (a1 + a2, a2)
            if 'twocol' in fam: assert lam == (2,) * b1 + (1,) * b2
            if 'hook' in fam: assert lam == (a1 + 1,) + (1,) * b2
            if 'F2' in fam: assert lam == (a2 + 1,) * b1 + (a2,)
            if 'F1' in fam: assert lam == (a2 + 1,) + (a2,) * b2
            if 'F4' in fam: assert lam == (a1 + 1,) * b1 + (1,)
            # membership in C
            assert in_C(lam) == bool(fam & {'tworow', 'twocol', 'hook'}), lam
    print(f"Lemma 2 (r=2 & A_1+B_1>=2 <=> one of six families; C = tworow|twocol|hook): OK for n<={N}")

# ---------------------------------------------------------------- Lemma 3: d_3 formulas and the d_3 bound
def d3_tworow(l1, l2):
    delta = l1 - l2
    assert delta >= 1 and l2 >= 1
    if delta >= 2 and l2 >= 2: return 8 - (delta == 2) - (l2 == 2)
    if delta == 1 and l2 >= 2: return 6 - (l2 == 2)
    if delta >= 2 and l2 == 1: return 4 - (delta == 2)
    return 2

def d3_hook(a, b):
    assert a >= 2 and b >= 1
    if a >= 3 and b >= 2: return 8 - (a == 3) - (b == 2)
    if a == 2 and b >= 2: return 4 - (b == 2)
    if a >= 3 and b == 1: return 4 - (a == 3)
    return 2

def check_lemma3(N):
    # all parameter pairs with n <= N
    for l2 in range(1, N):
        for l1 in range(l2 + 1, N - l2 + 1):
            assert d_vector((l1, l2))[3] == d3_tworow(l1, l2), (l1, l2)
    for a in range(2, N + 1):
        for b in range(1, N - a + 1):
            assert d_vector((a,) + (1,) * b)[3] == d3_hook(a, b), (a, b)
    for c in range(2, N + 1):
        for d in range(2, N + 1):
            if c + 1 + c * d > N: break
            assert d_vector((c + 1,) + (c,) * d)[3] == 10 - (c == 2) - (d == 2), ('F1', c, d)
    for b in range(2, N + 1):
        for c in range(2, N + 1):
            if (c + 1) * b + c > N: break
            assert d_vector((c + 1,) * b + (c,))[3] == 10 - (b == 2) - (c == 2), ('F2', b, c)
    for a in range(2, N + 1):
        for b in range(2, N + 1):
            if (a + 1) * b + 1 > N: break
            assert d_vector((a + 1,) * b + (1,))[3] == 10 - (a == 2) - (b == 2), ('F4', a, b)
    print(f"Lemma 3 (d_3 formulas for two-row, hook, F1, F2, F4): OK for all members with n <= {N}")

def check_lemma3b(N):
    # r=2 members of C have d_3 <= 8; r=2 non-C members of the six families have d_3 >= 9 except three shapes
    exc = set()
    for n in range(1, N + 1):
        for lam in partitions(n):
            if r_of(lam) != 2: continue
            fam = family_members(lam)
            if not fam: continue
            d3 = d_vector(lam)[3]
            if in_C(lam):
                assert d3 <= 8, lam
            else:
                if d3 <= 8: exc.add(lam)
    assert exc == {(3, 2, 2), (3, 3, 1), (3, 3, 2)}, exc
    for lam in exc: assert d_vector(lam)[3] == 8
    # C-members with d_3 = 8 and n in {7,8}
    c8 = [lam for n in (7, 8) for lam in partitions(n) if in_C(lam) and r_of(lam) == 2 and d_vector(lam)[3] == 8]
    assert set(c8) == {(4, 1, 1, 1), (5, 1, 1, 1), (4, 1, 1, 1, 1)}, c8
    assert f_hook((4, 1, 1, 1)) == 20 and f_hook((3, 2, 2)) == 21 and f_hook((3, 3, 1)) == 21
    assert f_hook((5, 1, 1, 1)) == 35 and f_hook((4, 1, 1, 1, 1)) == 35 and f_hook((3, 3, 2)) == 42
    print(f"Lemma 3b (d_3<=8 on C, >=9 off C except (3,2,2),(3,3,1),(3,3,2); exceptional f-values): OK for n<={N}")

# ---------------------------------------------------------------- Lemma 4: explicit d_j for two-row shapes (reflection) and hooks
def d_tworow(l1, l2, j):
    delta = l1 - l2
    k0 = max(0, -((delta - j) // 2))  # ceil((j-delta)/2)
    return sum(binom(j, k) - binom(j, j - delta - 1 - k) for k in range(k0, min(j, l2) + 1))

def d_tworow_sgn(l1, l2, j):  # the briefing/hint form
    delta = l1 - l2
    def sgn(x): return (x > 0) - (x < 0)
    return sum(comb(j, k) * sgn(delta + 1 - j + 2 * k) * (j - k <= l1 + 1) for k in range(0, min(j, l2) + 1))

def d_hook(p, q, j):  # hook (p+1, 1^q), j < n = p+q+1
    return sum(comb(j, i) for i in range(max(0, j - q), min(j, p) + 1))

def check_lemma4(N):
    for l2 in range(0, N + 1):
        for l1 in range(l2, N + 1):
            if l1 + l2 > N or l1 == 0: continue
            lam = (l1, l2) if l2 > 0 else (l1,)
            dv = d_vector(lam)
            for j in range(0, l1 + l2 + 1):
                assert dv[j] == d_tworow(l1, l2, j), (lam, j)
                assert dv[j] == d_tworow_sgn(l1, l2, j), (lam, j, 'sgn form')
    for p in range(0, N):
        for q in range(0, N - p):
            lam = (p + 1,) + (1,) * q
            n = p + q + 1
            dv = d_vector(lam)
            for j in range(0, n):
                assert dv[j] == d_hook(p, q, j), (lam, j)
            assert dv[n] == comb(n - 1, q) == f_hook(lam)
    print(f"Lemma 4 (two-row d_j formula, both forms; hook d_j formula, d_n=C(n-1,q)): OK for n<={N}")

# ---------------------------------------------------------------- Lemma 5: the invariant M and the values d_{M+1}, d_{M+2}
def M_of(dv):
    j = 0
    while j + 1 < len(dv) and dv[j + 1] == 2 ** (j + 1):
        j += 1
    return j

def check_lemma5(N):
    # two-row, l1 > l2 >= 1
    for l2 in range(1, N):
        for l1 in range(l2 + 1, N + 1 - l2):
            n = l1 + l2; delta = l1 - l2; M = min(delta, l2)
            dv = d_vector((l1, l2))
            assert M_of(dv) == M and M + 2 <= n, (l1, l2)
            if delta == l2:
                assert dv[M + 1] == 2 ** (M + 1) - 2 and dv[M + 2] == 2 ** (M + 2) - M - 5
            elif delta == M:  # delta < l2
                assert dv[M + 1] == 2 ** (M + 1) - 1
                assert dv[M + 2] == 2 ** (M + 2) - 2 - (l2 == M + 1)
            else:  # l2 == M < delta
                assert dv[M + 1] == 2 ** (M + 1) - 1
                assert dv[M + 2] == 2 ** (M + 2) - M - 3 - (delta == M + 1)
            assert dv[n] == comb(n, l2) - comb(n, l2 - 1) == comb(n - 1, l2) - binom(n - 1, l2 - 2)
    # hooks (p+1, 1^q), p >= q >= 1
    for q in range(1, N):
        for p in range(q, N - q):
            n = p + q + 1; M = q
            dv = d_vector((p + 1,) + (1,) * q)
            assert M_of(dv) == M and M + 1 <= n - 1
            assert dv[M + 1] == 2 ** (M + 1) - 1 - (p == q)
            if p >= 2:
                assert M + 2 <= n - 1
                if p >= q + 2: assert dv[M + 2] == 2 ** (M + 2) - M - 3
                elif p == q + 1: assert dv[M + 2] == 2 ** (M + 2) - M - 4
                else: assert dv[M + 2] == 2 ** (M + 2) - 2 * M - 6
    print(f"Lemma 5 (M = min(delta,l2) resp. min(p,q); values of d_(M+1), d_(M+2); f formulas): OK for n<={N}")

# ---------------------------------------------------------------- Lemma 6: rectangles
def check_lemma6(N):
    from young import f_hook as fh
    I = {}
    for j in range(0, N + 1):
        I[j] = sum(fh(rho) for rho in partitions(j))
    for a in range(1, N + 1):
        for b in range(1, N // a + 1):
            lam = (a,) * b; n = a * b
            dv = d_vector(lam)
            for j in range(0, n + 1):
                assert dv[j] == sum(fh(rho) for rho in partitions(j) if rho == () or (rho[0] <= a and len(rho) <= b)), (lam, j)
            m = min(a, b)
            assert all(dv[j] == I[j] for j in range(m + 1))
            if n >= m + 1: assert dv[m + 1] < I[m + 1]
            else: assert (a, b) == (1, 1)
    print(f"Lemma 6 (rectangles: d_j = sum of f^rho over rho in the box; m = min(a,b) read off): OK for ab<={N}")

# ---------------------------------------------------------------- Theorem: full statement by brute force
def check_theorem(N):
    for n in range(0, N + 1):
        table = {}
        for lam in partitions(n):
            table.setdefault(d_vector(lam), []).append(lam)
        for dv, group in table.items():
            for lam in group:
                if in_C(lam):
                    for mu in group:
                        assert mu == lam or mu == conjugate(lam), (lam, mu)
        print(f"  n={n}: theorem OK ({len(table)} distinct d-vectors among {sum(len(g) for g in table.values())} partitions)")
    print(f"THEOREM (lam in C, d(mu)=d(lam) => mu in {{lam, lam^t}}): OK for all n<={N}")

if __name__ == '__main__':
    NF = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    check_lemma1(22)
    check_lemma2(30)
    check_lemma3(40)
    check_lemma3b(30)
    check_lemma4(40)
    check_lemma5(40)
    check_lemma6(36)
    check_theorem(NF)
