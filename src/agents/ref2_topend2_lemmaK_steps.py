"""Referee #2: verify the interval-level claims (a),(b),(c) in the proof of Lemma K (topend.md T4.1) directly.
1-indexed rows/cols as in the proof. For every nu |- m <= 16 and every addable corner (i, lam_i+1):
  - build intervals of constant lam'_j on 1<=j<=lam_i;
  - (a) (L, j2) is a removable corner with content j2-L and h_{i j2} = x - y;
  - (b) (L+1, j1) is an addable corner with content j1-L-1 and h_{i j1}+1 = x - x';
  - (c) the set of corners from (a) == removable corners in rows >= i == {y in Y : y < x};
        the set from (b) == addable corners in rows > i == {x' in X : x' < x};
  - column part: the boxes (i', lam_i+1), i' < i, are exactly those whose hook changes besides row i;
    the product over them equals prod_{x'>x}(x'-x)/prod_{y>x}(y-x)."""
import os
import sys
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from ref2_topend2_check import parts, conj, rem_corners, add_corners, hooks, ad

bad = 0; cnt = 0
for m in range(0, 17):
    for nu in parts(m):
        ell = len(nu); lamc = conj(nu)
        lam = lambda r: nu[r-1] if 1 <= r <= ell else 0          # 1-indexed row lengths
        lamp = lambda c: lamc[c-1] if 1 <= c <= len(lamc) else 0  # 1-indexed column lengths
        h = {(r+1, c+1): v for (r, c), v in hooks(nu).items()}
        remc = [(r+1, nu[r]) for r in rem_corners(nu)]               # (row, col) 1-indexed
        addc = [(r+1, lam(r+1)+1) for r in add_corners(nu)]
        Y = {c-r: (r, c) for (r, c) in remc}; X = {c-r: (r, c) for (r, c) in addc}
        for (i, jx) in addc:
            cnt += 1
            x = jx - i
            assert lam(i-1) >= lam(i)+1 or i == 1
            # intervals of constant lam'_j on 1..lam_i
            intervals = []
            j = 1
            while j <= lam(i):
                L = lamp(j); j1 = j
                while j+1 <= lam(i) and lamp(j+1) == L: j += 1
                intervals.append((L, j1, j)); j += 1
            Ls = [t[0] for t in intervals]
            if Ls != sorted(Ls, reverse=True) or len(set(Ls)) != len(Ls): bad += 1; print("interval order FAIL", nu, i)
            got_a = set(); got_b = set()
            for (L, j1, j2) in intervals:
                if L < i: bad += 1; print("L<i FAIL", nu, i)
                # telescoping check
                prod = Fraction(1)
                for j in range(j1, j2+1): prod *= Fraction(h[(i, j)]+1, h[(i, j)])
                if prod != Fraction(h[(i, j1)]+1, h[(i, j2)]): bad += 1; print("telescope FAIL", nu, i, L)
                # (a)
                if (L, j2) not in remc: bad += 1; print("(a) removable FAIL", nu, i, (L, j2))
                y = j2 - L
                if h[(i, j2)] != x - y: bad += 1; print("(a) hook FAIL", nu, i)
                got_a.add((L, j2))
                # (b)
                if (L+1, j1) not in addc: bad += 1; print("(b) addable FAIL", nu, i, (L+1, j1))
                xp = j1 - L - 1
                if h[(i, j1)] + 1 != x - xp: bad += 1; print("(b) hook FAIL", nu, i)
                got_b.add((L+1, j1))
            # (c)
            exp_a = {(r, c) for (r, c) in remc if r >= i}
            exp_a2 = {Y[y] for y in Y if y < x}
            exp_b = {(r, c) for (r, c) in addc if r > i}
            exp_b2 = {X[xp] for xp in X if xp < x}
            if not (got_a == exp_a == exp_a2): bad += 1; print("(c) removable set FAIL", nu, i, got_a, exp_a, exp_a2)
            if not (got_b == exp_b == exp_b2): bad += 1; print("(c) addable set FAIL", nu, i, got_b, exp_b, exp_b2)
            if any(y == x for y in Y) or any(y > x for y in Y if Y[y][0] >= i): bad += 1; print("(c) content FAIL", nu, i)
            # column part
            nu2 = ad(nu, i-1); h2 = {(r+1, c+1): v for (r, c), v in hooks(nu2).items()}
            changed = {b for b in h if h2[b] != h[b]}
            exp_changed = {(i, j) for j in range(1, lam(i)+1)} | {(ip, jx) for ip in range(1, i)}
            if changed - exp_changed or (exp_changed - changed - {b for b in exp_changed if h2[b] == h[b]}): pass
            if changed != exp_changed: bad += 1; print("changed-set FAIL", nu, i, changed, exp_changed)
            if h2[(i, jx)] != 1: bad += 1; print("new box hook FAIL", nu, i)
            Cc = Fraction(1)
            for ip in range(1, i): Cc *= Fraction(h[(ip, jx)]+1, h[(ip, jx)])
            Cp = Fraction(1)
            for xp in X:
                if xp > x: Cp *= (xp - x)
            for y in Y:
                if y > x: Cp /= (y - x)
            if Cc != Cp: bad += 1; print("column product FAIL", nu, i, Cc, Cp)
print("Lemma K proof steps (a),(b),(c), telescoping, column part checked for", cnt, "addable corners, nu |- m<=16; failures:", bad)
