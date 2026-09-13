"""Verify every explicit formula used in the proof of Theorem 4.9 (class C):
 first-deviation index m and deficit for two-row shapes and hooks, and the d_{m+2} values in
 cases T1,T2,T3 (two-row) and H1,H2,H3 (hooks), plus the j = delta+1 comparison of (l1, m) vs (l1, 1^m)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from census import d_vector

def check():
    bad = 0
    N = 34
    # two-row shapes (l1 > l2 >= 1)
    for l2 in range(1, N):
        for l1 in range(l2 + 1, N - l2 + 1):
            n = l1 + l2; dv = d_vector((l1, l2)); delta = l1 - l2; m = min(delta, l2)
            for j in range(0, m + 1):
                if dv[j] != 2 ** j: bad += 1; print("T: 2^j fails", (l1, l2), j)
            if m + 1 <= n and dv[m + 1] != 2 ** (m + 1) - 1 - (delta == l2): bad += 1; print("T deficit fails", (l1, l2))
            if m + 2 <= n:
                if l2 == m and delta > m: exp = 2 ** (m + 2) - m - 3 - (delta == m + 1)      # T1
                elif delta == m and l2 > m: exp = 2 ** (m + 2) - 2 - (l2 == m + 1)           # T2
                else: exp = 2 ** (m + 2) - m - 5                                              # T3
                if dv[m + 2] != exp: bad += 1; print("T d_{m+2} fails", (l1, l2), dv[m + 2], exp)
    # hooks (a, 1^b), a >= 2, b >= 1
    for a in range(2, N):
        for b in range(1, N - a + 1):
            n = a + b; dv = d_vector((a,) + (1,) * b); m = min(a - 1, b)
            for j in range(0, m + 1):
                if dv[j] != 2 ** j: bad += 1; print("H: 2^j fails", (a, b), j)
            if m + 1 < n and dv[m + 1] != 2 ** (m + 1) - 1 - (a - 1 == b): bad += 1; print("H deficit fails", (a, b))
            if m + 2 < n:
                if a - 1 == m and b > m: exp = 2 ** (m + 2) - (m + 2) - 1 - (b == m + 1)      # H1
                elif b == m and a - 1 > m: exp = 2 ** (m + 2) - 1 - (m + 2) - (a - 1 == m + 1)  # H2
                else: exp = 2 ** (m + 2) - 2 * m - 6                                          # H3
                if dv[m + 2] != exp: bad += 1; print("H d_{m+2} fails", (a, b), dv[m + 2], exp)
    # (l1, m) vs hook (l1, 1^m), m >= 2: differ by exactly 1 at j = l1 - m + 1
    for m in range(2, 16):
        for l1 in range(m + 1, 34 - m):
            j = l1 - m + 1
            if d_vector((l1,) + (1,) * m)[j] - d_vector((l1, m))[j] != 1: bad += 1; print("delta+1 comparison fails", l1, m)
    print("class-C case formulas:", "OK" if bad == 0 else f"{bad} failures")
check()
