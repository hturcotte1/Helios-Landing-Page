"""Hand-verified small cases for src/shifted.py (shifted Young lattice)."""
import os
import sys
from fractions import Fraction
from math import factorial

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from shifted import (E_series, P_series, Pw_series, add_box, addable_boxes,
                     boxes, chain_to_tableau, commutator_DU_UD, d_vector,
                     down_set, down_set_weighted, dw_vector, egf, g_hook,
                     g_rec, g_skew, is_shifted_syt, is_strict_partition,
                     removable_boxes, remove_box, s_up, s_up_by_g, s_vector,
                     series_div, series_exp_poly, series_log, series_mul,
                     series_scale, shifted_chains, strict_partitions,
                     total_shifted_syt, up_set)


def test_strict_partitions():
    assert list(strict_partitions(6)) == [(6,), (5, 1), (4, 2), (3, 2, 1)]
    assert list(strict_partitions(0)) == [()]
    assert list(strict_partitions(1)) == [(1,)]
    assert list(strict_partitions(2)) == [(2,)]
    assert list(strict_partitions(3)) == [(3,), (2, 1)]
    counts = [len(list(strict_partitions(n))) for n in range(11)]
    assert counts == [1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10]
    for n in range(12):
        for lam in strict_partitions(n):
            assert is_strict_partition(lam) and sum(lam) == n
    assert not is_strict_partition((2, 2))


def test_boxes_and_covers():
    assert boxes((3, 1)) == [(0, 0), (0, 1), (0, 2), (1, 1)]
    assert boxes(()) == []
    # (1): only row 1 is addable (a new row needs last part >= 2)
    assert addable_boxes((1,)) == [(0, (0, 1), False)]
    assert addable_boxes(()) == [(0, (0, 0), True)]
    assert addable_boxes((2,)) == [(0, (0, 2), False), (1, (1, 1), True)]
    assert addable_boxes((2, 1)) == [(0, (0, 2), False)]
    assert addable_boxes((3, 1)) == [(0, (0, 3), False), (1, (1, 2), False)]
    assert addable_boxes((4, 2)) == [(0, (0, 4), False), (1, (1, 3), False), (2, (2, 2), True)]
    assert removable_boxes((3, 1)) == [(0, (0, 2), False), (1, (1, 1), True)]
    assert removable_boxes((2, 1)) == [(1, (1, 1), True)]
    assert removable_boxes((3, 2)) == [(1, (1, 2), False)]
    assert removable_boxes((1,)) == [(0, (0, 0), True)]
    assert removable_boxes(()) == []
    assert up_set((2, 1)) == [(3, 1)]
    assert up_set((3,)) == [(4,), (3, 1)]
    assert down_set((3, 2, 1)) == [(3, 2)]
    assert down_set((4, 2)) == [(3, 2), (4, 1)]
    assert down_set_weighted((3, 1)) == [((2, 1), 2), ((3,), 1)]
    assert add_box((2,), 1) == (2, 1) and remove_box((1,), 0) == ()
    for n in range(10):
        for lam in strict_partitions(n):
            for i, (r, c), diag in addable_boxes(lam):
                mu = add_box(lam, i)
                assert is_strict_partition(mu) and sum(mu) == n + 1
                assert set(boxes(mu)) - set(boxes(lam)) == {(r, c)}
                assert diag == (r == c)
                assert lam in down_set(mu)
            for i, (r, c), diag in removable_boxes(lam):
                mu = remove_box(lam, i)
                assert is_strict_partition(mu) and sum(mu) == n - 1
                assert set(boxes(lam)) - set(boxes(mu)) == {(r, c)}
                assert diag == (r == c)
                assert lam in up_set(mu)
            # #addable = #removable + 1 - [last part == 1]
            last1 = 1 if (lam and lam[-1] == 1) else 0
            assert len(addable_boxes(lam)) == len(removable_boxes(lam)) + 1 - last1
            # exactly one addable box iff lam is a staircase (l, l-1, ..., 1) or empty
            stair = lam == tuple(range(len(lam), 0, -1))
            assert (len(addable_boxes(lam)) == 1) == stair


def test_g_hand_values():
    # (3,1): chains (1)(2)(3)(3,1) and (1)(2)(2,1)(3,1)
    assert g_hook((3, 1)) == 2 and g_rec((3, 1)) == 2
    assert sorted(shifted_chains((3, 1))) == sorted([
        ((), (1,), (2,), (3,), (3, 1)), ((), (1,), (2,), (2, 1), (3, 1))])
    # (2,1): the unique chain (1)(2)(2,1)
    assert g_hook((2, 1)) == 1 and g_rec((2, 1)) == 1
    # (3,2,1): chains through (3,2) = (3,1)->(3,2) with g^(3,1) = 2, so 2
    assert g_hook((3, 2, 1)) == 2 and g_rec((3, 2, 1)) == 2
    # (4,2): g = g^(3,2) + g^(4,1) = 2 + 3 = 5  [g^(3,2)=g^(3,1)=2; g^(4,1)=g^(3,1)+g^(4)=3]
    assert g_hook((4, 2)) == 5 and g_rec((4, 2)) == 5
    assert g_hook((4, 1)) == 3 and g_hook((3, 2)) == 2 and g_hook((5, 1)) == 4
    assert g_hook(()) == 1 and g_hook((7,)) == 1
    assert g_hook((4, 3, 2, 1)) == 12   # 10!/(4!3!2!1!) * (1/7)(2/6)(3/5)(1/5)(2/4)(1/3)


def test_g_recursion_vs_formula_vs_chains():
    for n in range(11):
        for lam in strict_partitions(n):
            assert g_hook(lam) == g_rec(lam) == g_skew(lam, ())
            if n <= 8:
                chains = shifted_chains(lam)
                assert len(chains) == g_hook(lam)
                for ch in chains:
                    assert is_shifted_syt(chain_to_tableau(ch), lam)


def test_sum_2_pow_g_squared_is_factorial():
    for n in range(9):
        tot = sum(2 ** (n - len(lam)) * g_hook(lam) ** 2 for lam in strict_partitions(n))
        assert tot == factorial(n), (n, tot)


def test_total_shifted_syt():
    assert [total_shifted_syt(k) for k in range(9)] == [1, 1, 1, 2, 3, 6, 12, 27, 63]
    for k in range(10):
        assert s_up((), k) == total_shifted_syt(k)


def test_s_up_hand_values():
    # from (1): (2); then (3),(2,1); then (4),(3,1) from (3) and (3,1) from (2,1)
    assert s_vector((1,), 3) == (1, 1, 2, 3)
    assert s_vector((2, 1), 2) == (1, 1, 2)      # (3,1) then (4,1),(3,2)
    assert s_up((3, 2), 1) == 2                   # (4,2), (3,2,1)
    assert s_up((2,), 2) == 3                     # (3)->(4),(3)->(3,1),(2,1)->(3,1)
    for n in range(7):
        for lam in strict_partitions(n):
            for k in range(6):
                assert s_up(lam, k) == s_up_by_g(lam, k)


def test_down_censuses():
    assert d_vector(()) == (1,)
    assert d_vector((1,)) == (1, 1)
    assert d_vector((2, 1)) == (1, 1, 1, 1)
    assert d_vector((3, 1)) == (1, 2, 2, 2, 2)    # (3,1)->(2,1)|(3)->(2)->(1)->()
    assert dw_vector((3, 1)) == (1, 3, 4, 8, 8)   # weights: (2,1) via 2, (3) via 1; ...
    assert dw_vector((1,)) == (1, 1)
    assert dw_vector((2,)) == (1, 2, 2)
    assert dw_vector((2, 1)) == (1, 1, 2, 2)
    for n in range(10):
        for lam in strict_partitions(n):
            d = d_vector(lam)
            dw = dw_vector(lam)
            assert d[0] == dw[0] == 1
            if n:
                assert d[1] == len(removable_boxes(lam))
                assert d[n] == d[n - 1] == g_hook(lam)
                assert dw[1] == sum(w for _, w in down_set_weighted(lam))
                # every chain to () removes exactly l(lam) diagonal boxes
                assert dw[n] == 2 ** (n - len(lam)) * g_hook(lam)
                assert all(dw[j] >= d[j] for j in range(n + 1))
            # d_j via explicit chains from every nu contained in lam
            for j in range(n + 1):
                assert d[j] == sum(g_skew(lam, nu) for nu in strict_partitions(n - j))


def test_fomin_relations():
    for n in range(11):
        for lam in strict_partitions(n):
            assert commutator_DU_UD(lam, True) == {lam: 1}
            unw = commutator_DU_UD(lam, False)
            if lam and lam[-1] == 1:
                assert unw == {}
            else:
                assert unw == {lam: 1}


def test_series_helpers():
    K = 6
    one = [Fraction(1)] + [Fraction(0)] * K
    e = series_exp_poly(1, Fraction(1, 2), K)   # exp(t + t^2/2) = EGF of involutions
    assert [e[k] * factorial(k) for k in range(K + 1)] == [1, 1, 2, 4, 10, 26, 76]
    assert series_mul(e, series_div(one, e, K), K) == one
    assert series_log(e, K) == [0, 1, Fraction(1, 2), 0, 0, 0, 0]
    assert series_scale(e, 2, K)[2] == 4 * e[2]
    assert E_series((), K) == egf([total_shifted_syt(k) for k in range(K + 1)], K)
    assert P_series((1,), K) == [1, 1, 0, 0, 0, 0, 0]
    assert Pw_series((2,), K) == [1, 2, 1, 0, 0, 0, 0]


def test_truncated_census_dp_matches_brute_force():
    from shifted_scan import truncated_census, c_up
    cen = truncated_census(7, 6)
    for n in range(8):
        for lam in strict_partitions(n):
            assert cen[lam] == s_vector(lam, 6)
    assert c_up((), 3) == 1                      # shapes of size 3 ending in 1: (2,1) only, g=1
    assert c_up((2,), 1) == 1 and c_up((3, 1), 1) == 1


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
