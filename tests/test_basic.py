"""Hand-verified small cases for every function in young.py / census.py."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from young import (add_box, addable_corners, conjugate, corner_runs, f_hook,
                   f_rec, f_skew, from_corner_runs, involutions, is_symmetric,
                   partitions, remove_box, removable_corners, shape_of,
                   syt_of_shape, is_syt, transpose_tableau, hook_lengths)
from census import (d_vector, d_vector_by_skew, s_up, s_up_by_skew, s_from_d,
                    check_identity)


def test_partitions_count():
    counts = [len(list(partitions(n))) for n in range(11)]
    assert counts == [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42]
    assert list(partitions(4)) == [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]


def test_conjugate():
    assert conjugate((3, 1)) == (2, 1, 1)
    assert conjugate((2, 2)) == (2, 2)
    assert conjugate(()) == ()
    assert conjugate((4, 2, 1)) == (3, 2, 1, 1)
    assert is_symmetric((3, 2, 1)) and is_symmetric((2, 1)) and not is_symmetric((3, 1))
    for n in range(8):
        for lam in partitions(n):
            assert conjugate(conjugate(lam)) == lam


def test_corners():
    assert removable_corners((3, 1)) == [0, 1]
    assert removable_corners((2, 2)) == [1]
    assert removable_corners((3, 2, 1)) == [0, 1, 2]
    assert addable_corners((3, 2, 1)) == [0, 1, 2, 3]
    assert addable_corners(()) == [0]
    assert addable_corners((2, 2)) == [0, 2]
    assert remove_box((2, 2), 1) == (2, 1)
    assert remove_box((1,), 0) == ()
    assert add_box((2, 2), 2) == (2, 2, 1)
    assert add_box((), 0) == (1,)
    # number of addable corners = number of removable corners + 1
    for n in range(9):
        for lam in partitions(n):
            assert len(addable_corners(lam)) == len(removable_corners(lam)) + 1


def test_corner_runs():
    assert corner_runs((3, 2, 1)) == [(1, 1), (1, 1), (1, 1)]
    assert corner_runs((3, 3, 2, 2, 2, 2)) == [(1, 2), (2, 4)]
    assert corner_runs((4, 4)) == [(4, 2)]
    for n in range(9):
        for lam in partitions(n):
            runs = corner_runs(lam)
            assert from_corner_runs(runs) == lam
            assert from_corner_runs([(b, a) for (a, b) in reversed(runs)]) == conjugate(lam)


def test_hooks_and_f():
    assert hook_lengths((3, 2)) == {(0, 0): 4, (0, 1): 3, (0, 2): 1, (1, 0): 2, (1, 1): 1}
    assert f_hook((3, 2)) == 5 and f_rec((3, 2)) == 5
    assert f_hook((2, 2)) == 2 and f_rec((2, 2)) == 2
    assert f_hook((2, 1)) == 2
    assert f_hook((3, 2, 1)) == 16
    assert f_hook((4, 4)) == 14
    assert f_hook(()) == 1
    for n in range(10):
        total = 0
        for lam in partitions(n):
            assert f_hook(lam) == f_rec(lam) == len(syt_of_shape(lam))
            assert f_hook(lam) == f_hook(conjugate(lam))
            total += f_hook(lam)
        assert total == involutions(n)


def test_syt_explicit():
    for T in syt_of_shape((3, 2)):
        assert is_syt(T) and shape_of(T) == (3, 2)
        assert is_syt(transpose_tableau(T)) and shape_of(transpose_tableau(T)) == (2, 2, 1)
    assert involutions(0) == 1 and involutions(4) == 10 and involutions(6) == 76


def test_f_skew():
    assert f_skew((3, 2), (1,)) == 5
    assert f_skew((3, 2), (2,)) == 3      # chains (2)->(3)->(3,1)->(3,2), (2)->(2,1)->(3,1)->(3,2), (2)->(2,1)->(2,2)->(3,2)
    assert f_skew((3, 2), (1, 1)) == 2    # (1,1)->(2,1)->(3,1)->(3,2), (1,1)->(2,1)->(2,2)->(3,2)
    assert f_skew((2, 2), (2,)) == 1
    assert f_skew((2, 2), (1, 1)) == 1
    assert f_skew((2, 2), (3,)) == 0
    assert f_skew((3, 3), (2, 1)) == 2
    assert f_skew((3, 3), (3, 3)) == 1
    for n in range(7):
        for lam in partitions(n):
            assert f_skew(lam, ()) == f_hook(lam)
            for m in range(n + 1):
                # f^lam = sum_{nu |- m} f^nu f^{lam/nu}
                assert sum(f_hook(nu) * f_skew(lam, nu) for nu in partitions(m)) == f_hook(lam)


def test_d_vectors_size4():
    expected = {
        (4,): (1, 1, 1, 1, 1),
        (3, 1): (1, 2, 3, 3, 3),
        (2, 2): (1, 1, 2, 2, 2),
        (2, 1, 1): (1, 2, 3, 3, 3),
        (1, 1, 1, 1): (1, 1, 1, 1, 1),
    }
    for lam, dv in expected.items():
        assert d_vector(lam) == dv, (lam, d_vector(lam))
        assert d_vector_by_skew(lam) == dv
    assert d_vector(()) == (1,)
    assert d_vector((1,)) == (1, 1)
    assert d_vector((2, 1)) == (1, 2, 2, 2)
    for n in range(9):
        for lam in partitions(n):
            dv = d_vector(lam)
            assert dv == d_vector_by_skew(lam)
            assert dv == d_vector(conjugate(lam))
            assert dv[0] == 1
            if n >= 1:
                assert dv[1] == len(removable_corners(lam))
                assert dv[n] == f_hook(lam) and dv[n - 1] == f_hook(lam)
            if n >= 2:
                assert dv[n - 2] == f_hook(lam)


def test_s_up_small():
    assert s_up((2, 1), 1) == 3
    assert s_up((), 1) == 1 and s_up((), 2) == 2 and s_up((), 3) == 4 and s_up((), 4) == 10
    assert s_up((1,), 1) == 2
    assert s_up((1,), 2) == 4   # (2)->(3),(2,1); (1,1)->(2,1),(1,1,1)
    assert s_up((2, 2), 1) == 2
    for n in range(6):
        for lam in partitions(n):
            for k in range(5):
                assert s_up(lam, k) == s_up_by_skew(lam, k)
                assert s_up(lam, k) == s_up(conjugate(lam), k)


def test_identity_small():
    assert s_from_d((), 5) == involutions(5) == 26
    assert s_from_d((1,), 2) == 4
    assert check_identity(6, 6)


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
