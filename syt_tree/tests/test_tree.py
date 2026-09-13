"""Tests for tree.py: the truncated tree SYT, tau/theta automorphisms, parity invariants."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from tree import (build_tree, children, shape_of_vertex, transpose_vertex, tau, theta,
                  mirror_pairs, symmetric_vertices, permutation_of, check_edge_preserving,
                  parity_E, compose)
from young import involutions, is_symmetric, f_skew


def test_build_tree_counts():
    levels, index = build_tree(6)
    assert [len(l) for l in levels] == [involutions(k) for k in range(7)]  # 1,1,2,4,10,26,76
    assert len(index) == sum(involutions(k) for k in range(7))
    # children have distinct shapes, one per addable corner
    T = ((0, 0), (0, 1), (1, 0))  # shape (2,1)
    ch = children(T)
    assert len(ch) == 3 and len({shape_of_vertex(c) for c in ch}) == 3
    assert shape_of_vertex(T) == (2, 1)
    assert transpose_vertex(T) == ((0, 0), (1, 0), (0, 1))


def test_tau_theta_are_automorphisms():
    N = 6
    levels, index = build_tree(N)
    for T in symmetric_vertices(levels, N - 1):
        perm = permutation_of(lambda V, T=T: tau(T, V), levels, index)
        assert check_edge_preserving(perm, levels, index)
        assert compose(perm, perm) == list(range(len(index)))  # involution
        for p in mirror_pairs(shape_of_vertex(T)):
            q = permutation_of(lambda V, T=T, p=p: theta(T, p, V), levels, index)
            assert check_edge_preserving(q, levels, index)
            assert compose(q, q) == list(range(len(index)))


def test_mirror_pairs():
    assert mirror_pairs((3, 2, 1)) == [frozenset({(0, 3), (3, 0)}), frozenset({(1, 2), (2, 1)})]
    assert mirror_pairs((2, 1)) == [frozenset({(0, 2), (2, 0)})]   # (1,1) is diagonal
    assert mirror_pairs(()) == []                                   # only addable corner (0,0) is diagonal
    assert mirror_pairs((1,)) == [frozenset({(0, 1), (1, 0)})]


def test_parity_invariant_321():
    N = 7
    levels, index = build_tree(N)
    lam = (3, 2, 1)
    p1, p2 = mirror_pairs(lam)
    # every tau_T has equal parity for the two pairs, equal to f^{lam/lam0} mod 2
    for T in symmetric_vertices(levels, N - 1):
        perm = permutation_of(lambda V, T=T: tau(T, V), levels, index)
        e1, e2 = parity_E(perm, levels, index, lam, p1), parity_E(perm, levels, index, lam, p2)
        assert e1 == e2 == f_skew(lam, shape_of_vertex(T)) % 2
    # a single pair transpose has different parities
    T = [T for T in levels[6] if shape_of_vertex(T) == lam][0]
    q = permutation_of(lambda V: theta(T, p1, V), levels, index)
    assert parity_E(q, levels, index, lam, p1) == 1
    assert parity_E(q, levels, index, lam, p2) == 0


def test_tau_is_product_of_thetas_and_diagonal_tau():
    # tau_T = prod_p theta_{T,p} * tau_{T+c} for the diagonal corner c, on the truncated tree
    N = 6
    levels, index = build_tree(N)
    T = ((0, 0), (0, 1), (1, 0))  # shape (2,1): pair {(0,2),(2,0)} and diagonal corner (1,1)
    tau_perm = permutation_of(lambda V: tau(T, V), levels, index)
    prod = list(range(len(index)))
    for p in mirror_pairs((2, 1)):
        prod = compose(prod, permutation_of(lambda V, p=p: theta(T, p, V), levels, index))
    Td = T + ((1, 1),)
    assert is_symmetric(shape_of_vertex(Td))
    prod = compose(prod, permutation_of(lambda V: tau(Td, V), levels, index))
    assert prod == tau_perm


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
