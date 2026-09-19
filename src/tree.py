"""The tree SYT (truncated at depth N) and its partial-transpose automorphisms.

A vertex (a standard Young tableau) is encoded as the tuple of its boxes in
order of entry:  T = ((r_1,c_1), ..., (r_n,c_n))  meaning entry i sits in box
(r_i, c_i) (0-indexed).  The parent of T is T[:-1]; the children of T are
T + (box,) for the addable corners of shape(T).

tau_T   : whole-subtree partial transpose at a symmetric-shape vertex T.
theta_Tp: transpose only on the subtrees of the two children T+c, T+c^t
          (c an off-diagonal addable corner of the symmetric shape of T).
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from young import (Partition, addable_corners, conjugate, corner_boxes_addable,
                   is_symmetric)

Vertex = tuple  # tuple of (row, col) boxes


def shape_of_vertex(T: Vertex) -> Partition:
    rows: dict[int, int] = {}
    for (r, c) in T:
        rows[r] = rows.get(r, 0) + 1
    return tuple(rows[i] for i in range(len(rows)))


def children(T: Vertex) -> list[Vertex]:
    lam = shape_of_vertex(T)
    return [T + (box,) for box in corner_boxes_addable(lam)]


def transpose_vertex(T: Vertex) -> Vertex:
    return tuple((c, r) for (r, c) in T)


def build_tree(N: int):
    """All vertices of depth <= N, depth by depth; returns (levels, index)."""
    levels = [[()]]
    for d in range(N):
        nxt = []
        for T in levels[d]:
            nxt.extend(children(T))
        levels.append(nxt)
    index = {}
    for lev in levels:
        for T in lev:
            index[T] = len(index)
    return levels, index


def in_subtree(T: Vertex, V: Vertex) -> bool:
    """True iff V lies in the subtree rooted at T (V extends T)."""
    return len(V) >= len(T) and V[: len(T)] == T


def tau(T: Vertex, V: Vertex) -> Vertex:
    """tau_T applied to V."""
    if not in_subtree(T, V):
        return V
    return T + transpose_vertex(V[len(T):])


def theta(T: Vertex, pair: frozenset, V: Vertex) -> Vertex:
    """theta_{T,pair} applied to V, pair = {c, c^t} an off-diagonal mirror pair
    of addable corners (boxes) of the symmetric shape of T."""
    if not in_subtree(T, V) or len(V) == len(T) or V[len(T)] not in pair:
        return V
    return T + transpose_vertex(V[len(T):])


def mirror_pairs(lam: Partition) -> list[frozenset]:
    """Off-diagonal mirror pairs {c, c^t} of addable corners of a symmetric lam."""
    assert is_symmetric(lam)
    boxes = corner_boxes_addable(lam)
    pairs = set()
    for (r, c) in boxes:
        if r != c:
            pairs.add(frozenset([(r, c), (c, r)]))
    for p in pairs:
        for b in p:
            assert b in boxes
    return sorted(pairs, key=lambda p: sorted(p))


def symmetric_vertices(levels, max_depth):
    return [T for d in range(max_depth + 1) for T in levels[d] if is_symmetric(shape_of_vertex(T))]


def permutation_of(fn, levels, index):
    """Return the permutation (as a list) of all vertices of depth <= N induced by fn."""
    perm = [None] * len(index)
    for lev in levels:
        for T in lev:
            perm[index[T]] = index[fn(T)]
    assert sorted(perm) == list(range(len(index)))
    return perm


def is_automorphism(perm, levels, index):
    """Check perm preserves the parent relation (hence all edges) of the truncated tree."""
    for lev in levels[1:]:
        for T in lev:
            if perm[index[T[:-1]]] != index[levels_vertex_parent(perm, index, T, levels)]:
                return False
    return True


def levels_vertex_parent(perm, index, T, levels):
    # helper: parent of the image of T
    inv = _inverse_index(levels)
    return inv[perm[index[T]]][:-1]


_inv_cache = {}


def _inverse_index(levels):
    key = id(levels)
    if key not in _inv_cache:
        inv = []
        for lev in levels:
            inv.extend(lev)
        _inv_cache[key] = inv
    return _inv_cache[key]


def check_edge_preserving(perm, levels, index) -> bool:
    inv = _inverse_index(levels)
    for lev in levels[1:]:
        for T in lev:
            img = inv[perm[index[T]]]
            if img[:-1] != inv[perm[index[T[:-1]]]]:
                return False
    return True


# ---------------------------------------------------------------------------
# The parity homomorphisms E_{lam,p}
# ---------------------------------------------------------------------------

def parity_E(perm, levels, index, lam: Partition, pair: frozenset) -> int:
    """E_{lam,pair}(g) = #{T of shape lam : g(T + c) = g(T) + c^t} mod 2,
    for c in pair.  Requires depth(lam)+1 <= N."""
    inv = _inverse_index(levels)
    d = sum(lam)
    c = sorted(pair)[0]
    ct = (c[1], c[0])
    total = 0
    for T in levels[d]:
        if shape_of_vertex(T) != lam:
            continue
        gT = inv[perm[index[T]]]
        assert shape_of_vertex(gT) == lam
        gTc = inv[perm[index[T + (c,)]]]
        assert gTc[:-1] == gT
        if gTc[-1] == ct:
            total += 1
        else:
            assert gTc[-1] == c
    return total % 2


def compose(p, q):
    """(p o q)[i] = p[q[i]]"""
    return [p[q[i]] for i in range(len(q))]
