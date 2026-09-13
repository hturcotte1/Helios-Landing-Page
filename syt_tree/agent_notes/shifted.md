# Shifted case (Worley Problem 2): tree SSYT of shifted standard Young tableaux

All paths are under `/home/user/Helios-Landing-Page/syt_tree/`. Exact integer / `fractions.Fraction`
arithmetic everywhere. No existing file was modified.

## Files

| file | content |
|---|---|
| `src/shifted.py` | library: strict partitions, shifted covers up/down (box coordinates + diagonal flag), `g_rec`/`g_hook` (recursion vs shifted hook formula), explicit chains/tableaux, `s_up` (brute-force up-census), `d_vector` (down-census), `dw_vector` (weighted down-census d'), EGF helpers over `Fraction`, `commutator_DU_UD` |
| `tests/test_shifted.py` | 11 pytest tests, all hand-verified values (see Task 1); `python3 -m pytest tests/test_shifted.py` -> 11 passed (together with the old suite: 21 passed) |
| `src/shifted_scan.py` | modes `identity`, `fomin`, `census N K`, `siblings N K`; writes the `data/shifted_*.txt` files below |
| `data/shifted_identity.txt` | Task 2 output (`identity 8 10`) |
| `data/shifted_fomin.txt` | Fomin relations, all strict partitions of n <= 10 |
| `data/shifted_census_N40_K14.txt` | Task 3 (n <= 40, K = 14) |
| `data/shifted_census_N60_K21.txt` | Task 3 extension (n <= 60, K = 21) |
| `data/shifted_siblings_N39_K14.txt` | Task 4 (n <= 39, K = 14) |

Conventions (0-indexed in code, 1-indexed in the task): row i of the shifted diagram occupies boxes
(i, i), ..., (i, i+lam_i-1); a box is diagonal iff column == row. Adding a box to an existing row is
never diagonal; the new-row box (l, l) is always diagonal. Removing the last box of row i is diagonal iff
lam_i = 1 (only possible for the last row). Cover rules exactly as in the task statement; the tests check
that every cover changes the box set by exactly one box and that the diagonal flag agrees with (r == c).

## Task 1 — implementation and hand-checked values (tests/test_shifted.py)

* Strict partitions of 6: (6), (5,1), (4,2), (3,2,1). Counts for n = 0..10: 1,1,1,2,2,3,4,5,6,8,10.
* g^(3,1) = 2 with the two chains ()(1)(2)(3)(3,1) and ()(1)(2)(2,1)(3,1) (explicit enumeration
  `shifted_chains` is compared with the list). g^(2,1) = 1, g^(3,2,1) = 2, g^(4,1) = 3, g^(3,2) = 2,
  **g^(4,2) = 5** (= g^(3,2) + g^(4,1) = 2 + 3; the formula gives 6!/(4!2!) * 2/6 = 5, so the "5?" in the
  task is right), g^(5,1) = 4, g^(4,3,2,1) = 12.
* Recursion == shifted hook formula == chain count (g_skew from ()) for all strict lam, n <= 10; explicit
  chains converted to fillings are shifted SYT (rows/columns increasing) for n <= 8.
* sum_{strict lam |- n} 2^{n-l(lam)} (g^lam)^2 = n! for n <= 8.
* s_k(()) = total number of shifted SYT of size k = 1,1,1,2,3,6,12,27,63,154,398,1055 (k = 0..11);
  brute-force s_k agrees with the independent computation sum_{|mu| = n+k} g^{mu/lam} (n <= 6, k <= 5).
* d_j: d_1 = #removable, d_n = d_{n-1} = g^lam; d_j = sum_{nu |- n-j} g^{lam/nu} (cross-check).
  d'_j (weighted, off-diagonal removal = 2, diagonal = 1): d'_n = 2^{n-l(lam)} g^lam.
  Hand values: d(3,1) = (1,2,2,2,2), d'(3,1) = (1,3,4,8,8), d'(2) = (1,2,2), d'(2,1) = (1,1,2,2).
* Structural facts verified for n <= 9: #addable(lam) = #removable(lam) + 1 - [lam_l = 1];
  lam has exactly ONE addable box iff lam is empty or a staircase (l, l-1, ..., 1).

## Task 2 — identity search (data/shifted_identity.txt; all strict lam with n <= 8, series mod t^11)

Notation: E_lam(t) = sum_k s_k t^k/k!, P_lam = sum_j d_j t^j/j!, P'_lam = sum_j d'_j t^j/j!.

**Nothing of the requested form holds.** Every candidate fails already for lam = (1) at the coefficient of t^1:

| candidate | result | first failure |
|---|---|---|
| E_lam / P_lam independent of lam | NO (24/25 shapes differ from lam=()) | (1): t^1 — E_()/P_() = 1 + t + t^2/2 + t^3/3 + ..., E_(1)/P_(1) = 1 + 0 t + t^2 - t^3/2 + ... |
| E_lam / P'_lam independent of lam | NO (24/25) | (1): t^1 |
| E_lam = E_()(t) P_lam(t) | NO (24/25) | (1), t^1: LHS 1, RHS 2 |
| E_lam = E_()(t) P_lam(t/2) | NO | (1), t^1: LHS 1, RHS 3/2 |
| E_lam = E_()(t) P_lam(2t) | NO | (1), t^1: LHS 1, RHS 3 |
| E_lam = E_()(t) P'_lam(t) | NO | (1), t^1: LHS 1, RHS 2 |
| E_lam = E_()(t) P'_lam(t/2) | NO | (1), t^1: LHS 1, RHS 3/2 |
| E_lam = E_()(t) P'_lam(2t) | NO | (1), t^1: LHS 1, RHS 3 |
| same with an *arbitrary* universal Phi (ratio E_lam/Q_lam(ct) lam-independent), all 6 combinations | NO | (1), t^1 in every case |

(The universal factor is forced: Q_() = 1 for both P and P', so Phi = E_().)

**Structural obstruction (so no identity of this shape exists at all):** since () has the unique cover (1),
s_k((1)) = s_{k+1}(()), i.e. E_(1) = E_()'. Hence E_(1)/E_() = (log E_())' =
1 + t^2/2 - t^3/3 + 5t^4/24 - 19t^5/120 + 89t^6/720 - 227t^7/2520 + ... , which is not a polynomial
(nonzero coefficients through t^18 checked). Therefore there is NO universal Phi(t) and NO family of
polynomials Q_lam with deg Q_lam <= |lam| (whatever Q_lam is — P, P', rescaled, anything) with
E_lam = Phi(t) Q_lam(t): Phi would have to be E_() and Q_(1) would have to be the non-polynomial above.
Consequently there is no analogue of the unshifted reduction "census <-> finite d-vector", and the exact
down-census scan to n <= 90 was NOT run (it would not be an invariant of the census).

Identification of E_(): E_()(t) = sum_k (#shifted SYT of size k) t^k/k! with k! coefficients
1,1,1,2,3,6,12,27,63,154,398; log E_() = t + t^3/6 - t^4/12 + t^5/24 - 19t^6/720 + ... so E_() is not
exp(a t + b t^2) (first difference from exp(t + t^2/2) at t^2: 1 vs 2).

**What does hold** (verified for all |lam| <= 8, k < 10, `shifted_scan.py identity`):
with omega = all-ones functional, omega(U mu) = #addable = omega(D' mu)/2 + 1 - chi(mu)/2 where
chi(mu) = [mu_l = 1], and [D', U] = I gives [D', e^{tU}] = t e^{tU}; hence

    A_lam' = (1/2) A_{D' lam} + (1 + t/2) A_lam - (1/2) C_lam,      A_lam = E_lam,
    A_{D' lam} = sum_{mu covered by lam} w(mu) A_mu,   C_lam(t) = sum_k c_k(lam) t^k/k!,

c_k(lam) = number of chains of length k up from lam ending at a partition with last part 1; coefficientwise
2 s_{k+1}(lam) = sum_w w s_k(mu) + 2 s_k(lam) + k s_{k-1}(lam) - c_k(lam). Without the correction term
-c_k(lam) it fails for every one of the 25 shapes (first at k = 1 for lam = (), 2 s_2 = 2 vs 3). The extra
term C_lam is not expressible through the lam-data D'^j lam, which is why no e^{...} P_lam(ct) identity exists.
(Side observation: for one-row shapes (k) the ratio E_(k)/P_(k) agrees with E_() up to t^{k-1} and first
differs at t^k.)

**Fomin relations** (data/shifted_fomin.txt, all 43 strict partitions of n <= 10):
* D'U - UD' = I with U unweighted and D' weighted (off-diagonal 2, diagonal 1): HOLDS on every lam.
* DU - UD with unweighted D: FAILS; exactly (DU - UD) lam = lam if lam_l >= 2 (or lam = ()), and
  (DU - UD) lam = 0 if lam_l = 1. I.e. DU - UD = I - Pi, Pi = projection onto strict partitions with last part
  1 (all cross terms cancel; the diagonal coefficient is #addable - #removable = 1 - [lam_l = 1]).
  Examples: (DU-UD)(1) = 0, (DU-UD)(2) = (2), (DU-UD)(2,1) = 0, (DU-UD)(3,1) = 0, (DU-UD)(3,2) = (3,2).

## Task 3 — census injectivity (data/shifted_census_N40_K14.txt, data/shifted_census_N60_K21.txt)

`truncated_census(N, K)` computes (s_0, ..., s_K)(lam) for ALL strict partitions of n <= N by the level
recursion s_k(lam) = sum_{mu covers lam} s_{k-1}(mu), using all strict partitions up to N + K (54 for the
main run). Cross-checked against brute-force `s_up` for n <= 7, K = 6 (test).

**Result: no two distinct strict partitions of the same n <= 40 share (s_0..s_14)**; extension: no two
distinct strict partitions of the same n <= 60 share (s_0..s_21) (28.7 s). Since the truncated census
is a function of the full census, the full up-census determines the shape for every n <= 60 (rigorous).

Per n (n, q(n) = #strict partitions = #distinct truncated censuses, K_n = smallest K separating all strict partitions of n):

```
n : 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
q : 1 1 1 2 2 3 4 5 6 8 10 12 15 18 22 27 32 38 46 54 64 76 89 104 122 142 165 192 222 256 296
Kn: 0 0 0 1 5 4 2 4 4 4  4  4  5  6  6  6  5  6  5  6  7  6  7  8  7  8  9  8  9  10  9
n : 31 32 33 34 35 36 37 38 39 40 | 41 .. 60 (K=21 run)
q : 340 390 448 512 585 668 760 864 982 1113 | ... 10880 at n=60
Kn: 10 11 10 11 12 11 12 13 12 13 | 14 13 14 15 14 15 16 15 16 17 16 17 18 17 18 19 18 19 20 19
```
(q(n) and K_n for all n <= 60 are in the data files.) K_n grows like n/3: the hardest pairs are two-row
shapes, e.g. n=20: (14,6) vs (13,7) need K=7; n=30: (22,8) vs (19,11) need K=9; n=40: (28,12) vs
(26,14) need K=13; n=45: (32,13) vs (29,16) need K=14 — a chain of length k only "sees" a shape to depth
about k, and two-row shapes differ only when the row gap or the second row is exhausted. So a fixed K
cannot work for all n; K must be about n/3 (K=14 is just enough for n=40, where K_40 = 13).

Caveat on cross-size comparison: the truncated census with fixed K does NOT determine n — e.g. (15), (16),
..., (40) all have the same (s_0..s_14) (158 such cross-size groups for n <= 40, all of this truncation
type). This is irrelevant for the automorphism argument (depth is preserved) but means that "the census
determines the node" across different sizes has only been verified in the form: within each n, and the
full census distinguishes long rows once K exceeds their length.

Reading of the statement "each node of SSYT is uniquely determined by the census": the subtree below a
node depends only on its shape, so tableaux of the same shape have identical censuses; the meaningful
statement (and the one needed for automorphisms) is that the census determines the shape among nodes of
the same depth, equivalently among siblings. That is what was verified.

## Task 4 — sibling rigidity and triviality of Aut(SSYT) (data/shifted_siblings_N39_K14.txt)

For every strict lam with |lam| <= 39 and every pair of distinct covers mu != mu' of lam, the truncated
censuses (s_0..s_14)(mu), (s_0..s_14)(mu') differ: 43806 sibling pairs checked, 0 equal. The smallest K
separating all sibling pairs at level n (K_sib(n)) is <= 13 for n <= 39 (table in the data file, e.g.
K_sib(39) = 12, K_sib(37) = 13). This also follows from Task 3 (distinct shapes of size n+1 <= 40 have
distinct censuses), and Task 3's extension gives sibling rigidity for |lam| <= 59.

**Degrees.** Every strict partition has an addable box (row 1 is always addable), so a vertex T of shape lam
has degree #addable(lam) + 1 (its parent) if T is not the root, and the root has degree #addable(()) = 1.
Exactly one addable box occurs iff lam = () or lam is a staircase (l, l-1, ..., 1), l >= 1 (proved: #addable
= #removable + [lam_l >= 2], #removable >= 1 for lam != (), and #removable = 1 with lam_l = 1 forces
lam_i = lam_{i+1} + 1 for all i, i.e. a staircase; verified for n <= 9). Hence the staircase vertices
(including the vertex (1)) have degree 2, every other non-root vertex has degree >= 3, and the ROOT IS THE
UNIQUE VERTEX OF DEGREE 1. Any graph automorphism phi of SSYT therefore fixes the root, hence preserves
depth (= distance to the root) and maps the children of any vertex v onto the children of phi(v).

**Sibling rigidity => Aut = 1.** Induction on depth: suppose phi fixes every vertex of depth <= n (true
for n = 0). For a vertex v of depth n with shape lam, phi permutes the children of v, and phi restricted to
the subtree below a child c is an isomorphism of rooted trees onto the subtree below phi(c). The census
(s_k = number of vertices at depth k in the subtree) is an invariant of the rooted subtree, so
census(phi(c)) = census(c). The children of v have the distinct shapes lam + box, and by sibling rigidity
distinct children have distinct censuses, so phi(c) = c for every child c; thus phi fixes every vertex of
depth n+1. Hence sibling rigidity for all n implies Aut(SSYT) = {id}. Verified here for all vertices of
depth <= 39 (and, via the n <= 60 census run, depth <= 59): every automorphism of SSYT fixes every vertex
of depth <= 60 pointwise. For the full statement one needs sibling rigidity for all n; the growth K_n ~ n/3
means a proof cannot come from a bounded-depth invariant.

## Reproduction

```
cd /home/user/Helios-Landing-Page/syt_tree
python3 -m pytest tests/test_shifted.py -q          # 11 passed
python3 src/shifted_scan.py identity 8 10           # data/shifted_identity.txt
python3 src/shifted_scan.py fomin 10                # data/shifted_fomin.txt
python3 src/shifted_scan.py census 40 14            # data/shifted_census_N40_K14.txt  (0.9 s)
python3 src/shifted_scan.py census 60 21            # data/shifted_census_N60_K21.txt  (29 s)
python3 src/shifted_scan.py siblings 39 14          # data/shifted_siblings_N39_K14.txt (0.9 s)
```
