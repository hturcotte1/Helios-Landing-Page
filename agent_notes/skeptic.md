# Skeptic angle: trying to disprove Conjecture B

Author: skeptic agent. All scripts are in `src/agents/skeptic_*.py`, logs in `src/agents/skeptic_*.log`
(and `skeptic_near.json`, `skeptic_near12.json`). Exact integer arithmetic throughout; no code is shared with
`census.py`/`young.py`. Each claim is labelled PROVED / VERIFIED (with range and script) / CONJECTURAL.

Notation. `alpha|beta` = the partition whose rows are the rows of alpha followed by the rows of beta
(requires alpha_last >= beta_1). `(c^k)` = k rows of length c. `nu^t` = conjugate. `d(lam) = (d_0,...,d_n)`,
Hamming distance `ham(lam,mu) = #{j : d_j(lam) != d_j(mu)}`; "prefix agreement" a = first differing index,
"suffix agreement" b = n - (last differing index). C_k = sum of k-th powers of contents.

## 0. Summary

1. **Independent confirmation (VERIFIED).** A from-scratch implementation (`skeptic_dvec.py`, level-by-level recursion
   d_j(lam) = sum_c d_{j-1}(lam-c), two levels kept in memory) finds **zero collisions for all n <= 45**, checks
   d(lam) = d(lam^t) for every lam, and its per-n digests sum_{lam,j} d_j(lam) 7^j mod (2^61-1) **match
   `data/collisions_python.txt` for every n = 0..45** (`skeptic_digestdiff.py`, log `skeptic_digestdiff.log`, 46 lines all `MATCH`).
2. **Closest pairs (VERIFIED, n <= 40, Hamming radius 12).** For every 24 <= n <= 40 there are exactly 21 pairs of transpose classes
   at Hamming distance <= 12, and exactly 5 at distance <= 8; the minimum distance is 4 for every 8 <= n <= 40, attained only by
   (2,2,1^{n-4}) vs (3,1^{n-3}). **All 21 pairs are "attachment pairs" (N)|nu vs (N)|nu^t** (up to transposing both) with
   nu in {(2),(3),(3,1),(3,2),(3,3),(4),(4,1),(4,1,1),(4,2),(4,2,1),(4,2,2),(4,3),(4,3,1),(4,4),(5),(5,1),(5,1,1),(5,2),(6)}.
   Their agreement on d_0..d_{N-max(nu_1,l(nu))} is explained by the attachment lemma (Lemma A, PROVED below), which generalises
   the F_a family of results.md 4.10; the disagreement is always at the top end and involves f.
3. **Both-ends families (new).** There are infinite families whose prefix agreement grows linearly in n while the top four entries
   (f and C_2) agree, e.g.
   **H_k: lam = (5,3^k,2,2,2), mu = (4,4,3^k,1,1,1), n = 3k+11**: d_j(lam) = d_j(mu) for j <= k+3 and for j >= n-3,
   d_{n-4} differs (C_1 differs by exactly 2), Hamming distance 2k+4. PROVED: f equality (hook multisets coincide), C_2 equality,
   C_1(lam)-C_1(mu) = 2, and d_j agreement for j <= k conditional on one reduced identity R_K (VERIFIED K <= 18, equivalent to a
   3-row chain symmetry VERIFIED on a large range). The disagreement window [k+4, n-4] has length (2n-4)/3: **it is not bounded**.
   Two more such families H'_k, H''_k are recorded. No pair at n <= 40 agrees on the top 5 entries (f, C_2, C_1^2) and on d_0..d_3
   simultaneously.
4. **Targeted search (VERIFIED, n <= 63).** Pairs of transpose classes with the same multiset of corner boxes {a_i,b_i}, same f,
   same (C_2, C_1^2, C_4): **none for any n <= 63** (`skeptic_targeted.py`). Counts of the weaker filters are tabulated in Section 4;
   they grow, and a crude extrapolation puts the first survivor of the full prefilter near n ~ 90-110.
5. **Heuristic (Section 5).** An honest conditional-probability heuristic says: coincidences of any fixed finite set of top/bottom
   invariants become plentiful (their number grows without bound), but all coincidences ever observed are *structural*
   (explained by a factorisation lemma, a hook-multiset identity or a content-moment identity), and every structural mechanism found
   agrees with the d-vector only on a bottom window of length <= ~n/3 or a top window of length <= 8. A counterexample would need a
   structural identity that pins down all ~n-4 informative coordinates at once; the data give no hint of such a mechanism.
   I could not disprove Conjecture B, and I explain in Section 6 why I now think a counterexample is unlikely but not excluded.

## 1. Independent confirmation

`skeptic_dvec.py` (independent code: own partition generator, own corner detection, recursion d_j(lam) = sum_{c} d_{j-1}(lam - c)
computed for all partitions of n from those of n-1, keeping only two levels). Checks per n: (i) d(lam) = d(lam^t) for all lam
(would abort otherwise); (ii) grouping by exact d-vector, every group is {lam} (lam symmetric) or {lam, lam^t}.

* `python3 skeptic_dvec.py 30`: zero collision groups for all n <= 30 (1.3 s).
* `python3 skeptic_digestdiff.py 45`: zero collision groups for all n <= 45; the tuple (p(n), #distinct vectors, #symmetric,
  #transpose pairs, #collision groups, digest mod 2^61-1) **equals the corresponding line of `data/collisions_python.txt` for every
  n = 0..45**. Example n = 45: (89134, 44601, 68, 44533, 0, 1568118061668171670).

Status: **VERIFIED, n <= 45**, independently of the project's code (Theorem 4.1 of results.md covers n <= 75 with the project's code).

## 2. Near-collisions in Hamming distance, n <= 40

`skeptic_near12.py` (= `skeptic_near.py` with radius K = 12). For each n it takes one representative per transpose class and finds
all pairs at Hamming distance <= K by pigeonhole (coordinates split into K+1 interleaved residue classes mod K+1; a pair with
<= K disagreements agrees on a whole class; group by class, compare inside groups; every reported pair is re-checked coordinate by
coordinate with exact integers). It also computes, for each n, the pair maximising a+b (prefix+suffix agreement) by grouping on
prefixes and taking longest common suffixes among neighbours in reversed-sorted order (exact), and the counts of pairs agreeing on
the top s and bottom a entries.

Results (`skeptic_near12.log`, `skeptic_near12_classify.log`):

| n | #pairs ham<=8 | #pairs ham<=12 | min ham | attained by |
|---|---|---|---|---|
| 8..14 | 66,52,35,24,7,7,7 | 780 (n=12), 541, 344 | 4 | (2,2,1^{n-4}) vs (3,1^{n-3}) |
| 15..23 | 5,5,6,5,5,5,5,5,5 | 187,88,55,31,26,23,22,22,22 | 4 | same |
| 24..40 | 5 | 21 | 4 | same |

For n >= 24 the 21 pairs at distance <= 12 are (writing (N)|nu with N = n - |nu|; the pair is {(N)|nu, (N)|nu^t}):

| nu | Hamming | agree on |
|---|---|---|
| (2) | 4 | d_0..d_{n-4} |
| (3) | 6 | d_0..d_{n-6} |
| (3,1) | 7 | d_0..d_{n-7} |
| (4), (3,2) | 8 | d_0..d_{n-8} |
| (4,1), (3,3), (3,3,1)* | 9 | d_0..d_{n-9} (* one extra coincidence inside the window) |
| (5), (4,2), (4,1,1) | 10 | d_0..d_{n-10} |
| (5,1), (4,3), (4,2,1), (4,3,1,1)* | 11 | (* two extra coincidences) |
| (6), (5,2), (4,4), (5,1,1), (4,2,2), (4,3,1) | 12 | |

For 19 of the 21 the set of differing indices is exactly {n-|nu|-max(nu_1,l(nu))+1, ..., n}. These are exactly the pairs
predicted by the following lemma (with kappa = (N)); the F_a of results.md 4.10 are the case nu = (a-1).

**Lemma A (attachment; PROVED).** Let kappa, nu be partitions with kappa_last >= nu_1 and lam = kappa|nu. Then for
0 <= j <= kappa_last - nu_1,
    d_j(lam) = sum_{i=0}^{j} C(j,i) d_i(kappa) d_{j-i}(nu)     (with d_i(kappa) = 0 for i > |kappa|, etc.).
**Corollary A.** If also kappa_last >= l(nu), then d_j(kappa|nu) = d_j(kappa|nu^t) for all j <= kappa_last - max(nu_1, l(nu)).

*Proof.* Use d_j(lam) = sum over order filters S of lam with |S| = j of e(S) (BRIEFING, definition of d_j; a filter is a subset closed
under moving right and down inside lam, e(S) its number of linear extensions = number of orders in which S can be removed box by box).
Let S be a filter with |S| = j <= kappa_last - nu_1, S_1 = S ∩ kappa, S_2 = S ∩ nu (nu occupies rows l(kappa)+1, ...).
(a) S_1 contains no box in a column <= nu_1: if (i,c) in S_1 with c <= nu_1 then, moving right, S_1 contains (i,c'),...,(i,kappa_i),
which are kappa_i - c + 1 >= kappa_last - nu_1 + 1 > j boxes, contradiction.
(b) S_1 is a filter of kappa and S_2 a filter of nu (closure under moving right/down inside lam restricts to each part; moving down
from row l(kappa) inside a column > nu_1 leaves lam, so imposes nothing).
(c) No box of S_1 is comparable to a box of S_2 in the poset of lam: (i,c) in S_1 has c > nu_1 >= c' for every (i',c') in S_2,
while i < i'; comparability would need c <= c'.
(d) Conversely, for any filter S_1 of kappa with |S_1| = i <= j (which, by the argument in (a), lies in columns > nu_1) and any filter
S_2 of nu, S_1 ∪ S_2 is a filter of lam (moving down from S_1 leaves lam; moving right/down from S_2 stays in nu).
(e) For a disjoint union of two incomparable posets, e(S_1 ⊔ S_2) = C(|S|,|S_1|) e(S_1) e(S_2) (interleave two linear extensions).
Summing over (S_1,S_2) gives the formula. For the corollary apply the lemma to nu and nu^t and use d(nu) = d(nu^t) (BRIEFING fact 2/
transpose invariance). □

Tests (`skeptic_attach.py`, log `skeptic_attach.log`): Lemma A checked exactly for all kappa with |kappa| <= 10, all nu with
|nu| <= 6 and all admissible j: 2778 instances, 0 failures. In all 180 tested non-symmetric instances the first differing index was
exactly kappa_last - max(nu_1,l(nu)) + 1 (sharpness holds there; it is *not* claimed in general — the two starred pairs above
show extra coincidences for larger nu).

*f never agrees in attachment pairs (VERIFIED).* `skeptic_attach.py` also searched f^{kappa|nu} = f^{kappa|nu^t} (hook products)
for all non-symmetric nu with |nu| <= 12 and kappa = (N) (N <= 400), (N,N_2) (N <= 60, all N_2), (N,N,N) (N <= 400): **no solution**.
So Lemma A pairs always differ at d_n, d_{n-1}, d_{n-2}: they are bottom-only near-collisions.

## 3. Both-ends families

### 3.1 Data

`skeptic_ends.py` (log `skeptic_ends.log`) counts, for each n <= 40, the pairs agreeing on the top s entries and the bottom a
entries (grid s in {3,4,5,7,8}, a in {3,...,14}), and lists the pair with the largest prefix agreement among those agreeing on the
top s >= 4 (f and C_2) and s >= 5 (f, C_2, C_1^2).

* s >= 5 (top five = (n, f, C_2, C_1^2)): the largest prefix agreement is a <= 4 for all n <= 40 (a = 4 once, at n = 37; a <= 3
  otherwise; count of pairs with s >= 5 and a >= 3 is 0 or 1 for every n). **No pair agrees on d_0..d_3 and d_{n-4}..d_n.**
* s >= 7 or 8: no pair at all for n <= 40 (consistent with results.md 4.10: first at n = 49 resp. 84).
* s >= 4 (f and C_2): the maximal prefix agreement grows linearly: a = 9,10,11,12,13 at n = 26,29,32,35,38, always attained by
  (5,3^k,2,2,2) vs (4,4,3^k,1,1,1) (k = 5..9), and a = 7,8,10,11 at n = 28,31,37,40 by (6,4,3^k,2^4,1) vs (5,4,4,3^k,2,1^4),
  a = 9,10 at n = 30,34 by (6,4^k,3,3,2) vs (5,5,4^k,2,1,1).

### 3.2 The family H_k

**Theorem H (mixed status, itemised).** For k >= 1 let lam(k) = (5,3^k,2,2,2), mu(k) = (4,4,3^k,1,1,1), n = 3k+11.
Then mu(k) is neither lam(k) nor lam(k)^t, and
 (i)  f^{lam(k)} = f^{mu(k)}  — PROVED (hook multisets coincide);
 (ii) C_2(lam(k)) = C_2(mu(k)) and C_1(lam(k)) - C_1(mu(k)) = 2 — PROVED; hence d_j(lam) = d_j(mu) for j = n-3,...,n
      (results.md Thm 4.6: d_{n-3} is a function of (n,f,C_2)) and d_{n-4}(lam) != d_{n-4}(mu) (Cor. 4.10.2: d_{n-4} determines
      C_1^2 given n,f,C_2, and C_1(lam)^2 != C_1(mu)^2 since C_1(mu) = C_1(lam) - 2 and C_1(lam) = -(3k^2+9k-2)/2 <= -5, see the proof of (ii));
 (iii) d_j(lam(k)) = d_j(mu(k)) for all j <= k — PROVED conditional on identity R_K below (Lemma B is proved; R_K is VERIFIED for K <= 18);
 (iv) VERIFIED for k <= 14 (`skeptic_family.py`, log `skeptic_family.log`): the set of differing indices is exactly
      {k+4, ..., n-4}, so ham = 2k+4 = (2n-14)/3, prefix agreement k+4 = (n+1)/3, suffix agreement 4.

*Proof of (i).* Hook lengths (arm + leg + 1), rows indexed from 1.
lam(k): row 1 = (5): boxes (1,1..5) have arms 4,3,2,1,0 and legs k+3,k+3,k,0,0: hooks k+8, k+7, k+3, 2, 1.
Rows i = 2..k+1 (length 3): arms 2,1,0; legs (k+4-i), (k+4-i), (k+1-i): hooks k+7-i, k+6-i, k+2-i, i.e. as i runs the three columns
give the sets {6..k+5}, {5..k+4}, {1..k}. Rows k+2,k+3,k+4 (length 2): hooks (4,3), (3,2), (2,1).
mu(k): row 1 = (4): arms 3,2,1,0; legs k+4, k+1, k+1, 1: hooks k+8, k+4, k+3, 2. Row 2 = (4): arms 3,2,1,0; legs k+3, k, k, 0:
hooks k+7, k+3, k+2, 1. Rows i = 3..k+2 (length 3): arms 2,1,0; legs k+5-i, k+2-i, k+2-i: hooks k+8-i, k+4-i, k+3-i, giving
{6..k+5}, {2..k+1}, {1..k}. Rows k+3,k+4,k+5 (length 1): hooks 3,2,1.
Multisets: lam: {k+8,k+7,k+3,2,1} ∪ {6..k+5} ∪ {5..k+4} ∪ {1..k} ∪ {4,3,3,2,2,1};
mu: {k+8,k+4,k+3,2,k+7,k+3,k+2,1} ∪ {6..k+5} ∪ {2..k+1} ∪ {1..k} ∪ {3,2,1}.
Remove the common parts {k+8,k+7,k+3,2,1}, {6..k+5}, {1..k}, {3,2,1}: lam leaves {5..k+4} ∪ {4,3,2} = {2..k+4};
mu leaves {k+4,k+3,k+2} ∪ {2..k+1} = {2..k+4}. Equal multisets, so equal hook products and equal f = n!/prod(hooks). □
(Checked: `hooks_equal=True` for k = 1..14.)

*Proof of (ii).* Contents c - i. lam(k): row 1 contributes sum 10, sum of squares 30; row i in 2..k+1 contributes 6-3i and
3i^2-12i+14; row i in k+2..k+4 contributes 3-2i and 2i^2-6i+5. mu(k): rows 1,2 contribute 6+2 = 8 and 14+6 = 20; row i in
3..k+2 contributes 6-3i and 3i^2-12i+14; row i in k+3..k+5 contributes 1-i and (i-1)^2.
C_1 difference: (10-8) + [sum_{i=2}^{k+1} - sum_{i=3}^{k+2}](6-3i) + sum_{i=k+2}^{k+4}(3-2i) - sum_{i=k+3}^{k+5}(1-i)
 = 2 + [(6-6) - (6-3k-6)] + (-6k-9) - (-3k-9) = 2 + 3k - 6k - 9 + 3k + 9 = 2.
C_2 difference: (30-20) + [q(2) - q(k+2)] + sum_{i=k+2}^{k+4}(2i^2-6i+5) - sum_{i=k+3}^{k+5}(i-1)^2 with q(i) = 3i^2-12i+14:
q(2) - q(k+2) = 2 - (3k^2+2) = -3k^2; the two sums are 6k^2+18k+19 and 3k^2+18k+29; total 10 - 3k^2 + 6k^2 + 18k + 19 - 3k^2 - 18k - 29 = 0. □
(Checked numerically for k <= 14: `C2eq=True`, `C1 = -5,-7`, `-14,-16`, ....) Closed form: C_1(lam(k)) = 10 + sum_{i=2}^{k+1}(6-3i) + sum_{i=k+2}^{k+4}(3-2i) = 10 + (3k-3k^2)/2 - 6k - 9 = -(3k^2+9k-2)/2 (= -5, -14, -26 for k = 1,2,3, matching the log), so C_1(lam(k)) <= -5 and C_1(mu(k)) = C_1(lam(k)) - 2; C_1(lam)^2 = C_1(mu)^2 would need C_1(lam) = 1. Hence C_1^2 differs and d_{n-4} differs.

**Lemma B (body insertion; PROVED).** Let alpha, beta be partitions (possibly empty), c >= 1, k >= 1, with alpha_last >= c >= beta_1,
and lam = alpha|(c^k)|beta. Let alpha' = alpha/(c^{l(alpha)}) be the skew diagram of the boxes of alpha in columns > c, and let
D_i(alpha') = number of ways to remove i boxes one at a time from alpha' (= sum over filters S of alpha' with |S| = i of e(S)).
Then for 0 <= j <= k,
    d_j(lam) = sum_{i=0}^{j} C(j,i) D_i(alpha') d_{j-i}((c^k)|beta).
*Proof.* Let S be a filter of lam with |S| = j <= k. If S contained a box (i,c') of alpha with c' <= c, then moving down it would
contain (i',c') for all i <= i' <= l(alpha)+k (all these rows have length >= c >= c'), which are >= k+1 > j boxes. So S ⊆ alpha' ∪ P
where P = (c^k)|beta is the set of rows below alpha. S_1 = S ∩ alpha' is a filter of alpha' (moving right/down from a box of alpha' inside
lam stays in alpha': rows <= l(alpha), columns > c), S_2 = S ∩ P is a filter of the Young diagram P. No box (i,c') of S_1 is comparable
with a box (i',c'') of S_2: c' > c >= c'' and i < i'. Conversely for any filters S_1 of alpha', S_2 of P the union is a filter of lam
(moving down from row l(alpha) in a column > c leaves lam). Hence, as in Lemma A(e), d_j(lam) = sum_i C(j,i) D_i(alpha') D_{j-i}(P)
and D_{j-i}(P) = d_{j-i}((c^k)|beta). □
Tests (`skeptic_lemmaB.py`, log `skeptic_lemmaB.log`): all alpha with |alpha| <= 8, beta with |beta| <= 5, c in 1..4, k in 1..5,
all j <= k: 19360 instances, 0 failures.

*Proof of (iii) given R_k.* lam(k) = (5)|(3^k)|(2,2,2): alpha' = two boxes in one row, D(alpha') = (1,1,1).
mu(k) = (4,4)|(3^k)|(1,1,1): alpha' = two boxes in one column, D(alpha') = (1,1,1). By Lemma B, for j <= k,
d_j(lam(k)) - d_j(mu(k)) = sum_{i=0}^{2} C(j,i) [ d_{j-i}((3^k,2,2,2)) - d_{j-i}((3^k,1,1,1)) ], which vanishes if
    (R_K)   d_j((3^K,2,2,2)) = d_j((3^K,1,1,1)) for all j <= K.   □

**Identity R_K (VERIFIED K <= 18, `skeptic_lemmaB_reduced.log`).** In fact the first differing index of the two vectors is exactly
K+2 for every K = 1..18 (so R_K holds with one index to spare). Reformulation (PROVED equivalence): rotating the skew shapes
lam/nu by 180 degrees inside the box (K+3) x 3 — which preserves f^{lam/nu} — and transposing, one gets for j <= K
    d_j((3^K,2,2,2)) = N_j((3)),   d_j((3^K,1,1,1)) = N_j((3,3)),
where N_j(nu) = number of chains nu = kappa^0 ⊂ kappa^1 ⊂ ... ⊂ kappa^j in Young's lattice with all l(kappa^i) <= 3
(= sum_{l(kappa)<=3, |kappa/nu|=j} f^{kappa/nu}). [Derivation: for lam = (K+3,K+3,K) = (3^K,2,2,2)^t and |lam/nu| = j <= K the
complements in the 3 x (K+3) box are lam^c = (3) and nu^c = kappa with kappa ⊇ (3), |kappa| = j+3, l(kappa) <= 3, and every such
kappa occurs; similarly (K+3,K,K)^c = (3,3).] N_j(nu) depends only on (u,v) = (nu_1-nu_2+1, nu_2-nu_3+1), and R_K for all K is
equivalent to N_j(4,1) = N_j(1,4) for all j. `skeptic_walks.py` (log `skeptic_walks.log`) verifies the general symmetry
**N_j(u,v) = N_j(v,u) for all 1 <= u,v <= 8 and j <= 24** (0 failures), and re-derives d_j((3^K,2,2,2)) = N_j(1,4) = N_j(4,1) = d_j((3^K,1,1,1))
for j <= K, K <= 12. I did not find a proof of the symmetry; in representation-theoretic terms it says that V_nu ⊗ V^{⊗j} and
V_{nu*} ⊗ V^{⊗j} (V the defining representation of SL_3, nu* the dual) have the same number of irreducible constituents. It is
labelled CONJECTURAL as a general statement, VERIFIED in the stated range; the family statement (iii) needs only N_j(4,1) = N_j(1,4).

**Companion families (VERIFIED).** `skeptic_family.py` / `skeptic_lemmaB_reduced.log`:
* H'_k: (6,4,3^k,2^4,1) vs (5,4,4,3^k,2,1^4), n = 3k+19, k = 1..14: hooks equal, C_2 equal, C_1 differs by 4, differing indices
  exactly {k+4,...,n-4}. Lemma B reduces the bottom agreement (j <= k) to d_j((3^K,2^4,1)) = d_j((3^K,2,1^4)) for j <= K, which
  holds with first difference at K+3 for K <= 18.
* H''_k: (6,4^k,3,3,2) vs (5,5,4^k,2,1,1), n = 4k+14, k = 2..10 (k = 1 is a transpose pair): f equal, differing indices exactly
  {k+5,...,n-4}; reduced identity d_j((4^K,3,3,2)) = d_j((4^K,2,1,1)) for j <= K+1, K <= 18.
Both alpha' are again a 2-box row vs a 2-box column. The general recipe is: alpha|(c^k)|beta vs alpha*|(c^k)|beta* with
D(alpha') = D(alpha*') and (c^K)|beta, (c^K)|beta* agreeing at the bottom; f, C_2 agreement is an extra Diophantine coincidence.

**What these families do and do not give.** They are the first pairs with *linear* prefix agreement (about n/3) *and* agreement of
the top four entries. But: (a) the disagreement window [k+4, n-4] has length (2n-4)/3 and is not bounded; (b) they disagree at
d_{n-4} because C_1 differs (C_1(lam) - C_1(mu) = 2, so C_1^2 differs unless C_1(lam) = 1, which never happens);
(c) `skeptic_ends.log` shows that no pair at n <= 40 combines top-5 agreement with bottom-4 agreement. To turn such a family into a
candidate counterexample one would need (1) C_1(lam) = -C_1(mu) (mirror-type), (2) all C_{2i} equal and all products of odd moments
equal — by the content-multiset argument (the multiset of contents determines lam: the number of boxes of content c is the length of
the c-th diagonal, and the diagonal lengths determine lam) this forces the content multiset of mu to be that of lam or its negative,
i.e. mu in {lam, lam^t} — **unless** the top levels stop determining new content invariants, which they never do in the verified range
(results.md Thm 4.10.4 lists a new invariant at each of levels 7..10). So a counterexample cannot come from "content coincidence";
it must come from two shapes whose top-level polynomial invariants (one per level) coincide accidentally at every level. This is the
same conclusion as results.md 4.7, reached from the other side.

## 4. Targeted search: same corner multiset, same f, same (C_2, C_1^2, C_4), n <= 63

`skeptic_targeted.py` (logs `skeptic_targeted.log` for n = 10..50, `skeptic_targeted_51_64.log` for n = 51..63; n = 64 was still
running when this report was written). One representative per transpose class. Keys (all transpose-invariant): K0 = (n, multiset of
unordered corner-box pairs {a_i,b_i}); K0+C = K0 + (C_2, C_1^2, C_4); K1 = K0 + f (f computed by hook lengths only inside K0 groups of
size >= 2); K2 = K1 + (C_2, C_1^2, C_4). By the local product theorem (results.md 4.2) K0 fixes d_0..d_{2m}, and by 4.4/4.10 the
content data fix d_{n-6}..d_n, so K2 survivors are the pairs agreeing on both ends *as far as the known structure theorems reach*.

| n | classes | pairs same K0 | same K0+C | same K0+f | same K2 |
|---|---|---|---|---|---|
| 10 | 22 | 1 | 0 | 0 | 0 |
| 20 | 317 | 153 | 0 | 1 | 0 |
| 30 | 2811 | 3799 | 0 | 4 | 0 |
| 40 | 18692 | 58544 | 1 | 21 | 0 |
| 50 | 102162 | 710578 | 13 | 111 | 0 |
| 55 | 225710 | 2295094 | 45 | 287 | 0 |
| 60 | 483338 | 7174908 | 123 | 363 | 0 |
| 63 | 752877 | 13858631 | 202 | 618 | 0 |

**Zero K2 survivors for every 10 <= n <= 63** (VERIFIED). For the K1 survivors at n <= 50 the script computed exact d-vectors
(own memoised recursion): every K1 pair (same corners, same f, different contents) is separated already at d_{n-3} (C_2) and at a
small bottom index d_3..d_6 (prefix agreement = 2m+1 or 2m+2 typically), e.g. n = 42: (13,10,4,2^4,1^8) vs (12,9,3,3,2^4,1^9),
ham = 36, differing indices 6..41.

*Crude extrapolation (heuristic, not a claim).* If K0+f and K0+C coincidences were independent inside K0 groups, the expected number
of K2 pairs would be (K0+f)(K0+C)/K0 = 618*202/13858631 ~ 0.009 at n = 63. Fitting the growth on n = 50..63 (K0+f ~ 5x, K0+C ~ 15x,
K0 ~ 20x over that range, p(n) ~ 4.7x) suggests the expected number reaches 1 around n ~ 90-110. Such a survivor would still have to
agree on the remaining ~n-10 middle coordinates.

## 5. A probabilistic heuristic for the number of collisions at size n

Informative coordinates. d_0 = 1, d_1 = r, d_n = d_{n-1} = d_{n-2} = f, and d_{n-6} is determined by (n, d_n, d_{n-3}, d_{n-4}, d_{n-5})
(results.md Thm 4.10.3). So at most n-4 coordinates (d_2..d_{n-3} minus d_{n-6}) carry information beyond (n, r, f).
More structure: (results.md 4.2) d_2..d_{2m+1} are functions of the corner data; (4.10) d_{n-3}, d_{n-4}, d_{n-5} are functions of
(n, f, C_2, C_1^2, C_4), and level n-i for i >= 7 adds one polynomial content invariant each.

(a) *Naive birthday.* Treating the ~n-4 coordinates as independent uniform draws makes the expected number of collisions
C(p(n)/2, 2) / prod_j V_j, with V_j the number of achievable values of d_j; since d_j is of order j! up to polynomial factors,
prod_j V_j is super-exponential in n^2 while p(n)^2 = exp(O(sqrt n)), giving an expectation ~ 0 at every n. This is worthless:
the coordinates are extremely correlated (d_j and d_{j+1} differ by a factor ~ r, the bottom ones are functions of few parameters,
the top ones of few invariants), and the near-collision census of Sections 2-3 shows agreement propagating over windows of
length ~n/3.

(b) *Empirical conditional probabilities at n = 40* (`skeptic_near12.json`, 18692 classes, 1.747e8 pairs).
Top: P(same f) = 3404/1.747e8 = 1.9e-5; P(same d_{n-3} | same f) = 49/3404 = 0.014; P(same d_{n-4} | ...) = 4/49 = 0.08;
P(same d_{n-5} | ...) = 0/4. Bottom: P(same d_1) = 0.29, then conditional probabilities of agreeing on the next coordinate given
agreement on the previous ones are 0.30, 0.21, 0.19, 0.23, 0.35, 0.45, 0.52, 0.57, 0.60, 0.64, 0.67 (a = 3..13). The bottom
conditional probabilities *increase* with the depth of agreement: the survivors of deep bottom agreement are structured
(Lemma A/B type pairs) and structure persists. The top conditional probabilities decrease (each level is a new polynomial invariant
and the group of shapes sharing the previous invariants is small).

(c) *Expected collisions, structured heuristic.* Write E_n = sum over pairs of P(agree everywhere). Condition on the top:
T_5(n) = #pairs agreeing on (n,f,C_2,C_1^2,C_4) is 0 for n <= 48 and 76 at n = 84 (results.md 4.10, agent scan); T_8(84) = 1.
Each further top level i contributes a new polynomial in the contents of weighted degree i (Thm 4.10.4); for two shapes with the
same lower invariants the new invariant is an integer of size ~n^{i}, and — if it behaved like a random integer in a range of size
n^{c} relative to the fixed data — the conditional agreement probability per level would be ~n^{-c}. Over the ~n/2 top levels that
gives P ~ n^{-c n/2} = exp(-(c/2) n log n), which beats the p(n)^2 = exp(O(sqrt n)) pairs by an enormous margin, so E_n -> 0 summably.
The bottom and middle levels can only help. **The honest content of this heuristic is small**: its only real input is the
assumption that successive levels are "generic" relative to each other, and Sections 2-3 exhibit families where ~n/3 successive
levels agree for a structural reason. What the heuristic does say correctly is that a counterexample cannot be an accident: with
p(n)^2 pairs and n-4 coordinates each carrying ~n log n bits, accidental agreement of even ~n/10 unstructured middle coordinates
has probability exp(-Omega(n^2 log n)) per pair. A counterexample must be produced by an identity that controls *all* coordinates
simultaneously (as transposition does).

(d) *What such an identity would have to do.* The three mechanisms found (Lemma A/B factorisation at the bottom; hook-multiset
equality for f; content-moment equalities for the top 5-8 levels) are independent of each other and each controls a window of
length <= ~n/3 (bottom) or <= 8 (top). A collision needs, in addition, equality of the middle coordinates d_j for n/3 < j < n-8,
i.e. of sum_{|nu| = n-j} f^{lam/nu} where nu ranges over shapes of about half the size — quantities with no known local or
polynomial description. I found no mechanism affecting them.

## 6. Assessment

* Independent confirmation: **VERIFIED, no collisions for n <= 45**, digests identical to the project's data for all n <= 45.
* Closest pairs: minimum Hamming distance is 4 for all 8 <= n <= 40 (only (2,2,1^{n-4}) vs (3,1^{n-3})); all pairs within distance 12
  for 24 <= n <= 40 are the 21 attachment pairs (N)|nu vs (N)|nu^t, |nu| <= 6, explained by Lemma A (PROVED); none has equal f.
* Best both-ends pairs: the families H_k, H'_k, H''_k (Theorem H) agree on d_0..d_{~n/3} and d_{n-3}..d_n; f, C_2 equal (PROVED),
  C_1 differs by a constant (PROVED), bottom agreement PROVED modulo the verified identity R_K (a 3-row chain symmetry). No pair at
  n <= 40 agrees on the top 5 and the bottom 4 simultaneously; no pair at n <= 63 shares corner multiset, f and (C_2, C_1^2, C_4).
* Disproof: **not achieved**; no candidate counterexample found. Why a counterexample is unlikely: every observed coincidence is
  explained by a local mechanism whose reach is bounded by a fraction of the vector, and accidental agreement of the rest is
  astronomically improbable; a counterexample requires a currently unknown global identity. Why it is not excluded: the number of
  pairs agreeing on any fixed finite set of top and bottom invariants grows without bound (Section 4 table, results.md 4.10 scan to
  n = 84), so a proof cannot come from finitely many invariants, and the "middle" coordinates are not understood at all.

Open items suggested by this angle: (1) prove the chain symmetry N_j(u,v) = N_j(v,u) (would make Theorem H(iii) unconditional);
(2) find the general form of "body insertion" pairs with equal f and C_2 and decide whether C_1^2 can also be matched (I expect not:
in all three families the C_1 difference is a nonzero constant determined by alpha, beta); (3) push the K2 prefilter to n ~ 100 in C++.

## 7. Scripts (all in `src/agents/`)

* `skeptic_dvec.py` — independent d-vector code; `python3 skeptic_dvec.py 30`.
* `skeptic_digestdiff.py` — comparison with `data/collisions_python.txt`, n <= 45 (`skeptic_digestdiff.log`).
* `skeptic_near.py`, `skeptic_near12.py` — Hamming near-collision census n <= 40, K = 8 / 12 (`skeptic_near.log`, `skeptic_near12.log`,
  `skeptic_near.json`, `skeptic_near12.json`, `skeptic_near12_classify.log`).
* `skeptic_ends.py` — both-ends grid n <= 40 (`skeptic_ends.log`).
* `skeptic_attach.py` — Lemma A tests and f-equality search for attachment pairs (`skeptic_attach.log`).
* `skeptic_family.py` — families H, H' exact for k <= 14 (`skeptic_family.log`); H'' and reduced pairs in `skeptic_lemmaB_reduced.log`.
* `skeptic_lemmaB.py` — Lemma B tests, 19360 instances (`skeptic_lemmaB.log`).
* `skeptic_walks.py` — 3-row chain symmetry N_j(u,v) = N_j(v,u), u,v <= 8, j <= 24 (`skeptic_walks.log`).
* `skeptic_targeted.py` — corner-multiset/f/content prefilter search, n <= 63 (`skeptic_targeted.log`, `skeptic_targeted_51_64.log`).
