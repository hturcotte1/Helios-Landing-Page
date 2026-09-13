# Angle classC — Conjecture B for the class C of (3,2,1)-avoiding partitions

Author: classC agent. Scripts: `src/agents/classC_verify.py` (all lemmas + final brute force),
`src/agents/classC_extra.py` (independent recheck of the case analysis by direct enumeration).
All arithmetic exact (Python integers).

## 0. Statement and status

**Theorem (PROVED; VERIFIED for all n ≤ 30 by `check_theorem` in `classC_verify.py`).**
Let C be the set of partitions containing no (3,2,1), i.e. λ with λ_3 = 0 (two-row shapes, including one row),
or λ_1 ≤ 2 (two-column shapes), or λ_2 ≤ 1 (hooks). If λ ∈ C, μ ⊢ |λ| and d(μ) = d(λ)
(equal down-census vectors (d_0,…,d_n)), then μ = λ or μ = λ^t.

Everything below is elementary: the only ingredients are the definition of d_j as the number of
length-j removal sequences, the corner-run bookkeeping of Lemma 1, the reflection principle
(Lemma 4) and a few binomial sums. No character theory is used (the C_2 invariant of the briefing
would give an alternative for the two-row-vs-two-row case; it is not needed).

Conventions. λ = (λ_1 ≥ λ_2 ≥ …). A *removal sequence of length j* from λ is a sequence
λ = λ^0 ⊃ λ^1 ⊃ … ⊃ λ^j of partitions with |λ^{i-1}/λ^i| = 1; d_j(λ) is their number,
d_j(λ) = Σ_c d_{j-1}(λ − c) (sum over removable corners c), d_0 = 1, d_j = 0 for j > n.
Transposition is an automorphism of Young's lattice, so d(λ^t) = d(λ).
Corner runs: distinct part sizes α_1 > … > α_r > 0, b_i = multiplicity of α_i,
a_i = α_i − α_{i+1} (α_{r+1} = 0); r = number of removable corners; the removable corner of run i
is the last box of the last row of the block of rows of length α_i.
A_1 = #{i : a_i = 1}, B_1 = #{i : b_i = 1}. Binomial convention: C(j,k) = 0 unless 0 ≤ k ≤ j.

## 1. Lemma 1 (corner bookkeeping; d_1, d_2)  — PROVED; VERIFIED all partitions n ≤ 22

**Lemma 1.** (a) d_1(λ) = r. (b) If c_i is the removable corner of run i, then
r(λ − c_i) = r + 1 − [a_i = 1] − [b_i = 1]. (c) d_2(λ) = r(r+1) − A_1 − B_1 (valid for every
partition, including n ≤ 1 where both sides are 0).

*Proof.* (a) is the definition. (b) Removing c_i changes exactly one row, from length α_i to α_i − 1.
The multiset of positive part sizes changes as follows: the multiplicity of α_i drops from b_i to
b_i − 1, and (if α_i − 1 > 0) the multiplicity of α_i − 1 rises by one. The number of distinct
positive part sizes is therefore r − [b_i = 1] + [α_i − 1 is positive and was not already a part size].
If i < r then α_i − 1 ≥ α_{i+1} ≥ 1 is positive and it is a new size iff α_i − 1 > α_{i+1}, i.e.
iff a_i ≥ 2. If i = r then α_r − 1 is a new positive size iff α_r ≥ 2, i.e. iff a_r = α_r ≥ 2.
In both cases the last indicator is 1 − [a_i = 1]. Since r counts distinct positive part sizes, (b) follows.
(c) d_2(λ) = Σ_i d_1(λ − c_i) = Σ_i r(λ − c_i) = Σ_i (r + 1 − [a_i=1] − [b_i=1]) = r(r+1) − A_1 − B_1. ∎

Verified: `check_lemma1(22)` tests (a),(b),(c) for every partition of every n ≤ 22.

## 2. Lemma 2 (the six 2-corner families)  — PROVED; VERIFIED n ≤ 30

**Lemma 2.** Let μ have r = 2 with corner runs [(a_1,b_1),(a_2,b_2)]. Then A_1 + B_1 ≥ 2 iff at
least two of the four numbers a_1,b_1,a_2,b_2 equal 1, and the six possible pairs give
| pair equal to 1 | shape |
|---|---|
| b_1 = b_2 = 1 | two-row (a_1+a_2, a_2) |
| a_1 = a_2 = 1 | two-column (2^{b_1}, 1^{b_2}) |
| b_1 = a_2 = 1 | hook (a_1+1, 1^{b_2}) |
| a_1 = b_2 = 1 | F2 = ((a_2+1)^{b_1}, a_2) |
| a_1 = b_1 = 1 | F1 = (a_2+1, a_2^{b_2}) |
| a_2 = b_2 = 1 | F4 = ((a_1+1)^{b_1}, 1) = F1^t |
Moreover, C ∩ {r = 2} = two-row (l_1 > l_2 ≥ 1) ∪ two-column ∪ hooks (a ≥ 2, b ≥ 1); and
F1 = (c+1, c^d) lies in C iff c = 1 (hook) or d = 1 (two-row); F2 = ((c+1)^b, c) lies in C iff b = 1
(two-row) or c = 1 (two-column); F4 = ((a+1)^b, 1) lies in C iff a = 1 or b = 1.

*Proof.* From the definition of corner runs, α_2 = a_2, α_1 = a_1 + a_2 and μ = (α_1^{b_1}, α_2^{b_2});
substituting the two unit entries gives the table. Membership statements: (c+1,c^d) has three rows
iff d ≥ 2, first part ≥ 3 iff c ≥ 2, second part ≥ 2 iff c ≥ 2, so it contains (3,2,1) iff c,d ≥ 2;
similarly for the others. A partition is in C iff it does not contain (3,2,1) iff (λ_3 = 0 or λ_1 ≤ 2
or λ_2 ≤ 1). ∎  Verified: `check_lemma2(30)`.

## 3. Lemma 3 (d_3 on the six families)  — PROVED; VERIFIED all members with n ≤ 40

Write [P] for an indicator. Using Lemma 1(c) on each λ − c_i (whose corner runs are listed):

**Lemma 3.**
(i) Two-row λ = (l_1,l_2), δ = l_1 − l_2 ≥ 1, l_2 ≥ 1:
d_3 = 8 − [δ=2] − [l_2=2] if δ ≥ 2, l_2 ≥ 2;  = 6 − [l_2=2] if δ = 1, l_2 ≥ 2;
= 4 − [δ=2] if δ ≥ 2, l_2 = 1;  = 2 if λ = (2,1).  In all cases d_3 ≤ 8, with d_3 = 8 iff δ ≥ 3 and l_2 ≥ 3.
(ii) Hook λ = (a,1^b), a ≥ 2, b ≥ 1:
d_3 = 8 − [a=3] − [b=2] if a ≥ 3, b ≥ 2;  = 4 − [b=2] if a = 2, b ≥ 2;  = 4 − [a=3] if a ≥ 3, b = 1;
= 2 if λ = (2,1). In all cases d_3 ≤ 8, with d_3 = 8 iff a ≥ 4 and b ≥ 3.
(iii) F1 = (c+1, c^d) with c,d ≥ 2: d_3 = 10 − [c=2] − [d=2].
(iv) F2 = ((c+1)^b, c) with b,c ≥ 2: d_3 = 10 − [b=2] − [c=2].
(v) F4 = ((a+1)^b, 1) with a,b ≥ 2: d_3 = 10 − [a=2] − [b=2] (transpose of (iii)).

*Proof.* d_3(λ) = Σ_i d_2(λ − c_i), and each d_2 is r(r+1) − A_1 − B_1 by Lemma 1(c).
(i) λ − c_1 = (l_1 − 1, l_2). If δ ≥ 2 its runs are [(δ−1,1),(l_2,1)]: r = 2, A_1 = [δ=2] + [l_2=1],
B_1 = 2, so d_2 = 4 − [δ=2] − [l_2=1]. If δ = 1 it is (l_2,l_2) with runs [(l_2,2)]: d_2 = 2 − [l_2=1].
λ − c_2 = (l_1, l_2 − 1). If l_2 ≥ 2 its runs are [(δ+1,1),(l_2−1,1)]: d_2 = 6 − [l_2=2] − 2 = 4 − [l_2=2].
If l_2 = 1 it is (l_1), runs [(l_1,1)], d_2 = 2 − [l_1=1] − 1 = 1 (l_1 ≥ 2). Adding gives the four cases.
(ii) λ − c_1 = (a−1,1^b): for a ≥ 3 runs [(a−2,1),(1,b)], A_1 = [a=3]+1, B_1 = 1+[b=1], d_2 = 4 − [a=3] − [b=1];
for a = 2 it is (1^{b+1}), runs [(1,b+1)], d_2 = 2 − 1 − 0 = 1. λ − c_2 = (a,1^{b−1}): for b ≥ 2 runs
[(a−1,1),(1,b−1)], d_2 = 6 − (1+[a=2]) − (1+[b=2]) = 4 − [a=2] − [b=2]; for b = 1 it is (a), d_2 = 1. Add.
(iii) Corners of (c+1,c^d) (c,d ≥ 2): row 1 and row d+1. λ − c_1 = (c^{d+1}), runs [(c,d+1)], r = 1,
d_2 = 2 − [c=1] − [d+1=1] = 2. λ − c_2 = (c+1, c^{d−1}, c−1), runs [(1,1),(1,d−1),(c−1,1)] (three runs
since d−1 ≥ 1 and c−1 ≥ 1), r = 3, A_1 = 2 + [c=2], B_1 = 2 + [d=2], d_2 = 12 − 4 − [c=2] − [d=2].
Sum: 10 − [c=2] − [d=2].
(iv) Corners of ((c+1)^b, c) (b,c ≥ 2): row b and row b+1. λ − c_1 = ((c+1)^{b−1}, c, c), runs
[(1,b−1),(c,2)], A_1 = 1, B_1 = [b=2], d_2 = 5 − [b=2]. λ − c_2 = ((c+1)^b, c−1), runs [(2,b),(c−1,1)],
A_1 = [c=2], B_1 = 1, d_2 = 5 − [c=2]. Sum: 10 − [b=2] − [c=2].
(v) d(F4) = d(F1^t) = d(F1). ∎

**Corollary 3b (VERIFIED n ≤ 30, `check_lemma3b`).** Every r = 2 member of C has d_3 ≤ 8. Every member
of F1 ∪ F2 ∪ F4 not in C has d_3 ≥ 9, except exactly (3,2,2), (3,3,1), (3,3,2), which have d_3 = 8.
The r = 2 members of C with d_3 = 8 and n ∈ {7,8} are exactly (4,1^3) (n=7) and (5,1^3), (4,1^4) (n=8)
(two-row shapes need δ, l_2 ≥ 3, hence n ≥ 9; two-column shapes are their transposes; hooks need a ≥ 4, b ≥ 3).
Their SYT counts, f = C(n−1,b) (Lemma 4(ii)), are 20, 35, 35, whereas f^{(3,2,2)} = f^{(3,3,1)} = 21 and
f^{(3,3,2)} = 42 (hook length formula: hooks of (3,2,2) are 5,4,1,3,2,2,1, product 240, 7!/240 = 21;
hooks of (3,3,2) are 5,4,2,4,3,1,2,1, product 960, 8!/960 = 42).

## 4. Lemma 4 (explicit d_j for two-row shapes and hooks)  — PROVED; VERIFIED n ≤ 40

**Lemma 4.** (i) Let λ = (l_1,l_2) with l_1 ≥ l_2 ≥ 0, δ = l_1 − l_2. For every j ≥ 0,
  d_j(λ) = Σ_{k = max(0, ⌈(j−δ)/2⌉)}^{min(j, l_2)} [ C(j,k) − C(j, j−δ−1−k) ].
Equivalently (the form in the task statement) d_j = Σ_{k=0}^{min(j,l_2)} C(j,k) sgn(δ+1−j+2k) [j−k ≤ l_1+1].
(ii) Let λ = (p+1, 1^q), p,q ≥ 0, n = p+q+1. For 0 ≤ j ≤ n−1,
  d_j(λ) = Σ_{i = max(0, j−q)}^{min(j, p)} C(j,i),  and d_n(λ) = f^λ = C(n−1, q) = C(n−1, p).

*Proof of (i).* A removal sequence of length j is encoded by the word w ∈ {1,2}^j recording the row of
the box removed at each step; after the prefix u the shape is (l_1 − #_1(u), l_2 − #_2(u)), and the
sequence is valid iff for every prefix u this is a partition: l_1 − #_1(u) ≥ l_2 − #_2(u) ≥ 0. Since #_2
only grows, the condition "l_2 − #_2(u) ≥ 0 for all prefixes" is equivalent to k := #_2(w) ≤ l_2. The
condition l_1 − #_1(u) ≥ l_2 − #_2(u) reads #_1(u) − #_2(u) ≤ δ. The quantity h(u) = #_1(u) − #_2(u)
starts at 0 and changes by ±1 at each letter, so it exceeds δ for some prefix iff it equals δ+1 for some
prefix. Fix k with 0 ≤ k ≤ min(j, l_2); the words with k letters 2 number C(j,k). Call w *bad* if h(u) = δ+1
for some prefix u. A valid word must also satisfy the end condition h(w) = j − 2k ≤ δ, i.e. 2k ≥ j − δ;
if 2k < j − δ there are no valid words with k twos. Assume 2k ≥ j − δ. Reflection: for a bad w let u be the
shortest prefix with h(u) = δ+1 and write w = uv; let w' = u v̄ where v̄ is v with 1 ↔ 2 swapped. Then
h(w') = (δ+1) − (h(w) − h(u)) = 2δ + 2 − j + 2k, so w' has k' := j − δ − 1 − k letters 2, and u is still the
shortest prefix of w' with h = δ+1, so the map is an involution and hence a bijection from the bad words with
k twos onto the bad words with k' twos. Every word with k' twos is bad: its final value is
h = j − 2k' = 2δ + 2 − j + 2k ≥ δ + 2 > δ+1 (using 2k ≥ j − δ), and h moves by ±1 from 0, so it passes
through δ+1. Hence the number of bad words with k twos equals the number of all words with k' twos, which is
C(j, k') = C(j, j−δ−1−k) (zero if k' < 0, i.e. if #_1 = j−k ≤ δ, in which case indeed no word is bad).
Therefore the number of valid words with k twos is C(j,k) − C(j, j−δ−1−k), and summing over
max(0,⌈(j−δ)/2⌉) ≤ k ≤ min(j,l_2) gives the first formula.
For the equivalence with the sgn form: the positive terms are the C(j,k) with 2k ≥ j−δ (⇔ δ+1−j+2k ≥ 1),
0 ≤ k ≤ min(j,l_2); for these j − k ≤ (j+δ)/2 ≤ (n+δ)/2 = l_1, so the bracket [j−k ≤ l_1+1] is automatic.
The negative terms −C(j,k') with k' = j−δ−1−k: as k runs over the range above, k' runs over the integers with
2k' ≤ j − δ − 2 (⇔ δ+1−j+2k' ≤ −1 ⇔ sgn = −1), k' ≤ j−δ−1 (automatic from the previous when j ≥ δ, and when
j < δ there are no such k' ≥ 0 at all), and k' ≥ j − δ − 1 − l_2 = j − l_1 − 1 (⇔ [j−k' ≤ l_1+1]); moreover
k' ≤ (j−δ−2)/2 ≤ (n−δ−2)/2 = l_2 − 1 < l_2 and k' < j automatically. Terms with sgn = 0 contribute nothing.
So the two expressions agree term by term.

*Proof of (ii).* The box (1,1) of a hook is a removable corner only when the shape is (1). During the first
n−1 removals the shape has ≥ 2 boxes, so (1,1) is never removed for j ≤ n−1. The remaining boxes form two
chains: the arm (1,2),…,(1,p+1) and the leg (2,1),…,(q+1,1). At any moment the last remaining arm box is
removable (row 1 is then longer than row 2, which has length ≤ 1) and the last remaining leg box is removable
(it is the last row, of length 1). So a removal sequence of length j ≤ n−1 is exactly a word in {arm, leg}^j
with i ≤ p arm-letters and j − i ≤ q leg-letters, and there are C(j,i) such words for each i. For j = n: the
first n−1 steps remove all of the arm and the leg (C(n−1,p) = C(n−1,q) words) and the last step removes (1,1). ∎

Verified: `check_lemma4(40)` tests (i) in both forms for every l_1 ≥ l_2 ≥ 0 with l_1 + l_2 ≤ 40 and every
0 ≤ j ≤ n, and (ii) for every p,q ≥ 0 with p+q+1 ≤ 40 and every j.

## 5. Lemma 5 (the invariant M and the next two entries)  — PROVED; VERIFIED n ≤ 40

For a d-vector define M(λ) := max{ j ≥ 0 : d_i(λ) = 2^i for all 0 ≤ i ≤ j } (an invariant of the d-vector).

**Lemma 5A (two-row).** Let λ = (l_1,l_2), l_1 > l_2 ≥ 1, δ = l_1 − l_2, M := min(δ, l_2) ≥ 1, n = 2l_2 + δ ≥ 3M.
(a) M(λ) = M and M + 2 ≤ n.
(b) d_{M+1} = 2^{M+1} − 2 if δ = l_2 (= M), and 2^{M+1} − 1 otherwise.
(c) d_{M+2} = 2^{M+2} − M − 5 if δ = l_2 = M;
    d_{M+2} = 2^{M+2} − 2 − [l_2 = M+1] if δ = M < l_2;
    d_{M+2} = 2^{M+2} − M − 3 − [δ = M+1] if l_2 = M < δ.
(d) d_n = f^λ = C(n,l_2) − C(n,l_2−1) = C(n−1,l_2) − C(n−1,l_2−2).

*Proof.* Use Lemma 4(i); write k_0 = max(0, ⌈(j−δ)/2⌉) and N_k = C(j, j−δ−1−k) for the subtracted term.
(a) For j ≤ M: j − δ ≤ 0 so k_0 = 0; j ≤ l_2 so the upper limit is j; and j−δ−1−k ≤ −1 so every N_k = 0.
Hence d_j = Σ_{k=0}^{j} C(j,k) = 2^j. For j = M+1 (≤ n since n ≥ 3M ≥ M+1): if M = δ < l_2 then k_0 = 1,
upper limit j, N_k = C(j,−k) = 0, so d_j = 2^j − 1; if M = l_2 < δ then k_0 = 0, upper limit l_2 = j−1,
N_k = 0 (j−δ−1−k ≤ −1), so d_j = 2^j − 1; if δ = l_2 = M then k_0 = 1, upper limit j − 1, N_k = 0, so
d_j = 2^j − 2. In all cases d_{M+1} < 2^{M+1}, proving (a) and (b). n ≥ 3M ≥ M+2 since M ≥ 1.
(c) j = M+2. Case δ = l_2 = M: k_0 = ⌈2/2⌉ = 1, upper limit M = j−2; N_k = C(j, 1−k) is nonzero only for
k = 1, where it equals 1. So d_j = Σ_{k=1}^{M} C(M+2,k) − 1 = (2^{M+2} − 1 − (M+2) − 1) − 1 = 2^{M+2} − M − 5.
Case δ = M < l_2: k_0 = 1, upper limit min(M+2, l_2); N_k nonzero only for k = 1 (value 1). So
d_j = Σ_{k=1}^{min(M+2,l_2)} C(M+2,k) − 1, which is 2^{M+2} − 1 − 1 if l_2 ≥ M+2 and one less (the term
C(M+2,M+2) = 1 is missing) if l_2 = M+1.
Case l_2 = M < δ: upper limit M; k_0 = 0 if δ ≥ M+2 and k_0 = 1 if δ = M+1; N_k = C(j, M+1−δ−k) with
M+1−δ−k ≤ −k, nonzero only if δ = M+1 and k = 0, which is excluded by k_0 = 1 in that case. So
d_j = Σ_{k=k_0}^{M} C(M+2,k) = 2^{M+2} − (M+2) − 1 − [δ = M+1].
(d) j = n = 2l_2 + δ: k_0 = ⌈(n−δ)/2⌉ = l_2 = upper limit, and n−δ−1−l_2 = l_2 − 1, so
d_n = C(n,l_2) − C(n,l_2−1); Pascal's rule gives the second form (C(n−1,−1) = 0 when l_2 = 1). ∎

**Lemma 5B (hooks).** Let μ = (p+1, 1^q) with p ≥ q ≥ 1, n = p+q+1, M := q.
(a) M(μ) = q.  (b) d_{q+1} = 2^{q+1} − 1 − [p = q].
(c) If p ≥ 2 (so q+2 ≤ n−1): d_{q+2} = 2^{q+2} − q − 3 if p ≥ q+2; = 2^{q+2} − q − 4 if p = q+1;
    = 2^{q+2} − 2q − 6 if p = q.  (d) d_n = C(n−1, q).

*Proof.* Lemma 4(ii), valid for j ≤ n−1 = p+q. (a) For j ≤ q ≤ p the sum runs over all 0 ≤ i ≤ j: 2^j.
For j = q+1 ≤ p+q the sum runs over 1 ≤ i ≤ min(q+1,p), which misses i = 0, so d_{q+1} < 2^{q+1}.
(b) The sum is Σ_{i=1}^{q+1} = 2^{q+1} − 1 if p ≥ q+1, and Σ_{i=1}^{q} = 2^{q+1} − 2 if p = q.
(c) j = q+2 ≤ p+q iff p ≥ 2. The sum is Σ_{i=2}^{min(q+2,p)} C(q+2,i) = 2^{q+2} − 1 − (q+2) minus the terms
with i > p: none if p ≥ q+2; C(q+2,q+2) = 1 if p = q+1; C(q+2,q+1) + C(q+2,q+2) = q+3 if p = q.
(d) is Lemma 4(ii). ∎

Verified: `check_lemma5(40)` tests 5A(a)–(d) for all l_1 > l_2 ≥ 1 with n ≤ 40 and 5B(a)–(d) for all p ≥ q ≥ 1
with n ≤ 40.

## 6. Lemma 6 (rectangles, r = 1)  — PROVED; VERIFIED ab ≤ 36

(Fact 4 of the briefing; a self-contained proof is included for completeness.)

**Lemma 6.** For the rectangle a^b (a,b ≥ 1) and every j, d_j(a^b) = Σ_{ρ ⊢ j, ρ_1 ≤ a, ℓ(ρ) ≤ b} f^ρ.
Let I_j := Σ_{ρ ⊢ j} f^ρ (a universal sequence). Then d_j(a^b) = I_j for j ≤ m := min(a,b), and
d_{m+1}(a^b) < I_{m+1} whenever ab ≥ m+1 (i.e. unless a = b = 1). Consequently, if a^b and a'^{b'} have the
same d-vector then {a,b} = {a',b'}, i.e. the two rectangles are equal or transposes.

*Proof.* A removal sequence of length j from a^b removes the boxes of a skew shape (a^b)/ν, |ν| = n − j, in
an order that is a linear extension of the reversed containment order (a box is removed only after the boxes
to its right and below it). Rotating the rectangle by 180° maps (a^b)/ν onto a straight shape ρ ⊢ j with
ρ_1 ≤ a and ℓ(ρ) ≤ b, and reverses the order, so removal sequences with underlying set (a^b)/ν correspond
bijectively to standard Young tableaux of shape ρ; conversely each such ρ arises from exactly one ν. This
gives the formula. If j ≤ min(a,b) every ρ ⊢ j satisfies ρ_1 ≤ j ≤ a and ℓ(ρ) ≤ j ≤ b, so d_j = I_j.
For j = m+1 ≤ ab: if m = a then the row ρ = (m+1) violates ρ_1 ≤ a; if m = b then the column (1^{m+1})
violates ℓ(ρ) ≤ b; in either case a partition with f^ρ = 1 is omitted, so d_{m+1} ≤ I_{m+1} − 1.
Hence m = max{ j : d_i = I_i for all i ≤ j } is read off the d-vector (when a = b = 1 the vector is (1,1)
and n = 1 forces the shape). With n = ab also known, {a,b} = {m, n/m}. ∎  Verified: `check_lemma6(36)`.

## 7. Proof of the Theorem

Let λ ∈ C, μ ⊢ n, d(μ) = d(λ). Since d(ν^t) = d(ν) for all ν and (λ^t)^t = λ, the conclusion
"μ ∈ {λ, λ^t}" is invariant under replacing λ by λ^t and/or μ by μ^t; we do so freely.

**Step 0: r.** d_1 = r (Lemma 1(a)), so r(μ) = r(λ). Every subpartition of a member of C is in C
(containment is transitive), and (3,2,1) is the unique minimal partition with three removable corners
(a partition with three distinct part sizes α_1 > α_2 > α_3 ≥ 1 contains (α_1,α_2,α_3) ⊇ (3,2,1)),
so r(λ) ≤ 2. If n = 0 there is nothing to prove.

**Step 1: r = 1.** Then λ and μ are rectangles (r = 1 means one distinct part size) and Lemma 6 gives
μ ∈ {λ, λ^t}.

**Step 2: r = 2, reduction to μ ∈ C.** By Lemma 2, λ is a two-row shape (l_1 > l_2 ≥ 1), a two-column
shape, or a hook (a ≥ 2, b ≥ 1); for these Lemma 1(c) gives d_2(λ) = 6 − A_1 − B_1 with
A_1 + B_1 ≥ 2 (two-row: b_1 = b_2 = 1; two-column: a_1 = a_2 = 1; hook: b_1 = a_2 = 1), so d_2(λ) ≤ 4.
Then d_2(μ) ≤ 4 forces A_1(μ) + B_1(μ) ≥ 2, so by Lemma 2, μ belongs to one of the six families.
If μ ∉ C then μ ∈ F1 ∪ F2 ∪ F4 with both parameters ≥ 2 (Lemma 2), and by Lemma 3(iii)–(v)
d_3(μ) ≥ 9 unless μ ∈ {(3,2,2), (3,3,1), (3,3,2)}. But d_3(λ) ≤ 8 by Lemma 3(i),(ii) (two-column shapes
are transposes of two-row shapes). So μ is one of the three exceptions, d_3(λ) = 8, and n ∈ {7,8}; by
Corollary 3b, λ or λ^t is (4,1^3) (n = 7) or one of (5,1^3), (4,1^4) (n = 8) with d_n(λ) = f^λ ∈ {20, 35},
while d_n(μ) = f^μ ∈ {21, 42} — contradiction. Hence μ ∈ C, and r(μ) = 2.

**Step 3: both in C with r = 2.** Transposing as needed, we may assume that each of λ, μ is either in
R := {(l_1,l_2) : l_1 > l_2 ≥ 1} or in H := {(a,1^b) : a ≥ 2, b ≥ 1} (two-column shapes are transposes of
members of R). Note R ∩ H = {(l_1,1) : l_1 ≥ 2}.

*Case RR (λ, μ ∈ R).* By Lemma 5A(a) both have the same M = min(δ,l_2) = min(δ',l_2') ≥ 1. If
d_{M+1} = 2^{M+1} − 2 then by 5A(b) δ = l_2 = M and δ' = l_2' = M, so λ = μ = (2M, M). Otherwise
δ ≠ l_2 and δ' ≠ l_2', and by 5A(c) d_{M+2} ∈ {2^{M+2} − 2, 2^{M+2} − 3} if the minimum M is attained by δ,
while d_{M+2} ∈ {2^{M+2} − M − 3, 2^{M+2} − M − 4} if it is attained by l_2. These two sets are disjoint
because M ≥ 1 ({2,3} ∩ {M+3, M+4} = ∅). Hence λ and μ are of the same kind: either δ = δ' = M, and then
l_2 = l_2' = (n − M)/2; or l_2 = l_2' = M, and then δ = δ' = n − 2M. In both cases λ = μ.

*Case HH (λ, μ ∈ H).* Write λ = (p+1,1^q), μ = (p'+1,1^{q'}); transposing (which maps (p+1,1^q) to
(q+1,1^p)) we may assume p ≥ q and p' ≥ q'. By Lemma 5B(a), q = M = q', and p = n − 1 − q = p'. So λ = μ.

*Case RH (λ ∈ R, μ ∈ H, the case λ ∈ H, μ ∈ R being symmetric).* Write λ = (l_1,l_2), μ = (p+1,1^q).
If l_2 = 1 then λ ∈ H and Case HH applies. If q = 1 then μ ∈ R and Case RR applies. If p = 1 then
μ = (2,1^q) = (q+1,1)^t with (q+1,1) ∈ R, and Case RR applied to λ and μ^t gives μ^t = λ, i.e. μ = λ^t.
So assume l_2 ≥ 2, p ≥ 2, q ≥ 2. If δ = 1 then A_1(λ) = [δ=1] + [l_2=1] = 1, B_1(λ) = 2, so
d_2(λ) = 3, whereas d_2(μ) = 6 − (1 + [p=1]) − (1 + [q=1]) = 4: contradiction. So δ ≥ 2 as well.
Transposing μ we may assume p ≥ q ≥ 2. Now Lemma 5A(a) and 5B(a) give M := min(δ,l_2) = q ≥ 2, and
M + 2 ≤ min(n, n−1), so all entries used below are covered by Lemmas 5A/5B (5B(c) needs p ≥ 2, which holds).
 - If d_{M+1} = 2^{M+1} − 2 then δ = l_2 = M (5A(b)) and p = q (5B(b)); then 5A(c), 5B(c) give
   d_{M+2}(λ) = 2^{M+2} − M − 5 and d_{M+2}(μ) = 2^{M+2} − 2M − 6, which differ since M ≠ −1. Contradiction.
 - Otherwise p > q = M and δ ≠ l_2. If δ = M < l_2 then d_{M+2}(λ) ∈ {2^{M+2} − 2, 2^{M+2} − 3} by 5A(c),
   while d_{M+2}(μ) ∈ {2^{M+2} − M − 3, 2^{M+2} − M − 4} by 5B(c); disjoint since M ≥ 2 (indeed M ≥ 1
   suffices). Contradiction.
 - Hence l_2 = M = q < δ. Then by 5A(d) and 5B(d), d_n(λ) = C(n−1,q) − C(n−1,q−2) and d_n(μ) = C(n−1,q).
   Since q ≥ 2 and 0 ≤ q−2 ≤ n−1, C(n−1,q−2) ≥ 1, so d_n(λ) < d_n(μ). Contradiction.
So Case RH with l_2, p, q ≥ 2 cannot occur, and in the remaining sub-cases μ ∈ {λ, λ^t} was shown.

All cases give μ ∈ {λ, λ^t}. ∎

## 8. Verification summary (all in `src/agents/classC_verify.py`; every assertion passed)

| item | range |
|---|---|
| Lemma 1 (a),(b),(c) | all partitions, n ≤ 22 |
| Lemma 2 (six families, shapes, C-membership) | all r=2 partitions, n ≤ 30 |
| Lemma 3 (i)–(v) d_3 formulas | all parameters ≤ 26 |
| Corollary 3b (d_3 bound / exceptions / f-values) | all r=2 partitions n ≤ 30 |
| Lemma 4 (i) both forms, (ii) | all two-row shapes n ≤ 40, all j; all hooks n ≤ 40, all j |
| Lemma 5A, 5B | all two-row shapes and hooks n ≤ 40 |
| Lemma 6 | all rectangles ab ≤ 36 |
| Theorem, brute force (every λ ∈ C vs every μ ⊢ n) | n ≤ 30 |

Independent recheck `src/agents/classC_extra.py`: recomputes d-vectors of all two-row shapes and hooks
by direct enumeration of removal sequences (no recursion) for n ≤ 12, and re-runs the Case RR / HH / RH
decision procedure of Step 3 as an explicit algorithm on the d-vectors, recovering the shape (up to
transpose) for all C-shapes with n ≤ 60.

## 9. Remarks / what remains open

Nothing in the stated goal remains open. Remarks:
* The proof uses only d_1, d_2, d_3, d_{M+1}, d_{M+2} and d_n, i.e. the bottom of the vector plus f^λ.
  The pair (3,3,2) vs (5,1^3)/(4,1^4) shows that d_0..d_3 alone do not suffice (they agree: 1,2,4,8 for all).
* The reflection formula of Lemma 4(i) is the ℓ = 2 case of the briefing's Pfaffian formula
  ([x^p y^q](x−y)/((1−xy)(1−x)(1−y)) = sgn(p−q)); here it is proved directly.
* Beyond C: the natural next class is r = 2 in general (families F1, F2, F4 and all other 2-corner
  shapes); Lemma 1(c) and the d_3 computation of Lemma 3 extend verbatim to any 2-corner shape, but the
  small-j analysis of Lemma 5 would need the general 2-run local structure (briefing fact 3).
