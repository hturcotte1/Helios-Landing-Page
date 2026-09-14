# Angle `topend` — what the TOP of the d-vector knows: characters, Jucys–Murphy elements, content moments

Author: topend agent (second pass; the earlier `topend_lib/verify/fit/G` scripts in `src/agents/` are from a
first pass and are not used below except as an independent cross-check of the Murnaghan–Nakayama code).
All scripts of this pass are `src/agents/topend_*.py` with the names listed in §9; every log quoted is next to
its script.  All arithmetic is exact (Python integers / `Fraction`, sympy `Rational`).

Conventions.  λ ⊢ n, boxes (i,j) 0-indexed, content c(i,j) = j − i,  C_k(λ) = Σ_{boxes} c^k  (C_0 = n),
u_i(λ) = d_{n−i}(λ) = Σ_{ν ⊢ i} f^{λ/ν}  (u_0 = u_1 = u_2 = f := f^λ),  r = number of removable corners.
For a partition ρ with no part 1 and |ρ| ≤ n, K_ρ = K_ρ(n) ∈ Q[S_n] is the sum of all permutations of cycle
type ρ ∪ 1^{n−|ρ|} (K_∅ = 1; K_ρ(n) := 0 if n < |ρ|); |K_ρ| = (n)_{|ρ|}/z_ρ.  For a central z ∈ Q[S_n],
ω_λ(z) is the scalar by which z acts on the irreducible module V^λ (Schur); ω_λ(K_ρ) = |K_ρ| χ^λ(ρ)/f^λ.
m(ρ) := |ρ| − ℓ(ρ) ("reduced size").  J_k = Σ_{i<k} (i k) (k = 1..n, J_1 = 0) are the Jucys–Murphy elements,
p_m(J) = Σ_k J_k^m, p_μ(J) = Π_t p_{μ_t}(J).  (n)_i = n(n−1)…(n−i+1).  I_k = number of involutions of [k].

---------------------------------------------------------------------------------------------------------
## 0. Results at a glance

| # | statement | status |
|---|-----------|--------|
| A | Explicit class expansions (I1)–(I6) in Q[S_n] (p_1..p_4 of the JM elements, K_2², K_3²), hence exact formulas ω_λ(K_3), ω_λ(K_{22}), ω_λ(K_5), ω_λ(K_{33}) in terms of n, C_1², C_2, C_4 | PROVED (§2–3); VERIFIED in Z[S_n] for n ≤ 8 and via MN for all λ ⊢ n ≤ 16 |
| B | u_i = Σ_{ρ⊢i} (σ(ρ)/z_ρ) χ^λ(ρ∪1^{n−i}), explicit coefficient tables for i ≤ 10; χ^λ(3-cycle) = 3u_3 − 2f, χ^λ((2,2)) = 4u_4 − 4u_3 + f, χ^λ(5-cycle) = 5u_5 − 5u_4 + f, χ^λ((3,3)) = (9/2)u_6 − (9/2)u_5 + (3/2)u_3 − f/2; hence the d-vector determines n, f, C_2, C_1², C_4 | PROVED (§4); VERIFIED n ≤ 20 using only the d-vector recursion and contents |
| P | Polynomiality with degree bound: for every ρ without 1's there is F_ρ ∈ Q[n, p_1, …, p_m], m = m(ρ), of weighted degree ≤ m (wt n = 1, wt p_k = k) with K_ρ(n) = F_ρ(n; p(J)) in Q[S_n] for all n ≥ 0, so ω_λ(K_ρ) = F_ρ(n; C_1(λ), …, C_m(λ)) | PROVED (§5) |
| Ω | The 16 polynomials Ω_ρ for all even-type ρ without 1's, |ρ| ≤ 10 (table in §5.4) | PROVED for all n (Theorem P + full-rank certificate on n ≤ 16 + exact solve; verified n ≤ 18) |
| C | Level-6 identity: (n)_6 u_6/f is an explicit polynomial in n, C_1², C_2, C_4, i.e. d_{n−6} is a function of (n, d_n, d_{n−3}, d_{n−4}, d_{n−5}) | PROVED (§6); VERIFIED n ≤ 20 |
| L | Levels 7–10 each carry exactly one new content polynomial: h_7 = C_6 − 16C_1C_3; at level 8 a combination of C_3² and (n)C_1C_3; at level 9 C_8 − 24C_1C_5 + (…)C_1C_3; at level 10 C_1(C_2C_3, C_3, C_5)-terms (explicit in §6.3). The first function-level witnesses: level 7 at n = 18, level 8 at n = 40, levels 9, 10 have no witness for n ≤ 40 | PROVED (polynomial identities) / VERIFIED (witnesses) |
| S | (n, f, C_2, C_1², C_4) = data of (n, d_n, d_{n−3}, d_{n−4}, d_{n−5}) separates all partitions up to transpose for n ≤ 48, fails at n = 49 (1 pair), 50 (3 pairs), …; adding u_7 separates all n ≤ 58; depth 3 fails from n = 14, depth 4 from n = 28 | VERIFIED |
| S′ | Within two-row shapes, within hooks, within rectangles, (n, C_2) alone determines the shape up to transpose | PROVED (§7.2) |
| M | "Mirror box moves": pairs λ ↦ μ = λ − box(content −c) + box(content c) with equal f and C_1(λ) = −c have all even C_{2k}, C_1², f, n equal; they exist (n = 50, 51, 53, 54, 57, …), so no finite set of the invariants (f, n, C_1², C_2, C_4, C_6, …) can prove Conjecture B; they are separated by C_1C_3 (level 7) | PROVED (mechanism) / VERIFIED (existence, n ≤ 62) |
| H | Lemma U (upper-interval isomorphism ⇒ equal d_j for j ≤ n − s) and the families F_a: (a, 1^{n−a}) vs (2^{a−1}, 1^{n−2a+2}) agree on d_0..d_{n−2a+2} for every n ≥ 2a − 1, a ≥ 3 | PROVED (§8) ; VERIFIED a ≤ 7, n ≤ 44 (and sharp there) |
| H′ | For 8 ≤ n ≤ 40 the ONLY pairs of transpose classes agreeing on d_0..d_{n−4} are F_3 = {(3,1^{n−3}), (2,2,1^{n−4})}; the only ones agreeing on d_0..d_{n−5} are again F_3; on d_0..d_{n−6}: F_3 and F_4 (n ≥ 11) | VERIFIED |

---------------------------------------------------------------------------------------------------------
## 1. Classical facts used (stated precisely)

(F1) **Schur.** A central element z of Q[S_n] acts on the irreducible module V^λ as the scalar ω_λ(z) = tr(z|V^λ)/f^λ; ω_λ is a Q-algebra homomorphism Z(Q[S_n]) → Q. For z = K_ρ this gives ω_λ(K_ρ) = |K_ρ| χ^λ(ρ)/f^λ.  An element Σ_g a_g g is central iff a_g depends only on the cycle type of g; the class sums form a basis of the center.

(F2) **Jucys–Murphy theorem** (Jucys 1974; Murphy 1981).  V^λ has a basis {v_T : T a standard Young tableau of shape λ} (Young's seminormal basis) with J_k v_T = c_T(k) v_T for every k, where c_T(k) is the content (column − row) of the box of T containing k.

(F3) **Murnaghan–Nakayama rule** (used only inside the computer verifications, implemented independently in `topend_mn.py` by beta-numbers; cross-checked against two other implementations and against column orthogonality for n ≤ 10, log `topend_algebra.log` (A)).

(F4) **Branching rule.** Res^{S_n}_{S_i} V^λ ≅ ⊕_{ν ⊢ i} (V^ν)^{⊕ f^{λ/ν}}; hence χ^λ(g) = Σ_{ν⊢i} f^{λ/ν} χ^ν(g) for g ∈ S_i ⊂ S_n.

(F5) **Frobenius–Schur for S_n.** Every irreducible representation of S_n is realisable over Q, so all Frobenius–Schur indicators are 1 and #{h ∈ S_i : h² = g} = Σ_{ν ⊢ i} χ^ν(g).  (Also checked for all classes of S_i, i ≤ 10, by MN: `topend_levels.log` item (6).)

(F6) Hook length formula f^λ = n!/Π hooks.

Nothing else from the literature is used in a proof.  (The Kerov–Olshanski theory of shifted symmetric functions is mentioned once in a remark, never in a proof.)

---------------------------------------------------------------------------------------------------------
## 2. Group-algebra identities   (PROVED; VERIFIED in Z[S_n], 2 ≤ n ≤ 8: `topend_algebra.py` check_B)

Products are composed right-to-left: (pq)(x) = p(q(x)).  (All identities below are invariant under the
anti-automorphism g ↦ g^{-1}, which fixes every J_k and every K_ρ, so the convention is immaterial.)

**Lemma 2.1 (star products).** For distinct i_1, …, i_μ all different from k:
  (i_1 k)(i_2 k) ⋯ (i_μ k) = (i_μ i_{μ−1} … i_1 k)   (a (μ+1)-cycle).
Consequently, for a fixed (μ+1)-set S with maximum k, the map (i_1,…,i_μ) ↦ (i_1 k)⋯(i_μ k) is a bijection from
orderings of S∖{k} onto the μ! cycles of length μ+1 on S.

*Proof.* Induction on μ; μ = 1 is trivial.  Let γ = (i_{μ−1} … i_1 k) = (i_1 k)⋯(i_{μ−1} k).  In γ·(i_μ k) the
transposition acts first: i_μ ↦ k ↦ γ(k) = i_{μ−1};  k ↦ i_μ ↦ i_μ (γ fixes i_μ);  i_t ↦ i_{t−1} for 2 ≤ t ≤ μ−1
and i_1 ↦ k as in γ.  This is the cycle (i_μ i_{μ−1} … i_1 k).  Every cycle on S can be written uniquely as
(x_1 … x_μ k) with k last, which gives the bijection. ∎

**Lemma 2.2 (commutation and centrality).** (a) J_k J_l = J_l J_k.  (b) For s_k = (k k+1):
s_k J_k s_k = J_{k+1} − s_k, s_k J_{k+1} s_k = J_k + s_k, and s_k J_i s_k = J_i for i ∉ {k, k+1}.
(c) For every symmetric polynomial F ∈ Q[x_1,…,x_n]^{S_n}, F(J_1,…,J_n) is central in Q[S_n].

*Proof.* (a) For k < l, J_k ∈ Q[S_{l−1}] and J_l = Σ_{i<l}(i l) is fixed by conjugation by S_{l−1}.
(b) Conjugation by s_k maps (i k) ↦ (i k+1) and (i k+1) ↦ (i k) for i < k, and (k k+1) ↦ itself; so
s_k J_k s_k = Σ_{i<k} (i k+1) = J_{k+1} − (k k+1), and s_k J_{k+1} s_k = Σ_{i<k}(i k) + (k k+1) = J_k + s_k.
For i < k, J_i ∈ Q[S_{i}] commutes with s_k (disjoint supports); for i > k+1, s_k ∈ S_{i−1} fixes J_i.
(c) Put e_1 = J_k + J_{k+1}, e_2 = J_k J_{k+1}.  By (b), s_k e_1 s_k = e_1, and
s_k e_2 s_k = (J_{k+1} − s_k)(J_k + s_k) = J_{k+1}J_k + (J_{k+1} s_k − s_k J_k) − 1.  From
s_k J_k s_k = J_{k+1} − s_k we get J_{k+1} s_k = s_k J_k + 1, so the middle bracket is 1 and
s_k e_2 s_k = J_{k+1}J_k = e_2 by (a).  Since the J's commute, F(J) is a polynomial in the J_i (i ≠ k, k+1) and in
e_1, e_2 (F is symmetric in x_k, x_{k+1}); all of these commute with s_k, hence so does F(J).  The s_k generate S_n. ∎

**Proposition 2.3 (six class expansions).** In Q[S_n], for every n ≥ 1 (classes that do not fit in S_n are 0):

 (I1) p_1(J) = K_2.
 (I2) p_2(J) = C(n,2)·1 + K_3.
 (I3) K_2² = C(n,2)·1 + 3K_3 + 2K_{22}.       (hence e_2(J) = (p_1² − p_2)/2 = K_3 + K_{22})
 (I4) p_3(J) = K_4 + (2n − 3)K_2.
 (I5) p_4(J) = K_5 + (3n − 4)K_3 + 4K_{22} + n(n−1)(4n−5)/6 · 1.
 (I6) K_3² = 2C(n,3)·1 + (3n − 8)K_3 + 8K_{22} + 5K_5 + 2K_{33}.

*Proof.* (I1) is the definition.

(I2) J_k² = Σ_{i,j<k} (ik)(jk) = (k−1)·1 + Σ_{i≠j<k} (j i k) by Lemma 2.1.  As (i,j) runs over ordered pairs of
distinct elements of [k−1], (j i k) runs exactly once over the 3-cycles on {i,j,k}, i.e. once over all 3-cycles whose
largest point is k.  Summing over k, every 3-cycle appears once; Σ_k (k−1) = C(n,2).

(I3) K_2² = Σ over ordered pairs (t,t′) of transpositions.  t = t′: C(n,2) copies of 1.  t,t′ share one point:
(ab)(ac) is the 3-cycle (a c b); a given 3-cycle γ on {x,y,z} arises exactly from the 3 choices of a ∈ {x,y,z} (then
c = γ(a), b = γ^{-1}(a) are forced).  t, t′ disjoint: (ab)(cd) = (cd)(ab), each (2,2)-element arises twice.

(I4) J_k³ = Σ_{i,j,l<k} (ik)(jk)(lk).  All distinct: 4-cycles (l j i k), each 4-cycle with maximum k exactly once
(Lemma 2.1).  Coincidences: i = j ≠ l gives (lk); i ≠ j = l gives (ik); i = l ≠ j gives (ik)(jk)(ik) = (ij)
(conjugation); i = j = l gives (ik).  Fix a transposition (x k) with x < k: it arises from the triples (i,i,x) with
i ≠ x (k−2 of them), from (x,j,j) with j ≠ x (k−2), and from (x,x,x) (1): 2k − 3 triples.  A transposition (xy)
with x,y < k arises from (x,y,x) and (y,x,y): 2 triples.  Summing over k: a transposition (a b), a < b, receives
2b − 3 from k = b and 2 from each k = b+1, …, n: total 2n − 3.

(I5) J_k⁴ = Σ_{i,j,l,m<k} (ik)(jk)(lk)(mk).  We classify the 15 equality patterns of (i,j,l,m), using
(ak)(bk) = (b a k) for a ≠ b and (ak)² = 1:
 • all distinct: the 5-cycle (m l j i k); each 5-cycle with maximum k once.
 • {i=j}, {l},{m} distinct: (lk)(mk) = (m l k), 3-cycle through k, i free: k−3 choices per 3-cycle;  likewise
   {l=m}: (j i k), l free: k−3;  {j=l}: (ik)(mk) = (m i k), j free: k−3.
 • {i=m}: (ik)(jk)(lk)(ik) = (ik)(l j k)(ik) = (l j i): a 3-cycle on three points < k; each such 3-cycle arises from
   its 3 rotations: 3 times.
 • {i=l}: (ik)(jk)(ik)(mk) = (ij)(mk): a (2,2)-element containing k, arising 2 times (order of i,j);  likewise
   {j=m}: (ik)(lj): 2 times.
 • {i=j},{l=m}, i≠l: 1, (k−1)(k−2) times;  {i=m},{j=l}, i ≠ j: 1, (k−1)(k−2) times;
   {i=l},{j=m}, i≠j: ((ik)(jk))² = (j i k)² = (i j k): once per 3-cycle through k.
 • three equal: {i=j=l}: (ik)(mk) = (m i k); {i=j=m}: (lk)(ik) = (i l k); {i=l=m}: (ik)(jk) = (j i k);
   {j=l=m}: (ik)(jk) = (j i k): each once per 3-cycle through k.
 • all equal: 1, k−1 times.
Totals for fixed k: 5-cycles with maximum k: 1;  3-cycles through k (k maximal): 3(k−3) + 1 + 4 = 3k − 4;
3-cycles on points < k: 3;  (2,2)-elements containing k with the other points < k: 4;  identity: 2(k−1)(k−2) + (k−1)
= (k−1)(2k−3).  Summing over k: a 3-cycle with maximum b gets (3b − 4) + 3(n − b) = 3n − 4; a (2,2)-element gets 4
(only from k = its maximum); Σ_{k=1}^n (k−1)(2k−3) = Σ_{m=0}^{n−1} m(2m−1) = n(n−1)(4n−5)/6.

(I6) The coefficient of g in K_3² is N(g) = #{a a 3-cycle : a^{-1}g is a 3-cycle}.  g = 1: N = |K_3| = 2C(n,3).
g = (1 2 3): let S = supp(a).  |S ∩ {1,2,3}| = 3: a = (132) works ((123)(123) = (132)), a = (123) gives 1: one.
|S ∩ {1,2,3}| = 2, S = {x,y,w}, w ≥ 4: a^{-1}g is an even permutation of {1,2,3,w}; e.g. for {x,y} = {1,2}, w = 4:
(142)(123) = (2 3 4) works and (124)(123) = (14)(23) does not; the stabiliser ⟨g⟩ of g permutes the three pairs
{x,y} transitively, so for each pair and each w exactly one of the two 3-cycles works: 3(n−3).
|S ∩ {1,2,3}| ≤ 1: every point of S ∪ {1,2,3} is moved by a^{-1}g (a point y ∈ {1,2,3}∖S goes to g(y) ∈ {1,2,3},
then to a^{-1}(g(y)) ∈ S if g(y) ∈ S, else stays at g(y) ≠ y; a point of S∖{1,2,3} is moved by a^{-1}), at least 5
points: impossible.  Total 3n − 8.
g = (12)(34): as above every point of the symmetric difference S △ {1,2,3,4} is moved by a^{-1}g, so
|S ∩ {1,2,3,4}| ≥ 2.  If |S ∩ {1,2,3,4}| = 2, S = {x,y,w}: a^{-1}g must fix x and y, i.e. a(x) = g(x), a(y) = g(y)
∈ S, forcing g(x) = y, g(y) = x and a(x) = y, a(y) = x, impossible for a 3-cycle.  If S ⊂ {1,2,3,4}: all 8 choices
work, e.g. (123)(12)(34) = (1 3 4), (132)(12)(34) = (2 3 4), and the other sets follow by the symmetry of g.  Total 8.
g = (1 2 3 4 5): the moved set of a^{-1}g contains {1..5}∖S and S∖{1..5}, so S ⊂ {1..5} and a^{-1}g is a 3-cycle on
({1..5}∖S) ∪ {w} fixing the two other points x, y of S: a(x) = g(x) ∈ S, a(y) = g(y) ∈ S.  For a 5-cycle this forces
S = {x, g(x), g²(x)} and a = (x g(x) g²(x)) (the case a(x) = g(x) = w, a(y) = g(y) = x is the same family started at y);
then a^{-1}g is an even permutation moving exactly 3 points, hence a 3-cycle.  Total 5 (choice of x).
g of type (3,3): a and a^{-1}g are 3-cycles with product g, so their supports are the two cycles of g: 2 ways.
Any other even type with ≤ 6 moved points does not occur (a product of two 3-cycles moves ≤ 6 points; the even
types on ≤ 6 points are 1, (3), (2,2), (5), (3,3), (4,2)? no — (4,2) is odd; (2,2,2) is odd). ∎

---------------------------------------------------------------------------------------------------------
## 3. The JM scalar lemma and Theorem A   (PROVED; VERIFIED n ≤ 16 by MN, n ≤ 20 by the d-vector)

**Lemma 3.1.** For a symmetric polynomial F, F(J) acts on V^λ as the scalar F(contents of λ); thus
ω_λ(F(J)) = F(c_1, …, c_n) where (c_1,…,c_n) is the content multiset of λ, and in particular ω_λ(p_m(J)) = C_m(λ).
If moreover F(J) = Σ_ρ a_ρ K_ρ then F(contents of λ) = Σ_ρ a_ρ ω_λ(K_ρ).

*Proof.* By (F2), F(J)v_T = F(c_T(1),…,c_T(n)) v_T, and {c_T(k)}_k is the content multiset of λ for every T. ∎

**Theorem A.** For every n ≥ 0 and λ ⊢ n (with ω_λ(K_ρ) := 0 when n < |ρ|):
 (A1) ω_λ(K_3)  = C_2 − C(n,2).
 (A2) ω_λ(K_{22}) = C_1²/2 − (3/2)C_2 + C(n,2).
 (A3) ω_λ(K_5)  = C_4 − (3n − 10)C_2 − 2C_1² + n(n−1)(5n−19)/6.
 (A4) ω_λ(K_{33}) = ½[ω(K_3)² − 2C(n,3) − (3n−8)ω(K_3) − 8ω(K_{22}) − 5ω(K_5)]
                = C_2²/2 − (5/2)C_4 + 3C_1² − ½n²C_2 + (13/2)nC_2 − 15C_2 + n⁴/8 − (7/4)n³ + (47/8)n² − (17/4)n.
Equivalently, with |K_3| = 2C(n,3), |K_{22}| = 3C(n,4), |K_5| = 24C(n,5), |K_{33}| = (n)_6/18:
 χ^λ(3-cycle) = f(3C_2 − 3n(n−1)/2)/(n)_3,  χ^λ((2,2)) = f(C_1² − 3C_2 + n(n−1))/(6C(n,4)), etc.

*Proof.* Apply Lemma 3.1 to (I2), (I3), (I5) and use ω_λ(1) = 1, ω_λ(p_m(J)) = C_m:
(I2) gives C_2 = C(n,2) + ω(K_3).  (I3) with ω(K_2) = C_1 (from (I1)) and ω(K_2²) = ω(K_2)² (ω_λ is a homomorphism)
gives C_1² = C(n,2) + 3ω(K_3) + 2ω(K_{22}), i.e. (A2).  (I5) gives C_4 = ω(K_5) + (3n−4)ω(K_3) + 4ω(K_{22}) +
n(n−1)(4n−5)/6; substituting (A1), (A2): ω(K_5) = C_4 − (3n−4)(C_2 − C(n,2)) − 2C_1² + 6C_2 − 4C(n,2) −
n(n−1)(4n−5)/6, whose n-part is n(n−1)[(3n−4)/2 − 2 − (4n−5)/6] = n(n−1)(5n−19)/6.  (I6) gives (A4) by ω(K_3²) =
ω(K_3)²; the expansion was done symbolically (`topend_final_checks.log` (d)).  For n < |ρ| both sides vanish because
(I2)–(I6) hold with the empty classes equal to 0. ∎

Verification: `topend_algebra.py` (E): (A1)–(A4) for all λ ⊢ n ≤ 16 by MN (3660 checks); (F): the corresponding
d-vector identities of Theorem B for all λ ⊢ n ≤ 20 (10856 checks) without any character computation.

---------------------------------------------------------------------------------------------------------
## 4. Theorem B: u_i as character combinations; what the d-vector determines   (PROVED)

**Lemma 4.1 (square roots).** For g ∈ S_i of cycle type ρ, σ(ρ) := #{h : h² = g} equals
Π_{k odd} Σ_j C(m_k, 2j)(2j−1)!! k^j · Π_{k even} [m_k even]·(m_k − 1)!! k^{m_k/2},  m_k = multiplicity of k in ρ.
In particular σ(ρ) = 0 unless ρ is of *even type* (every even part has even multiplicity), σ(ρ∪1^m) = σ(ρ)·I_m if
ρ has no 1's, and σ(ρ) = Σ_{ν⊢i} χ^ν(ρ) (F5).

*Proof.* The square of an odd cycle of length k is a k-cycle on the same support; the square of a 2k-cycle is a pair
of k-cycles.  Hence a square root h of g is obtained by choosing a partial matching of the cycles of g into pairs of
equal length (all even-length cycles must be matched), choosing for each matched pair (two k-cycles) one of the k
2k-cycles whose square is that pair, and taking on each unmatched odd cycle its unique square root inside the cycle
(the ((k+1)/2)-th power).  Counting gives the formula; for k = 1 the factor is Σ_j C(m,2j)(2j−1)!! = I_m. ∎
(Formula checked against brute force for i ≤ 8 by the first-pass script `topend_verify.py`; σ = Σχ checked by MN
for i ≤ 10, `topend_levels.log` (6).)

**Theorem B.** For 0 ≤ i ≤ n and λ ⊢ n,
  u_i(λ) = Σ_{ρ ⊢ i} (σ(ρ)/z_ρ) χ^λ(ρ ∪ 1^{n−i}),   only even-type ρ contributing.
Explicitly (all coefficients from `topend_levels.log` (1); classes written without the padding 1's):
  u_3 = (2/3)f + (1/3)χ(3)
  u_4 = (5/12)f + (1/3)χ(3) + (1/4)χ(22)
  u_5 = (13/60)f + (1/3)χ(3) + (1/4)χ(22) + (1/5)χ(5)
  u_6 = (19/180)f + (2/9)χ(3) + (1/4)χ(22) + (1/5)χ(5) + (2/9)χ(33)
  u_7 = (29/630)f + (5/36)χ(3) + (1/6)χ(22) + (1/5)χ(5) + (2/9)χ(33) + (1/12)χ(322) + (1/7)χ(7)
  u_8 = (191/10080)f + (13/180)χ(3) + (5/48)χ(22) + (2/15)χ(5) + (2/9)χ(33) + (1/12)χ(322) + (1/7)χ(7) + (1/32)χ(2222) + (1/8)χ(44) + (1/15)χ(53)
  u_9 = (131/18144)f + (19/540)χ(3) + (13/240)χ(22) + (1/12)χ(5) + (4/27)χ(33) + (1/12)χ(322) + (1/7)χ(7) + (1/32)χ(2222) + (1/8)χ(44) + (1/15)χ(53) + (1/9)χ(9) + (1/20)χ(522) + (5/81)χ(333)
  u_10 = (1187/453600)f + (29/1890)χ(3) + (19/720)χ(22) + (13/300)χ(5) + (5/54)χ(33) + (1/18)χ(322) + (2/21)χ(7) + (1/32)χ(2222) + (1/8)χ(44) + (1/15)χ(53) + (1/9)χ(9) + (1/20)χ(522) + (5/81)χ(333) + (1/21)χ(73) + (3/25)χ(55) + (1/18)χ(3322)
Inverting the first four:
  χ^λ(3)  = 3u_3 − 2f,   χ^λ(22) = 4u_4 − 4u_3 + f,   χ^λ(5) = 5u_5 − 5u_4 + f,
  χ^λ(33) = (9/2)u_6 − (9/2)u_5 + (3/2)u_3 − f/2.
**Corollary B′.** The d-vector determines n, f, ω_λ(K_3), ω_λ(K_{22}), ω_λ(K_5), ω_λ(K_{33}); by Theorem A it
therefore determines C_2, then C_1², then C_4 (triangularly):
  C_2 = C(n,2) + (n)_3(3u_3 − 2f)/(6f),   C_1² = 3C_2 − n(n−1) + 6C(n,4)(4u_4 − 4u_3 + f)/f,
  C_4 = 24C(n,5)(5u_5 − 5u_4 + f)/f + (3n−10)C_2 + 2C_1² − n(n−1)(5n−19)/6.

*Proof of Theorem B.* By (F4), χ^λ(ρ∪1^{n−i}) = Σ_{ν⊢i} f^{λ/ν} χ^ν(ρ).  Hence
Σ_{ρ⊢i} (σ(ρ)/z_ρ) χ^λ(ρ∪1^{n−i}) = Σ_ν f^{λ/ν} Σ_ρ σ(ρ)χ^ν(ρ)/z_ρ = Σ_ν f^{λ/ν} (1/i!) Σ_{g∈S_i} σ(g)χ^ν(g)
(the class of ρ has i!/z_ρ elements) = Σ_ν f^{λ/ν} Σ_{μ⊢i} ⟨χ^μ, χ^ν⟩ = Σ_ν f^{λ/ν} = u_i, using (F5) and the
orthonormality of the (real) irreducible characters.  The coefficient tables are σ(ρ)/z_ρ from Lemma 4.1.  The
inversions are linear algebra on the first four lines (e.g. 4u_4 − 4u_3 + f = (5/3 − 8/3 + 1)f + (4/3 − 4/3)χ(3) +
χ(22)).  Corollary B′ then follows from Theorem A since |K_ρ|/f ≠ 0. ∎

Verification independent of all character theory: `topend_algebra.py` (F) checks, for every λ ⊢ n ≤ 20, the three
identities (3u_3−2f)·2C(n,3) = f·(A1), (4u_4−4u_3+f)·3C(n,4) = f·(A2), (5u_5−5u_4+f)·24C(n,5) = f·(A3) with u_i
from the Young-lattice recursion; `topend_final_checks.py` (a) checks the χ(33) inversion against MN for n ≤ 16.

---------------------------------------------------------------------------------------------------------
## 5. Theorem P: class sums are polynomials in n and the JM power sums, with a degree bound   (PROVED)

**Lemma 5.1 (fixed points cost).** Let t_1 ⋯ t_m = g be a product of m transpositions of a finite set, V the set
of points moved by at least one t_i, c_V(g) the number of cycles of g on V (fixed points counted), and
φ = #{x ∈ V : g(x) = x}.  Then  m ≥ (|V| − c_V(g)) + φ.   Note |V| − c_V(g) = m(type of g) (the reduced size).

*Proof.* Induction on m; m = 0 is trivial (V = ∅).  Let t_m = (x y), g′ = t_1⋯t_{m−1} = g t_m, V′ the points used
by t_1..t_{m−1}, and by induction m − 1 ≥ |V′| − c_{V′}(g′) + φ′.  Write R = |V| − c_V(g) + φ and
R′ = |V′| − c_{V′}(g′) + φ′; we show R ≤ R′ + 1 in every case.
 (i) x, y ∈ V′ (so V = V′).  g = g′(x y) merges the two cycles of g′ through x and y or splits the one cycle
 containing both.  Merge: c = c′ − 1; φ = φ′ − (number of x,y fixed by g′) ∈ {φ′, φ′−1, φ′−2}, so R ≤ R′ + 1.
 Split: c = c′ + 1; the split can create 0, 1 or 2 fixed points (φ ≤ φ′ + 2), and φ = φ′ + 2 only when the cycle
 is (x y) itself; in all cases R = R′ − 1 + (φ − φ′) ≤ R′ + 1.
 (ii) x ∈ V′, y ∉ V′: |V| = |V′| + 1, y is a fixed point of g′ not counted in φ′; g attaches y to the cycle of x, so
 c_V(g) = c_{V′}(g′), and φ = φ′ − [x fixed by g′] ≤ φ′.  R = |V′| + 1 − c_{V′}(g′) + φ ≤ R′ + 1.
 (iii) x, y ∉ V′: |V| = |V′| + 2, c_V(g) = c_{V′}(g′) + 1, φ = φ′: R = R′ + 1. ∎
(VERIFIED for all 813 616 sequences with m ≤ 5 on ≤ 6 points: `topend_algebra.log` (D).)

**Lemma 5.2 (minimal case).** If m = |V| − c_V(g) then φ = 0, the transposition multigraph on V (edge {a,b} for each
t_i = (a b), with multiplicity) is a forest without multiple edges, and the product of the transpositions of each
component (in the given order) is a single cycle on the vertex set of that component.

*Proof.* Transpositions from different components commute, so g restricted to a component's vertex set V_i is the
product of that component's e_i transpositions.  Lemma 5.1 gives e_i ≥ |V_i| − c_i + φ_i ≥ |V_i| − c_i and
connectedness gives e_i ≥ |V_i| − 1.  Summing, m ≥ Σ(|V_i| − c_i) = |V| − c_V(g) with equality only if every
e_i = |V_i| − c_i and φ_i = 0; then c_i ≥ 1 forces e_i ≤ |V_i| − 1 ≤ e_i, so c_i = 1 and e_i = |V_i| − 1: a
connected multigraph with |V_i| − 1 edges is a tree (no repeated edge). ∎

**Theorem P.** Let μ be a partition, m = |μ|, ℓ = ℓ(μ).  In Q[S_n], for every n ≥ 0,
  p_μ(J) = Σ_{ρ′ (no 1's)} a_{μ,ρ′}(n) K_{ρ′}(n),
where (i) a_{μ,ρ′} is a polynomial in n; (ii) a_{μ,ρ′} = 0 unless m(ρ′) ≤ m, and deg_n a_{μ,ρ′} ≤ m − m(ρ′);
(iii) if m(ρ′) = m then a_{μ,ρ′} is a constant, equal to 0 unless ℓ(ρ′) ≤ ℓ, and a_{μ,μ+1} = Π_k m_k(μ)!
where μ+1 = (μ_1+1, …, μ_ℓ+1) and m_k(μ) is the multiplicity of k in μ.
Consequently, for every partition ρ without 1's there is a polynomial F_ρ ∈ Q[n, p_1, …, p_{m(ρ)}] of weighted
degree ≤ m(ρ) (wt n = 1, wt p_k = k) with K_ρ(n) = F_ρ(n; p_1(J), …, p_{m(ρ)}(J)) for all n ≥ 0, and hence
  ω_λ(K_ρ) = F_ρ(n; C_1(λ), …, C_{m(ρ)}(λ))   for every n ≥ 0 and λ ⊢ n.

*Proof.* p_μ(J) = Π_t (Σ_k J_k^{μ_t}) is the sum over all sequences of m transpositions arranged in ℓ consecutive
"stars": the t-th star consists of μ_t transpositions (i k_t) with i < k_t (a common centre k_t).  By Lemma 2.2(c)
p_μ(J) is central, so the coefficient of g depends only on the cycle type ρ′ of g; take g with support {1,…,s},
s = |ρ′|.  A sequence with product g uses a point set V ⊇ {1..s}; the points of W := V ∖ {1..s} ⊆ {s+1,…,n} are fixed
points of g.  For two sets W, W′ ⊂ {s+1..n} of the same size j, the order-preserving bijection W → W′, extended by
the identity on {1..s} (and arbitrarily elsewhere) to a permutation π of [n], maps sequences using exactly W to
sequences using exactly W′ (it preserves all the constraints i < k because every element of W, W′ exceeds every
element of {1..s}, and the product is conjugated by π, which fixes g).  Hence the coefficient of g is
a_{μ,ρ′}(n) = Σ_j N_j C(n − s, j), where N_j is the number of sequences using exactly the outside set {s+1..s+j}
(N_j does not depend on n), a polynomial in n; for n < s the class is empty, so the identity holds for all n with
this polynomial.  By Lemma 5.1 applied to a sequence using j outside points, m ≥ m(ρ′) + j, so N_j = 0 for
j > m − m(ρ′): this is (ii), and (i).  If m(ρ′) = m then j = 0 (constant coefficient) and Lemma 5.2 applies: every
component of the transposition graph is a tree carrying one cycle of g, every star has distinct points (a repeated
point would be a repeated edge) and is a star subtree, two stars never share an edge, so a component is an
edge-disjoint union of star subtrees and has 1 + Σ_{stars in it} μ_t vertices.  Thus the parts of ρ′ are
1 + Σ_{t∈B} μ_t for the blocks B of a set partition of the stars, whence ℓ(ρ′) ≤ ℓ with equality iff every block
is a singleton iff ρ′ = μ + 1.  In that case each star is a factorisation of one cycle of g of length μ_t + 1 with
centre the maximum of the cycle's support, which by Lemma 2.1 exists in exactly one way; the stars of equal size
can be assigned to the cycles of equal length in Π_k m_k(μ)! ways; since the cycles of g commute, every assignment
gives product g.  This is (iii).
For the consequence, order the partitions ρ (no 1's) by (m(ρ), ℓ(ρ)) lexicographically and induct: K_∅ = 1; for ρ
with m(ρ) = m put μ = ρ − 1 (each part minus 1, so μ+1 = ρ, |μ| = m).  By (i)–(iii),
  K_ρ(n) = (1/a_{μ,ρ}) [ p_μ(J) − Σ_{ρ′ ≠ ρ} a_{μ,ρ′}(n) K_{ρ′}(n) ],
where every ρ′ occurring has either m(ρ′) < m or (m(ρ′) = m and ℓ(ρ′) < ℓ(ρ)), so K_{ρ′}(n) = F_{ρ′}(n; p(J))
for all n by induction.  This defines F_ρ ∈ Q[n, p_1..p_m], and its weighted degree is ≤ m because
deg(p_μ) = |μ| = m and deg(a_{μ,ρ′} F_{ρ′}) ≤ (m − m(ρ′)) + m(ρ′) = m.  The identity holds in Q[S_n] for every
n ≥ 0 (including n < |ρ|, where both sides are 0 by induction).  Apply ω_λ (a homomorphism, Lemma 3.1). ∎

### 5.3 From Theorem P to PROVED explicit polynomials

Let V_m = span_Q{ n^j C_μ : j + |μ| ≤ m } (C_μ = Π C_{μ_t}), dim V_m = Σ_{k≤m} p(k)(m + 1 − k) (= 7, 26, 75, 187 for
m = 2, 4, 6, 8).  By Theorem P, the function λ ↦ ω_λ(K_ρ) is the restriction of an element of V_{m(ρ)}.
If the evaluation map V_m → Q^T, T = {λ ⊢ n : n ≤ N}, is injective, then the unique element of V_m that agrees
with ω(K_ρ) on T is F_ρ, and the resulting identity holds for ALL n.  `topend_omega.py` does exactly this for the
16 even-type ρ with |ρ| ≤ 10: it computes ω_λ(K_ρ) on T by MN (N = 16), certifies injectivity by computing the rank
of the integer evaluation matrix modulo the prime 2^61 − 1 (a maximal minor nonzero mod p is nonzero over Q; the
ranks equal 1, 7, 7, 26, 26, 75, 26, 75, 75, 26, 187, 75, 75, 187, 187, 75 = dim V_m in every case; N = 16 is
needed: for m = 8 the rank on n ≤ 15 is only 186, `topend_rankscan.py`), solves the square pivot system exactly
over Q, and re-verifies the resulting polynomial on all 1597 partitions with n ≤ 18 (`topend_omega.log`).

### 5.4 The polynomials Ω_ρ(n; C) = ω_λ(K_ρ)   (PROVED for all n ≥ 0, all λ ⊢ n)

```
Ω_∅ = 1
Ω_(3) = C2 - n^2/2 + n/2
Ω_(2,2) = C1^2/2 - 3C2/2 + n^2/2 - n/2
Ω_(5) = -2C1^2 - 3C2 n + 10C2 + C4 + 5n^3/6 - 4n^2 + 19n/6
Ω_(3,3) = 3C1^2 + C2^2/2 - C2 n^2/2 + 13C2 n/2 - 15C2 - 5C4/2 + n^4/8 - 7n^3/4 + 47n^2/8 - 17n/4
Ω_(7) = 14C1^2 n - 72C1^2 - 8C1C3 - 9C2^2/2 + 21C2 n^2/2 - 241C2 n/2 + 252C2 - 5C4 n + 105C4/2 + C6
        - 49n^4/24 + 329n^3/12 - 2111n^2/24 + 751n/12
Ω_(3,2,2) = C1^2C2/2 - C1^2 n^2/4 + 25C1^2 n/4 - 26C1^2 - 4C1C3 - 3C2^2/2 + 5C2 n^2/4 - 125C2 n/4 + 75C2 + 15C4
        - n^4/4 + 13n^3/2 - 97n^2/4 + 18n
Ω_(5,3) = -2C1^2C2 + C1^2 n^2 - 61C1^2 n + 260C1^2 + 40C1C3 - 3C2^2 n + 65C2^2/2 + C2C4 + 7C2 n^3/3 - 63C2 n^2
        + 1442C2 n/3 - 868C2 - C4 n^2/2 + 61C4 n/2 - 455C4/2 - 7C6 - 5n^5/12 + 283n^4/24 - 629n^3/6 + 7013n^2/24 - 795n/4
Ω_(4,4) = 2C1^2 n^2 - 38C1^2 n + 267C1^2/2 - 2C1C3 n + 23C1C3 + 9C2^2 - 24C2 n^2 + 219C2 n - 819C2/2 + C3^2/2
        + 15C4 n - 105C4 - 7C6/2 + 4n^4 - 136n^3/3 + 134n^2 - 278n/3
Ω_(2,2,2,2) = C1^4/24 - 3C1^2C2/4 + C1^2 n^2/4 - 17C1^2 n/4 + 14C1^2 + 10C1C3/3 + 9C2^2/8 - 3C2 n^2/4 + 63C2 n/4
        - 35C2 - 35C4/4 + n^4/8 - 35n^3/12 + 83n^2/8 - 91n/12
Ω_(9) = 54C1^2C2 - 81C1^2 n^2 + 1305C1^2 n - 4112C1^2 + 72C1C3 n - 908C1C3 - 12C1C5 + 81C2^2 n/2 - 480C2^2 - 15C2C4
        - 81C2 n^3/2 + 1071C2 n^2 - 14549C2 n/2 + 12176C2 - 8C3^2 + 45C4 n^2/2 - 670C4 n + 3738C4 - 7C6 n + 168C6 + C8
        + 243n^5/40 - 351n^4/2 + 11637n^3/8 - 7713n^2/2 + 25713n/10
Ω_(5,2,2) = -C1^4 - 3C1^2C2 n/2 + 38C1^2C2 + C1^2C4/2 + 5C1^2 n^3/12 - 28C1^2 n^2 + 4291C1^2 n/12 - 1044C1^2
        + 20C1C3 n - 270C1C3 - 6C1C5 + 9C2^2 n/2 - 105C2^2 - 3C2C4/2 - 11C2 n^3/4 + 355C2 n^2/2 - 6159C2 n/4 + 2772C2
        + C4 n^2/2 - 211C4 n/2 + 840C4 + 28C6 + 5n^5/12 - 329n^4/12 + 1161n^3/4 - 10039n^2/12 + 1720n/3
Ω_(3,3,3) = 3C1^2C2 - 3C1^2 n^2/2 + 111C1^2 n/2 - 208C1^2 - 40C1C3 + C2^3/6 - C2^2 n^2/4 + 25C2^2 n/4 - 39C2^2
        - 5C2C4/2 + C2 n^4/8 - 19C2 n^3/4 + 563C2 n^2/8 - 1663C2 n/4 + 672C2 + 5C4 n^2/4 - 145C4 n/4 + 210C4 + 28C6/3
        - n^6/48 + 13n^5/16 - 203n^4/16 + 4201n^3/48 - 5263n^2/24 + 431n/3
Ω_(7,3) = 14C1^2C2 n - 366C1^2C2 - 7C1^2 n^3 + 484C1^2 n^2 - 6077C1^2 n + 17136C1^2 - 8C1C2C3 + 4C1C3 n^2 - 452C1C3 n
        + 4564C1C3 + 84C1C5 - 9C2^3/2 + 51C2^2 n^2/4 - 1499C2^2 n/4 + 2688C2^2 - 5C2C4 n + 315C2C4/2 + C2C6
        - 175C2 n^4/24 + 4055C2 n^3/12 - 133841C2 n^2/24 + 380137C2 n/12 - 48816C2 + 56C3^2 + 5C4 n^3/2 - 745C4 n^2/4
        + 14735C4 n/4 - 17262C4 - C6 n^2/2 + 113C6 n/2 - 1008C6 - 9C8 + 49n^6/48 - 11767n^5/240 + 14055n^4/16
        - 292573n^3/48 + 359053n^2/24 - 96957n/10
Ω_(5,5) = 2C1^4 + 6C1^2C2 n - 170C1^2C2 - 2C1^2C4 - 5C1^2 n^3/3 + 208C1^2 n^2 - 7819C1^2 n/3 + 7272C1^2 - 200C1C3 n
        + 1970C1C3 + 42C1C5 + 9C2^2 n^2/2 - 285C2^2 n/2 + 1085C2^2 - 3C2C4 n + 55C2C4 - 5C2 n^4/2 + 797C2 n^3/6
        - 4599C2 n^2/2 + 79543C2 n/6 - 20520C2 + 20C3^2 + C4^2/2 + 5C4 n^3/6 - 79C4 n^2 + 9289C4 n/6 - 7245C4 + 28C6 n
        - 420C6 - 9C8/2 + 25n^6/72 - 455n^5/24 + 3227n^4/9 - 60859n^3/24 + 451087n^2/72 - 16277n/4
Ω_(3,3,2,2) = 3C1^4/2 + C1^2C2^2/4 - C1^2C2 n^2/4 + 37C1^2C2 n/4 - 110C1^2C2 - 5C1^2C4/4 + C1^2 n^4/16 - 31C1^2 n^3/8
        + 1623C1^2 n^2/16 - 8269C1^2 n/8 + 2742C1^2 - 4C1C2C3 + 2C1C3 n^2 - 86C1C3 n + 817C1C3 + 21C1C5 - 3C2^3/4
        + C2^2 n^2 - 40C2^2 n + 771C2^2/2 + 75C2C4/4 - 7C2 n^4/16 + 217C2 n^3/8 - 10609C2 n^2/16 + 35331C2 n/8 - 7182C2
        + 8C3^2 - 35C4 n^2/4 + 1715C4 n/4 - 2520C4 - 126C6 + n^6/16 - 63n^5/16 + 1597n^4/16 - 12993n^3/16 + 16961n^2/8 - 1404n
```
The first five agree with Theorem A (an independent hand proof).  As predicted by the transpose symmetry
(χ^{λ'} = χ^λ on even-type classes, C_k(λ') = (−1)^k C_k(λ)) every monomial has an even number of odd indices.

Remark (not used).  By Kerov–Olshanski the functions n, C_1, C_2, … are algebraically independent on the set of all
partitions (C_k has top shifted-symmetric part p*_{k+1}/(k+1)); the full-rank certificate above is a finite, checkable
substitute for this fact in the range needed.

---------------------------------------------------------------------------------------------------------
## 6. Task (2): the levels u_3, …, u_10 — which content polynomial each one contributes

Define the *level polynomials* g_i(λ) := (n)_i u_i(λ)/f^λ.  From Theorem B, writing ρ = ρ̄ ∪ 1^{i−|ρ̄|} with ρ̄
free of 1's, σ(ρ) = σ(ρ̄) I_{i−|ρ̄|}, z_ρ = z_ρ̄ (i−|ρ̄|)!, and χ^λ(ρ̄∪1^{n−|ρ̄|}) = f ω_λ(K_ρ̄) z_ρ̄ (n−|ρ̄|)!/n!:

  **g_i = Σ_{ρ̄ even type, no 1's, |ρ̄| ≤ i} σ(ρ̄) · I_{i−|ρ̄|} · C(n − |ρ̄|, i − |ρ̄|) · Ω_ρ̄(n; C).**      (6.1)

`topend_levels.py` builds G_i := the right-hand side symbolically (sympy, exact) for i ≤ 10 from the proved Ω's and
verifies G_i(λ) = (n)_i d_{n−i}(λ)/f^λ for every λ ⊢ n ≤ 20 and every i ≤ 10 (29 570 checks, using only the
Young-lattice recursion for d) — `topend_levels.log` item (3).  Hence (PROVED) g_i = G_i for all n.

### 6.1 Levels 3, 4, 5 (PROVED)
  G_3 = C_2 + n(n−1)(4n−11)/6
  G_4 = C_1² + (n − 6)C_2 + 5n⁴/12 − 3n³ + 91n²/12 − 5n
  G_5 = C_4 + (n − 6)C_1² + (n² − 13n + 34)C_2 + 13n⁵/60 − 8n⁴/3 + 161n³/12 − 88n²/3 + 551n/30
So A_5 := Q[n, g_3, g_4, g_5] = Q[n, C_2, C_1², C_4] (triangular with unit leading coefficients).

### 6.2 Theorem C — level 6 is redundant (PROVED; VERIFIED n ≤ 20)
  G_6 = (n² − 11n + 42)C_1² + 2C_2² + (2n³/3 − 16n² + 328n/3 − 210)C_2 + (n − 15)C_4
        + 19n⁶/180 − 23n⁵/12 + 563n⁴/36 − 823n³/12 + 12623n²/90 − 171n/2,
i.e. **d_{n−6}(λ) = f^λ · G_6(n; C_1², C_2, C_4)/(n)_6 is a function of (n, d_n, d_{n−3}, d_{n−4}, d_{n−5})**.
*Proof.* (6.1) for i = 6 involves ρ̄ ∈ {∅, (3), (2,2), (5), (3,3)} with σ = 1, 1, 2, 1, 4, and Ω_(3,3) is a
polynomial in n, C_1², C_2, C_4 (Theorem A (A4)); collect terms.  Equivalently: χ^λ(3,3,1^{n−6}) = (9/2)u_6 −
(9/2)u_5 + (3/2)u_3 − f/2 (Theorem B) and (A4). ∎  Computationally there is no pair λ, μ ⊢ n ≤ 40 with equal
u_3/f, u_4/f, u_5/f and different u_6/f (`topend_witness.py`, `topend_witness2.py`), as it must be.

### 6.3 Levels 7–10: the new invariants (PROVED as polynomial identities; witnesses VERIFIED)
Write G_i = N_i + (element of A_{i−1}) where A_i := Q[n, g_3, …, g_i].  Using the successive normal forms
C_6 → h_7 + 16C_1C_3, C_3² → (h_8 + (8n+164)C_1C_3)/2, C_8 → h_9 + …, one finds (`topend_levels.log` (4)):

  level 7:  G_7 = h_7 + P_7,  **h_7 := C_6 − 16 C_1 C_3**,
            P_7 = C_1²C_2 + (2n³/3 − 25n²/2 + 659n/6 − 336)C_1² + (2n − 39/2)C_2² + (5n⁴/12 − 29n³/2 + 2083n²/12 − 1779n/2 + 1452)C_2
                  + (n² − 26n + 345/2)C_4 + 29n⁷/630 − 47n⁶/40 + 1009n⁵/72 − 395n⁴/4 + 146957n³/360 − 32243n²/40 + 40633n/84 ∈ A_5.
            (It comes from σ(7)Ω_(7) + σ(3,2,2)Ω_(3,2,2) = Ω_(7) + 2Ω_(3,2,2): C_6 − 8C_1C_3 + 2(−4C_1C_3) + ….)
  level 8:  new part  **h_8 := 2C_3² − (8n + 164) C_1 C_3**  (G_8 contains (n − 28)C_6 + (284 − 24n)C_1C_3 + 2C_3² + …; subtracting
            (n − 28)h_7 leaves the C_1C_3 coefficient −24n + 284 + 16(n − 28) = −8n − 164).
  level 9:  new part  **h_9 := C_8 − 24 C_1 C_5 − (32n − 7720/3) C_1 C_3**  (from Ω_(9) + 2Ω_(5,2,2) + 10Ω_(3,3,3)).
  level 10: new part  **h_10 := −24 C_1C_2C_3 + (12n² + 1060n + 62664) C_1C_3 − 360 C_1C_5**
            (from Ω_(7,3) + 6Ω_(5,5) + 8Ω_(3,3,2,2), after eliminating C_6, C_3², C_8 by h_7, h_8, h_9).
So: level 7 brings C_6 − 16C_1C_3; level 8 brings C_3² (modulo C_1C_3); level 9 brings C_8 − 24C_1C_5; level 10 brings
C_1C_2C_3 (modulo C_1C_3, C_1C_5).  Every new invariant beyond level 6 involves *products of odd content moments*;
the even moments C_6, C_8 never appear alone.  As polynomials in the formal variables n, C_1, C_2, … these are
genuinely new at each level (the reduced forms involve monomials absent from A_{i−1}); at the level of functions on
partitions the first witnesses are:
  level 7: n = 18, λ = (5,5,3,3,1,1), μ = (5,5,5,1,1,1): equal u_3/f, …, u_6/f, but u_7/f = 197/4641 vs 263/6188;
           also (6,4,4,2,2) vs (6,3,3,3,3) (n = 18) have equal (C_2, C_1², C_4) but h_7 = 21513 vs 30153.
  level 8: n = 40, λ = (9,7,7,7,2,2,2,2,1,1), μ = (9,8,7,4,4,3,2,1,1,1): equal u_3/f..u_7/f, u_8/f = 5310275/287110824
           vs 5310191/287110824.
  levels 9, 10: no witness with n ≤ 40 (`topend_witness2.log`).

---------------------------------------------------------------------------------------------------------
## 7. Task (3): does the top determine λ?   (VERIFIED data; PROVED class statements)

### 7.1 Data (`topend_separation.log` n ≤ 40, `topend_collisions.log` n ≤ 50, `topend_collisions2.log` n ≤ 58)
Key K_i = (n, r, f, u_3, …, u_i) (equivalently (n, r, f, C_2, C_1², C_4, …)); K′_i = the same without r.
Number of colliding pairs of transpose classes:
  depth 3 (n, r, f, C_2):  0 for n ≤ 13; first collision n = 14: (4,4,3,1,1,1) ~ (5,3,2,2,2); counts n = 17: 2, 19: 1, 20: 1, 22: 3, 23: 3, 24: 1, 25: 3, 26: 5, 27: 4, 28: 8, …, 40: 28 (zero at the other n ≤ 21).
  depth 4 (n, r, f, C_2, C_1²): 0 for n ≤ 27; collisions at n = 28 ((7,7,5,3,1^6) ~ (8,8,4,3,1^5)), 33, 36, 37, 38 (3), 40, 41 (4), …; without r also n = 29, 32.
  depth 5 (n, f, C_2, C_1², C_4), even without r: **0 collisions for all n ≤ 48**; n = 49: 1 pair,
      (9,8,8,6,3,3,3,3,3,1,1,1) ~ (11,7,5,5,4,4,4,4,3,1,1)  [r = 5 vs 6, so K_5 with r still separates];
      n = 50: 3 pairs, one of them with equal r = 6:  (12,10,8,3,3,2,2,2,2,2,2,1,1) ~ (12,11,8,3,3,2,2,2,2,2,1,1,1);
      n = 51..58: 3, 0, 4, 3, 3, 3, 3, 4 pairs (with r: 1, 0, 2, 1, 1, 1, 1, 1).
  depth 7 (n, f, C_2, C_1², C_4, h_7) i.e. (n, f, u_3, u_4, u_5, u_7): 0 collisions for all n ≤ 58.
The d-vectors of the depth-5 collision pairs differ from the top exactly from u_7 on (and from the bottom from
d_1 or d_2 on): `topend_pairs.log`.
Interpretation: the "top separating depth" is not bounded by 5; it is 3 for n ≤ 13, 4 for n ≤ 27, 5 for n ≤ 48, ≥ 6
(in fact ≥ 7, since level 6 is redundant) from n = 49 on.  Since deeper levels bring only products of odd content
moments (§6.3), the growth is expected to continue; a proof of Conjecture B from the top alone would have to use
unboundedly many levels.

### 7.2 Proved positive statements (PROVED)
Let S(x) = x(x+1)(2x+1)/6 = Σ_{c=1}^x c² for integers x ≥ −1 (S(−1) = S(0) = 0).
**Proposition 7.1.** (a) *Two-row shapes.* If λ = (a, b), a ≥ b ≥ 1, then C_2(λ) = S(a−1) + S(b−2) + 1, and for fixed
n = a + b the map (a,b) ↦ C_2 is injective.  (b) *Hooks.* If λ = (c, 1^d), then C_2 = S(c−1) + S(d), and for fixed n
= c + d the map is injective up to transpose ((c,1^d)^t = (d+1, 1^{c−1})).  (c) *Rectangles.* If λ = (a^b),
C_2 = n(a² + b² − 2)/12, so (n, C_2) determines {a, b}.  Hence (n, d_n, d_{n−3}) — through C_2 — determines the shape
up to transpose inside each of the three classes.
*Proof.* (a) Row 1 has contents 0..a−1, row 2 has −1, 0, …, b−2; so C_2 = S(a−1) + 1 + S(b−2).  Put x = a − 1,
y = b − 2, x + y = n − 3, x ≥ y + 1.  φ(x) := S(x) + S(n−3−x) satisfies φ(x+1) − φ(x) = (x+1)² − (n−3−x)² > 0 iff
x + 1 > n − 3 − x iff x > (n−5)/2, which holds since x ≥ (n−2)/2.  So φ is strictly increasing on the admissible
range and injective.  (b) Contents are 0..c−1 and −1..−d, C_2 = S(c−1) + S(d); with x = c − 1, y = d, x + y = n − 1
and transposition swapping x, y, assume x ≥ y, i.e. x ≥ (n−1)/2; ψ(x) = S(x) + S(n−1−x) has ψ(x+1) − ψ(x) =
(x+1)² − (n−1−x)² > 0 iff x > (n−2)/2, true.  (c) Σ_{i<b, j<a}(j − i)² = Σ((j − (a−1)/2) − (i − (b−1)/2))² =
b·Σ_j (j−(a−1)/2)² + a·Σ_i (i−(b−1)/2)² (the cross term vanishes because both sums are centred) =
ba(a²−1)/12 + ab(b²−1)/12 = n(a² + b² − 2)/12; then a + b = √(a² + b² + 2n) and ab = n determine {a, b}. ∎
Cross-class remark: (n, C_2) does *not* separate hooks from two-row shapes ((8,4) and (6,1^6) both have C_2 = 146),
but (n, C_2, f) separates the whole class C (two-row ∪ two-column ∪ hooks) up to transpose for all n ≤ 200
(`topend_classC_top.py`, VERIFIED only; the classC agent proved class C with the full vector).

### 7.3 Mirror box moves: why even moments plus f can never suffice (PROVED mechanism, VERIFIED existence)
**Proposition 7.2.** Let μ be obtained from λ ⊢ n by removing a corner box of content −c and adding a box of content
+c (c ≠ 0).  Then C_k(μ) = C_k(λ) for all even k, C_k(μ) = C_k(λ) + 2c^k for odd k; C_1(μ)² = C_1(λ)² iff C_1(λ) = −c,
and then C_1C_3(μ) − C_1C_3(λ) = 2c(C_3(λ) + c³).  In particular, if also f^μ = f^λ, then λ and μ share n, f, C_1²
and every even content moment, hence u_3, u_4, u_5, u_6 (Theorems B, C), while u_7 separates them as soon as
C_3(λ) ≠ −c³.
*Proof.* Immediate from C_k(μ) = C_k(λ) − (−c)^k + c^k. ∎
Existence (`topend_boxmove.py`, `topend_boxmove62.log`): for n ≤ 49 there is no such pair with f^μ = f^λ; for
n = 50, 51, 53, 54, 57, 58, 59, 60, 61 there are such pairs, all of the form
  λ = (λ_1, 10, λ_3, …, λ_10, 2, 1, [1]) with λ_1 ∈ {12,13},  μ = λ − (row 11, col 2) + (row 2, col 11),
i.e. the box of content −9 moved to its mirror position of content +9 (so C_1(λ) = −9), e.g.
  n = 50: (12,10,8,3,3,2^6,1,1) ↦ (12,11,8,3,3,2^5,1,1,1), f = 81059447377611523143373050000 for both,
  C_1: −9 ↦ 9, C_1C_3: 20169 vs −7047.
These are exactly the equal-r depth-5 collisions of §7.1 (the n = 50 pair with r = 6 is this one).  Whether the
family is infinite is open (the frame rows 2/11 confine the inner shape to a box, so other frames would be needed);
what is proved is that any argument using only (n, f, C_1², C_2, C_4, C_6, …) — all even moments — cannot separate
these pairs, so **a proof of Conjecture B must use either odd-moment products (levels ≥ 7) or the bottom of the vector**.

---------------------------------------------------------------------------------------------------------
## 8. Task (4): the hard pairs and the bottom of the vector

### 8.1 Lemma U (PROVED)
For s ≥ 0 let Y_{≥s}(λ) = {ν ⊆ λ : |ν| ≥ s} with the order induced from Young's lattice (graded by size).
**Lemma U.** If there is a rank-preserving poset isomorphism Y_{≥s}(λ) → Y_{≥s}(μ) sending λ to μ (|λ| = |μ| = n),
then d_j(λ) = d_j(μ) for all 0 ≤ j ≤ n − s.
*Proof.* d_j(λ) is the number of chains λ = ν^0 ⊃ ν^1 ⊃ … ⊃ ν^j with |ν^{t−1}/ν^t| = 1; all ν^t have size ≥ n − j ≥ s,
and the relation |ν/ν′| = 1 with ν′ ⊂ ν is "cover in Y_{≥s}" (a rank-preserving isomorphism preserves covers and
sizes), so chains correspond bijectively. ∎

### 8.2 Theorem F_a (PROVED; VERIFIED a ≤ 7, n ≤ 44 with sharpness: `topend_family.py`)
For a ≥ 3 and n ≥ 2a − 1 let λ = (a, 1^{n−a}) and μ = (2^{a−1}, 1^{n−2a+2}).  Then d_j(λ) = d_j(μ) for all
j ≤ n − 2a + 2.  (For a = 3 this is the known pair (3,1^{n−3}) ~ (2,2,1^{n−4}) agreeing on d_0..d_{n−4}; for a = 4,
(4,1^{n−4}) ~ (2,2,2,1^{n−6}) agree on d_0..d_{n−6}; computationally the agreement is sharp: d_{n−2a+3} differ for
a ≤ 7, n ≤ 44.)
*Proof.* Put b = n − a and s = 2a − 2.  The subshapes of λ are ∅ and the hooks (c, 1^k), 1 ≤ c ≤ a, 0 ≤ k ≤ b;
the subshapes of μ are ∅ and the shapes (2^{c′}, 1^{k′}), 0 ≤ c′ ≤ a − 1, k′ ≥ 0, c′ + k′ ≤ b + 1.  Hence
  Y_{≥s}(λ) = {(c,1^k) : 1 ≤ c ≤ a, 0 ≤ k ≤ b, c + k ≥ s},
  Y_{≥s}(μ) = {(2^{c′},1^{k′}) : 0 ≤ c′ ≤ a−1, k′ ≥ 0, c′ + k′ ≤ b+1, 2c′ + k′ ≥ s}.
Define φ(c,1^k) = (2^{c−1}, 1^{k−c+2}).  (Well defined on Y_{≥s}(λ): c′ = c−1 ∈ [0, a−1]; k′ = k − c + 2 ≥ 0
because c + k ≥ 2a − 2 and c ≤ a give k ≥ 2a − 2 − c ≥ c − 2; c′ + k′ = k + 1 ≤ b + 1; the size c + k = 2c′ + k′ is
preserved.)  φ is injective, and surjective: (2^{c′},1^{k′}) ∈ Y_{≥s}(μ) is the image of (c′+1, 1^{k′+c′−1}), where
k′ + c′ − 1 ≥ 0 because 2c′ + k′ ≥ s ≥ 4, k′ + c′ − 1 ≤ b, and 1 ≤ c′+1 ≤ a.  So φ is a size-preserving bijection.
Covers.  Inside λ, the shapes of size |ν| + 1 containing ν = (c,1^k) are (c+1, 1^k) (present iff c + 1 ≤ a) and
(c, 1^{k+1}) (iff k + 1 ≤ b); these are the covers of ν in Y_{≥s}(λ) (they automatically have size ≥ s).  Inside μ,
the shapes of size +1 containing ν′ = (2^{c′},1^{k′}) are (2^{c′+1}, 1^{k′−1}) (iff k′ ≥ 1 and c′+1 ≤ a−1; a box can
only be added to the first row of length 1) and (2^{c′}, 1^{k′+1}) (iff c′ + k′ + 1 ≤ b + 1).  Now
φ(c+1, 1^k) = (2^{c}, 1^{k−c+1}) = (2^{c′+1}, 1^{k′−1}) and φ(c, 1^{k+1}) = (2^{c′}, 1^{k′+1}), and the existence
conditions correspond: c + 1 ≤ a ⟺ c′ + 1 ≤ a − 1, and in that case k′ = k − c + 2 ≥ 1 holds automatically in
Y_{≥s}(λ) (c ≤ a − 1 and c + k ≥ 2a − 2 give k ≥ a − 1 ≥ c); k + 1 ≤ b ⟺ c′ + (k′+1) ≤ b + 1.  Thus φ maps the set of
cover relations of Y_{≥s}(λ) bijectively onto that of Y_{≥s}(μ).  Both posets are graded by size and their order is
the transitive closure of their cover relations, so φ is a rank-preserving poset isomorphism with φ(λ) = μ (take
c = a, k = b: φ(a,1^b) = (2^{a−1}, 1^{b−a+2}) = μ).  Lemma U gives d_j(λ) = d_j(μ) for j ≤ n − s = n − 2a + 2. ∎

### 8.3 Classification of pairs agreeing on the bottom (VERIFIED, `topend_hardpairs.py` n ≤ 30, `topend_hardpairs2.py` n ≤ 40)
Among all pairs of transpose classes of partitions of n:
 • agree on d_0..d_{n−4}: for n ≥ 8 **only** F_3 = {(3,1^{n−3}), (2,2,1^{n−4})}; sporadic small cases n = 4, 5
   (three-element groups) and n = 7: (4,1,1,1) ~ (3,2,2) (this one is *not* explained by Lemma U: the two shapes
   have 4 resp. 3 subshapes of size 4, yet d_3 = 8 for both).
 • agree on d_0..d_{n−5}: for n ≥ 9 only F_3 (extra cases: n = 8: (4,1^4) ~ (3,3,2) and (2,2,2,1,1) joining F_3).
 • agree on d_0..d_{n−6}: for n ≥ 11 exactly F_3 and F_4; n = 10 additionally (3,3,2,1,1) ~ (5,2,1,1,1).
Thus F_3 is the unique infinite family with "only the top three values differ" in the range examined, and it is
explained (and proved infinite) by Theorem F_a; more generally F_a agrees up to depth n − 2a + 2.

### 8.4 What this says about a proof of Conjecture B
 • The pair F_3 has f^{(3,1^{n−3})} − f^{(2,2,1^{n−4})} = C(n−1,2) − n(n−3)/2 = 1 and C_2 gap 4 (both trivial from
   the hook formula / contents; checked n ≤ 40, `topend_final_checks.log` (c)).  Every proof must therefore use
   d_n = f or d_{n−3} (= (2f + χ(3))/3): the bottom n − 3 entries d_0, …, d_{n−4} — which by the local product
   theorem already encode r, m, the multiset of {a_i, b_i} mod t^{2m+1}, etc. — cannot separate F_3.
 • Conversely (§7.1) the top alone, at any fixed depth, fails for large n: depth 3 from n = 14, depth 4 from n = 28,
   depth 5 from n = 49; and the mirror-box-move mechanism shows that even f together with *all* even content
   moments and C_1² cannot suffice in general.
 • Hence a proof must combine both ends (or use unboundedly many levels), consistent with the "separating depth"
   observations of the other angles.  What the top supplies in closed form is the list of Corollary B′ / §6.3:
   n, f, C_2, C_1², C_4, C_6 − 16C_1C_3, 2C_3² − (8n+164)C_1C_3, C_8 − 24C_1C_5 − …, ….

---------------------------------------------------------------------------------------------------------
## 9. Scripts of this pass (all in `src/agents/`, logs alongside)

| script | what it checks / produces | range |
|---|---|---|
| `topend_mn.py` | MN characters by beta-numbers, contents, σ, z, class sizes, ω, u-table DP | library |
| `topend_algebra.py` → `topend_algebra.log` | (A) MN cross-checks (3 implementations, orthogonality, transpose rule) n ≤ 10; (B) identities (I1)–(I6) in Z[S_n] 2 ≤ n ≤ 8; (C) JM trace identity for p_1..p_4 on every V^λ, n ≤ 8; (D) Lemma 5.1 on all 813 616 sequences (m ≤ 5, 6 points); (E) Theorem A via MN, all λ ⊢ n ≤ 16; (F) Theorem B/A identities and the level-6 identity using only the d-vector, all λ ⊢ n ≤ 20 | as stated |
| `topend_rankscan.py` | rank of V_m evaluation matrix vs N: full at N = 12 (m = 6), N = 16 (m = 8) | |
| `topend_omega.py` → `topend_omega.log`, `topend_omega_results.pkl` | the 16 Ω_ρ: full-rank certificate on n ≤ 16, exact solve, verification on all λ ⊢ n ≤ 18 | |
| `topend_levels.py` → `topend_levels.log` | σ = Σχ (i ≤ 10); u_i coefficient tables; G_i for i ≤ 10 and their verification against d-vectors n ≤ 20; reductions §6.3; level-7 witness | |
| `topend_witness.py`, `topend_witness2.py` | level witnesses (equal u_3/f..u_{i−1}/f, different u_i/f) n ≤ 30 / n ≤ 40 | |
| `topend_separation.py` → `topend_separation.log` | collisions of (n,r,f,u_3..u_i) and without r, i ≤ 10, n ≤ 40 | |
| `topend_collisions.py`, `topend_collisions2.py` → logs | explicit collision pairs, keys of depth 3, 4, 5 (n ≤ 50) and 5, 7 (n ≤ 58) | |
| `topend_pairs.py` → `topend_pairs.log` | anatomy of the depth-5 collisions at n = 49, 50 | |
| `topend_boxmove.py` → `topend_boxmove.log`, `topend_boxmove62.log` | mirror-box-move pairs with equal f and C_1² | n ≤ 62 |
| `topend_classC_top.py` → log | (n, C_2) vs (n, C_2, f) on class C | n ≤ 200 |
| `topend_hardpairs.py`, `topend_hardpairs2.py` → logs | pairs agreeing on d_0..d_{n−k}, k = 4..8 (n ≤ 30), k = 4,5,6 (n ≤ 40) | |
| `topend_family.py` → log | Theorem F_a and its sharpness, a = 3..7, n ≤ 44 | |
| `topend_final_checks.py` → log | χ(3,3) inversion (n ≤ 16); coefficient tables; F_3 gaps; symbolic expansion of (A4) | |

---------------------------------------------------------------------------------------------------------
## 10. Open problems / obstacles
1. Is the top separating depth unbounded?  (Data: 3 → 4 → 5 → ≥ 7 at n = 14, 28, 49.)  A construction of pairs
   agreeing on u_0..u_i for arbitrary i would show that the top alone can never prove B.
2. Are the mirror-box-move pairs with f^λ = f^μ (§7.3) infinite in number?  The examples all have the frame
   λ_2 = 10, λ_11 = 2; a hook-length identity behind the f-equality was not found.
3. Prove the classification of §8.3 (only F_3 agrees on d_0..d_{n−4} for n ≥ 8).  The local product theorem gives
   d_j for j ≤ 2m only; the pairs F_a suggest studying "upper interval isomorphisms" systematically.
4. Levels 9 and 10 have no function-level witness for n ≤ 40 although they are new as polynomials; find the first
   witnesses (probably n > 40, as at level 8 the first is n = 40).
