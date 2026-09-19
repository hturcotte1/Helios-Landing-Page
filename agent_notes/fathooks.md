# Angle "fathooks": Conjecture B for two-corner shapes λ(a,b,c,d) = ((a+c)^b, c^d)

Notation. For a,b,c,d ≥ 1 let λ = λ(a,b,c,d) = ((a+c)^b, c^d); its corner runs are ((a,b),(c,d)), n = ab+bc+cd,
and λ^t = λ(d,c,b,a). Every partition with exactly two removable corners is of this form (r = 2 ⇔ two distinct
part sizes), so "Conjecture B for two-corner shapes" is: **if μ ⊢ n and d(μ) = d(λ(a,b,c,d)) then μ ∈ {λ, λ^t}.**
Since d_1(μ) = r(μ) (Briefing fact 2), such a μ is itself two-corner, μ = λ(a',b',c',d'); so the conjecture is
exactly the statement that the d-vector determines (a,b,c,d) up to reversal (a,b,c,d) ↦ (d,c,b,a) — no
"among two-corner shapes" restriction is needed. Series: I(t) = e^{t+t²/2} = Σ_k I_k t^k/k! (I_k = involution
numbers), G_{x,y}(t) = P_{x^y}(t) = Σ_{ρ ⊆ y×x box} f^ρ t^{|ρ|}/|ρ|! (Lemma 4.1 of results.md; G_{x,y} = G_{y,x}),
D_x(t) = Σ_{ρ : ρ_1 > x} f^ρ t^{|ρ|}/|ρ|!, E_{x,y}(t) = Σ_{ρ : ρ_1 > x, ℓ(ρ) > y} f^ρ t^{|ρ|}/|ρ|! (sums over all
partitions ρ of all sizes). "Filter" = subset S ⊆ λ closed under moving right/down inside λ, e(S) = number of
removal orders, d_j = Σ_{|S|=j} e(S) (Briefing). N_k := #{slots of (a,b,c,d) equal to k}; x1 ≤ x2 ≤ x3 ≤ x4 the
sorted parameters (so x1 = m of Theorem 4.2 and N_{x1} = A_m + B_m).

Scripts (all exact integer / Fraction arithmetic; all in `src/agents/`):
* `fathooks2_lemmas.py N NB` — verifies Lemmas 1, 2, Theorem 3 (reading algorithm), the formulas (S1)–(S4) and
  Theorem F (see the log `fathooks2_lemmas.log`, run with N = 40, NB = 30).
* `fathooks2_invariants.py N` — separation power of the closed-form invariants (n, f, C_2, C_1², C_4, local data)
  on all two-corner shapes with n ≤ N (logs: N = 150 in the terminal, N = 300 in `fathooks2_invariants300.log`).
* `fathooks2_c2.py` — closed form of n, C_1, C_2 (sympy) checked against direct summation for n ≤ 60.
* `fathooks2_arrangement.py M` — the "arrangement" statement for all multisets with values ≤ M.

Summary of status.
* PROVED: Lemma 1 (two-window locality up to degree a+d, with the exact first correction), Lemma 2 (box/up-census
  window up to degree min(b,c), with the exact first correction), Theorem 3 (the d-vector determines x1, x2 and
  N_k for every k ≤ x1+x2−1; in particular the whole multiset {a,b,c,d} when x4 ≤ x1+x2−1, and it detects whether
  this is the case), Theorem F (Conjecture B, against **all** partitions μ, for every two-corner shape with at least
  two parameters equal to 1 — these are class C together with the three families F_1, F_2, F_4 of Theorem 4.9).
* VERIFIED (not proved): (n, f^λ, C_2) separates all two-corner shapes up to transpose for n ≤ 300 (604 513 shapes);
  since n, f^λ, C_2 are PROVED to be determined by the d-vector (Thm 4.2, fact 2, Thm 4.6), Conjecture B for all
  two-corner shapes would follow from this purely Diophantine statement (Conjecture D below). The weaker
  "arrangement" statement with (n, C_2) alone is FALSE: (6,7,8,10) and (7,10,6,8) (n = 178) have the same multiset,
  the same n and the same C_2 = 9385 (they have different f).
* Conjecture B itself is VERIFIED for all two-corner shapes with n ≤ 40 against each other (`fathooks2_lemmas.py`),
  and for all two-corner shapes with ≥ 2 unit parameters against all partitions of n ≤ 30.

---------------------------------------------------------------------------------------------------------------

## 1. Two exact windows (PROVED)

Write the boxes of λ as (i,j), 1 ≤ i ≤ b+d, 1 ≤ j ≤ λ_i, and split λ into
R1 = {i ≤ b, c < j ≤ a+c} (a copy of the rectangle a^b), R2 = {b < i ≤ b+d, j ≤ c} (a copy of c^d) and
R0 = {i ≤ b, j ≤ c}.

**Lemma 1 (PROVED).** (i) For 0 ≤ j ≤ a+d: d_j(λ) = j! [t^j] G_{a,b}(t) G_{c,d}(t), i.e.
d_j(λ) = Σ_{j1+j2=j} C(j,j1) d_{j1}(a^b) d_{j2}(c^d).
(ii) If a+d+1 ≤ n: d_{a+d+1}(λ) = (a+d+1)! [t^{a+d+1}] G_{a,b} G_{c,d} + C(a+d, a).

*Proof.* Let S be a filter of λ containing a box (i,j) ∈ R0. Then S contains the whole part of row i to the right
of (i,j), namely the a+c−j+1 ≥ a+1 boxes (i,j'), j ≤ j' ≤ a+c (all lie in λ since i ≤ b), and the whole part of
column j below (i,j), namely the b+d−i+1 ≥ d+1 boxes (i',j), i ≤ i' ≤ b+d (all lie in λ since j ≤ c). These two
sets share only (i,j), so |S| ≥ (a+1)+(d+1)−1 = a+d+1, with equality only if i = b, j = c and S is exactly the
hook H := {(b,j') : c ≤ j' ≤ a+c} ∪ {(i',c) : b ≤ i' ≤ b+d}.

Hence a filter S with |S| ≤ a+d satisfies S ⊆ R1 ∪ R2. Put S1 = S ∩ R1, S2 = S ∩ R2. S1 is a filter of the
rectangle R1 (if (i,j) ∈ S1 and (i',j') ∈ R1 with i' ≥ i, j' ≥ j then (i',j') ∈ S because S is a filter of λ,
and it lies in R1), and likewise S2 is a filter of R2. Conversely, for any filters S1 of R1 and S2 of R2 the union
is a filter of λ: if (i,j) ∈ S1 and (i',j') ∈ λ with i' ≥ i, j' ≥ j, then j' ≥ j > c forces i' ≤ b, so
(i',j') ∈ R1 and therefore (i',j') ∈ S1; if (i,j) ∈ S2 then i' ≥ i > b forces j' ≤ c, so (i',j') ∈ R2 and
(i',j') ∈ S2. No box of R1 lies weakly right-and-below a box of R2 or vice versa (rows of R1 are < b+1 ≤ rows of
R2, columns of R1 are > c ≥ columns of R2), so the removal poset of S is the disjoint union of the posets of S1
and S2 with no relations between them, and e(S) = C(|S|,|S1|) e(S1) e(S2) (interleave two removal orders).
Summing over all pairs (S1,S2) with |S1|+|S2| = j gives d_j(λ) = Σ_{j1+j2=j} C(j,j1) d_{j1}(a^b) d_{j2}(c^d), which
is (i), because P_{a^b} = G_{a,b} and P_{c^d} = G_{c,d}.

For |S| = a+d+1 the same count applies to the filters S ⊆ R1 ∪ R2, and the only other filter is H (it is a filter:
the boxes of λ weakly right-and-below (b,j') with j' > c are the (b,j''), j'' ≥ j', which lie in H; those below
(i',c) with i' > b are the (i'',c), i'' ≥ i', in H; those of (b,c) are row b from column c on, column c from row b
on, and nothing else because λ has no box (i'',j'') with i'' > b, j'' > c). In the removal order of H, the box
(b,c) is removed last (every other box of H is to its right or below it), the arm (b,c+1),…,(b,a+c) must be removed
from right to left and the leg (b+1,c),…,(b+d,c) from bottom to top; so the removal orders are the interleavings of
a chain of length a with a chain of length d: e(H) = C(a+d,a). This gives (ii). ∎

**Lemma 2 (PROVED).** Let B be the (b+d)×(a+c) box and m = min(b,c).
(i) (box formulation) d_j(λ) = number of saturated chains a^d = ρ^0 ⊂ ρ^1 ⊂ … ⊂ ρ^j in Young's lattice with
ρ^j ⊆ B (j boxes added one at a time to the rectangle a^d, never leaving the box).
(ii) For 0 ≤ j ≤ m: d_j(λ) = s_j(a^d) = j! [t^j] I(t) G_{a,d}(t) (up-census of the rectangle a^d).
(iii) If m+1 ≤ n: d_{m+1}(λ) = (m+1)! [t^{m+1}] I G_{a,d} − [b=m] − [c=m].

*Proof.* (i) Let ν ⊆ λ. Its complement B∖ν in the box is closed under moving right/down in B (ν is closed under
moving left/up), so its image under the rotation (i,j) ↦ (b+d+1−i, a+c+1−j) of B is closed under moving left/up,
i.e. it is the diagram of a partition ρ(ν) ⊆ B, and ρ(ν) ⊇ ρ(λ) = rotation of B∖λ = {b<i≤b+d, c<j≤a+c}, which is
the rectangle a^d (d rows of length a). The map ν ↦ ρ(ν) is a bijection from {ν ⊆ λ} onto {ρ : a^d ⊆ ρ ⊆ B}
(inverse: ν = B ∖ rotation of ρ, which is closed under moving left/up and contained in λ), with |ρ(ν)| = ad+|λ/ν|,
and it carries the skew diagram λ/ν onto ρ(ν)/a^d reversing the order; reading a standard filling of λ/ν with the
entries k ↦ |λ/ν|+1−k gives a standard filling of ρ(ν)/a^d, so f^{λ/ν} = f^{ρ(ν)/a^d}. Since f^{ρ/a^d} is the number
of saturated chains from a^d to ρ, d_j(λ) = Σ_{|λ/ν|=j} f^{λ/ν} = Σ_{a^d ⊆ ρ ⊆ B, |ρ|=ad+j} f^{ρ/a^d} is the number
of j-step chains from a^d whose top lies in B (a chain lies in B iff its top does).
(ii) Adding j ≤ m boxes to a^d (d rows, a columns) produces ρ with ℓ(ρ) ≤ d+j ≤ d+b and ρ_1 ≤ a+j ≤ a+c, so every
j-step chain from a^d stays in B and d_j(λ) = Σ_{ρ ⊇ a^d, |ρ/a^d|=j} f^{ρ/a^d} = s_j(a^d), and the up-census
identity Σ_k s_k(ν) t^k/k! = e^{t+t²/2} P_ν(t) (Briefing, PROVED) gives the series form.
(iii) A chain of m+1 steps from a^d leaves B iff its top ρ has ℓ(ρ) = d+m+1 or ρ_1 = a+m+1. The first needs
b = m and forces each added box to lie in a new row, so ρ = (a^d,1^{m+1}) and ρ/a^d is a single column
(one chain); the second needs c = m and ρ = (a+m+1, a^{d−1}), a single row (one chain). Both cannot happen for the
same ρ (m+1 ≥ 2 boxes are not simultaneously a column and a row). Subtracting these chains from s_{m+1}(a^d) gives
(iii). ∎

*Computational verification* (`fathooks2_lemmas.py 40 30`): Lemma 1 (i),(ii) and Lemma 2 (ii),(iii) for all
two-corner shapes with n ≤ 40 (every j in the stated ranges), Lemma 2 (i) as the full vector for n ≤ 24.
See §6 for the log line.

## 2. What the bottom of the vector reads off: local multiset data (PROVED)

**Theorem 3 (PROVED).** Let λ = λ(a,b,c,d), x1 ≤ x2 ≤ x3 ≤ x4 its sorted parameters, N_k = #{slots equal to k}.
Then, with Q(t) := I(t) − P_λ(t)/I(t),
  Q(t) ≡ Σ_k N_k D_k(t)  (mod t^{x1+x2+1}).                                                    (3.1)
Consequently the d-vector determines x1, x2 and N_k for every k ≤ x1+x2−1, by the following algorithm which uses
only the d-vector: let q_j = j! [t^j] Q. Then q_j = 0 for j ≤ x1 and q_{x1+1} = N_{x1} ≥ 1, which gives x1 and
N_{x1}. If N_{x1} ≥ 2 then x2 = x1. Otherwise x2 is the least k > x1 with N_k ≠ 0, found by peeling: for
k = x1+1, x1+2, … compute N_k = q_{k+1} − Σ_{x<k} N_x D_x^{(k+1)}, where D_x^{(k+1)} := (k+1)! [t^{k+1}] D_x =
Σ_{ρ ⊢ k+1, ρ_1 > x} f^ρ is a universal constant; all these reads are valid while k ≤ x1+x2−1, and x2 itself is
≤ x1+x2−1, so x2 is reached; then continue peeling up to k = x1+x2−1.

*Proof.* By Lemma 4.1 and inclusion–exclusion on the two conditions ρ_1 ≤ x, ℓ(ρ) ≤ y (using f^ρ = f^{ρ^t}, so that
Σ_{ℓ(ρ)>y} f^ρ t^{|ρ|}/|ρ|! = D_y), G_{x,y} = I − D_x − D_y + E_{x,y}. Orders of vanishing: D_x = t^{x+1}/(x+1)! +
O(t^{x+2}) (the only ρ ⊢ x+1 with ρ_1 > x is (x+1), f = 1); E_{x,y} = O(t^{x+y+1}) (ρ_1 ≥ x+1 and ℓ(ρ) ≥ y+1 give
|ρ| ≥ ρ_1+ℓ(ρ)−1 ≥ x+y+1). By Lemma 1(i), P_λ ≡ G_{a,b}G_{c,d} (mod t^{a+d+1}), and a+d ≥ x1+x2 because a and d
are two distinct slots (any two distinct slots have sum ≥ x1+x2). Expanding the product modulo t^{x1+x2+1}:
I·E_{a,b} and I·E_{c,d} vanish (orders a+b+1, c+d+1 ≥ x1+x2+1); (D_a+D_b)(D_c+D_d) vanishes (order ≥
min(a,b)+min(c,d)+2 ≥ x1+x2+2, as min(a,b), min(c,d) are distinct slots); all products involving an E and a D or two
E's have order ≥ x1+x2+1+x1+1. Hence P_λ ≡ I² − I·(D_a+D_b+D_c+D_d), and dividing by the invertible series I
(I(0) = 1) gives (3.1). For the reading: (k+1)![t^{k+1}] D_x = D_x^{(k+1)} equals 0 for x > k and 1 for x = k,
so for k+1 ≤ x1+x2, q_{k+1} = N_k + Σ_{x<k} N_x D_x^{(k+1)}: q_j = 0 for j ≤ x1 (no slot is < x1), q_{x1+1} =
N_{x1}, and the displayed peeling recursion holds for all k ≤ x1+x2−1. The reader knows the validity range only
after finding x2, but every index used before that is ≤ x2+1 ≤ x1+x2 (as x1 ≥ 1), so the algorithm never reads
outside the valid range. ∎

**Corollary 3.1 (PROVED).** (a) The d-vector determines the multiset of those parameters that are ≤ x1+x2−1, and
the number 4 − Σ_{k ≤ x1+x2−1} N_k of parameters that are ≥ x1+x2. (b) In particular, if x4 ≤ x1+x2−1
("balanced" shapes) the d-vector determines the whole multiset {a,b,c,d}, and it recognises balanced shapes. (c) For
x1 = 1 it determines N_1 (this is Theorem 4.2(4): d_2 = 6 − N_1) and, when N_1 = 1, also x2 and N_{x2}.

*Computational verification*: the reading algorithm (implemented exactly as stated, `read_local` in
`fathooks2_lemmas.py`) returns the correct (x1, x2, (N_k)_{k ≤ x1+x2−1}) for every two-corner shape with n ≤ 40.

*Remark (what the local data cannot do).* Even together with n, the local data are far from sufficient: e.g.
(n, x1, x2, (N_k)) collides already at n = 5 (λ(1,1,1,3) vs λ(1,1,2,1)); see §5.

## 3. Closed formulas for the F-families and rectangles (PROVED)

Let F_1(c,d) = (c+1, c^d) = λ(1,1,c,d), F_2(b,c) = ((c+1)^b, c) = λ(1,b,c,1), F_4(a,b) = ((a+1)^b, 1) = λ(a,b,1,1)
= F_1(b,a)^t (Theorem 4.9 notation), and note F_2(b,c)^t = F_2(c,b). Put x = min of the two parameters.

**(S4) Rectangles** (from Lemma 4.1: d_i(c^d) = Σ_{ρ ⊢ i, ρ_1 ≤ c, ℓ(ρ) ≤ d} f^ρ), x = min(c,d):
d_i(c^d) = I_i for i ≤ x; d_{x+1}(c^d) = I_{x+1} − [c=x] − [d=x] (if x+1 ≤ cd);
d_{x+2}(c^d) = I_{x+2} − (x+2)([c=x]+[d=x]) − [c=x+1] − [d=x+1] (if x+2 ≤ cd).
*Proof.* For i ≤ x every ρ ⊢ i fits in the box. For i = x+1 the only ρ ⊢ x+1 with ρ_1 > c is (x+1) when c = x, the
only one with ℓ(ρ) > d is 1^{x+1} when d = x; they are distinct (x+1 ≥ 2). For i = x+2: the ρ ⊢ x+2 with ρ_1 ≥ c+1
are, if c = x, (x+2) and (x+1,1) with f = 1 and x+1 (total x+2), and if c = x+1, only (x+2); likewise for ℓ(ρ) ≥ d+1;
a ρ ⊢ x+2 with both ρ_1 ≥ c+1 and ℓ(ρ) ≥ d+1 would need c+d+2 ≤ ρ_1+ℓ(ρ) ≤ x+3, i.e. c+d ≤ x+1, which forces
c = d = x = 1, excluded by x+2 ≤ cd. ∎

**(S1),(S2) F_2.** The rectangle (c+1)^{b+1} has a single removable corner and removing it gives F_2(b,c); hence by
the recursion d_{j+1}((c+1)^{b+1}) = d_j(F_2(b,c)) for all j ≥ 0 (S1). With y = min(b,c), by (S4) applied to the
rectangle (c+1)^{b+1} (whose minimum side is y+1):
d_j(F_2) = I_{j+1} (j ≤ y); d_{y+1}(F_2) = I_{y+2} − [b=y] − [c=y];
d_{y+2}(F_2) = I_{y+3} − (y+3)[b=y] − [b=y+1] − (y+3)[c=y] − [c=y+1]     (S2)
(the last two need y+2 ≤ (b+1)(c+1)−1 = n, true as n ≥ y²+2y ≥ y+2).

**(S3) F_1.** F_1(c,d) = λ(1,1,c,d), so Lemma 1 applies with a+d = d+1 and G_{1,1} = 1+t:
d_j(F_1) = d_j(c^d) + j·d_{j−1}(c^d) for j ≤ d+1, and d_{d+2}(F_1) = d_{d+2}(c^d) + (d+2) d_{d+1}(c^d) + (d+1)
(if d+2 ≤ n). With x = min(c,d) and the involution recursion I_{k+1} = I_k + k I_{k−1}, (S4) gives:
d_j(F_1) = I_{j+1} (j ≤ x); d_{x+1}(F_1) = I_{x+2} − [c=x] − [d=x];
if c = x < d: d_{x+2}(F_1) = I_{x+3} − 2(x+2) − [d=x+1];
if d = x ≤ c and c,d ≥ 2: d_{x+2}(F_1) = I_{x+3} − (x+3) − 2(x+2)[c=x] − [c=x+1].     (S3)
*Proof.* j ≤ x ≤ d+1: I_j + jI_{j−1} = I_{j+1}. j = x+1 ≤ d+1: (I_{x+1} − [c=x] − [d=x]) + (x+1)I_x. If c = x < d then
x+2 ≤ d+1 and x+2 ≤ cd, so d_{x+2}(F_1) = [I_{x+2} − (x+2) − [d=x+1]] + (x+2)[I_{x+1} − 1] = I_{x+3} − 2(x+2) −
[d=x+1]. If d = x ≤ c with c,d ≥ 2, then x+2 = d+2 = a+d+1 and x+2 ≤ cd (cd ≥ 2x ≥ x+2), so by the corrected formula
d_{x+2}(F_1) = [I_{x+2} − (x+2)(1+[c=x]) − [c=x+1]] + (x+2)[I_{x+1} − 1 − [c=x]] + (x+1) = I_{x+3} − (x+3) −
2(x+2)[c=x] − [c=x+1]. ∎

*Computational verification* (`fathooks2_lemmas.py`, part "S1–S4"): all four formula groups for all parameter
values with n ≤ 60 (S3 last line for c,d ≥ 2), 0 failures; see §6.

## 4. Theorem F: Conjecture B for two-corner shapes with at least two unit parameters (PROVED)

**Theorem F (PROVED).** Let λ have exactly two removable corners, corner runs ((a,b),(c,d)), and suppose at least
two of a,b,c,d equal 1. If μ ⊢ n = |λ| and d(μ) = d(λ), then μ ∈ {λ, λ^t}.

The shapes in question are: b = d = 1 (two-row), a = c = 1 (two-column), b = c = 1 (hook (a+1,1^d)) — these are
class C of Theorem 4.9 — and a = b = 1 (F_1(c,d)), a = d = 1 (F_2(b,c)), c = d = 1 (F_4(a,b)). If an F-shape has a
parameter equal to 1 it lies in C (F_1(1,d) = (2,1^d), F_1(c,1) = (c+1,c), F_2(1,c) = (c+1,c), F_2(b,1) = (2^b,1),
F_4 = F_1^t), so the new cases are F_1(c,d), F_2(b,c), F_4(a,b) with both parameters ≥ 2.

*Proof.* If λ ∈ C, this is Theorem 4.9. So let λ be F_1, F_2 or F_4 with both parameters ≥ 2; then N_1(λ) = 2
exactly. Since d_1(μ) = d_1(λ) = 2, μ = λ(a',b',c',d') is two-corner, and by Theorem 4.2(4), d_2 = 6 − N_1 for
two-corner shapes, so N_1(μ) = 2 and μ ∈ C or μ is an F-shape with both parameters ≥ 2. If μ ∈ C, Theorem 4.9
applied to the pair (μ, λ) gives λ ∈ {μ, μ^t} ⊆ C, a contradiction. Hence both λ and μ are F-shapes with parameters
≥ 2. The hypothesis and the conclusion are invariant under replacing λ by λ^t and/or μ by μ^t (d(ν) = d(ν^t)), and
F_4(a,b) = F_1(b,a)^t, F_2(b,c)^t = F_2(c,b); so we may assume each of λ, μ is F_1(c,d) with c,d ≥ 2 or F_2(b,c) with
2 ≤ b ≤ c.

*Step 1 (the invariants x and ε).* For both families with parameters ≥ 2 and x the minimum parameter, (S2)/(S3)
give d_j = I_{j+1} for j ≤ x and d_{x+1} = I_{x+2} − ε with ε = [c=x]+[d=x] ∈ {1,2} for F_1 and ε = [b=x]+[c=x] =
1+[b=c] for F_2 (x+1 ≤ n in all cases). Hence x = min{j : d_j ≠ I_{j+1}} − 1 and ε = I_{x+2} − d_{x+1} are
determined by the d-vector, and they coincide for λ and μ.

*Step 2 (F_2 vs F_2).* If λ = F_2(b,c), μ = F_2(b',c'), then by (S1) the rectangles (c+1)^{b+1} and (c'+1)^{b'+1}
(both of size n+1) have d_{j+1} equal for all j ≥ 0, and d_0 = 1 for both, so they have equal d-vectors; by Theorem
4.4, {b+1,c+1} = {b'+1,c'+1}, so μ ∈ {λ, λ^t}.

*Step 3 (F_1 vs F_1).* λ = F_1(c,d), μ = F_1(c',d'), c,d,c',d' ≥ 2, n = cd+c+1 = c'd'+c'+1. If ε = 2 then
c = d = x = c' = d'. If ε = 1, each shape is in exactly one of the sub-cases "c = x < d" or "d = x < c". If both are
in the first, n = x(d+1)+1 = x(d'+1)+1 gives d = d'; if both are in the second, n = c(x+1)+1 = c'(x+1)+1 gives
c = c'. If λ is in the first and μ in the second, (S3) gives d_{x+2}(λ) = I_{x+3} − 2(x+2) − [d=x+1] and
d_{x+2}(μ) = I_{x+3} − (x+3) − [c'=x+1] (x+2 ≤ n), whose difference is −(x+1) − [d=x+1] + [c'=x+1] ≤ −x < 0, a
contradiction. So μ = λ.

*Step 4 (F_1 vs F_2).* λ = F_1(c,d) with c,d ≥ 2, μ = F_2(b',c') with 2 ≤ b' ≤ c', n = cd+c+1 = (b'+1)(c'+1)−1.
If ε = 2: c = d = x and b' = c' = x, so x²+x+1 = x²+2x, i.e. x = 1, contradiction. If ε = 1: b' = x < c', and
d_{x+2}(μ) = I_{x+3} − (x+3) − [c'=x+1] by (S2). If c = x < d, then d_{x+2}(λ) = I_{x+3} − 2(x+2) − [d=x+1] and
d_{x+2}(λ) − d_{x+2}(μ) = −(x+1) − [d=x+1] + [c'=x+1] < 0, contradiction. If d = x < c, then
n = c(x+1)+1 = (x+1)(c'+1)−1 gives (x+1)(c'+1−c) = 2, so x+1 divides 2 and x = 1, contradiction. ∎

*Computational verification* (`fathooks2_lemmas.py 40 30`): for every partition λ of n ≤ 30 with r = 2 and at least
two unit parameters, no partition μ ∉ {λ,λ^t} of n has d(μ) = d(λ) (the d-vectors of all partitions of n ≤ 30 were
computed and grouped); and no two two-corner shapes with n ≤ 40 outside a transpose pair share a d-vector.

## 5. The remaining cases: reduction to a Diophantine statement (VERIFIED, not proved)

Theorem 3 shows that the bottom of the vector is *local* (it sees only the parameters ≤ x1+x2−1), so any proof for
shapes with a large parameter must use global invariants. The PROVED-determined global invariants are f^λ
(Briefing fact 2), C_2 (Theorem 4.6), C_1² and C_4 (Corollary 4.10.2). Closed forms (`fathooks2_c2.py`, checked by
direct summation for n ≤ 60): n = ab+bc+cd,
 2C_1 = a²b − ab² + 2abc − b²c + bc² − 2bcd + c²d − cd²  (antisymmetric under reversal),
 6C_2 = 2a³b − 3a²b² + 6a²bc + 2ab³ − 6ab²c + 6abc² − ab + 2b³c − 3b²c² + 6b²cd + 2bc³ − 6bc²d + 6bcd² − bc
        + 2c³d − 3c²d² + 2cd³ − cd  (symmetric under reversal),
and f^λ = n!/H with H = Π_{i<b,j<a}(a+b−1−i−j) · Π_{i<d,j<c}(c+d−1−i−j) · Π_{i<b,j<c}(a+b+c+d−1−i−j).

`fathooks2_invariants.py 300` (604 513 two-corner shapes, n ≤ 300, closed forms only) gives, for the number of
colliding keys among transpose classes (first collision):
 (n,H): 57 (n = 7: λ(1,1,3,1) vs λ(1,2,1,3));   (n,C_2): 17 217 (n = 9);   **(n,H,C_2): 0**;
 (n,C_2,C_1²,C_4): 10 (n = 50: λ(4,5,3,5) vs λ(5,1,5,8));   (n,H,C_2,C_1²,C_4): 0;
 (n,LOC): 14 740 (n = 5);  (n,LOC,H): 33 (n = 14);  (n,LOC,C_2): 1274 (n = 12);  (n,LOC,H,C_2): 0;
 (n,multiset): 805 (n = 13);  (n,multiset,H): 0;  (n,multiset,C_2): 1 (n = 178: (6,7,8,10) vs (7,10,6,8)).
Here LOC = (x1, x2, (N_k)_{k≤x1+x2−1}) is the local data of Theorem 3.

**Conjecture D (VERIFIED for n ≤ 300).** (n, f^λ, C_2) determines a two-corner shape up to transpose. Since all three
are determined by the d-vector (PROVED), Conjecture D implies Conjecture B for all two-corner shapes. Equivalently
(given n): the hook product H(a,b,c,d) and the cubic polynomial C_2(a,b,c,d) together determine (a,b,c,d) up to
reversal. Note that neither invariant alone suffices, and that the "arrangement" version with (n, C_2) in place of
(n, H, C_2) is false (the n = 178 pair above; `fathooks2_arrangement.py` reproduces it), while (n, multiset, H) is
collision-free up to n ≤ 300.

I did not find a proof of Conjecture D; the natural route (compare the largest hook lengths and the 2-adic or
prime structure of H, then C_2) is Diophantine and was not completed in the time available. Everything in §§1–4 is
independent of it.

## 6. Verification log

`python3 src/agents/fathooks2_lemmas.py 40 30` (log `src/agents/fathooks2_lemmas.log`):
```
[W1,W1',W2,W2',READ: n<=40; BOX: n<=24] 3917 shapes, failures = 0
Conjecture B among 2-corner shapes, n <= 40: non-transpose collisions = 0
[S1-S4: all parameter values with n <= 60] 190 (c,d) pairs, failures = 0
[Theorem F vs all partitions, n <= 30] 861 shapes, failures = 0
```
(W1 = Lemma 1(i), W1' = Lemma 1(ii), W2 = Lemma 2(ii), W2' = Lemma 2(iii), BOX = Lemma 2(i), READ = the algorithm of
Theorem 3; S1–S4 as in §3, with the F_1 formula for d = x ≤ c tested for c,d ≥ 2 and the one for c = x < d for all
c,d; "Theorem F vs all partitions" groups the d-vectors of all partitions of each n ≤ 30.)

`python3 src/agents/fathooks2_invariants.py 300` (log `src/agents/fathooks2_invariants300.log`): the collision
counts quoted in §5 (604 513 shapes). `python3 src/agents/fathooks2_c2.py`: C_2 formula, 11 289 shapes with n ≤ 60,
0 failures, and the symmetry/antisymmetry of C_2/C_1 under reversal checked symbolically.
`python3 src/agents/fathooks2_arrangement.py 14`: over all 2380 multisets with values ≤ 14, the 12 arrangements are
separated by (n,H) (0 collisions) but not by (n,C_2) (exactly one collision: (178, 9385) for (6,7,8,10) vs
(7,10,6,8)).

## 7. Assessment

Complete, airtight: the two exact windows (Lemmas 1–2), the self-certifying local reading (Theorem 3: the d-vector
determines x1, x2 and N_k for k ≤ x1+x2−1, hence the full multiset for balanced shapes), and Theorem F, which
extends Theorem 4.9 from class C to all two-corner shapes with at least two unit parameters, against arbitrary μ.
Not achieved: the general two-corner case. The obstruction is structural — the bottom of the vector is local and
cannot see a parameter ≥ x1+x2, so the general case needs the global invariants f, C_2 (and possibly C_1², C_4);
the data strongly suggest that (n, f, C_2) already suffices (Conjecture D, verified to n = 300), but I have no proof,
and the same-multiset version with (n, C_2) alone is false at n = 178, so the argument must genuinely use f.
