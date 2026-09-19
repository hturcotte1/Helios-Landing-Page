# Skeptic angle, second independent pass: trying to disprove Conjecture B

Author: skeptic agent (second pass). `agent_notes/skeptic.md` (first pass) already exists and, by the rule "do not modify
existing files", this report is written as `skeptic2.md`. All scripts are `src/agents/skeptic2_*.py`, logs `src/agents/skeptic2_*.log`
(`skeptic2_window.json` holds the raw near-collision data). Exact integer arithmetic throughout. No code is shared with
`census.py`, `young.py` or the first-pass `skeptic_*.py` (only my own `skeptic2_dvec.py` / `skeptic2_targeted.py` are imported by
the other `skeptic2_*` scripts). Every claim is labelled PROVED / VERIFIED (range + script) / CONJECTURAL.

Notation as in BRIEFING/results.md: `d(lam) = (d_0,...,d_n)`; for a pair of transpose classes `a` = first differing index,
`b = n - (last differing index)`, `w = n + 1 - a - b` = length of the disagreement window; `C_k` = k-th content power sum;
`(c^k)` = k rows of length c.

## 0. Summary

1. **Independent confirmation (VERIFIED).** From-scratch d-vector code (`skeptic2_dvec.py`): zero collisions for all n <= 42, transpose
   invariance checked for every shape, and the tuple (p(n), #distinct vectors, #symmetric, #transpose pairs, #collision groups,
   digest sum d_j 7^j mod 2^61-1) equals the corresponding line of `data/collisions_python.txt` for every n = 0..42
   (`skeptic2_dvec.log` n <= 30, `skeptic2_dvec42.log` n <= 42; 73 lines, all `MATCH`). The recursion was cross-checked against a
   definitionally different computation (sum over nu subset lam of f^{lam/nu} by Aitken's determinant) on all partitions of n <= 12
   and 20 random partitions of 20.
2. **Closest pairs / windows (VERIFIED, n <= 42).** Over all pairs the minimal window is w = 4 for every 8 <= n <= 42, attained at
   the top end by (2,2,1^{n-4}) vs (3,1^{n-3}) (the first pass showed by a Hamming-distance census that this is the only pair at distance 4). Restricted to pairs with equal f the minimal window grows linearly:
   wmin_f = 18, 22, 21, 20, 25, 24, 22, 26, 26, 24, 27 for n = 32..42 (roughly 2n/3, attained by the families H_k, H'_k, H''_k of
   the first pass). No pair with n <= 42 agrees on the top six entries; pairs agreeing on the top five have prefix agreement a <= 4;
   pairs agreeing on the top four (f and C_2) have a <= (n+1)/3 (H_k). **No family with a bounded disagreement window away from the
   top end exists for n <= 42, and the data show the window of equal-f pairs growing linearly.**
3. **New PROVED result (Theorem S3 and Corollary C).** The 3-row chain symmetry N_j(u,v) = N_j(v,u), left CONJECTURAL in the first pass,
   is proved here: for GL_3, the number of irreducible constituents of V_lambda (x) W equals that of V_lambda (x) W^* for every
   finite-dimensional representation W. The proof reduces, via the Racah-Speiser form of Weyl's character formula, to an identity
   S3 about signs of sorted 3-vectors, which is piecewise constant on the cells of a central hyperplane arrangement with
   {0,+-1}-coefficients and is therefore proved by a finite check (`skeptic2_chain.py`). Consequences: the identities R_K
   (d_j((3^K,2,2,2)) = d_j((3^K,1,1,1)) for j <= K) hold for all K, so Theorem H(iii) of `skeptic.md` (bottom agreement of the family
   (5,3^k,2,2,2) vs (4,4,3^k,1,1,1) for j <= k) is now unconditional given Lemma B there. So the family H_k is a PROVED infinite family
   agreeing on d_0..d_k and on d_{n-3}..d_n with f and C_2 equal; its disagreement window has length (2n-14)/3 (VERIFIED k <= 12 here).
4. **Targeted search (VERIFIED, n <= 52).** Pairs of transpose classes with the same multiset of corner boxes {a_i,b_i}, the same f
   and the same (C_2, C_1^2, C_4): none for any n <= 52 (`skeptic2_targeted.py`); the counts of the weaker filters agree exactly with
   the first pass's independent code at n = 30, 40, 50, confirming both. Exact d-vectors were computed for every survivor of the
   weaker filters that reached the K2 stage (there were none), so no collision candidate exists in this range.
5. **Heuristic (Section 5).** With the empirical conditional agreement probabilities measured at n = 30, 40, 42, the expected number of
   collisions is astronomically small *provided* successive top levels behave like independent polynomial invariants; the only
   place where the data show "stickiness" (conditional agreement probabilities rising to 0.75) is the bottom of the vector, and all
   such sticky pairs are structural (attachment/insertion pairs) and provably disagree at the top. I could not disprove B; Section 6
   states what a counterexample would have to look like.

## 1. Independent confirmation

`skeptic2_dvec.py`: own partition generator (decreasing tuples), own corner detection, recursion
`d_j(lam) = sum_{c removable} d_{j-1}(lam - c)`, computed level by level (two levels in memory). Per level: assert
`d(lam) = d(lam^t)` for all lam; group by exact vector; every group must be {lam} (symmetric) or {lam, lam^t}; digest
`sum_{lam ⊢ n} sum_j d_j(lam) 7^j mod (2^61-1)`. Independent cross-check of the recursion: `dvector_by_definition(lam)` enumerates all
nu subset lam and sums `f^{lam/nu}` computed by Aitken's determinant `f^{lam/nu} = |lam/nu|! det[1/(lam_i - nu_j - i + j)!]`
(exact rationals) — agreement on all 272 partitions of n <= 12 and on 20 random partitions of 20.

Result (`skeptic2_dvec42.log`): for every n <= 42, zero collision groups, and all six statistics match `data/collisions_python.txt`
(e.g. n = 42: p(42) = 53174, 26613 distinct vectors, 52 symmetric, 26561 transpose pairs, digest 577640239776089963).
Status: **VERIFIED, n <= 42**, with code independent of everything else in the project.

## 2. Near-collisions: windows, both-ends agreement, n <= 42

`skeptic2_window.py` (one representative per transpose class; exact). For each n it computes
* `wmin_all`: min over all pairs of w (exact: for each prefix length a0, group by d_0..d_{a0-1}, and inside a group the maximal common
  suffix is found by sorting the reversed vectors — the maximizing pair is adjacent in sorted order);
* for all pairs with equal f (grouped by d_n, all pairs inside a group compared): `wmin_f`, and `amax(b>=B)` = largest prefix agreement
  among pairs agreeing on the top B entries, B = 4 (f, C_2), 5 (f, C_2, C_1^2), 6;
* the numbers of pairs agreeing on the bottom a0 / top b0 entries and V_j = number of distinct values of d_j (Section 5).

| n | classes | wmin_all (pair) | #pairs equal f | wmin_f (pair, a, b) | amax(b>=4) | amax(b>=5) | amax(b>=6) |
|---|---|---|---|---|---|---|---|
| 20 | 317 | 4 (F_3) | 44 | 10: (4,4,3,3,3,1^3) vs (5,3,3,3,2,2,2), a=7, b=4 | 7 | - | - |
| 26 | 1224 | 4 | 258 | 14: (4,4,3^5,1^3) vs (5,3^5,2^3) [H_5], a=9, b=4 | 9 | - | - |
| 30 | 2811 | 4 | 470 | 18: (5,5,4^4,2,1,1) vs (6,4^4,3,3,2) [H''_4], a=9, b=4 | 9 | - | - |
| 34 | 6168 | 4 | 1577 | 21: (5,5,4^5,2,1,1) vs (6,4^5,3,3,2) [H''_5], a=10, b=4 | 10 | - | - |
| 36 | 9005 | 4 | 1328 | 25: (5,5,4,4,3^4,1^6) vs (7,5,3^4,2^6), a=8, b=4 | 8 | 2 | - |
| 37 | 10836 | 4 | 2721 | 24: (5,4,4,3^6,2,1^4) vs (6,4,3^6,2^4,1) [H'_6], a=10, b=4 | 10 | 4 | - |
| 38 | 13026 | 4 | 2851 | 22: (4,4,3^9,1^3) vs (5,3^9,2^3) [H_9], a=13, b=4 | 13 | 2 | - |
| 40 | 18692 | 4 | 3404 | 26: (5,4,4,3^7,2,1^4) vs (6,4,3^7,2^4,1) [H'_7], a=11, b=4 | 11 | 3 | - |
| 41 | 22316 | 4 | 6383 | 24: (4,4,3^10,1^3) vs (5,3^10,2^3) [H_10], a=14, b=4 | 14 | 3 | - |
| 42 | 26613 | 4 | 4834 | 27: (5,5,4^7,2,1,1) vs (6,4^7,3,3,2) [H''_7], a=12, b=4 | 12 | 1 | - |

(Full table for 6 <= n <= 42 in `skeptic2_window.log`; "-" means no such pair.) Observations (all VERIFIED, n <= 42):
* `wmin_all = 4` for all 8 <= n <= 42, attained by the pair (2,2,1^{n-4}) vs (3,1^{n-3}) (my script records one attaining pair; uniqueness is the first pass's result) (window = the top four entries; this is
  the pair of results.md Prop. 4.6.1 with m = 2).
* Among pairs with equal f the minimal window grows linearly (wmin_f >= 16 for n >= 29, >= 20 for n >= 33) and is attained, for every n >= 26 in the table,
  by a member of the "body insertion" families H_k, H'_k, H''_k of the first pass, whose windows (differing indices
  {k+4..n-4} for H_k, H'_k and {k+5..n-4} for H''_k, VERIFIED there) have length H_k (n = 3k+11): w = 2k+4 = (2n-14)/3;
  H'_k (n = 3k+19): w = 2k+12 = (2n-2)/3; H''_k (n = 4k+14): w = 3k+6 = (3n-18)/4 — e.g. n = 40: H'_7 gives 26, n = 42: H''_7 gives 27.
* Every pair agreeing on the top five entries (f, C_2, C_1^2) has a <= 4, i.e. differs already at d_1, d_2, d_3 or d_4.
* No pair with n <= 42 agrees on the top six entries (consistent with results.md 4.10: the first such pair is at n = 49).

**Lemma W (PROVED, trivial).** For non-transpose lam, mu ⊢ n >= 6, b ∉ {1, 2, 6}: if the last differing index is >= n-2 it is n
(d_n = d_{n-1} = d_{n-2} = f), and it is never n-6 because d_{n-6} is a function of (n, d_n, d_{n-3}, d_{n-4}, d_{n-5})
(results.md Thm 4.10.3). In particular a pair with equal f has b >= 3, and b = 6 never occurs. (Checked implicitly by the tables:
b takes the values 0, 3, 4, 5 only for n <= 42.)

Conclusion of this section: the most promising route to a counterexample named in the task — a family with the disagreement confined
to a bounded window while both ends agree — does not exist for n <= 42 (the window of every equal-f pair is >= 16 for n >= 29, >= 20 for n >= 33, and
grows linearly along the only families that come close).

## 3. The 3-row chain symmetry is a theorem (PROVED)

The first pass reduced the bottom agreement of the family H_k to the identity
R_K: d_j((3^K,2,2,2)) = d_j((3^K,1,1,1)) for j <= K, and further to a symmetry N_j(u,v) = N_j(v,u) of chain counts in Young's lattice
restricted to three rows, which it could only verify numerically. Here is a proof.

**Definitions.** For nu = (nu_1 >= nu_2 >= nu_3 >= 0) and j >= 0 let N_j(nu) = number of chains nu = kappa^0 ⊂ kappa^1 ⊂ ... ⊂ kappa^j,
|kappa^i / kappa^{i-1}| = 1, all kappa^i partitions with at most 3 rows. Put nu* = (nu_1 - nu_3, nu_1 - nu_2, 0).
For v in Z^3 let sgn(v) = prod_{i<j} sgn(v_i - v_j) in {-1, 0, 1} (0 iff two entries coincide; otherwise the sign of the permutation
that sorts v into strictly decreasing order, since prod_{i<j} sgn(v_i - v_j) = (-1)^{#{i<j: v_i<v_j}}).

**Theorem S3 (PROVED by finite check, `skeptic2_chain.py`).** For every strictly decreasing a in Z^3 and every b in Z^3,
    sum_{sigma in S_3} sgn(a + sigma b) = sum_{sigma in S_3} sgn(a - sigma b),
where (sigma b)_i = b_{sigma(i)}.

*Proof.* Both sides are unchanged when a or b is shifted by a constant vector (all pairwise differences are unchanged), and the
right-hand side is symmetric under permuting b. So we may take a = (A+B, B, 0) with integers A, B >= 1 and b = (p, q, 0) with
p >= q >= 0. Each factor sgn((a +- sigma b)_i - (a +- sigma b)_j) (i < j) is the sign of a linear form in x = (A, B, p, q) with
coefficients in {0, 1, -1}: (a +- sigma b)_i - (a +- sigma b)_j = (A, B or A+B) +- (b_{sigma(i)} - b_{sigma(j)}) and
b_{sigma(i)} - b_{sigma(j)} in {+-p, +-q, +-(p-q)}. Let H be the finite set of these forms together with A, B, p, q, p-q, A+B,
and consider the central hyperplane arrangement {l = 0 : l in H} in R^4. Both sides of S3 are constant on each cell of the
arrangement (a cell = set of points with a prescribed sign vector (sgn l(x))_{l in H}), and our region {A >= 1, B >= 1, p >= q >= 0}
is a union of cells (its defining forms are in H). Since A, B, p, q in H, the only point where all forms vanish is 0, so the closure
of every cell F is a pointed polyhedral cone. Its extreme rays are 1-dimensional faces, i.e. kernels of rank-3 systems of forms from H;
choosing three independent forms, the kernel is spanned by the vector of signed 3x3 minors of a 3x4 matrix with entries in {0, +-1},
whose absolute values are at most 4 (the determinant is multilinear in the rows, so its maximum over entries in [-1,1] is attained at a
+-1 matrix, and the determinant of a 3x3 +-1 matrix is divisible by 4 and bounded by Hadamard's 3^{3/2} < 8, hence is 0 or +-4).
Dividing by the gcd, every extreme ray has a primitive integer direction r with |r_i| <= 4. By Caratheodory's theorem for cones,
any x in F is a positive combination of a linearly independent set R of at most 4 extreme rays of cl(F) (drop rays with coefficient 0).
Then y = sum_{r in R} r lies in F: for l in H with l = 0 on F, l(r) = 0 for r in cl(F); for l > 0 on F, l(r) >= 0 for all r in R and
l(x) = sum c_r l(r) > 0 forces some l(r) > 0, so l(y) > 0; similarly for l < 0. Thus every cell meeting the region contains an integer
point y with 1 <= A, B <= 16 and 0 <= q <= p <= 16 (|y_i| <= 4*4). `skeptic2_chain.py` evaluates both sides of S3 at all 39168 such
points (A, B in [1,16], 0 <= q <= p <= 16) and finds equality at every one of them. Since both sides are constant on cells, S3 holds on
every cell of the region, i.e. for all integer (and real) parameters. □

**Theorem C1 (PROVED).** For GL_3(C), every dominant weight lambda and every finite-dimensional rational representation W: the
number of irreducible constituents of V_lambda (x) W, counted with multiplicity, equals that of V_lambda (x) W^*.

*Proof.* Weyl's character formula: ch V_lambda = A_{lambda+rho} / A_rho with A_gamma = sum_{w in S_3} sgn(w) x^{w gamma},
rho = (2,1,0). For any gamma in Z^3, A_gamma = 0 if gamma has two equal entries, and otherwise A_gamma = sgn(w) A_{w gamma} where w gamma
is strictly decreasing, so A_gamma / A_rho = sgn(gamma) ch V_{(gamma)^+ - rho} with (gamma)^+ the decreasing rearrangement (this is the
Racah-Speiser / Brauer-Klimyk form). Writing ch W = sum_beta m_W(beta) x^beta,
    ch V_lambda * ch W = sum_beta m_W(beta) A_{lambda+rho+beta} / A_rho = sum_beta m_W(beta) sgn(lambda+rho+beta) ch V_{(lambda+rho+beta)^+ - rho},
so the total number of constituents is sum_beta m_W(beta) sgn(lambda+rho+beta). W^* has the weights -beta with the same multiplicities,
so its count is sum_beta m_W(beta) sgn(lambda+rho-beta). The multiplicity function m_W is S_3-invariant; grouping beta into S_3-orbits O
with representative b and stabilizer of order s_b, sum_{beta in O} sgn(a +- beta) = (1/s_b) sum_sigma sgn(a +- sigma b) with
a = lambda+rho strictly decreasing, and Theorem S3 says the two sums agree orbit by orbit. □

**Corollary C2 (PROVED).** N_j(nu) = N_j(nu*) for all nu and j.

*Proof.* By Pieri's rule in three variables, ch V_kappa * ch V = sum ch V_{kappa'} over kappa' = kappa + one box with at most 3 rows
(V = C^3). Iterating, the multiplicity of V_kappa in V_nu (x) V^{(x) j} is the number of chains nu -> kappa of length j through partitions
with <= 3 rows, so N_j(nu) = #const(V_nu (x) V^{(x) j}). By Theorem C1 with W = V^{(x) j}, this equals #const(V_nu (x) (V^*)^{(x) j})
= #const((V_{nu}^* (x) V^{(x) j})^*) = #const(V_nu^* (x) V^{(x) j}) (dualizing permutes the irreducibles). Finally V_nu^* = V_{(-nu_3,-nu_2,-nu_1)}
= V_{nu*} (x) det^{-nu_1}, and tensoring with a power of det permutes the irreducibles, so #const(V_nu^* (x) V^{(x) j}) = N_j(nu*). □
(In the notation of the first pass, (u,v) = (nu_1 - nu_2 + 1, nu_2 - nu_3 + 1) and nu* has (v,u).)
Checked numerically: N_j(nu) = N_j(nu*) for all nu with nu_1 <= 8 and all j <= 20 (3465 cases, `skeptic2_chain.log`).

**Corollary C3 (PROVED): R_K holds for every K >= 1**, i.e. d_j((3^K,2,2,2)) = d_j((3^K,1,1,1)) for 0 <= j <= K.

*Proof.* Let lam = (3^K,2,2,2), n = 3K+6, and j <= K. d_j(lam) = sum_{nu ⊆ lam, |lam/nu| = j} f^{lam/nu} = d_j(lam^t) with
lam^t = (K+3, K+3, K). For nu^t ⊆ lam^t, rotate the 3 x (K+3) box by 180 degrees: lam^t becomes the complement of (3) (rows of the box
minus rows of lam^t: 0, 0, 3, rotated), nu^t becomes the complement of a partition kappa ⊇ (3) with |kappa| = j + 3 and at most 3 rows,
and the skew shape lam^t/nu^t becomes kappa/(3) rotated by 180 degrees, which has the same number of standard fillings
(rotation reverses the order of the poset, and reversing a linear extension of a reversed poset is a bijection). Conversely every
kappa ⊇ (3) with |kappa| = j+3, l(kappa) <= 3 fits in the box because kappa_1 <= 3 + j <= K + 3. Hence
d_j((3^K,2,2,2)) = sum_{kappa ⊇ (3), |kappa| = j+3, l(kappa) <= 3} f^{kappa/(3)} = N_j((3,0,0)). The same argument for
(3^K,1,1,1)^t = (K+3, K, K), whose complement in the box is (3,3), gives d_j((3^K,1,1,1)) = N_j((3,3,0)). Since (3,0,0)* = (3,3,0),
Corollary C2 finishes the proof. □
Numerically (`skeptic2_chain.log`): R_K holds for K <= 22, the reductions d_j = N_j((3)) resp. N_j((3,3)) were checked for j <= K,
and the first differing index of the two vectors is exactly K+2 for every K <= 22 (so R_K is sharp up to one index; sharpness is not
claimed in general).

**Consequence for the family H_k.** The first pass proved (its Lemma B, body insertion) that for lam(k) = (5,3^k,2,2,2),
mu(k) = (4,4,3^k,1,1,1), n = 3k+11, d_j(lam(k)) - d_j(mu(k)) = sum_{i=0}^{2} C(j,i) [d_{j-i}((3^k,2,2,2)) - d_{j-i}((3^k,1,1,1))] for
j <= k; with Corollary C3 the bracket vanishes, so **d_j(lam(k)) = d_j(mu(k)) for all j <= k, for every k >= 1 (PROVED)**, in addition to
f, C_2 equal (PROVED there) and hence d_{n-3},...,d_n equal, and C_1 different by 2, hence d_{n-4} different. I re-verified the whole
statement with my own code: for k <= 12 the set of differing indices is exactly {k+4, ..., n-4} (`skeptic2_chain.log`). So H_k is a
PROVED infinite family with prefix agreement >= k+1 = (n-8)/3, suffix agreement exactly 4, and window (2n-14)/3 (VERIFIED k <= 12 for the
exact window; PROVED that the window is contained in [k+1, n-4]).

## 4. Targeted search: same corner multiset, same f, same (C_2, C_1^2, C_4), n <= 52

`skeptic2_targeted.py` (own corner runs, hook lengths, contents; log `skeptic2_targeted.log`). One representative per transpose class.
Keys: K0 = (n, sorted multiset of unordered pairs {a_i,b_i}) — fixes d_0..d_{2m} by the local product theorem; K0+f; K0+C with
C = (C_2, C_1^2, C_4) — with f fixes d_{n-6}..d_n (results.md 4.4, 4.10); K2 = K0 + f + C; K3 = K2 + (C_6 - 16 C_1 C_3) (the level-7
invariant of results.md Thm 4.10.4; not used for any claim). For every K2 group of size >= 2 the exact d-vectors would be computed and
compared (memoised recursion); no such group occurred.

| n | classes | pairs same K0 | same K0+f | same K0+C | same K2 | collisions |
|---|---|---|---|---|---|---|
| 20 | 317 | 153 | 1 | 0 | 0 | 0 |
| 30 | 2811 | 3799 | 4 | 0 | 0 | 0 |
| 40 | 18692 | 58544 | 21 | 1 | 0 | 0 |
| 45 | 44601 | 210472 | 64 | 3 | 0 | 0 |
| 48 | 73680 | 439170 | 72 | 8 | 0 | 0 |
| 50 | 102162 | 710578 | 111 | 13 | 0 | 0 |
| 51 | 120025 | 901087 | 145 | 20 | 0 | 0 |
| 52 | 140853 | 1143677 | 124 | 25 | 0 | 0 |

**Zero K2 (and K3) survivors for every 8 <= n <= 52 (VERIFIED).** The counts at n = 30, 40, 50 coincide exactly with the first pass's
table (3799/4/0; 58544/21/1; 710578/111/13), an independent confirmation of both implementations. The K0+f and K0+C counts grow
(roughly x5 and x13 from n = 40 to 52 while K0 grows x20); a survivor of K2 would still have to agree on the ~n-10 middle coordinates.

`skeptic2_level7.py` (own code) re-verified the top-end facts used as prefilters on all partitions of n <= 28: d_{n-3} is a function of
(n, f, C_2); d_{n-4}, d_{n-5}, d_{n-6} are functions of (n, f, C_2, C_1^2, C_4); d_{n-7} is a function of these and C_6 - 16 C_1 C_3.
(For n <= 28 the last statement is vacuous in the sense that (n,f,C_2,C_1^2,C_4) already separates everything, so it does not test the
form of the level-7 invariant; that is why K3 is reported but not used.)

## 5. A probabilistic heuristic for the expected number of collisions

*Bookkeeping.* d_0 = 1, d_1 = r, d_n = d_{n-1} = d_{n-2} = f, and d_{n-6} is determined by d_n, d_{n-3}, d_{n-4}, d_{n-5}; so at most
n - 4 coordinates (d_2, ..., d_{n-3} without d_{n-6}) carry information beyond (n, r, f). Structure: d_2..d_{2m} are functions of the
multiset of corner boxes (local product theorem; for the typical shape m = 1, so this pins only d_2); d_{n-3}, d_{n-4}, d_{n-5} are
functions of (n, f, C_2, C_1^2, C_4); level n - i for 7 <= i adds one polynomial content invariant (results.md Thm 4.10.4).

*Empirical numbers (`skeptic2_window.json`, n = 42, 26613 classes, 3.54e8 pairs).* Number of distinct values of d_j among classes:
V_1 = 8, V_2 = 45, V_3 = 209, V_4 = 1126, V_5 = 4459, V_6 = 10280, V_7 = 15392, ..., V_j > 26400 for 23 <= j <= 39, V_40 = 26524,
V_{41} = V_{42} = 23153 (= number of distinct f). Bottom: P(agree d_1) = 0.288, then conditional probabilities of agreeing on d_{a0} given
agreement on d_0..d_{a0-1} for a0 = 2..15: 0.30, 0.20, 0.18, 0.22, 0.34, 0.44, 0.51, 0.56, 0.60, 0.64, 0.67, 0.71, 0.73, 0.76 (the same
numbers within 0.03 at n = 30 and 40). Top: P(same f) = 1.4e-5 (4834 pairs); P(same d_{n-3} | same f) = 81/4834 = 0.017;
P(same d_{n-4} | ...) = 1/81 = 0.012; P(same d_{n-5} | ...) = 0/1.

*(a) The naive birthday estimate is worthless.* Treating the n - 4 informative coordinates as independent, E[#collisions] ≈
C(#classes, 2) / prod_j V_j; with V_j ≈ #classes for the middle coordinates this is ≈ #classes^{-(n-8)} · #classes^2/2 ≈ 0 for every
n >= 12. It is worthless because the coordinates are extremely correlated: the bottom conditional probabilities above are not 1/V_j
(≈ 4e-5) but 0.2-0.76, and they *increase* with depth.

*(b) Why the bottom is sticky, and why that does not produce collisions.* The pairs that survive deep bottom agreement are the attachment
pairs (N)|nu vs (N)|nu^t of the first pass (Lemma A there) and the body-insertion families of Section 3: their agreement on
d_0..d_{~n/3} is forced by a factorisation of the filter poset, not by chance. All attachment pairs have different f (VERIFIED in the
first pass for |nu| <= 12 over large N), and all body-insertion pairs found so far with equal f and C_2 have C_1 differing by a nonzero
constant (PROVED for H_k), hence differ at d_{n-4}. In the data, every pair agreeing on the top five entries has a <= 4 (Section 2).
So the correct model is not "bottom stickiness × top stickiness" but: bottom-sticky pairs are separated at the top with (empirical)
probability 1 - O(1/#pairs).

*(c) Top-down estimate.* Let T_5(n) be the number of pairs agreeing on (n, f, C_2, C_1^2, C_4). Data: T_5(n) = 0 for n <= 48, T_5(49) >= 1,
T_5(84) = 76 (results.md 4.10, agent scan); so T_5 grows, but slowly (sub-polynomially in p(n), which is ≈ 10^{10} at n = 84). Given
agreement on levels n, n-3, n-4, n-5 (hence n-6), each further level n-i (i = 7, 8, ...) is a new integer-valued polynomial invariant of
weighted degree i in the contents. If, conditionally on the previous levels, agreement at each new level had probability at most c < 1
(the measured conditional probabilities at the visible top levels are 0.017, 0.012 and < 1), then P(all top levels agree) <= c^{n/2 - 6},
and E[#collisions at n] <= T_5(n) c^{n/2-6}. With T_5(n) <= p(n)^2 = exp(O(sqrt n)) this sums to a finite total over n and is < 10^{-10}
for every n >= 60 if c <= 0.1. The bottom and middle coordinates can only decrease this.

*(d) What the heuristic is worth.* Its only real input is the assumption that the top levels are "generic" relative to each other, i.e.
that no structural mechanism forces many consecutive top levels to agree. The known mechanisms at the top (mirror box moves, results.md
Prop. 4.10.5: levels n..n-6; Prop. M0: levels n..n-10) have bounded reach; the known mechanisms at the bottom have reach ~n/3 and
provably break at the top. A collision requires simultaneous control of all ~n-4 coordinates by a single identity, as transposition
provides. The heuristic says a counterexample cannot be an accident; it says nothing about whether such an identity exists. Honest
verdict: the heuristic supports Conjecture B only conditionally, and the growth of T_5(n) (and of the K0+f, K0+C counts in Section 4)
shows that no finite set of invariants can replace a proof.

## 6. Assessment

* Disproof: **not achieved**; no candidate counterexample for n <= 52 by the targeted prefilter, none for n <= 42 by exhaustive census
  (independent code, digests matching the project's data for all n <= 42).
* Closest pairs: window 4 at the top end (F_3) for every n <= 42; among equal-f pairs the window is >= 16 for n >= 29 (>= 20 for n >= 33) and grows like
  ~2n/3 along the body-insertion families; no pair with n <= 42 agrees on the top six entries.
* The family H_k is now a PROVED infinite family of pairs agreeing on d_0..d_k and d_{n-3}..d_n (Theorem S3, Corollaries C1-C3 plus
  Lemma B of the first pass). It is not a route to a counterexample: C_1(lam) - C_1(mu) = 2 is a constant, so d_{n-4} always differs,
  and the window grows linearly.
* What a counterexample would need: two non-conjugate shapes with (i) equal f, (ii) equal C_2, C_1^2, C_4 and equal values of one new
  content polynomial per level down to level ~n/2, (iii) the same multiset of corner boxes or at least equal d_2..d_{2m+1}, and (iv)
  equal middle coordinates. Items (i)-(iii) are satisfied by nobody up to n = 52 (Section 4), and every mechanism known for (i)-(iii)
  separately has bounded reach. I therefore regard a counterexample as unlikely but not excluded: the middle coordinates are not
  understood, and the numbers of pairs satisfying any fixed finite subset of (i)-(iii) grow without bound.

Open items: (1) a direct combinatorial proof of Theorem C1 (the finite-check proof of S3 is rigorous but opaque; the statement
"#const(V_lambda (x) W) = #const(V_lambda (x) W^*)" was also VERIFIED for GL_4 and GL_5 orbit-wise on small boxes in `skeptic2_signs.py`
— n = 4: a_1 <= 7, b_1 <= 6, 5880 cases; n = 5: a_1 <= 6, b_1 <= 4, 1470 cases; 0 failures — and is CONJECTURAL for GL_n, n >= 4);
(2) whether any body-insertion family can have C_1(lam) = -C_1(mu) (mirror type), which would push the agreement to d_{n-4}; the data
to n = 42 say no; (3) pushing the K2 prefilter beyond n ~ 60 needs C++.

## 7. Scripts (all in `src/agents/`)

* `skeptic2_dvec.py` — independent d-vector code, Aitken cross-check, census with digest comparison (`skeptic2_dvec.log` n <= 30,
  `skeptic2_dvec42.log` n <= 42).
* `skeptic2_window.py` — windows / both-ends agreement / conditional agreement counts, n <= 42 (`skeptic2_window.log`, `skeptic2_window.json`).
* `skeptic2_targeted.py` — corner-multiset / f / content prefilters with exact d-vectors for survivors, n <= 52 (`skeptic2_targeted.log`).
* `skeptic2_level7.py` — functional dependences at the top, n <= 28 (`skeptic2_level7.log`).
* `skeptic2_signs.py` — sign identity S(a,b) tested for n = 2..5 on boxes (`skeptic2_signs.log`).
* `skeptic2_chain.py` — finite-check proof of S3 (box 16), chain symmetry N_j(nu) = N_j(nu*), R_K for K <= 22, family H_k for k <= 12
  (`skeptic2_chain.log`).
