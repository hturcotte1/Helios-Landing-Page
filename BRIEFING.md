# Briefing for research agents — Conjecture B (census determines shape up to transpose)

All code is in `src/` (run scripts from that directory or add it to `sys.path`).
Exact integer arithmetic only. Do not modify existing files; create new files with your own name prefix.

## Definitions
* Partition λ ⊢ n, tuple of decreasing positive ints. λ^t = conjugate. Corner-run parametrisation:
  `corner_runs(λ) = [(a_1,b_1),...,(a_r,b_r)]` (distinct part sizes α_1>…>α_r, a_i = α_i − α_{i+1}, b_i = multiplicity of α_i);
  r = number of removable corners; transposition reverses the list and swaps entries.
* f^λ = number of SYT (`f_hook`), f^{λ/ν} = number of saturated chains ν → λ in Young's lattice (`f_skew(lam, nu)`).
* **Down-census** d_j(λ) = Σ_{ν ⊆ λ, |λ/ν| = j} f^{λ/ν} = number of ways to remove j boxes one at a time
  (`d_vector(lam)` returns (d_0,…,d_n)). Recursion d_j(λ) = Σ_c d_{j−1}(λ − c). P_λ(t) := Σ_j d_j t^j/j!.
  Equivalently d_j(λ) = Σ_{S} e(S) over order filters S ⊆ λ (subsets closed under moving right/down inside λ)
  with |S| = j, e(S) = number of linear extensions of S.
* **Up-census** s_k(λ) = Σ_{μ ⊇ λ, |μ/λ| = k} f^{μ/λ}. PROVED: Σ_k s_k t^k/k! = exp(t + t²/2) P_λ(t).
  Hence **Conjecture B** ⟺ **the vector (d_0(λ),…,d_n(λ)) determines λ up to transpose.**
* u_i(λ) := d_{n−i}(λ) = Σ_{ν ⊢ i} f^{λ/ν}. Facts: u_0 = u_1 = u_2 = f^λ; u_3 = f^λ − f^{λ/(2,1)};
  u_4 = f^λ − 2 f^{λ/(2,1)} + f^{λ/(2,2)}.

## Established facts (proved, and computationally verified)
1. COMPUTATION: B holds for all n ≤ 60 exactly (and mod a 61-bit prime to n ≥ 60, run in progress to 75): no two
   non-transpose partitions of the same n ≤ 60 share a d-vector.
2. d_1 = r; d_n = d_{n−1} = d_{n−2} = f^λ.
3. LOCAL PRODUCT THEOREM (proved). Let m = min_i min(a_i, b_i). For j ≤ 2m,
   d_j(λ) = j! [t^j] Π_{i=1}^r P_{(a_i^{b_i})}(t), where for a rectangle a^b, d_j(a^b) = Σ_{ν ⊢ j, ν ⊆ b×a box} f^ν.
   In particular d_j(λ) = j! [t^j] exp(r(t + t²/2)) for j ≤ m, and
   d_{m+1}(λ) = (m+1)! [t^{m+1}] exp(r(t+t²/2)) − (#{i : a_i = m} + #{i : b_i = m}).
   So the d-vector determines n, r, m, A_m + B_m, and Π_i P_{a_i^{b_i}} mod t^{2m+1} (which depends only on the
   multiset of unordered pairs {a_i,b_i}).
4. RECTANGLES (r = 1): B holds (m and n determine {a,b}).
5. C_2 := Σ_{boxes} (col − row)² is determined: 3 d_{n−3} − 2 f^λ = χ^λ(3-cycle) = f^λ (3C_2 − 3n(n−1)/2)/(n(n−1)(n−2)).
   (Jucys–Murphy: Σ_k J_k² = C(n,2) + (sum of all 3-cycles).)
   Also χ^λ((12)(34)) is determined (from u_3, u_4). In general u_i = Σ_{ρ ⊢ i} (σ(ρ)/z_ρ) χ^λ(ρ ∪ 1^{n−i}) with
   σ(ρ) = number of square roots of a permutation of cycle type ρ; only even-type ρ (each even part with even
   multiplicity) contribute. Character values at odd permutations are NOT determined (they change sign under transpose).
6. PFAFFIAN FORMULA (verified numerically, proof via Schur's Pfaffian identity): for even ℓ ≥ ℓ(λ), α = λ + (ℓ−1,…,0),
   P_λ(t) = Pf( B ), B_{ij} = [x^{α_i} y^{α_j}] (x − y) e^{t(x+y)} / ((1−xy)(1−x)(1−y)).
7. Symmetric functions: d_j(λ) = ⟨ s_λ , p_1^j F ⟩ with F = Σ_μ s_μ; Σ_λ P_λ(t) s_λ = e^{t p_1} F.
   Reversed polynomial: n! P_λ(t)/f^λ = Σ_i g_i(λ) t^{n−i}, where g_i(λ) = (n)_i u_i/f^λ = Σ_{ν⊢i} s*_ν(λ)
   (Okounkov–Olshanski shifted Schur functions, s*_ν(λ) = (n)_{|ν|} f^{λ/ν}/f^λ).
8. Hard pairs: (2,2,1^{n−4}) and (3,1^{n−3}) have identical d_j for all j ≤ n−4 and differ only at j ≥ n−3
   (different f^λ). So the "top" of the vector (f^λ, u_3, u_4, …) is essential.
9. NOT sufficient: (f^λ, C_2) alone has collisions from n=14 on; (r, f^λ) has many.

## Useful API
```python
import sys; sys.path.insert(0, 'src')
from young import partitions, conjugate, corner_runs, from_corner_runs, f_hook, f_skew, boxes, removable_corners
from census import d_vector, s_up
```
`partitions(n)` yields all partitions of n; `d_vector(lam)` is memoised (fine up to n≈40 for individual shapes).
