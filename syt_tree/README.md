# The tree of standard Young tableaux: automorphisms and census (Worley, Problems 2–3)

**What was accomplished.** Worley's Problem 3 concerns the infinite tree SYT whose vertices are all standard Young tableaux, joined by adding one box. It makes two conjectures: (A) every automorphism of the tree is a composition of "partial transposes" (transpose everything below a vertex of symmetric shape), and (B) the census of the subtree below a vertex — the number of vertices at each depth — determines the shape of the vertex up to transpose. We first proved (and found in the literature: Stanley 2003) that the infinite census is equivalent to a finite vector, the *down-census* (d_0, …, d_n) counting the ways to peel j boxes off the shape, via the identity Σ s_k t^k/k! = exp(t + t²/2) Σ d_j t^j/j!. Using this we verified Conjecture B for **every shape with at most 75 boxes** (8.1 million partitions at n = 75; exact arithmetic to 60, modulo two 61/62-bit primes to 75 — a rigorous negative), far beyond the unspecified range of Worley's computation. For Conjecture A we found that **the conjecture as literally stated is false**: at a vertex of shape (3,2,1) one can transpose the subtrees of just one mirror pair of children, and a parity invariant shows this automorphism is not even a limit of compositions of Worley's whole-subtree transposes (confirmed by an exact group computation). The correct generators are these *pair transposes*, and we proved that the corrected Conjecture A is *equivalent* to a weak form of Conjecture B ("sibling subtrees are non-isomorphic unless mirror images"); hence, unconditionally, every automorphism agrees on the first 75 levels with a finite product of pair transposes, and under the weak form the automorphism group is a pro-2 group of continuum cardinality in which every element is a unique convergent product of pair transposes.

**What remains open.** Conjecture B itself (equivalently: the down-census vector determines the shape up to transpose) is proved here only for special classes — rectangles, and the classes listed in `results.md` §4.9 — and we proved that several statistics are recoverable from the vector (number of corners, the smallest corner box, a "local product" of rectangle polynomials, the number of tableaux f^λ, the sum of squared contents). We also found families of non-conjugate shapes whose down-census vectors agree in all but the last three coordinates, which shows that any proof must use the top of the vector (essentially f^λ and character values at 3-cycles), not just local corner structure. A Pfaffian formula for the generating polynomial and a reformulation via shifted Schur functions are recorded as possible starting points. The shifted analogue (Problem 2) was studied computationally; see `results.md` §5. Every claim in `results.md` is labelled PROVED, COMPUTATIONALLY VERIFIED (with range and code), or CONJECTURAL, and the last section gives a candid assessment of what a human referee should check.

## Layout

* `results.md` — the mathematical write-up (statements, complete proofs, computational ranges, open questions, assessment).
* `NOTES.md` — running log including dead ends.
* `BRIEFING.md` — the summary of established facts handed to the parallel research agents.
* `src/` — Python (exact integer arithmetic) and one C++ program:
  `young.py` (partitions, corners, hook lengths, f^λ two ways, skew chain counts), `census.py` (up/down census, the identity), `tree.py` (the tree, τ and θ automorphisms, parity invariants), `dcensus_scan.py` (exact level-by-level scan), `dcensus_modp.cpp` (mod-prime scan to n = 75), `automorphism_check.py` (sympy group computations), `verify_identity.py`, `check_formulas.py`, `check_characters.py`, `check_schur_pfaffian.py`, `check_classC_cases.py`, `crosscheck_skew.py`, `independent_check.py` (from-scratch reimplementation of the scan), `separating_depth.py`, `explore_*.py`, `shifted*.py` (Problem 2), and `agents/` (scripts written by the research agents).
* `tests/` — pytest suites with hand-verified values (run `python3 -m pytest tests`).
* `data/` — d-vectors for n ≤ 30 (gzip), per-n collision reports and digests (`collisions_python.txt`, `scan_modp_*.txt`), identity check output, literature search synthesis.
* `agent_notes/` — reports of the background research agents (referee-checked claims are incorporated into `results.md`).

## Reproducing

```
cd syt_tree
python3 -m pytest tests -q
python3 src/verify_identity.py 10 10          # identity check
python3 src/dcensus_scan.py 40 30             # exact scan to n=40 (seconds); 60 takes ~6 min
g++ -O2 -std=c++17 -o src/dcensus_modp src/dcensus_modp.cpp && src/dcensus_modp 60 0
python3 src/automorphism_check.py 5 6 7      # group computations (sympy)
```
