# The tree of standard Young tableaux: automorphisms and the rank census

Code, data and write-up for a study of Problems 2 and 3 of D. R. Worley's notebook of open problems on the combinatorics of tableaux (arXiv:2509.25446). Problem 3 concerns the tree SYT of all standard Young tableaux ordered by containment: Conjecture A says that every automorphism of the tree is a composition of "partial transposes" at vertices of symmetric shape, and Conjecture B says that the rank census of the subtree below a vertex determines the shape of the vertex up to transpose. Problem 2 is the shifted analogue.

Prepared by an AI research session (Claude) directed and reviewed by [NAMES]. Every claim is labelled PROVED, COMPUTATIONALLY VERIFIED (with range and code), or CONJECTURAL.

## Main results

* Conjecture A is false as stated: at a vertex of shape (3,2,1) the subtrees of a single mirror pair of children can be transposed, and a parity invariant shows that this automorphism is not a limit of compositions of the whole-subtree transposes. The corrected generators are these pair transposes; the corrected statement is equivalent to a weak form of Conjecture B ("non-mirror siblings have non-isomorphic subtrees"), which holds to depth 75.
* The infinite census is equivalent to a finite vector (d_0, …, d_n) counting the ways to peel j boxes off the shape, by an identity that is Theorem 2.2 of R. P. Stanley, *On the enumeration of skew Young tableaux* (2003). Conjecture B is verified for every shape with at most 75 boxes and proved for rectangles, for all shapes avoiding (3,2,1), for two-corner shapes with two unit parameters, and for pairs of shapes each having at most four rows.
* The statistics that the vector provably encodes are determined, and proved families show that neither end of the vector suffices by itself.
* Shifted tree (Problem 2): both halves verified to depth 60; there is no shifted analogue of the generating-function identity with a universal factor.

Conjecture B in general, and its weak form, remain open.

## Where to read

* `for_worley/manuscript.pdf` (source `for_worley/manuscript.md`) — the manuscript: statements, complete proofs, computational ranges, open questions, and an assessment listing what a referee should check first.
* `results.md` — the working write-up from which the manuscript was generated (`for_worley/build_manuscript.py`); the same mathematics, with notes on how each part was obtained and checked.
* `agent_notes/` — supplementary reports with the full proofs and data behind the results summarised in Sections 4.10–4.13 and 5 of the manuscript.
* `NOTES.md` — running log of the work, including dead ends. `BRIEFING.md` — the internal summary of established facts used during the work.

## Code and data

`src/` — Python with exact integer arithmetic (sympy for the group computations) and one C++ program:

* `young.py` (partitions, corners, hook lengths, f^λ two ways, skew chain counts), `census.py` (up- and down-census, the identity), `tree.py` (the truncated tree, the τ and θ automorphisms, parity invariants).
* `dcensus_scan.py` (exact level-by-level scan of Conjecture B), `dcensus_modp.cpp` (the same scan modulo a 61- or 62-bit prime, used to n = 75), `independent_check.py` (from-scratch reimplementation of the scan).
* `automorphism_check.py` (sympy group computations), `verify_identity.py`, `check_formulas.py`, `check_characters.py`, `check_schur_pfaffian.py`, `check_classC_cases.py`, `check_topend.py`, `check_agents2.py`, `crosscheck_skew.py`, `separating_depth.py`, `explore_*.py` (checks and explorations behind individual theorems), `shifted.py`, `shifted_scan.py` (Problem 2), `render_manuscript.py` (renders `results.md` to HTML).
* `agents/` — scripts and logs of the separate computations cited in the manuscript.

`tests/` — pytest suites with hand-verified values.

`data/` — d-vectors for n ≤ 30 (gzip), per-n collision reports and digests of the scans (`collisions_python.txt`, `scan_modp_*.txt`, `independent_check_n40.txt`), the identity check output, the shifted census data (`shifted_census_N40_K14.txt`, `shifted_census_N60_K21.txt`, `shifted_siblings_N39_K14.txt`), and the literature-search synthesis.

## Running the tests and scans

```
python3 -m pytest tests -q
python3 src/verify_identity.py 10 10          # identity check
python3 src/dcensus_scan.py 40 30             # exact scan to n=40 (seconds); 60 takes ~6 min
g++ -O2 -std=c++17 -o src/dcensus_modp src/dcensus_modp.cpp && src/dcensus_modp 60 0
python3 src/automorphism_check.py 5 6 7      # group computations (sympy)
```

Requirements: Python 3 with `sympy` and `pytest`; a C++17 compiler for the mod-p scan. `dcensus_scan.py N D` scans all n ≤ N exactly and writes the d-vectors for n ≤ D; `dcensus_modp N p` scans to n = N modulo the prime selected by p (0 for 2^61 − 1, 1 for 2^62 − 57); n = 75 takes a few minutes and about 9 GB of memory. Each scan reports, for every n, the number of distinct vectors and a digest of the whole set, and lists every group of partitions with equal vectors other than {λ} (λ symmetric) or {λ, λ^t}: the Python scan writes them to `data/collisions_python.txt`, the C++ scan prints `COLLISION` lines. No such group occurs for n ≤ 75. The shifted-tree computations are `python3 src/shifted_scan.py census 40 14` (truncated census of all strict partitions of n ≤ 40 to depth 14; `census 60 21` for the n ≤ 60 run), `python3 src/shifted_scan.py identity` and `python3 src/shifted_scan.py fomin` (see the docstring).

## Test suite status

`python3 -m pytest tests -q` (28 tests, run from the repository root): `28 passed in 0.65s`
