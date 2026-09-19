"""Build for_worley/manuscript.md from ../results.md.

The mathematics is frozen.  This script applies (1) an explicit, asserted list of textual
substitutions that remove internal artifacts and neutralise provenance wording, (2) mechanical
formatting: heading identifiers and hyperlinked section cross-references, (3) new front matter
(title, provenance block, abstract, guide to the reader) and a references list.  Every substitution
is recorded in REPHRASINGS (printed with --list) so that the diff against results.md is auditable.
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "results.md")
OUT = os.path.join(HERE, "manuscript.md")

src = open(SRC, encoding="utf-8").read()

# ----------------------------------------------------------------------------------------------
# (1) explicit substitutions: (old, new).  Each `old` must occur exactly once.
# ----------------------------------------------------------------------------------------------
REPHRASINGS = [
 # §1 pointer to an internal file
 ("**The conjectures** (Worley, Problem 3, verbatim in `data/literature_search.md`):",
  "**The conjectures** (our formulation of the two halves of Worley's Problem 3):"),
 # §0 item 5
 ("(Section 4.10, from a background agent's work that was re-verified here)",
  "(Section 4.10; obtained in a separate computation and independently re-verified)"),
 # §4.9 Remark (ii)
 ("(ii) An independent proof of the same theorem, produced by a background research agent and checked by two adversarial referees, is in `agent_notes/classC.md` (scripts `src/agents/classC_verify.py`, `classC_extra.py`); it proves",
  "(ii) An independent proof of the same theorem was obtained separately and checked by two independent reviewers (scripts `src/agents/classC_verify.py`, `src/agents/classC_extra.py`); it proves"),
 # §4.10 heading and intro
 ("### 4.10 What the top of the vector determines (top-end angle, re-verified)",
  "### 4.10 What the top of the vector determines (re-verified)"),
 ("The background \"top-end\" agent (`agent_notes/topend.md`, scripts `src/agents/topend_*.py`) extended Theorem 4.6 systematically. Its proofs are elementary — explicit class expansions in $\\mathbb Q[S_n]$ of power sums of the Jucys–Murphy elements, plus the fact that these act by contents — and everything below that is marked PROVED has a complete proof in that report which I read line by line; the formulas were **re-verified with my own independent code** in `src/check_topend.py`",
  "A separate computation (scripts `src/agents/topend_*.py`) extended Theorem 4.6 systematically. Its proofs are elementary — explicit class expansions in $\\mathbb Q[S_n]$ of power sums of the Jucys–Murphy elements, plus the fact that these act by contents — and everything below that is marked PROVED has a complete proof that was checked line by line; the formulas were independently re-verified with separately written code in `src/check_topend.py`"),
 # §4.10 Theorem 4.10.3 parenthesis
 ("(The polynomial as printed in the agent report omits the pure-$n$ term; I recovered it by interpolation and verified the corrected identity on all partitions of $6\\le n\\le22$, and the functional dependence on all partitions of $n\\le30$.)",
  "(The polynomial as originally computed omitted the pure-$n$ term; it was recovered by interpolation, and the corrected identity was verified on all partitions of $6\\le n\\le22$, the functional dependence on all partitions of $n\\le30$.)"),
 # §4.10 Theorem 4.10.4 header and provenance sentence
 ("**Theorem 4.10.4 (polynomiality; PROVED in the report, proof read and found complete).**",
  "**Theorem 4.10.4 (polynomiality; PROVED — proof obtained separately and checked line by line).**"),
 ("The report also determines (by full-rank certificates on $n\\le16$ and exact solving, re-verified by the agent for $n\\le20$; not independently recomputed here)",
  "The same computation also determines (by full-rank certificates on $n\\le16$ and exact solving, re-verified there for $n\\le20$; not independently recomputed)"),
 # §4.10 referee paragraph
 ("**Referee verdicts.** Six independent adversarial referee agents (two each) checked Theorem 4.10.1 with its group-algebra identities, Theorem 4.10.4 with the rank-certificate argument, and the appendix results (closed forms for hooks and two-row shapes, sharpness of $F_a$, Kerov's transition formula), re-deriving every formula with their own code (identities in $\\mathbb Z[S_n]$ for $n\\le9$, Theorem 4.10.1 for all $\\lambda\\vdash n\\le24$, etc.). All six verdicts were HOLDS.",
  "**Independent review.** Six independent adversarial re-derivations (two each) checked Theorem 4.10.1 with its group-algebra identities, Theorem 4.10.4 with the rank-certificate argument, and further results (closed forms for hooks and two-row shapes, sharpness of $F_a$, Kerov's transition formula), re-deriving every formula with separately written code (identities in $\\mathbb Z[S_n]$ for $n\\le9$, Theorem 4.10.1 for all $\\lambda\\vdash n\\le24$, etc.). All six confirmed the statements they checked."),
 # §4.10 third-pass paragraph
 ("**Third pass of the top-end angle (appendix of `agent_notes/topend.md`; key data re-verified here).**",
  "**Further results of the same computation (key data independently re-verified).**"),
 ("proved there by the two formulas and in the appendix again via Aitken's determinant",
  "proved there by the two formulas and again, separately, via Aitken's determinant"),
 ("(PROVED, elementary hook-length proof in the appendix; verified for all $\\nu\\vdash m\\le26$)",
  "(PROVED, elementary hook-length proof; verified for all $\\nu\\vdash m\\le26$)"),
 ("(VERIFIED by the agent with a C++ sieve plus exact re-checks; the key pair re-verified here from hook lengths and contents)",
  "(VERIFIED with a C++ sieve plus exact re-checks; the key pair independently re-verified from hook lengths and contents)"),
 (" (e) The printed inversion formula for $C_2$ in the agent's Corollary B′ has a slip ($3f$, not $6f$, in the denominator); Theorem 4.6 above is the correct form.",
  ""),
 # §0 item 5 (second sentence), §4.10 separation data, §4.10 Lemma U paragraph
 ("Additional results from the referee-checked background work are collected in Section 4.10.",
  "Additional results, obtained in separate computations and independently checked, are collected in Section 4.10."),
 ("(agent data $n\\le58$; the $n\\le49$ statement re-verified here)",
  "(separate computation to $n\\le58$; the $n\\le49$ statement independently re-verified)"),
 ("The agent's exhaustive check for $n\\le40$ shows that",
  "An exhaustive check for $n\\le40$ shows that"),
 # §4.11 heading and intro, proof-sketch pointer
 ("### 4.11 Two-corner shapes (fat-hook angle, re-verified)",
  "### 4.11 Two-corner shapes (re-verified)"),
 ("The background \"fat-hook\" agent (`agent_notes/fathooks.md`) proved the following; its Lemma 1 was checked by two referees (HOLDS), and I re-verified the closed forms and Theorem 4.11.2 independently (`src/check_agents2.py`).",
  "The following was obtained in a separate computation; Lemma 4.11.1 was checked by two independent reviewers, and the closed forms and Theorem 4.11.2 were independently re-verified (`src/check_agents2.py`)."),
 ("*Proof sketch (full proof in the report).*",
  "*Proof sketch (the full case analysis is included with the code).*"),
 # §4.12 heading and intro
 ("### 4.12 Structural results (structure angle, re-verified)",
  "### 4.12 Structural results (re-verified)"),
 ("The \"structure\" agent (`agent_notes/structure.md`) found the following; the character-sum formula is proved here directly from Lemma 4.5, and the two row-bounded theorems were checked line by line and their sign thresholds re-verified (`src/check_agents2.py`).",
  "The following results were obtained in a separate computation; the character-sum formula is proved here directly from Lemma 4.5, and the two row-bounded theorems were checked line by line and their sign thresholds independently re-verified (`src/check_agents2.py`)."),
 ("(the agent stated $n\\ge21$; its referees and my check give $17$)",
  "(the original derivation stated $n\\ge21$; the exhaustive check gives $17$)"),
 ("**Theorem 4.12.4 (Pfaffian expansion; PROVED by the agent, three-row corollary re-verified here for all three-row shapes with $n\\le18$).**",
  "**Theorem 4.12.4 (Pfaffian expansion; PROVED separately, three-row corollary independently re-verified for all three-row shapes with $n\\le18$).**"),
 ("(A negative result from the same angle, VERIFIED for $n\\le10$:",
  "(A negative result from the same computation, VERIFIED for $n\\le10$:"),
 # §4.13 heading and intro
 ("### 4.13 The skeptic angle: independent confirmation, a both-ends family, and a heuristic",
  "### 4.13 Attempts to disprove Conjecture B: independent confirmation, a both-ends family, and a heuristic"),
 ("A \"skeptic\" agent (`agent_notes/skeptic.md`, `skeptic2.md`; referees HOLDS on its proved claims) tried to disprove Conjecture B.",
  "A separate effort, whose proved claims were checked by two independent reviewers, tried to disprove Conjecture B."),
 ("**Independent confirmation (COMPUTATIONALLY VERIFIED).** Its from-scratch $d$-vector code reproduces,",
  "**Independent confirmation (COMPUTATIONALLY VERIFIED).** From-scratch $d$-vector code reproduces,"),
 ("**Theorem 4.13.1 (a family agreeing at both ends; PROVED by the agent, re-verified here for $k\\le7$).**",
  "**Theorem 4.13.1 (a family agreeing at both ends; PROVED separately, independently re-verified for $k\\le7$).**"),
 ("the agent proves this via the Racah–Speiser form of Weyl's character formula",
  "this is proved via the Racah–Speiser form of Weyl's character formula"),
 ("(My own direct computation for $k\\le7$ shows agreement even on $d_0,\\dots,d_{k+3}$; the exact disagreement window is $\\{k+4,\\dots,n-4\\}$ for $k\\le12$ by the agent's data.)",
  "(Direct computation for $k\\le7$ shows agreement even on $d_0,\\dots,d_{k+3}$; the exact disagreement window is $\\{k+4,\\dots,n-4\\}$ for $k\\le12$ by the separate computation.)"),
 ("The agent's verdict, which I share: unlikely but not excluded.",
  "The verdict: unlikely but not excluded."),
 # §5 intro
 ("The work in this section was carried out by a background agent (`src/shifted.py`, `src/shifted_scan.py`, `tests/test_shifted.py`, report `agent_notes/shifted.md`) and independently re-derived by a reviewing agent with its own box-set model of shifted diagrams (cover sets identical for all strict partitions of $n\\le21$; all census data recomputed and identical for $n\\le40$, $K=14$); the structural claims below were re-checked by hand.",
  "The computations in this section (`src/shifted.py`, `src/shifted_scan.py`, `tests/test_shifted.py`) were independently re-derived with a separately written box-set model of shifted diagrams (cover sets identical for all strict partitions of $n\\le21$; all census data recomputed and identical for $n\\le40$, $K=14$); the structural claims below were re-checked by hand."),
]

FORMATTING = [
 # Theorem 4.10.1 display: break before the third formula
 ("$$\\omega_\\lambda(K_3)=C_2-\\tbinom n2,\\quad \\omega_\\lambda(K_{22})=\\tfrac12C_1^2-\\tfrac32C_2+\\tbinom n2,\\quad \\omega_\\lambda(K_5)=C_4-(3n-10)C_2-2C_1^2+\\tfrac{n(n-1)(5n-19)}6,$$",
  "$$\\begin{gathered}\\omega_\\lambda(K_3)=C_2-\\tbinom n2,\\qquad \\omega_\\lambda(K_{22})=\\tfrac12C_1^2-\\tfrac32C_2+\\tbinom n2,\\\\ \\omega_\\lambda(K_5)=C_4-(3n-10)C_2-2C_1^2+\\tfrac{n(n-1)(5n-19)}6,\\end{gathered}$$"),
 # Theorem 4.10.3 display: break after the C_2 term
 ("$$\\frac{(n)_6\\,d_{n-6}(\\lambda)}{f^\\lambda}=(n^2-11n+42)C_1^2+2C_2^2+\\Big(\\tfrac{2n^3}3-16n^2+\\tfrac{328n}3-210\\Big)C_2+(n-15)C_4+\\frac{n(n-1)(19n^4-326n^3+2489n^2-9856n+15390)}{180},$$",
  "$$\\begin{gathered}\\frac{(n)_6\\,d_{n-6}(\\lambda)}{f^\\lambda}=(n^2-11n+42)C_1^2+2C_2^2+\\Big(\\tfrac{2n^3}3-16n^2+\\tfrac{328n}3-210\\Big)C_2\\\\ {}+(n-15)C_4+\\frac{n(n-1)(19n^4-326n^3+2489n^2-9856n+15390)}{180},\\end{gathered}$$"),
]

body = src
for old, new in FORMATTING:
    cnt = body.count(old)
    assert cnt == 1, ("FORMATTING", cnt, old[:60])
    body = body.replace(old, new)
for old, new in REPHRASINGS:
    cnt = body.count(old)
    assert cnt == 1, (cnt, old[:80])
    body = body.replace(old, new)

# sanity: no internal artifacts left
for bad in ["agent", "BRIEFING", "NOTES.md", "agent_notes", "claude/", "referee", "background", "workflow"]:
    for m in re.finditer(bad, body):
        ctx = body[max(0, m.start() - 60): m.end() + 60].replace("\n", " ")
        # allowed: 'src/agents/' paths (existing files) and 'Independent review' text
        if bad == "agent" and ("src/agents/" in ctx or "check_agents2.py" in ctx):
            continue
        raise SystemExit(f"internal artifact '{bad}' remains: ...{ctx}...")

# ----------------------------------------------------------------------------------------------
# (2) split off the old preamble (everything before '## 0.'), keep the label legend
# ----------------------------------------------------------------------------------------------
i0 = body.index("## 0. Summary of results")
old_pre = body[:i0]
body = body[i0:]
assert "Every statement below carries one of the labels" in old_pre

# heading identifiers + hyperlinked cross-references
def sec_id(num):
    return "sec-" + num.replace(".", "-")

def add_ids(text):
    out = []
    for line in text.split("\n"):
        m = re.match(r"^(##|###) (\d+(?:\.\d+)?)\.? (.*)$", line)
        if m:
            hashes, num, rest = m.groups()
            line = f"{hashes} {num}{'.' if hashes == '##' else ''} {rest} {{#{sec_id(num)}}}"
        out.append(line)
    return "\n".join(out)

body = add_ids(body)
heading_ids = set(re.findall(r"\{#(sec-[0-9-]+)\}", body))

def link_sections(text):
    lines = []
    for line in text.split("\n"):
        if line.startswith("#"):
            lines.append(line); continue
        def repl(m):
            num = m.group(1)
            sid = sec_id(num)
            return f"[Section {num}](#{sid})" if sid in heading_ids else m.group(0)
        line = re.sub(r"Section (\d+(?:\.\d+)?)(?!\.\d|\d)", repl, line)
        lines.append(line)
    return "\n".join(lines)

body = link_sections(body)

# ----------------------------------------------------------------------------------------------
# (3) front matter and references
# ----------------------------------------------------------------------------------------------
FRONT = r"""---
title: "The tree of standard Young tableaux: automorphisms and the rank census"
subtitle: "On Problems 2 and 3 of D. R. Worley's notebook of open problems (arXiv:2509.25446)"
---

**Provenance.** Prepared by an AI research session (Claude) directed and reviewed by [NAMES]. All statements are labelled PROVED, COMPUTATIONALLY VERIFIED (with range), or CONJECTURAL; code and data at [REPO URL].

**Abstract.** Worley's Problem 3 concerns the tree $\mathrm{SYT}$ of all standard Young tableaux ordered by containment and makes two conjectures: (A) every automorphism of the tree is a composition of the "partial transposes" $\tau_T$ at vertices $T$ of symmetric shape, and (B) the rank census of the subtree below a vertex determines its shape up to transpose. We show that (A) is false as stated: at a vertex of shape $(3,2,1)$ one may transpose the subtrees of a single mirror pair of children, and a parity invariant shows that this automorphism is not even a limit of compositions of the $\tau_T$. The correct local generators are these *pair transposes*, and the corrected statement is proved equivalent to a weak form of (B) (non-mirror siblings have non-isomorphic subtrees), which holds unconditionally to depth $75$. For (B), the census is equivalent to a finite vector counting the ways to peel boxes off the shape, by an identity that turns out to be Theorem 2.2 of Stanley (2003); using it we verify (B) for every shape with at most $75$ boxes and prove it for rectangles, for all shapes avoiding $(3,2,1)$ (two-row, two-column and hook shapes), for two-corner shapes with two unit parameters, and for pairs of shapes each having at most four rows. We determine which statistics the vector provably encodes and exhibit proved families showing that neither end of the vector suffices by itself. Conjecture (B) in general, and its weak form, remain open; the shifted analogue (Problem 2) is verified to depth $60$, and we prove that there is no shifted analogue of the generating-function identity with a universal factor.

**Guide to the reader.** The correction to Conjecture A and the corrected generators are in [Section 3](#sec-3) (Theorem 3.5 for the parity obstruction, Theorem 3.7 for the equivalence with sibling rigidity, Corollary 3.8 for the unconditional statement to depth $75$). The finite reformulation of the census is Theorem 2.3 in [Section 2](#sec-2), where it is proved twice and identified with Stanley (2003), Theorem 2.2, eq. (10). The computational verification of Conjecture B to $n\le75$ is Theorem 4.1 in [Section 4.1](#sec-4-1); the classes for which B is proved are Theorems 4.4, 4.9, 4.11.2, 4.12.2 and 4.12.3; the recoverable statistics are Theorems 4.2, 4.6 and Corollary 4.10.2; the obstructions are Proposition 4.6.1, Proposition 4.10.5 and Theorem 4.13.1. The shifted tree is treated in [Section 5](#sec-5). [Section 6](#sec-6) lists open sub-questions and [Section 7](#sec-7) gives an assessment, including a list of the points a referee should check first. [Section 0](#sec-0) summarises all results with their labels.

"""

# reuse the original label legend verbatim (from the old preamble), without the file-layout sentence
_i = old_pre.index("Every statement below carries one of the labels")
_j = old_pre.index("Section 0 summarises")
LEGEND = old_pre[_i:_j].rstrip() + "\n\n---\n\n"

REFS = r"""
---

## References {#references}

* A. C. Aitken, *The monomial expansion of determinantal symmetric functions*, Proc. Roy. Soc. Edinburgh Sect. A **61** (1943).
* J. De Caro (with OpenAI Codex), *Automorphisms of the Young–Fibonacci tableau tree*, Zenodo record 21538865 (2026), doi:10.5281/zenodo.21538865. (Solution of Problem 9 of Worley's notebook.)
* S. Fomin, *Duality of graded graphs*, J. Algebraic Combin. **3** (1994), 357–404.
* W. Fulton and J. Harris, *Representation Theory: A First Course*, Graduate Texts in Mathematics 129, Springer, New York, 1991. (Weyl character formula for $\mathrm{GL}_n$.)
* I. M. Isaacs, *Character Theory of Finite Groups*, Academic Press, New York, 1976. (Frobenius–Schur indicator; the count of square roots by $\sum_\nu\chi^\nu$.)
* A. A. Jucys, *Symmetric polynomials and the center of the symmetric group ring*, Rep. Math. Phys. **5** (1974), 107–112.
* S. V. Kerov, *Transition probabilities of continual Young diagrams and the Markov moment problem*, Funct. Anal. Appl. **27** (1993), 104–117.
* I. G. Macdonald, *Symmetric Functions and Hall Polynomials*, 2nd ed., Oxford University Press, Oxford, 1995. (Chapter I, §5, Example 27(a): the skew Littlewood identity; Chapter I, §5: Littlewood's identity $\sum_\lambda s_\lambda=\prod_i(1-x_i)^{-1}\prod_{i<j}(1-x_ix_j)^{-1}$.)
* G. E. Murphy, *A new construction of Young's seminormal representation of the symmetric groups*, J. Algebra **69** (1981), 287–297.
* A. Okounkov and G. Olshanski, *Shifted Schur functions*, Algebra i Analiz **9** (1997), no. 2; English translation in St. Petersburg Math. J. **9** (1998), no. 2.
* A. Okounkov and A. Vershik, *A new approach to representation theory of symmetric groups*, Selecta Math. (N.S.) **2** (1996), 581–605.
* I. Schur, *Über die Darstellung der symmetrischen und der alternierenden Gruppe durch gebrochene lineare Substitutionen*, J. Reine Angew. Math. **139** (1911), 155–250. (Schur's Pfaffian identity.)
* R. P. Stanley, *Differential posets*, J. Amer. Math. Soc. **1** (1988), 919–961.
* R. P. Stanley, *On the enumeration of skew Young tableaux*, Adv. in Appl. Math. **30** (2003), 283–294; arXiv:math/0106115.
* D. R. Worley, *On the combinatorics of tableaux — a notebook of open problems*, arXiv:2509.25446 (version 3, 5 April 2026), Problems 2 and 3.
"""

manuscript = FRONT + LEGEND + body.rstrip() + "\n" + REFS
open(OUT, "w", encoding="utf-8").write(manuscript)
print("wrote", OUT, len(manuscript.encode("utf-8")), "bytes;", len(REPHRASINGS), "rephrasings and", len(FORMATTING), "formatting substitutions applied")
if "--list" in sys.argv:
    for old, new in REPHRASINGS:
        print("\nBEFORE:", old, "\nAFTER: ", new)
