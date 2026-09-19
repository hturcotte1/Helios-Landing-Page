# Review flags for the Worley package

**Mathematical errors noticed while reformatting: none.** Nothing in the mathematical body of
`results.md` (Sections 0–7) was changed. Audit (`build_manuscript.py` is the only path from
`results.md` to `manuscript.md`; the diff is in `body_diff.txt`): 403 body lines compared after
undoing the heading identifiers and cross-reference links; 27 lines differ; every difference is one
of the 35 listed substitutions (33 provenance rephrasings and 2 line breaks inserted inside
displayed formulas). No statement, proof, label, theorem number or computational range was touched.

The items below are **not errors**. They are decisions made while packaging, or points a reader
might raise, recorded here so that nothing is decided silently. Each should be looked at by the
humans named in the provenance block before the package is sent.

## Decisions made while packaging

1. **Erratum item (e) dropped.** The paragraph "Further results of the same computation" in
   Section 4.10 (results.md: "Third pass of the top-end angle") ends with five lettered items
   (a)–(e). Item (e) — "The printed inversion formula for $C_2$ in the agent's Corollary B′ has a
   slip ($3f$, not $6f$, in the denominator); Theorem 4.6 above is the correct form." — corrects a
   formula in an internal report (`agent_notes/topend.md`) that does not appear in the manuscript,
   so it was removed rather than rephrased. Items (a)–(d) are unchanged. Theorem 4.6 in the
   manuscript is the correct form, as results.md says.

2. **"Our formulation" wording.** results.md introduces Conjectures A and B with "(Worley,
   Problem 3, verbatim in `data/literature_search.md`)", pointing to the file that holds Worley's
   verbatim text; the two bullets that follow are the write-up's own formulation, not Worley's
   wording. The manuscript says "(our formulation of the two halves of Worley's Problem 3)". Check
   that this reading is the intended one.

3. **Section 0 kept.** The manuscript keeps results.md's Section 0 ("Summary of results") after
   the new abstract and guide to the reader, so the numbering of Sections 1–7 and all theorem
   numbers are unchanged. There is some overlap between the abstract and Section 0.

4. **`src/agents/` paths kept.** Script paths cited in results.md were kept whenever the files
   exist (`src/agents/topend_*.py`, `src/agents/classC_verify.py`, `src/agents/classC_extra.py`,
   `src/check_agents2.py`). The directory name `agents` therefore remains visible in the
   manuscript. Renaming would require changing the paths in the manuscript and the repository
   together.

5. **Pointers to `agent_notes/` removed.** Sentences that pointed to `agent_notes/*.md` for full
   proofs now read "obtained separately", "checked line by line", "the full case analysis is
   included with the code", etc. (all listed in the build script). For Lemma 4.11.1 / Theorem
   4.11.2 the manuscript says the full case analysis is included with the code; the full text is in
   `agent_notes/fathooks.md`, which the repository `README.md` points to. If the repository is published
   without `agent_notes/`, these sentences need changing.

6. **New text.** The abstract, the guide to the reader and the reference list are new. The abstract
   was checked against Section 0 and against the theorem statements it summarises (in particular:
   "pairs of shapes each having at most four rows" = Theorem 4.12.3 with Theorem 4.1 below 17;
   "no shifted analogue of the generating-function identity with a universal factor" = Proposition
   5.2). Please read both paragraphs once more against Section 0.

7. **Placeholders.** `[NAMES]` and `[REPO URL]` in `manuscript.md` (provenance block); `[NAME]`,
   `[REPO URL]`, `[EMAIL]` in `cover_email.txt`; `[NAMES]` in the repository `README.md`. The PDF must be
   rebuilt after filling them in (the pandoc command is in the repository `README.md`, section
   "Rebuilding the manuscript PDF").

## Points a reader may raise (no change made)

8. **Theorem 4.12.3 threshold.** The statement uses $n\ge17$, with the parenthetical "(the
   original derivation stated $n\ge21$; the exhaustive check gives 17)". As results.md explains,
   the sign lemma gives $C_1>0$ for $n\ge21$ and an exhaustive check of the $\le4$-row shapes with
   $n\le20$ extends this to $n\ge17$; so the PROVED label for $17\le n\le20$ rests on that finite
   check. This is consistent, but a referee may ask for the check to be stated as a computation.

9. **Theorem 4.13.1: proved prefix vs. observed prefix.** The theorem proves agreement of
   $d_j(\lambda_k)=d_j(\mu_k)$ for $0\le j\le k$; the following parenthesis reports agreement on
   $d_0,\dots,d_{k+3}$ for $k\le7$ by direct computation and a disagreement window
   $\{k+4,\dots,n-4\}$ for $k\le12$. The proved statement is weaker than what is observed. This is
   as in results.md; nothing was changed.

10. **References.** Every entry corresponds to a work named in the text. Two entries support
    standard results that the text names without a bibliographic citation and could be dropped:
    Fulton–Harris (Weyl's character formula, used in the proof of Theorem 4.13.1 via its
    Racah–Speiser form) and Isaacs (the Frobenius–Schur count of square roots, Lemma 4.5).
    Page numbers were omitted where they were not verified: Aitken (1943) and
    Okounkov–Olshanski (1997/1998). Worley's notebook is cited as arXiv version 3 (5 April 2026),
    as recorded in `data/literature_search.md`; that file also mentions an author-hosted draft
    dated August 2026 in which Problems 2 and 3 are unchanged. Check for a newer arXiv version
    before sending.

11. **Headings "(re-verified)".** The headings of Sections 4.10–4.12 keep results.md's
    "(re-verified)" qualifier, now without the internal angle names. Harmless, but unusual in a
    manuscript heading; the humans may prefer to delete the qualifier (a formatting change only).
