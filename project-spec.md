# Project spec: LIONESS edges are quadratic functions of expression: closed forms, an interferon artifact, and equivalence with first-order features

| | |
|---|---|
| Proposal | `proposal.md`, final 2026-11-16 (re-scope) |
| Authors | as in `DESCRIPTION`, set by `tools/authors.R` before the final render |
| Written | 2026-11-16 |
| Status | draft; section 3 completed after the verdict. The pregnancy spec is the record in git (`git show c677b11:project-spec.md`) |

## 1. Question and why it matters

**Question (one sentence):** Over repeated blood samples of the same patient, how stable are her LIONESS edges, can they identify her, and do they predict her next disease activity, beyond her interferon level, her measured blood-cell counts, her gene levels, and the edges' own first-order formula?

**Why it matters:** A personal network is a biomarker only if it is reproducible and carries something the cheap measures do not. If its stability is the person's stable interferon level and her own gene levels, a user would measure those directly; if a corrected network keeps a stable personal part, that part is the candidate biomarker.

## 2. Mechanism we suspect

A LIONESS edge equals, to first order, $r + z_{qi} z_{qj} - \tfrac{r}{2}(z_{qi}^2 + z_{qj}^2)$, so it carries the square of any factor that moves both genes. For a factor with within-person intraclass correlation $\kappa$, the intraclass correlation of $F^2$ is $\kappa^2$ (proposal section 4), so edges of genes sharing a stable factor look personal with no personal co-expression; a factor that changes between visits adds within-person noise. Each person's own stable gene levels enter the product term too, so every edge is partly person-specific through levels alone. The gap is zero for a factor with $\kappa = 0$ (the null cell) and grows with $\kappa$. If the mechanism is wrong, partial-correlation edges would keep the factor-driven stability, or edges and their formula would differ.

## 3. Pilot

**Code:** `pilot-d/pilot.R` (`Rscript pilot-d/pilot.R --reps <n>`), then `Rscript pilot-d/verdict.R`; data description `Rscript pilot-d/describe_data.R`.

**Pass rule:** `pilot-d/preregistration.md`, section 2.

**Verdict file:** `pilot-d/results/verdict.rds`, after the run.

**Plot:** `pilot-d/results/stability.png`.

**Status:** `piloted-close-call` (2026-11-16), 200 replicates per cell, 6 cells, three cores, 2,714 seconds; verdict file `pilot-d/results/verdict.rds`. Numbers rounded from `pilot-d/results/cell_summaries.rds`.

1. Does the effect exist? Edges equal their formula in stability (largest difference under 0.0001; 1(b) pass). Stable interferon raised the interferon pairs' stability (from about 0.07 to 0.14 at the measured level, 0.22 when very stable), but cross pairs rose too, and the pre-registered excess was 0.064 at the stress level, below 0.15 (1(a) fail).
2. Large enough at the measured level? No: excess 0.016 (2(a) fail); the fix changed identification on interferon pairs by 0.004 (2(b) fail).
3. Is there a fix? Yes: the fixed excess was 0.002 to 0.007 (3(a) pass), and the fix kept a real personal module (3(b) pass). Edges did not predict the next visit's activity (concordance near 0.48; 3(c) fail in cell 6). Single-sample edges were unreliable over visits (intraclass correlation 0.07 to 0.22).

**Pilot-e (network metrics, proposal section 12):** `piloted-fail` (2026-11-16), 120 replicates per cell, 4 cells, 3,442 seconds; verdict file `pilot-e/results/verdict.rds`; numbers rounded from `pilot-e/results/cell_summaries.rds`.

1. Closed forms: a module's connectivity correlated about 0.97 with its squared expression score in every cell; raw and first-order edges gave the same stability for all 19 metrics (largest difference under 0.0004).
2. Stability: IFN connectivity was more stable than IFN edges by about 0.10, and as stable as the squared IFN score (difference about -0.01).
3. Correction: with the five-gene score, the corrected IFN connectivity stayed more stable than the cross connectivity by about 0.09 and 0.15 (fail); the corrected module connectivity exceeded the cross by about 0.056 with a real personal module (close call).

**Pilot-f (own-module correction, proposal section 13):** `piloted-fail` (2026-11-16), 120 replicates per cell, 4 cells, 3,136 seconds; verdict file `pilot-f/results/verdict.rds`. Clauses 1 and 2 passed as in pilot-e. The own-module correction left the IFN set's connectivity as stable as raw (about 0.25 at the measured level, 0.62 when interferon was very stable; 3(a) fail), and the personal module's excess was about 0.036 (3(b) fail): a set's connectivity is its squared expression score up to scale, which no within-set adjustment changes (`pilot-f/decisions.md`, D2).

## 4. ADEMP skeleton

**Aims.** (1) How much of the stability of LIONESS edges over repeated samples comes from a stable shared factor? (2) Does a person's network identify her, and is that identification the factor? (3) Do partial-correlation edges remove the factor's part and keep a real personal network? (4) Do edges predict the next visit's activity beyond the genes, and do they equal their formula throughout?

**Data-generating mechanisms.**

- Factors varied: the IFN factor's within-person intraclass correlation $\kappa_F$ (0, 0.55 as measured, 0.9); a real personal module ($\tau$ 0 or 0.25); visits per person (6, the GSE65391 median, or 3); whether activity changes the module.
- Factors fixed, and why: 131 people (the GSE65391 children with two or more visits); 50 genes and the loadings of the first project; the cell-mix factor's intraclass correlation 0.45 (measured 0.41 and 0.48) and activity's 0.39; person-specific gene levels of variance 0.5 (people differ in level); normal noise.
- A cell where the answer is known: the null cell, with no stable factor and no personal module.

**Estimands.** The intraclass correlation of each edge over visits; the excess of the factor pairs' intraclass correlation over that of cross pairs (one factor gene, one module gene), which share the person-level gene effects but not the squared factor; fingerprinting accuracy; the concordance of next-visit activity predictions.

**Methods compared.** Raw LIONESS edges, first-order edges, partial-correlation LIONESS edges given the IFN score and the observed cell-mix factor; for prediction, $R_0$, $R_1$, E1, F1, N1 as in the first project, Gaussian elastic net, folds by person.

**Performance measures.** Mean intraclass correlation by pair group; fingerprinting accuracy; concordance; each with its MCSE.

**Replicates and Monte Carlo SE target.** 200 per cell if the time box allows; the MCSE of a mean intraclass correlation over 105 pairs is then well under 0.01.

## 5. Application

**Dataset:** GSE65391 (@banchereau2016personalized), GEO; the series matrix and the GPL10558 table in `data/raw/`, fetched by `analysis/00_download.R`, not committed.

**Comparison:** within each of the 131 children with two or more visits; 897 visits; the authors' training (91 children) and test (40) split as a second check.

**What to show:** how much of a child's network stability and identifiability is her interferon level, her cell mix, and her gene levels, and whether corrected networks predict her next disease activity.

## 6. Target journal

**Target journal:** BMC Bioinformatics

**Fallback:** Statistical Applications in Genetics and Molecular Biology

**The target's word limit and AI policy:** [TO CONFIRM: copied from the journal profile at stage 3]{.tc}

## 7. Milestones

| Stage | Ends with | Due |
|---|---|---|
| 2. Pilot | `pilot-d/pilot-report.html`, the verdict file, `pilot-d/decisions.md` | 2026-11-16 |
| 3. Manuscript, design stop | `protocol.md` from the approved design | [TO CONFIRM: date]{.tc} |
| 3. Manuscript, first draft | `manuscript/manuscript.html` in draft mode | [TO CONFIRM: date]{.tc} |
| Submission | the final render and the checker in final mode | [TO CONFIRM: date]{.tc} |

## 8. Out of scope

- SLE against RA (the user's decision, 2026-11-16).
- Lupus against healthy as a classification (too easy, the supervisor).
- Other single-sample methods; drug effects; ancestry and sex differences; subtypes (each set aside in `log/01-idea-review.md`, D-R3).
- Merging GSE65391 with any other study.

## 9. Reading

1. @kuijjer2019lioness. After reading, a reader can answer: what reproducibility of LIONESS networks did the authors test, and against what?
2. @petri2019stability. After reading, a reader can answer: which lupus gene signature is most stable within patients, and does its change track disease activity?
3. @finn2015fingerprinting. After reading, a reader can answer: how is fingerprinting accuracy defined, and what made individuals identifiable?
