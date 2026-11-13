# Project spec: LIONESS edges across pregnancy in SLE and RA

| | |
|---|---|
| Proposal | `proposal.md`, final 2026-11-13 |
| Authors | as in `DESCRIPTION`, set by `tools/authors.R` before the final render |
| Written | 2026-11-13 |
| Status | draft; section 3 completed 2026-11-13 after the verdict |

## 1. Question and why it matters

**Question (one sentence):** Do LIONESS edges separate pregnant from non-pregnant samples of the same women beyond the IFN score, cell composition, their squares, the gene levels, and their own first-order features, and does the answer differ between SLE and RA?

**Why it matters:** A change in LIONESS edges over pregnancy is read as rewiring of immune co-expression. If the mechanism holds, much of that change is the squared change of the factors pregnancy moves (composition in every woman, IFN in SLE), and a user would fit the base with the squared score and the gene levels before claiming that edges add anything, and would check every network result against the first-order formula.

## 2. Mechanism we suspect

A LIONESS edge equals, to first order, $r + z_{qi} z_{qj} - \tfrac{r}{2}(z_{qi}^2 + z_{qj}^2)$. Under a shared factor $F$ with loadings $\lambda$, its conditional mean is linear in $F_q^2$ with coefficient $\lambda_i \lambda_j - \tfrac{r}{2}(\lambda_i^2 + \lambda_j^2)$, so the difference in mean edge between pregnancy and the postpartum period is that coefficient times the difference in $\mathrm{E}[F^2]$, with the period-specific correlation unchanged. The gap is zero when pregnancy moves neither the mean nor the spread of $F$ (the null cell), and grows with the shift and the scale of $F$ (the factor cells). If the mechanism is wrong, edges would gain over E1 in factor-only cells, or edges and the first-order formula would classify differently.

## 3. Pilot

**Code:** `pilot/pilot.R` (`Rscript pilot/pilot.R --reps <n>`), then `Rscript pilot/verdict.R`; data description `Rscript pilot/describe_data.R`.

**Pass rule:** `pilot/preregistration.md`, section 2, copied from proposal section 7.

**Verdict file:** `pilot/results/verdict.rds`, after the run.

**Plot:** `pilot/results/auc_gains.png`; pilot-b: `pilot-b/results/shares.png`, `pilot-b/results/model_gains.png`; pilot-c: `pilot-c/results/shares.png`, `pilot-c/results/model_gains.png`.

**Second pilot (`pilot-b/`, proposal section 12, the user's three aims):** `piloted-close-call` (2026-11-13), 200 replicates per cell, 10 cells, three cores, 965 seconds; verdict file `pilot-b/results/verdict.rds`. Numbers rounded from `pilot-b/results/cell_summaries.rds`.

1. Does the fix work? Yes: raw edges called 0.17 of unchanged factor pairs changed at the stress level; partial-correlation edges gave 0.044 to 0.049 in every factor-only cell, within and between diseases, and kept a real module change (0.22 against 0.23).
2. Does it matter at the study's levels? Modestly: no false rewiring at the SLE level, 0.107 at the RA level (close call), a false SLE-against-RA difference of 0.079 (below 0.10, fail), and a real SLE-only change found in 0.157 of module pairs (close call).
3. Prediction: edges equalled their formula (N1 against F1 within 0.0007); the base with the squared score removed the activity artifact; edges did not detect activity-linked co-expression across women (fail).

**Third pilot (`pilot-c/`, proposal section 13: every visit, marker-gene cell mix, activity within woman; pass rule unchanged):** `piloted-close-call` (2026-11-13), 400 replicates per cell, 10 cells, three cores, 2,233 seconds; verdict file `pilot-c/results/verdict.rds`. Numbers rounded from `pilot-c/results/cell_summaries.rds`.

1. Does the fix work? Yes, again: partial-correlation edges gave 0.047 to 0.050 in every factor-only cell (raw 0.18 at the stress level), and kept real change within a disease (0.28 against 0.28) and between diseases.
2. Does it matter at the study's levels? Little: with the cell mix measured from marker genes, raw edges called 0.046 (SLE level) and 0.059 (RA level) of factor pairs changed, and 0.059 differently changed between diseases (2(a), 2(b) fail); a real SLE-only change was found in 0.18 of module pairs (2(c) pass).
3. Prediction: edges equalled their formula (within 0.0003); the squared-score base removed the activity artifact within woman; edges found no activity-linked co-expression within woman (3(c) fail).

**Status:** `piloted-close-call` (2026-11-13), 200 replicates per cell, 8 cells, three cores, 1,622 seconds.

**The three pilot answers** (`pilot/results/verdict.rds`; numbers rounded from `pilot/results/cell_summaries.rds`):

1. Does the effect exist? Edges equal their formula: N1 and F1 differed by at most 0.0011 in AUC in any cell (1(b) pass). With a module change at 40 women, N1 gained 0.043 (MCSE 0.0039) over E1, within 3 MCSE of 0.05 (1(a) close call).
2. Is it large enough at realistic sizes? Not by classification at 20 women: the artifact was 0.005 (MCSE 0.0035), below 0.02 (2(a) fail), and real rewiring added 0.009 (MCSE 0.0034) over E1, below 0.03 (2(b) fail). By the paired edge test, yes: with no change in co-expression, 0.45 of factor-gene pairs were called rewired at 20 women under a large factor change, against 0.047 in the null cell.
3. Is there a fix? For classification, yes: with the E1 base edges added between -0.0004 and 0.0015 in the factor-only cells (3(a) pass). For the paired edge test, none was tested.

## 4. ADEMP skeleton

**Aims.** (1) Do LIONESS edges and their first-order formula classify pregnancy stage alike? (2) How much do edges gain from a factor that changes over pregnancy, and does the E1 base remove it? (3) When co-expression truly changes, do edges detect it beyond E1 at the study's size?

**Data-generating mechanisms.**

- Factors varied, with levels: what pregnancy changes (nothing; the factor $F$, moderate or large; the module's correlation by 0.2 or 0.4; both); the number of women (20, 40).
- Factors fixed, and why: four samples per woman (two in each period, the median the study provides per disease); 50 genes (as in the first project, so results are comparable); loadings 0.8 and 0.4; woman-level gene offsets of variance 0.5 (the between-person part of whole-blood expression); normal noise (heavy tails were tested in the first project).
- A cell where the answer is known: the null cell, every AUC 0.5.

**Estimands.** The out-of-sample AUC of each model for the period of a new woman's sample; the paired AUC differences $\Delta_1$, $\Delta_2$, and the artifact $G$.

**Methods compared.** $R_0$, $R_1$ (logistic); E1, E5, F1, N1, N1 with $R_0$ (elastic net, `glmnet` `cv.glmnet`, mixing 0.5, 5 inner folds grouped by woman, binomial deviance, `lambda.min`, base terms unpenalized). Outer 5-fold cross-validation by woman.

**Performance measures.** Mean out-of-fold AUC and paired differences with MCSE; the share of gene pairs called different between periods by a paired test (apparent rewiring).

**Replicates and Monte Carlo SE target.** 200 per cell in the pilot if the time box allows, an MCSE of about 0.003 for a difference with standard deviation 0.04.

## 5. Application

**Dataset:** GSE235508 (Lien et al. 2024), GEO; the series matrix in `data/raw/` for the pilot's description, the counts at stage 3 after approval.

**Comparison:** within SLE and within RA separately, pregnant (trimesters 1 to 3) against non-pregnant (pre-pregnancy and postpartum) samples of the same women, folds by woman; counts per group from `pilot/results/data_summary.rds`.

**What to show:** whether the apparent rewiring over pregnancy, and any gain of edges over the base, disappears once the squared IFN and composition scores and the gene levels are in the base, and whether that differs between SLE and RA.

## 6. Target journal

**Target journal:** BMC Bioinformatics

**Fallback:** Statistical Applications in Genetics and Molecular Biology

**The target's word limit and AI policy:** [TO CONFIRM: copied from the journal profile at stage 3]{.tc}

## 7. Milestones

| Stage | Ends with | Due |
|---|---|---|
| 2. Pilot | `pilot/pilot-report.html`, the verdict file, `pilot/decisions.md` | 2026-11-13 |
| 3. Manuscript, design stop | `protocol.md` from the approved design | [TO CONFIRM: date]{.tc} |
| 3. Manuscript, first draft | `manuscript/manuscript.html` in draft mode | [TO CONFIRM: date]{.tc} |
| 3. Manuscript, settled | the draft after the revision rounds | [TO CONFIRM: date]{.tc} |
| 4. Background report | `background/background.html` | [TO CONFIRM: date]{.tc} |
| 5. Mock review | the last `review/round-N/` | [TO CONFIRM: date]{.tc} |
| Submission | the final render and the checker in final mode | [TO CONFIRM: date]{.tc} |

## 8. Out of scope

- A second dataset merged with GSE235508 (the user's rule: no studies combined).
- Pregnancy outcomes or flare prediction (a predictive project of its own).
- Other single-sample network methods (SSN, iENA, SWEET, BONOBO): covered in the first project.
- Partial-correlation edges and regularized discriminant analysis (the first project's N2 and E4).
- Mixed models of the edges over time: a different estimand, parked as its own idea.

## 9. Reading

1. @kuijjer2019lioness. After reading, a reader can answer: how is a LIONESS sample network defined, and what use across time do the authors propose?
2. @lien2024pregnancy. After reading, a reader can answer: what happens to the IFN signature and to cell composition over pregnancy in SLE and in RA?
3. @devlin1975robust. After reading, a reader can answer: what is the influence function of the correlation coefficient?
