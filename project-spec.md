# Project spec: LIONESS edges across pregnancy in SLE and RA

| | |
|---|---|
| Proposal | `proposal.md`, final 2026-11-13 |
| Authors | as in `DESCRIPTION`, set by `tools/authors.R` before the final render |
| Written | 2026-11-13 |
| Status | draft, section 3 completed after the verdict |

## 1. Question and why it matters

**Question (one sentence):** Do LIONESS edges separate pregnant from non-pregnant samples of the same women beyond the IFN score, cell composition, their squares, the gene levels, and their own first-order features, and does the answer differ between SLE and RA?

**Why it matters:** A change in LIONESS edges over pregnancy is read as rewiring of immune co-expression. If the mechanism holds, much of that change is the squared change of the factors pregnancy moves (composition in every woman, IFN in SLE), and a user would fit the base with the squared score and the gene levels before claiming that edges add anything, and would check every network result against the first-order formula.

## 2. Mechanism we suspect

A LIONESS edge equals, to first order, $r + z_{qi} z_{qj} - \tfrac{r}{2}(z_{qi}^2 + z_{qj}^2)$. Under a shared factor $F$ with loadings $\lambda$, its conditional mean is linear in $F_q^2$ with coefficient $\lambda_i \lambda_j - \tfrac{r}{2}(\lambda_i^2 + \lambda_j^2)$, so the difference in mean edge between pregnancy and the postpartum period is that coefficient times the difference in $\mathrm{E}[F^2]$, with the period-specific correlation unchanged. The gap is zero when pregnancy moves neither the mean nor the spread of $F$ (the null cell), and grows with the shift and the scale of $F$ (the factor cells). If the mechanism is wrong, edges would gain over E1 in factor-only cells, or edges and the first-order formula would classify differently.

## 3. Pilot

**Code:** `pilot/pilot.R` (`Rscript pilot/pilot.R --reps <n>`), then `Rscript pilot/verdict.R`; data description `Rscript pilot/describe_data.R`.

**Pass rule:** `pilot/preregistration.md`, section 2, copied from proposal section 7.

**Verdict file:** `pilot/results/verdict.rds`, after the run.

**Plot:** `pilot/results/auc_gains.png`.

**The three pilot answers:** to be filled after the verdict.

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
