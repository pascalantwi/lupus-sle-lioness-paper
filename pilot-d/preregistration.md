# Pre-registration: `pilot-d/`

| | |
|---|---|
| Proposal | `proposal.md`, final 2026-11-16 (re-scope: personal-network stability), section 7, with the change of section 5 below |
| Written | 2026-11-16 |
| Committed | in the commit "Stage 2, step 2 (pilot-d): pre-registration", before any pilot-d code or run |

## 1. The three questions

1. Does the effect exist? When a shared factor is stable within persons, do LIONESS edges of the genes it moves look more personal (stable over visits) than other edges, with no personal co-expression; and do LIONESS edges and their first-order formula show the same stability?
2. Is it large enough at the levels measured in GSE65391? At the measured interferon stability and 131 people with six visits, is the excess stability visible, and does it make people identifiable?
3. Is there a fix? Do partial-correlation edges given the IFN score and the observed cell mix remove the factor's excess, keep a real personal module, and do edges predict the next visit's activity beyond the genes when activity changes co-expression, and not otherwise?

## 2. The pass rule

For each replicate and edge type $T$ (raw LIONESS, first-order F1, partial), $I_T(g)$ is the mean intraclass correlation over visits of the edges in pair group $g$: IFN pairs (within genes 6 to 20, 105 pairs), module pairs (within genes 21 to 30, 45), cell-mix pairs (within genes 31 to 50, 190), and cross pairs (one gene of 6 to 20 and one of 21 to 30, 150). The excess is $X_T(g) = I_T(g) - I_T(\text{cross})$, computed within the replicate. $P_T(g)$ is the fingerprinting accuracy on group $g$: the share of people whose visit-2 edge vector over the group has its highest Pearson correlation with their own visit-1 vector. $C(M)$ is the pooled out-of-fold concordance of model $M$ for the next visit's activity. Means over replicates, MCSE as the replicate standard deviation over the square root of the replicates; differences paired within a replicate.

| Q | Clause | Quantity | Threshold | Cells | Within its MCSE |
|---|---|---|---|---|---|
| 1 | (a) | $X_{\text{raw}}(\text{IFN})$ | $\ge 0.15$ | 3 | $\lvert v - 0.15\rvert \le 3$ MCSE: close call |
| 1 | (b) | $\lvert I_{\text{raw}}(g) - I_{\text{F1}}(g)\rvert$ for every group | $\le \max(0.02, 3\,\text{MCSE})$ | 1 to 6 | $\max(0.02, 3\,\text{MCSE}) < v \le 0.02 + 3\,\text{MCSE}$: close call |
| 2 | (a) | $X_{\text{raw}}(\text{IFN})$ | $\ge 0.05$ | 2 | as 1(a) |
| 2 | (b) | $P_{\text{raw}}(\text{IFN}) - P_{\text{partial}}(\text{IFN})$ | $\ge 0.05$ | 2 | as 1(a) |
| 3 | (a) | $X_{\text{partial}}(\text{IFN})$ | $\lvert v\rvert \le \max(0.03, 3\,\text{MCSE})$ | 2, 3 | $\max(0.03, 3\,\text{MCSE}) < \lvert v\rvert \le 0.03 + 3\,\text{MCSE}$: close call |
| 3 | (b) | $X_{\text{partial}}(\text{module}) - 0.75\, X_{\text{raw}}(\text{module})$ | $\ge 0$ | 4 | $\lvert v\rvert \le 3$ MCSE: close call |
| 3 | (c) | $C(\text{N1}) - C(\text{E1})$ | $\ge 0.03$ in cell 6; $\lvert v\rvert \le \max(0.01, 3\,\text{MCSE})$ in cells 1 and 2 | 1, 2, 6 | as 1(a) in cell 6, as 3(a) with 0.01 in cells 1 and 2 |

**How the answers combine.** Pass when every clause passes. Fail when 1(b), 3(a), or 3(b) fails. Close call otherwise: when 1(a), 2(a), 2(b), or 3(c) fails, or when any clause is a close call for itself. A close call is stated both ways at the verdict.

## 3. Design

- **Data-generating mechanism** (`SimulateRepeatedCohort()`). $W = 131$ people with $V$ visits each. Fifty genes. The IFN factor $F_{wv} = a_w + b_{wv}$ with $a_w \sim N(0, \kappa_F)$ and $b_{wv} \sim N(0, 1 - \kappa_F)$; genes 1 to 20 load 0.8 on it. The cell-mix factor $C$, built the same way with $\kappa_C = 0.22$, is observed; genes 31 to 50 load 0.6 on it. Every gene carries the person's own level $u_{wg} \sim N(0, 0.5)$, the same at every visit, and unit noise. The noise of genes 21 to 30 (the module) has correlation $\rho_{wv} = 0.3 + \tau\, s_w + \gamma\,(\Phi(\alpha_{wv}) - 0.5)$, with $s_w$ uniform on $(-1, 1)$ fixed for the person (a real personal network when $\tau > 0$), activity $\alpha_{wv}$ = the person's level ($N(0, 0.39)$) plus visit noise ($N(0, 0.61)$), and $\gamma$ = 0.4 in cell 6 and 0 elsewhere.
- **Cells (6).**

| Cell | $\kappa_F$ | $\tau$ | Visits | $\gamma$ | Purpose |
|---|---|---|---|---|---|
| 1 | 0 | 0 | 6 | 0 | null: no stable factor, no personal module |
| 2 | 0.55 | 0 | 6 | 0 | measured IFN stability, no personal module |
| 3 | 0.9 | 0 | 6 | 0 | stress: very stable IFN |
| 4 | 0.55 | 0.25 | 6 | 0 | real personal module |
| 5 | 0.55 | 0.25 | 3 | 0 | three visits (reported, judged only in 1(b)) |
| 6 | 0.55 | 0 | 6 | 0.4 | activity changes the module |

- **Edges.** On all $W V$ samples of a replicate (the cohort as reference, as users do): raw LIONESS edges (`LionessEdges()`), first-order edges (`LionessEdgesFirstOrder()`), partial-correlation LIONESS edges (`LionessPartialEdges()`) given the IFN score (mean of genes 1 to 5 standardized over the samples) and the observed $C$, over genes 6 to 50.
- **Intraclass correlation.** One-way, for balanced visits, per edge: $(\text{MSB} - \text{MSW}) / (\text{MSB} + (V - 1)\,\text{MSW})$.
- **Fingerprinting.** Visit 2 against visit 1, as defined in section 2.
- **Prediction** (cells 1, 2, 6). Each visit with a next visit is a sample; its features predict the next visit's activity. Models and features as in `pilot/` (`PredictModelsB()`, Gaussian): $R_0$ (observed $C$, the IFN score's residual on it), $R_1$ (with squares), E1, F1, N1, N1 with $R_0$; outer and inner folds by person; the concordance index.
- **Replicates per cell.** From a timing run of ten replicates per cell with the workers loaded: 200 if the projected run fits one hour on three cores, otherwise the largest multiple of 20 that fits, not below 60.
- **Seeds.** Cell $c$ has seed $21061116 + 1000c$; replicate $r$ uses the cell seed plus $r$.
- **Cells with a known answer.** Cell 1 ($X$ near 0 for every edge type; fingerprinting by group as the person levels allow; $C(\text{N1}) - C(\text{E1})$ near 0). Tests: the intraclass correlation against a closed form on a balanced design; the intraclass correlation of $F^2$ against $\kappa^2$ in a large simulation; fingerprinting accuracy 1 for identical visits and near $1/W$ for independent ones; the partial-correlation edges against the inverse-covariance formula (first project).
- **Time box.** One hour of compute on three cores.

## 4. What is reported whatever the outcome

$I_T(g)$, $X_T(g)$, and $P_T(g)$ for every type, group, and cell; the concordance of every model in cells 1, 2, 6; the verdict table; the known-answer table; one figure of the excess stability and one of fingerprinting by cell, in `pilot-d/results/`; the data description of GSE65391 (`pilot-d/results/data_summary.rds`).

## 5. Readings of open cases and the change from the proposal

- 2026-11-16, before any code. **Change from the proposal's draft rule, and why.** The proposal judged absolute intraclass correlations and fingerprinting on all pairs. With person-specific gene levels in the model (people differ in level, as in real blood), every edge is partly person-specific through the product of the person's own levels, whatever the factors do; absolute thresholds would judge that, not the factor. The rule therefore judges the excess of the IFN pairs over cross pairs within the same replicate, and fingerprinting on the IFN pairs raw against partial. The thresholds were re-set on that scale (0.15, 0.05, 0.03) before any run. Under the workflow (stage 2, step 1.4) this change is reported to the user before anything runs.
- 2026-11-16. Clause 3(a) is two-sided as written.

## 6. Amendments before the run
- 2026-11-16, after the data description and before any run. **Cell-mix stability corrected from 0.22 to 0.45.** The proposal's value came from the agent's scratch check (`log/01-rescope-data-check/g65_icc.R`), whose intraclass-correlation helper misaligned values and groups whenever a child had a missing count. The package's tested `EdgeICC()`, in `pilot-d/describe_data.R` (`pilot-d/results/data_summary.rds`), gives 0.41 for neutrophil and 0.48 for lymphocyte percentage; the IFN score (0.55) and SLEDAI (0.39) have no missing values and are unchanged. $\kappa_C$ is set to 0.45, between the two, in `PilotDDesign()`.
