# Pre-registration: `pilot-b/`

| | |
|---|---|
| Proposal | `proposal.md`, section 12 (revision 1, the user's three aims), approved by the user 2026-11-13 |
| Written | 2026-11-13 |
| Committed | in the commit "Stage 2, step 2 (pilot-b): pre-registration", before any pilot-b code or run |

## 1. The three questions

1. **Does the fix work?** When pregnancy moves only the shared factors, does the paired test of raw LIONESS edges call rewiring, and do partial-correlation LIONESS edges (given the IFN score and composition) bring it back to its nominal level, for a change within a disease (aim 1) and for a difference between SLE and RA (aim 3), while keeping a real change?
2. **Does it matter at the study's levels and sizes?** At the factor changes measured in GSE235508 (`results/factor_shift.rds`) and with 18 SLE and 27 RA women, do raw edges call false rewiring, and a false SLE-against-RA difference, and can the fixed edges detect a real SLE-only change?
3. **Prediction (aims 1 and 2).** Do LIONESS edges equal their first-order formula in predicting pregnancy stage and disease activity, and does the E1 base remove the artifact from activity prediction while edges still detect a real activity-linked change in co-expression?

## 2. The pass rule

A share is the proportion of gene pairs with $p < 0.05$ in a replicate, averaged over replicates. Factor pairs are the pairs within genes 6 to 20 (IFN factor, 105 pairs) and within genes 31 to 50 (composition factor, 190 pairs), 295 in all; module pairs are the 45 pairs within genes 21 to 30. "Raw" means LIONESS Pearson edges; "partial" means LIONESS edges of the partial correlation given the IFN score and composition. $A(M)$ is the pooled out-of-fold AUC (pregnancy stage) or concordance index (activity) of model $M$. Paired differences and their MCSE as in `pilot/preregistration.md`.

| Q | Clause | Quantity | Threshold | Cells | Within its MCSE |
|---|---|---|---|---|---|
| 1 | (a) | raw share, factor pairs (paired test) | $\ge 0.10$ | 4 | $\lvert v - 0.10\rvert \le 3$ MCSE: close call |
| 1 | (b) | $P_c$ = partial share, factor pairs (paired test in 2, 3, 4; between-disease test in 8, 10) | $P_c \le 0.05 + \max(0.01, 3\,\text{MCSE}_c)$ in each cell | 2, 3, 4, 8, 10 | $0.05 + \max(0.01, 3\,\text{MCSE}_c) < P_c \le 0.06 + 3\,\text{MCSE}_c$: close call |
| 1 | (c) | partial share minus 0.75 times raw share, module pairs (paired test), per replicate | $\ge 0$ | 5 | $\lvert v \rvert \le 3$ MCSE: close call |
| 2 | (a) | raw share, factor pairs (paired test) | $\ge 0.10$ in cell 2 or in cell 3 | 2, 3 | the larger of the two; $\lvert v - 0.10\rvert \le 3$ MCSE: close call |
| 2 | (b) | raw share, factor pairs (between-disease test) | $\ge 0.10$ | 8 | as 1(a) |
| 2 | (c) | partial share, module pairs (between-disease test) | $\ge 0.15$ | 9 | $\lvert v - 0.15\rvert \le 3$ MCSE: close call |
| 3 | (a) | $D_c = A(\text{N1}) - A(\text{F1})$ | $\lvert D_c\rvert \le \max(0.01, 3\,\text{MCSE}_c)$ in every cell | 1 to 7 | $\max(0.01, 3\,\text{MCSE}_c) < \lvert D_c\rvert \le 0.01 + 3\,\text{MCSE}_c$: close call |
| 3 | (b) | $A(\text{N1}) - A(\text{E1})$, activity | as 3(a) | 6 | as 3(a) |
| 3 | (c) | $A(\text{N1}) - A(\text{E1})$, activity | $\ge 0.03$ | 7 | $\lvert v - 0.03\rvert \le 3$ MCSE: close call |

**How the answers combine.** Pass when every clause passes. Fail when 1(b), 3(a), or 3(b) fails. Close call otherwise: when 1(a), 1(c), 2(a), 2(b), 2(c), or 3(c) fails, or when any clause is a close call for itself. A close call is stated both ways at the verdict.

## 3. Design

- **Data-generating mechanism** (`SimulateCohortB()`). $W$ women, four samples each. Fifty genes. Two shared factors, each the woman's level (variance $\kappa$) plus sample noise (variance $1 - \kappa$): the IFN factor $F$, on which genes 1 to 20 load 0.8, and the composition factor $C$, on which genes 31 to 50 load 0.6, observed without error. In pregnant samples each factor becomes $\mu + \sigma \times$ (its non-pregnant draw). Every gene carries the woman's offset (variance 0.5) and unit noise; the noise of genes 21 to 30 has correlation $\rho$, 0.3 by default.
- **Factor levels** (from `results/factor_shift.rds`, rounded): SLE: $F$ $(\mu, \sigma, \kappa) = (0.2, 1.1, 0.9)$, $C$ $(-0.1, 1.0, 0.55)$; RA: $F$ $(-0.3, 0.65, 0.4)$, $C$ $(0.45, 1.45, 0.05)$; none: $(0, 1, 0.5)$ for both; stress: $F$ $(1, 2.25, 0.5)$, $C$ none.
- **Pregnancy (aims 1 and 3).** Two pregnant and two non-pregnant samples per woman. In cells with a module change, $\rho = 0.7$ in pregnant samples.
- **Activity (aim 2).** Four visits per woman, the first two pregnant; activity $a$ is the woman's level (variance 0.5) plus visit noise (variance 0.5). In cell 6 the IFN factor of each visit is multiplied by $\exp(0.5\,a)$ after the pregnancy change (activity raises its spread, not its mean). In cell 7 the module's correlation at a visit is $0.3 + 0.4\,\Phi(a)$, with $\Phi$ the standard normal distribution function.
- **Cells (10).**

| Cell | Aim | Factor levels | Module | Activity | Women |
|---|---|---|---|---|---|
| 1 | 1 | none | none | none | 18 |
| 2 | 1 | SLE | none | none | 18 |
| 3 | 1 | RA | none | none | 27 |
| 4 | 1 | stress | none | none | 18 |
| 5 | 1 | SLE | +0.4 in pregnancy | none | 18 |
| 6 | 2 | RA | none | raises the spread of $F$ | 36 |
| 7 | 2 | RA | rises with activity | linked to the module | 36 |
| 8 | 3 | SLE against RA | none | none | 18 and 27 |
| 9 | 3 | SLE against RA | +0.4 in SLE pregnancy | none | 18 and 27 |
| 10 | 3 | stress (SLE) against none (RA) | none | none | 18 and 27 |

- **Tests (aims 1 and 3).** Edges are computed on all samples of one disease (raw: `LionessEdges()`; partial: `LionessPartialEdges()` given the IFN score, the mean of genes 1 to 5 standardized over the disease's samples, and composition). For each woman and pair, her mean edge in pregnancy minus her mean edge outside it. Aim 1: a one-sample $t$-test of the changes across women. Aim 3: a Welch two-sample $t$-test of the SLE women's changes against the RA women's, each disease's edges from its own reference.
- **Prediction (cells 1 to 7).** As in `pilot/`: outer 5-fold cross-validation by woman, inner folds by woman, features inside the training fold. Models $R_0$, $R_1$, E1, F1, N1, N1 with $R_0$ (E5 is dropped). Pregnancy stage: logistic, AUC. Activity: linear regression for $R_0$ and $R_1$, Gaussian elastic net (mixing 0.5, `lambda.min` on mean squared error) for the others; the concordance index (the share of pairs of visits with different activity in which the visit with higher activity has the higher prediction, ties in prediction counted one half), with the Pearson correlation, the out-of-fold $R^2$, and the root mean squared error reported beside it.
- **Measures.** Every quantity in section 2 in every cell where it is defined; the shares for raw and partial edges, for factor and module pairs; the AUC or concordance of every model; $R_1 - R_0$, E1 $- R_1$, N1 with $R_0$ $- R_0$.
- **Replicates.** From a timing run of ten replicates per cell with the workers loaded: 200 if the projected run fits one hour on three cores, otherwise the largest multiple of 20 that fits, not below 60; entered in section 6 and committed before the run.
- **Seeds.** Cell $c$ has seed $20861113 + 1000c$; replicate $r$ uses the cell seed plus $r$.
- **Cells with a known answer.** Cell 1 (shares near 0.05 for raw and partial edges; AUC near 0.5). Tests: the partial-correlation edges against the inverse-covariance formula and against Pearson edges with no conditioning; the concordance index on known cases; the folds by woman.
- **Time box.** One hour of compute on three cores.

## 4. What is reported whatever the outcome

Every measure of section 3 in every cell with its MCSE; the verdict table clause by clause; the known-answer table; one figure of the shares (raw against partial) by cell and one of the paired model differences by cell, in `pilot-b/results/`.

## 5. Readings of open cases

- 2026-11-13, before the run. The main activity measure is the concordance index, not the correlation named in proposal section 12: the agent suggested it as the continuous counterpart of the AUC, and the user's next message was "APPROVE"; the correlation is reported beside it.
- 2026-11-13, before the run. Clause 3(a) "in every cell" includes the null and factor-only cells.

## 6. Amendments before the run

- 2026-11-13, after the timing run and before the run. Replicates per cell: 200. The timing run (10 replicates per cell, three workers loaded, the run's own scheduling; `pilot-b/results/timing.rds`) took 52 seconds, which projects 200 replicates to about 1,040 seconds, within the one-hour box.
