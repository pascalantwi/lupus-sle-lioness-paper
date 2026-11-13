# Pre-registration: `pilot/`

| | |
|---|---|
| Proposal | `proposal.md`, section 7 |
| Written | 2026-11-13 |
| Committed | in the commit "Stage 2, step 2: pre-registration", before any pilot code or run |

## 1. The three questions

1. Does the effect exist, and do LIONESS edges equal their formula over time? When co-expression changes in pregnancy, do edges gain over E1, and do LIONESS edges and first-order edges added to E1 give the same out-of-sample AUC in every cell?
2. Is it large enough to matter at the study's size? At 20 women, do edges gain over a base without squares when only the shared factor changes (the artifact), and do they still detect a real change in co-expression?
3. Is there a fix? When only the shared factor changes, does the E1 base (the score and composition, their squares, and the genes) leave edges nothing to add?

## 2. The pass rule

$A(M)$ is the pooled out-of-fold AUC of model $M$ in one replicate. A difference is paired within a replicate and averaged over replicates; its MCSE is the replicate standard deviation of the difference divided by the square root of the number of replicates. $\Delta_1 = A(\text{N1}) - A(\text{E1})$, $\Delta_2 = A(\text{N1}) - A(\text{F1})$, $G = A(\text{N1 with } R_0) - A(R_0)$.

| Question | Clause | Quantity | Threshold | Cells | Within its MCSE |
|---|---|---|---|---|---|
| 1 | (a) | mean of $\Delta_1$ | $\ge 0.05$ and $> 3$ MCSE | 7 | $\lvert v - 0.05 \rvert \le 3$ MCSE: close call |
| 1 | (b) | $D_c$ = mean of $\Delta_2$ | $\lvert D_c \rvert \le \max(0.01, 3\,\text{MCSE}_c)$ in every cell | 1 to 8 | $\max(0.01, 3\,\text{MCSE}_c) < \lvert D_c \rvert \le 0.01 + 3\,\text{MCSE}_c$: close call |
| 2 | (a) | mean of $G$ | $\ge 0.02$ | 2 | $\lvert v - 0.02 \rvert \le 3$ MCSE: close call |
| 2 | (b) | mean of $\Delta_1$ | $\ge 0.03$ | 5 | $\lvert v - 0.03 \rvert \le 3$ MCSE: close call |
| 3 | (a) | $E_c$ = mean of $\Delta_1$ | $\lvert E_c \rvert \le \max(0.01, 3\,\text{MCSE}_c)$ in each cell | 2, 3, 4 | as clause 1(b) |

**How the answers combine.** Pass when every clause passes. Fail when 1(a), 1(b), or 3(a) fails. Close call otherwise: when 2(a) or 2(b) fails, or when any clause is a close call for itself. A close call is stated both ways at the verdict.

## 3. Design

- **Data-generating mechanism** (`SimulatePregnancyCohort()`). $W$ women, each with four samples: two pregnant ($t = 1$) and two non-pregnant ($t = 0$). Fifty genes. For sample $s$ of woman $w$:
  - the shared factor: $F^*_{ws} = a_w + e_{ws}$ with $a_w \sim N(0, 0.5)$ and $e_{ws} \sim N(0, 0.5)$ (unit variance, within-woman correlation 0.5); $F_{ws} = F^*_{ws}$ when $t = 0$, and $F_{ws} = \mu + \sigma F^*_{ws}$ when $t = 1$, with $(\mu, \sigma)$ = (0, 1) unless the factor changes, (0.5, 1.5) moderate, (1, 2.25) large;
  - composition: $C_{ws} = b_w + c_{ws}$ with $b_w \sim N(0, 0.8)$ and $c_{ws} \sim N(0, 0.2)$, the same in both periods, observed without error;
  - genes 1 to 20 load 0.8 on $F$; genes 31 to 50 load 0.4 on $C$; every gene has the woman's offset $u_{wg} \sim N(0, 0.5)$ and independent noise $N(0, 1)$, except that the noise of genes 21 to 30 (the module) has correlation 0.3, and $0.3 + \Delta$ in pregnant samples when the module changes.
- **Cells (8).** $W$ women, $4W$ samples.

| Cell | What pregnancy changes | Factor | Module change $\Delta$ | $W$ |
|---|---|---|---|---|
| 1 | nothing (null) | none | 0 | 20 |
| 2 | factor | large | 0 | 20 |
| 3 | factor | moderate | 0 | 20 |
| 4 | factor | large | 0 | 40 |
| 5 | module | none | 0.4 | 20 |
| 6 | module | none | 0.2 | 20 |
| 7 | module | none | 0.4 | 40 |
| 8 | factor and module | large | 0.4 | 40 |

- **Derived variables, inside each training fold.** As in the user's first project: the IFN score is the mean of genes 1 to 5 standardized by the training fold, and its residual on composition uses the training regression; the gene set $G$ is genes 6 to 50 (990 pairs); products, squares, and first-order edges use the training fold's correlation and standardization for training and test samples; LIONESS edges of training samples use the training fold as reference (exact closed form); test samples are scored by the add-one form.
- **Methods.** Logistic regression: $R_0$ (composition, residual IFN score), $R_1$ ($R_0$ with both squares). Elastic net (`cv.glmnet`, mixing 0.5, 5 inner folds with each woman's samples in one fold, binomial deviance, `lambda.min`, the base terms with penalty factor 0): E1 (the 45 genes, $R_1$ base), E5 (the genes and their 1,035 products and squares, $R_1$ base), F1 (the genes and the 990 first-order edges, $R_1$ base), N1 (the genes and the 990 LIONESS edges, $R_1$ base), N1 with $R_0$ (the 990 LIONESS edges, $R_0$ base). Outer 5-fold cross-validation with each woman's four samples in one fold.
- **Measures.** $A(M)$ for the seven models; $\Delta_1$, $\Delta_2$, $G$; and, judged in no clause: $A(R_1) - A(R_0)$, $A(\text{E5}) - A(\text{E1})$, $A(\text{N1}) - A(\text{E5})$; apparent rewiring: LIONESS edges on all $4W$ samples, each woman's mean edge in pregnancy minus her mean edge outside it, a one-sample $t$-test across women for each pair, and the share of pairs with $p < 0.05$ among the 105 pairs of genes 6 to 20 (factor genes outside the score) and among the 45 module pairs.
- **Replicates per cell.** From a timing run of ten replicates per cell with the workers loaded: 200 if the projected run fits one hour on three cores, otherwise the largest multiple of 20 that fits, not below 60; entered in section 6 and committed before the run.
- **Seeds.** Cell $c$ has seed $20761113 + 1000c$; replicate $r$ uses the cell seed plus $r$.
- **Cells with a known answer.** Cell 1 (every AUC near 0.5; apparent rewiring near 0.05). Tests: the exact training edges against literal leave-one-out; the add-one test edges against literal addition; the first-order remainder shrinking with $N$; the folds keeping each woman's samples together.
- **Time box.** One hour of compute on three cores.

## 4. What is reported whatever the outcome

The AUC of every model in every cell with its MCSE; every quantity of section 3, "Measures", in every cell; the verdict table clause by clause; the known-answer table; one figure of the differences by cell (`pilot/results/auc_gains.png`); the data description of GSE235508 (`pilot/results/data_summary.rds`).

## 5. Readings of open cases

- 2026-11-13, before the run. Clause 1(b) "in every cell" includes the null and factor-only cells, where N1 and F1 may both add nothing; clause 1(a) guards against passing on two nulls.
- 2026-11-13, before the run. Clause 3(a) is two-sided as written.

## 6. Amendments before the run

- 2026-11-13, after the timing run and before the run. Replicates per cell: 200. The timing run (10 replicates per cell, three workers loaded, the run's own scheduling; `pilot/results/timing.rds`) took 84 seconds, which projects 200 replicates to about 1,700 seconds, within the one-hour box.
