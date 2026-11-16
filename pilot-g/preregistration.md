# Pre-registration: `pilot-g/`

| | |
|---|---|
| Proposal | `proposal.md`, section 14 (revision 3, a pathway approach), the user's "do it" (2026-11-16) |
| Written | 2026-11-16 |
| Committed | in the commit "Stage 2, step 2 (pilot-g): pre-registration", before any pilot-g code or run |

## 1. The three questions

1. Do the closed forms hold? Is each sample's LIONESS pathway connectivity and crosstalk a linear function of its expression terms (proposal section 14.1), and is connectivity close to the squared pathway score?
2. Are pathway metrics more stable than single edges, and only as stable as the expression score? At the measured interferon stability, does the interferon pathway's connectivity match the squared score's stability, and exceed the stability of its single edges?
3. Is anything network-specific left? After the expression terms, does any pathway's connectivity keep stability over visits, including the pathway with real personal co-expression; and is that personal co-expression visible in connectivity through the terms?

## 2. The pass rule

For each replicate, LIONESS edges (and first-order edges) over the 90 pathway genes are computed on all samples; each sample's within-pathway connectivity (6 pathways) and crosstalk (15 pairs) are mean edges. $R^2_{\text{terms}}$ is the $R^2$ over samples of the regression of a pathway metric on its three expression terms (section 14.1); the residual is that regression's residual. $I(\cdot)$ is the intraclass correlation over visits (`EdgeICC()`). Means over replicates, MCSE as before.

| Q | Clause | Quantity | Threshold | Cells | Within its MCSE |
|---|---|---|---|---|---|
| 1 | (a) | smallest $R^2_{\text{terms}}$ of within-pathway connectivity over the 6 pathways | $\ge 0.99$ | 1 to 3, each | $\lvert v - 0.99\rvert \le 3$ MCSE: close call |
| 1 | (b) | smallest $R^2_{\text{terms}}$ of crosstalk over the 15 pairs | $\ge 0.99$ | 1 to 3, each | as 1(a) |
| 1 | (c) | correlation over samples of interferon-pathway connectivity with its squared score | $\ge 0.9$ | 1 to 3, each | $\lvert v - 0.9\rvert \le 3$ MCSE: close call |
| 2 | (a) | $I(\text{IFN connectivity}) - I(\text{IFN score}^2)$ | $\lvert v\rvert \le \max(0.05, 3\,\text{MCSE})$ | 2 | $\max(0.05, 3\,\text{MCSE}) < \lvert v\rvert \le 0.05 + 3\,\text{MCSE}$: close call |
| 2 | (b) | $I(\text{IFN connectivity}) - \bar I(\text{IFN edges})$ | $\ge 0.05$ | 2 | $\lvert v - 0.05\rvert \le 3$ MCSE: close call |
| 3 | (a) | $I(\text{residual connectivity})$ of every pathway | $\lvert v\rvert \le \max(0.05, 3\,\text{MCSE})$ | 1 to 3 | as 2(a) |
| 3 | (b) | $I(\text{personal pathway connectivity}) - I(\text{null pathway connectivity})$ | $\ge 0.05$ | 2 | as 2(b) |

**How the answers combine.** Pass when every clause passes. Fail when 1(a), 1(b), or 3(a) fails. Close call otherwise.

## 3. Design

- **Mechanism** (`SimulatePathwayCohort()`). 131 people, six visits. Ninety genes in six pathways of 15 genes: P1 interferon (factor loading 0.8, within-person stability $\kappa_1$), P2 neutrophil-like (loading 0.6, stability 0.41), P3 lymphocyte-like (0.6, 0.48), P4 plasmablast-like (0.6, stability 0.39, the activity's), P5 personal co-expression (no factor; noise correlation $0.3 + 0.25\,s_w$, with $s_w$ uniform on $(-1, 1)$ fixed for the person), P6 null (no factor; noise correlation 0.3). Every gene carries the person's own level (variance 0.5) and unit noise. The factors are independent.
- **Cells (3).** 1: every factor stability 0 and P5 not personal (null); 2: the measured stabilities ($\kappa_1$ = 0.55); 3: as 2 with $\kappa_1$ = 0.9.
- **Metrics.** Within-pathway connectivity and crosstalk (mean edges); the three expression terms of each; the squared pathway score (mean standardized expression); the mean intraclass correlation of the single edges within each pathway.
- **Replicates.** From a timing run of ten replicates per cell with the workers loaded: 200 if the projected run fits one hour on three cores, otherwise the largest multiple of 20 that fits, not below 60.
- **Seeds.** Cell $c$ has seed $21361116 + 1000c$.
- **Known answers.** Test: first-order connectivity and crosstalk equal their closed forms exactly on random data. Cell 1: residual stabilities near 0.
- **Not tested here.** Prediction; real gene sets; GSE65391 edges (stage 3, after the protocol).
- **Time box.** One hour on three cores.

## 4. What is reported whatever the outcome

Every pathway's and pair's $R^2_{\text{terms}}$, stabilities of connectivity, score, and residual, the verdict table, the known-answer table, and one figure, in `pilot-g/results/`.

## 5. Readings of open cases

- 2026-11-16, before the run. Clause 3(a) holds the personal pathway P5 to the same bound: the theory predicts its personal co-expression appears in the expression terms, not in the residual.

## 6. Amendments before the run

