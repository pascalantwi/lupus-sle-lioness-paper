# Pre-registration: `pilot-f/`

| | |
|---|---|
| Proposal | `proposal.md`, section 13 (revision 2, the own-module correction), the user's choice ("option 1", 2026-11-16) |
| Written | 2026-11-16 |
| Committed | in the commit "Stage 2, step 2 (pilot-f): pre-registration", before any pilot-f code or run |

## 1. The three questions

As `pilot-e/preregistration.md`, section 1, with question 3 asked of the own-module correction: do edges built from each gene's residual on its own set's leave-one-gene-out score and the cell mix remove the interferon set's excess connectivity stability, and keep a real personal module's?

## 2. The pass rule

As `pilot-e/preregistration.md`, section 2, unchanged in every clause, quantity, threshold, cell, close-call reading, and combination, except that "partial" in clauses 3(a) and 3(b) means the own-module correction (edge type `ownset`).

| Q | Clause | Quantity | Threshold | Cells |
|---|---|---|---|---|
| 1 | (a) | correlation over samples of raw IFN connectivity with the squared IFN score | $\ge 0.9$ | 1 to 4, each |
| 1 | (b) | $\lvert I_{\text{raw}}(m) - I_{\text{F1}}(m)\rvert$, every metric | $\le \max(0.02, 3\,\text{MCSE})$ | 1 to 4 |
| 2 | (a) | $I_{\text{raw}}(\text{IFN connectivity}) - \bar I_{\text{raw}}(\text{IFN edges})$ | $\ge 0.05$ | 2 |
| 2 | (b) | $I_{\text{raw}}(\text{IFN connectivity}) - I(\text{IFN score}^2)$ | $\lvert v\rvert \le \max(0.05, 3\,\text{MCSE})$ | 2 |
| 3 | (a) | $I_{\text{ownset}}(\text{IFN connectivity}) - I_{\text{ownset}}(\text{cross connectivity})$ | $\lvert v\rvert \le \max(0.05, 3\,\text{MCSE})$ | 2, 3 |
| 3 | (b) | $I_{\text{ownset}}(\text{module connectivity}) - I_{\text{ownset}}(\text{cross connectivity})$ | $\ge 0.05$ | 4 |

Pass when every clause passes; fail when 1(a), 1(b), or 3(a) fails; close call otherwise.

## 3. Design

As `pilot-e/preregistration.md`, section 3, with these changes:

- **Edge types.** Raw LIONESS, first-order, and the own-module correction (`OwnSetResiduals()` then `LionessEdges()`): for each gene of genes 6 to 50, the residual of its regression, over all samples, on an intercept, the mean of the other genes of its set (genes standardized over the samples; sets IFN 6 to 20, module 21 to 30, cell mix 31 to 50), and the observed cell mix. The partial-correlation edges of pilot-e are not recomputed (their result is pilot-e's).
- **Seeds.** Cell $c$ has seed $21261116 + 1000c$.
- **Replicates.** From a timing run of ten replicates per cell with the workers loaded: 200 if the projected run fits one hour on three cores, otherwise the largest multiple of 20 that fits, not below 60.
- **Known answers.** As pilot-e, plus a test that the own-module residuals are uncorrelated over samples with their own regressors.

## 4. What is reported whatever the outcome

As pilot-e, with pilot-e's partial-correlation values beside the own-module ones in the report.

## 5. Readings of open cases

- 2026-11-16, before the run. The regression that removes a set's shared factor may also remove a real personal module's shared component; clause 3(b) is the test of that, and a fail there is reported as such.

## 6. Amendments before the run

- 2026-11-16, after the timing run and before the run. Replicates per cell: 120. The timing run (10 replicates per cell, three workers loaded; `pilot-f/results/timing.rds`) took 262 seconds, which projects 200 replicates to about 5,240 seconds, over the hour, and 120 to about 3,150 seconds, the largest multiple of 20 that fits.
