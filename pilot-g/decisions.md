# Decisions log: pilot-g (the pathway approach)

## D0. The literature pass

As pilot-e. Nothing read since. The gene-set collections for the application (BloodGen3, MSigDB Hallmark, Reactome) and the eigengene literature (Horvath and Dong 2008) are still to be opened at stage 3; no claim here rests on them.

## D1. The status call: `piloted-fail`

One deciding clause failed in one cell: clause 3(a) in cell 3, where interferon is very stable (0.9). There, the residual of the interferon pathway's connectivity, after its three expression terms, had stability 0.23 against the bound of 0.05 (`pilot-g/results/verdict.rds`). The other 17 residual checks passed, including all six pathways at the measured interferon stability (cell 2, largest 0.033) and the null cell. Every other clause passed:

- 1(a) and 1(b): every pathway's connectivity and every pair's crosstalk were linear in their expression terms, with $R^2$ above 0.99999 in every cell (the share left over was about 2 parts per million).
- 1(c): interferon connectivity correlated 0.965 with the squared pathway score in every cell.
- 2(a): at the measured level, interferon connectivity was as stable as the squared score (difference -0.012).
- 2(b): it was more stable than the pathway's single edges by 0.114.
- 3(b): the personal pathway's connectivity was more stable than the null pathway's by 0.084.

Under the combination rule a failed 3(a) is a fail, so the status is `piloted-fail`. There is no reading of the rule under which it passes.

## D2. Why the one clause failed: the residual is the next-order expression term

After the run, not pre-registered, the agent checked five replicates of cell 3 (`log/02-pilot-g-check/residual_order.R` and its output). The residual is about $2 \times 10^{-6}$ of connectivity's variance. Its stability, 0.18 on average in those replicates, fell to 0.009 when fourth-order terms of the pathway's own standardized expression were added to the regression ($S^4$, $S^2 Q$, $Q^2$, and $\sum z^4$). The section 14.1 closed form is the first-order expansion of the edge. The exact LIONESS edge also carries higher powers of the sample's own expression, and when the factor is very stable those powers are stable too.

So the failed clause does not show a network part beyond expression. It shows that the pre-registered terms stop at the first order, and that a stability measure can pick up a remainder of negligible size. The size is no guide: the intraclass correlation does not depend on scale, so a remainder of 2 parts per million can still be stable. The pre-registered rule did not foresee this. The agent's reading in section 5 expected the personal pathway to be the only risk.

The check was run after the result was known. It is a reading, not a rescue: the status stays `piloted-fail`. It makes a testable claim for stage 3's protocol: connectivity is a function of the pathway's own expression to every order, and its residual on second- and fourth-order terms has stability near 0 at any factor stability.

## D3. What the pathway view adds

- **Pathway metrics are expression.** Connectivity was the squared pathway score (correlation 0.965), and its stability matched the score's in every pathway at the measured level (for example neutrophil-like 0.139 against 0.142, lymphocyte-like 0.184 against 0.190). Crosstalk was linear in the product of the two pathway scores and their weighted squares.
- **Averaging over a pathway steadies the network.** Single interferon edges had stability 0.146; the pathway's connectivity, 0.260.
- **Personal co-expression enters through expression.** The personal pathway's connectivity was steadier than the null pathway's (0.103 against 0.019), but its residual on the terms was not (0.011). Its squared score was also steadier (0.060 against 0.008). A person whose genes co-vary more has a more variable pathway score, so the squared score carries the personal correlation. The network quantity sees the personal co-expression only through the expression terms.
- The null cell's residual stabilities were small but above 0 (0.001 to 0.006, about 1 to 3.5 MCSE). The person's own gene levels are stable in every cell, and their higher powers enter the remainder as in D2.

## D4. The run

200 replicates per cell (the amendment before the run), 98 seconds on three workers, one pass.

## D5. What was not done

- Prediction of next-visit activity from pathway metrics (pre-registration, section 3).
- Real gene sets.
- Any edge on GSE65391.
- A pathway approach with partial correlations or with an external correction. Pilots e and f showed that within-set adjustment cannot change connectivity (pilot-f, D2).
