# Decisions log: pilot-d (personal-network stability)

## D0. The literature pass

The re-scope proposal's prior-work section (`proposal.md`, section 5) and `log/01-idea-review.md` (re-scope section, 19 queries). Nothing read since.

## D1. The status call: `piloted-close-call`, with both readings

Every clause the rule names as deciding passed: edges equal their first-order formula in stability (1(b)), the partial-correlation edges remove the IFN pairs' excess stability (3(a)), and they keep a real personal module (3(b)). Clauses 1(a), 2(a), 2(b), and 3(c) in cell 6 failed. Values in `pilot-d/results/verdict.rds`.

- **The reading that passes.** The paper's three claims held in every cell with 200 replicates: LIONESS edges and their formula have the same stability to within $10^{-4}$; the fix brings the IFN excess to zero; the fix keeps, and sharpens, a real personal module. The failed clauses are about size: the IFN excess was smaller than the thresholds set from the sketch, and edges carried no next-visit activity signal.
- **The reading that fails.** At the IFN stability measured in GSE65391, IFN adds little stability that the cross pairs do not also carry (excess about 0.016), it barely changes identification (raw against fixed, about 0.004), and edges do not predict the next visit's activity even when activity changes co-expression. The application would then show little IFN-driven "personal network" to remove, and no clinical prediction.

The agent's reading: a pass for the method's claims, a fail for a large IFN-driven stability artifact. What the pilot found instead (D2, D3) is a different and possibly more useful message: single-sample edges are unreliable over visits, and what stability they have is mostly the person's own gene levels.

## D2. Why the excess was small: noise and the cross pairs

Two reasons, both visible in `cell_summaries.rds`.

- **Edges are noisy.** Every edge's intraclass correlation was between about 0.07 and 0.22: a single-visit edge is dominated by the product of the visit's noise terms, so the part that comes from $F^2$ is a small share of its variance. The sketch's $\kappa^2$ multiplies that share; the thresholds (0.15, 0.05) were set without computing it, which was the agent's error.
- **The cross pairs carry IFN too.** An edge between an IFN gene and a module gene contains the product of the IFN factor and the person's own level of the module gene, which is stable with intraclass correlation $\kappa_F$. The cross pairs' stability rose from about 0.075 (no stable IFN) to 0.128 (measured IFN stability), so subtracting them removes part of the IFN effect. Not judged, and reported: the IFN pairs' own stability about doubled, from about 0.069 to 0.144, between those two cells.

The rule is not changed after the run; both points go to stage 3's design if the user goes on.

## D3. Edges are barely personal, and mostly through gene levels

Fingerprinting on all 990 pairs identified about a quarter of the people from one other visit (chance under 1%), with or without the fix, and slightly more after it; on the IFN or module pairs alone, about 5% to 8%. The person's own gene levels (variance 0.5 in the simulation, as people differ in level) give every edge a stable part, which is why cross and module pairs are as stable as IFN pairs. A real personal module with person-to-person spread 0.25 raised its pairs' stability by only about 0.01. So an edge's "personal" content is mostly the person's gene levels, which the genes themselves carry.

## D4. Next-visit activity: no signal

Every model's concordance sat near 0.48 in cells 1, 2, and 6, below 0.5 for the reason of pilot-b's D3 (activity has a person-level part that grouped folds cannot learn). Edges added nothing over the genes (about -0.0006 in cell 6). A single-visit edge cannot carry a change of at most 0.2 in correlation at that visit.

## D5. The run

200 replicates per cell, 2,714 seconds on three workers, one pass (`run_info.rds`, `run.log`). The pre-registration, the change of rule, the user's acceptance, and both amendments were committed and pushed before the run.

## D6. What was not done

- No edge was computed on GSE65391 (only its description).
- One reference choice (the whole cohort); one noise level; six visits.
- The within-person version of the prediction.
