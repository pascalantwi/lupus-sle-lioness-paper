# Decisions log: pilot-e (network metrics)

## D0. The literature pass

Proposal section 12 and the two searches of 2026-11-16 on LIONESS network metrics and the test-retest reliability of graph metrics. Nothing read since.

## D1. The status call: `piloted-fail`, and what failed

The rule (`pilot-e/preregistration.md`) gives a fail because clause 3(a) failed in both cells: with partial-correlation edges, the interferon set's connectivity stayed more stable than the cross connectivity (`pilot-e/results/verdict.rds`). Every other deciding clause passed: the closed form held (module connectivity correlated about 0.97 with the squared IFN score in every cell), and raw and first-order edges gave the same stability for all 19 metrics in every cell (largest difference under 0.0004). Clause 2(a) and 2(b) passed: connectivity was more stable than single edges, and exactly as stable as the squared IFN score. Clause 3(b) was a close call: the real personal module was visible in connectivity, at about the threshold.

There is no reading under which the rule passes: 3(a) is a deciding clause and missed by far more than its MCSE. What the fail says is specific (D2): the correction, conditioning on a five-gene IFN score, does not remove the IFN set's shared signal once the edges are averaged into a connectivity. What it does not say: the closed forms, the equivalence, and the finding that connectivity is the squared expression score all held.

## D2. Why the correction failed for connectivity

Two mechanisms, both visible in `cell_summaries.rds`.

- **The score measures the factor with error.** The score averages five genes (1 to 5) that carry their own noise and their own person-specific levels. Residuals of genes 6 to 20 on that score keep a shared part of the factor; one edge carries it weakly, but connectivity averages 105 edges and keeps its square. The partial IFN connectivity's stability rose with the factor's stability (about 0.07, 0.14, 0.20 in cells 1 to 3).
- **The score carries person-specific levels.** Even with no stable factor (cell 1), the partial IFN connectivity was more stable than the raw one: regressing every IFN gene on a score that holds the person's own levels of genes 1 to 5 puts a shared person-stable term into all their residuals.

The theory section of the proposal anticipated the first ("the score's measurement error leaves a residual"); its size at the module level was not computed. A correction using the module's own score (each gene adjusted for the mean of the other genes of its set) would address both, but is a post hoc idea and would need its own pre-registration.

## D3. What the metrics showed (reported, not judged)

- Metrics that sum or average edges (strength, eigenvector centrality, efficiency, density, path length) became more stable as the IFN factor became more stable, and the correction brought most of them back to the null cell's level (density, path length, strength).
- Betweenness, local clustering, and hub overlap were barely stable in any cell (intraclass correlations under about 0.08; hub overlap under 0.03).
- The real personal module raised its connectivity's stability from about 0.03 to 0.09, for raw and corrected edges alike: averaging made a personal network visible that single edges hardly showed (pilot-d).

## D4. The run

120 replicates per cell (the largest multiple of 20 that fitted the hour; amendment before the run), 3,442 seconds on three workers, one pass.

## D5. What was not done

- No metric was computed on GSE65391.
- One threshold for the unweighted graph (absolute edge above 1); one distance (inverse absolute edge).
- No prediction of the next visit's activity.
