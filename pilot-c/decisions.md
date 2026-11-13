# Decisions log: pilot-c

## D0. The literature pass

As for pilot-b (`pilot-b/decisions.md`, D0). New since: the marker-gene composition score (`analysis/02_factor_shift_markers.R`); its marker list carries a `[TO CONFIRM]` marker in proposal section 13 until a published source is opened.

## D1. The status call: `piloted-close-call`, with both readings

Every deciding clause passed again (1(b), 3(a), 3(b)), and with 400 replicates clause 2(c) moved from a close call to a pass. Clauses 2(a), 2(b), and 3(c) failed. Values, with pilot-b's beside them, in `pilot-c/results/verdict.rds`.

- **The reading that passes.** The paper's claims held in every cell, with twice the replicates of pilot-b: partial-correlation edges kept the share of factor pairs called changed at the test's level in every factor-only cell, within and between diseases; they kept a real change in co-expression within a disease and a real SLE-only change between diseases (2(c) now passes); edges equalled their first-order formula in every prediction cell; and the base with the squared score removed the artifact from activity prediction. The failures of 2(a) and 2(b) say that, at the factor changes measured in GSE235508 with marker genes, raw edges call little false rewiring: a property of this dataset, which the theory explains (D2), not a failure of the method.
- **The reading that fails.** The real-data application of aims 1 and 3 will show little or no false rewiring to remove in GSE235508, and aim 2 cannot be answered: edges did not detect an activity-linked change in co-expression, across women (pilot-b) or within women (pilot-c). The pregnancy data then illustrate the fix weakly.

The agent's reading: a pass for the method claims and for detecting real SLE-only change; a fail for a strong real-data illustration of the artifact in this dataset and for activity prediction. The user decides.

## D2. Why the false rewiring is small at the measured levels: the unequal-groups law

With the marker score, pregnancy shifts the composition axis by about 0.6 (SLE) and 1.4 (RA) standard deviations with a spread ratio of about 0.85 (`results/factor_shift_markers.rds`). To first order, a woman's change in edge for two genes loading on a factor $C$ has mean proportional to $(1 - 2\pi)\Delta^2 + (\sigma_1^2 - \sigma_0^2)$, where $\pi$ is the share of pregnant samples, $\Delta$ the shift, and $\sigma_1, \sigma_0$ the spreads: the first project's unequal-groups law (its supplement, S1.5) plus the spread term. A pure mean shift cancels in a balanced design; at the RA level ($\pi = 0.4$) the shift term ($0.2 \times 1.96$) is partly offset by the narrower spread ($0.85^2 - 1$), so the false rewiring is small. Pilot-b's principal-component proxy gave a spread ratio of 1.45 at the RA level, which is where its larger false rewiring came from. This is a prediction the real data can check at stage 3: the artifact should follow the spread change and the imbalance, not the shift alone. The derivation here is a sketch; its numerical check is the comparison of cells 2, 3, 4 across pilots b and c, and a dedicated cell at stage 3.

## D3. Activity within woman

Centering activity and features on each woman removed the downward bias of pilot-b's across-woman analysis (concordances now near 0.5 rather than below it in cell 7), but edges still found no activity-linked change in co-expression (`cell_summaries.rds`, cell 7). A sample's edge is a single noisy product; a change in correlation of at most 0.4 spread over a woman's visits is too small to see in single samples at 36 women. Aim 2 by single-sample edges is not supported at this size.

## D4. Seeds and the run

400 replicates per cell, 2,233 seconds on three workers, one pass (`pilot-c/results/run_info.rds`, `pilot-c/run.log`). The pre-registration and the amendment were committed and pushed before the run.

## D5. What was not done

- No edge was computed on GSE235508.
- The marker list has no opened published source yet.
- Composition observed without error; two factors; normal noise.
- No dedicated cell for the spread term of D2.
