# Decisions log: pilot-b

## D0. The literature pass

As for the first pilot (`pilot/decisions.md`, D0). Read since: the GSE235508 counts and sample records for the factor shift (`results/factor_shift.rds`), and one search for other pregnancy datasets with SLE and RA (logged in `log/02-pilot.md`); none was found.

## D1. The status call: `piloted-close-call`, with both readings

Every clause the rule names as deciding passed: the fix (1(b)), the equivalence of edges and their formula (3(a)), and the removal of the activity artifact by the E1 base (3(b)). Two clauses failed and two were close calls; values in `pilot-b/results/verdict.rds`.

- **The reading that passes.** The paper's claims held in every cell: partial-correlation edges brought the share of factor pairs called changed back to its nominal 5% in every cell where only the shared factors moved, within a disease and between diseases, at the real levels and at the stress level, while keeping nearly all of a real change (1(b), 1(c)). Edges equalled their first-order formula for pregnancy stage and for activity (3(a)). The failures are about size at the study's levels: the false SLE-against-RA difference was above its nominal level but below the 0.10 threshold (2(b)), and the false rewiring at the RA level was within 3 MCSE of 0.10 (2(a), close). An artifact at about twice the nominal rate, removed by the fix, is a real-data illustration a methods paper can use.
- **The reading that fails.** At the levels measured in GSE235508 the artifact is modest: none at the SLE level, where the IFN score barely moves and is stable within women (`cell_summaries.rds`, cell 2), about twice the nominal rate at the RA level (cell 3), and about one and a half times between diseases (cell 8). And aim 2 cannot be answered by prediction at 36 women: edges did not detect an activity-linked change in co-expression (3(c) fail), and every model's concordance in that cell was below 0.5.

The agent's reading: the close call is a pass for the paper's methodological claims and a fail for clinical prediction of activity in this dataset. The user decides.

## D2. Why the SLE level shows no false rewiring

At the SLE level the IFN factor's shift (0.2) and spread ratio (1.1) change $\mathrm{E}[F^2]$ by about 0.25 in units of the non-pregnant variance, and its within-woman correlation of 0.9 leaves little within-woman change for the paired test to see; the composition factor barely moves. The RA level moves composition by a shift of 0.45 and a spread ratio of 1.45, a change in $\mathrm{E}[C^2]$ of about 1.3, which is where the false rewiring comes from. This is the theory's prediction (the artifact scales with the change in the squared factor), visible across cells 2, 3, and 4.

## D3. The activity cells: concordance below 0.5

In cell 7 every model's pooled out-of-fold concordance was below 0.5, and the out-of-fold $R^2$ was negative in cells 6 and 7 (`cell_summaries.rds`). Activity has a woman-level part (half its variance); with folds by woman, each held-out woman's level is unknown to the model, and the training mean moves opposite to the held-out women's mean, which biases pooled out-of-fold predictions downward against the outcome. The paired differences the clauses judge share this bias. For stage 3, activity should be analysed within woman (each visit's activity relative to the woman's own mean), which removes the between-woman part; this is a reading for the design stop, not a change to this verdict.

## D4. The artifact in activity prediction

In cell 6 the edges with $R_0$ gained over $R_0$, and $R_1$ (with the squared score) gained more; with the E1 base edges added nothing (3(b) pass). The artifact the first project found in classification appears in a continuous outcome, and the same base removes it.

## D5. Seeds and restart

The machine restarted once during this session, before the pilot-b run started; the run itself completed in one pass (`pilot-b/results/run_info.rds`, `pilot-b/run.log`). The pre-registration and the amendment were committed and pushed before the run.

## D6. What was not done

- No edge was computed on GSE235508; only the factor shift was measured.
- One loading per factor, normal noise, four samples per woman.
- No within-woman activity analysis (D3), no plasmode, no healthy pregnancies; these are stage 3 options.
