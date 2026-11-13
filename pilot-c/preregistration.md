# Pre-registration: `pilot-c/`

| | |
|---|---|
| Proposal | `proposal.md`, section 13 (revision 2), on section 12; the user asked to go straight on ("go to C", 2026-11-13) |
| Written | 2026-11-13 |
| Committed | in the commit "Stage 2, step 2 (pilot-c): pre-registration", before any pilot-c code or run |

## 1. The three questions

As `pilot-b/preregistration.md`, section 1, unchanged.

## 2. The pass rule

As `pilot-b/preregistration.md`, section 2, unchanged: the same nine clauses, thresholds, cells, close-call readings, and combination (fail when 1(b), 3(a), or 3(b) fails; pass when every clause passes; close call otherwise). In clauses 3(b) and 3(c) the concordance is that of the within-woman analysis (section 3).

## 3. Design

As `pilot-b/preregistration.md`, section 3, with these changes:

- **Visits.** Five samples per woman: three pregnant and two non-pregnant at the SLE levels (and in the null and stress cells), two pregnant and three non-pregnant at the RA levels; in aim 3 cells, the SLE cohort has the SLE pattern and the RA cohort the RA pattern. Activity cells (6, 7) use the RA pattern.
- **Factor levels** (`results/factor_shift_markers.rds` for composition, `results/factor_shift.rds` for IFN, rounded): SLE: $F$ $(0.2, 1.1, 0.9)$, $C$ $(0.6, 0.85, 0.8)$; RA: $F$ $(-0.3, 0.65, 0.4)$, $C$ $(1.4, 0.85, 0.1)$; none and stress as in pilot-b.
- **Activity within woman** (cells 6 and 7). Inside each fold, the activity and every feature of every model ($R_0$ and $R_1$ terms, genes, first-order edges, LIONESS edges) are centered on the woman's own mean over her samples, training and test women alike (a woman's features use only her own samples; the test women's activity is used only to score). Models are fitted to the centered activity and scored by the concordance index against it, with the correlation, $R^2$, and root mean squared error reported.
- **Replicates.** From a timing run of ten replicates per cell with the workers loaded: 400 if the projected run fits one hour on three cores, otherwise 200 if it fits, otherwise the largest multiple of 20 that fits, not below 60; entered in section 6 and committed before the run.
- **Seeds.** Cell $c$ has seed $20961113 + 1000c$; replicate $r$ uses the cell seed plus $r$.

## 4. What is reported whatever the outcome

As pilot-b, in `pilot-c/results/`, with pilot-b's values beside pilot-c's for every clause.

## 5. Readings of open cases

- 2026-11-13, before the run. Clause 2(a) "cell 2 or cell 3" is judged on the larger of the two, as in pilot-b.

## 6. Amendments before the run

- 2026-11-13, after the timing run and before the run. Replicates per cell: 400. The timing run (10 replicates per cell, three workers loaded, the run's own scheduling; `pilot-c/results/timing.rds`) took 60 seconds, which projects 400 replicates to about 2,400 seconds, within the one-hour box.
