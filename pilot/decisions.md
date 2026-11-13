# Decisions log: pilot

## D0. The literature pass

The proposal's prior-work section (`proposal.md`, section 5) and `log/01-idea-review.md` hold the literature read at stage 1: the mechanism's sources carried from the first project, the dataset paper (Lien and colleagues 2024), and three new searches, none of which found a single-sample network analysis of pregnancy. Nothing was read since.

## D1. The status call: `piloted-close-call`, with both readings

The rule (`pilot/preregistration.md`, section 2) gives a close call: no deciding clause failed (1(b) and 3(a) passed, 1(a) is within 3 MCSE of its threshold), and both clauses of question 2 failed. The values are in `pilot/results/verdict.rds`.

- **The reading that passes.** The claims the critique rests on held: LIONESS edges and their first-order formula gave the same AUC in every cell (1(b)), and with the E1 base the edges added nothing when only the shared factor changed (3(a)). Clause 1(a) fell short of 0.05 by about 1.7 MCSE, a value the rule reads as a close call, not a fail. The failures of question 2 concern AUC at 20 women; the artifact the critique is about is large at 20 women in the analysis that longitudinal users run, the paired test of edges between periods: with no change in co-expression, the share of factor-gene pairs called rewired was about nine times its nominal 5% in the large factor-change cell (`pilot/results/cell_summaries.rds`, `rewire_factor`, cell 2), and about three times in the moderate cell. Read this way, the project goes on with the paired edge test as the primary measure of the artifact.
- **The reading that fails.** At the study's size (20 women), neither the artifact nor a real change in co-expression moved the AUC enough to matter: the artifact was about a quarter of its 0.02 threshold, and N1 gained under a third of 0.03 over E1 in the rewiring cell. The classification analysis that the proposal planned for GSE235508 (18 SLE women, 27 RA women with samples in both periods) would then be uninformative in both directions: a null result would not separate "no rewiring" from "too few women". Read this way, the classification aim is dropped or the application is re-scoped before stage 3.

The agent's reading: the deciding claims held, and the result that matters most for the paper is outside the pass rule (the apparent rewiring), so the close call is closer to a pass for a critique built on the paired edge test, and closer to a fail for the classification aim. The user decides.

## D2. Folds by woman, inner and outer

Every woman's four samples are in one outer fold, and the inner folds of `cv.glmnet` are also by woman (`GroupFolds()`). Without this, a woman's offsets in the training fold would score her own test samples. Alternative: stratified folds by sample, as in the first project; rejected, since it leaks the woman.

## D3. The artifact grows with the number of women

The artifact $A(\text{N1 with } R_0) - A(R_0)$ was small at 20 women and about four times larger at 40 in the same large factor-change setting (cells 2 and 4, `cell_summaries.rds`). The elastic net needs enough samples to use 990 edges; at 80 samples it uses them little. So the AUC artifact is a function of size as well as of the factor change, and question 2's threshold, set at 20 women, judged it where it is weakest.

## D4. Real rewiring is hard to see by classification at 20 women

In the rewiring cell at 20 women (cell 5) N1 gained under 0.01 over E1, and at 40 women (cell 7) about 0.04; with a change of 0.2 at 20 women (cell 6) nothing. The paired edge test found the module change in about a quarter of the module pairs at 20 women and in half at 40. The woman offsets (variance 0.5 on every gene) dilute the module's correlation in the pooled reference, which is the realistic case for repeated samples.

## D5. The paired edge test is fooled in proportion to the factor change, and the fix for it was not tested

The share of factor-gene pairs called rewired with no change in co-expression grew from about the nominal rate (null) to about 0.14 (moderate), 0.45 (large, 20 women), and 0.82 (large, 40 women). No fix for the paired test (an edge adjusted for the squared score, or partial-correlation edges) was in the pre-registration, so none was run; that is the obvious content of a re-pilot if the user re-scopes to it.

## D6. The null cell

Every model's AUC in the null cell sat between about 0.49 and 0.50, within about 2 MCSE of 0.5 except E1 and R1 (about 1.9 MCSE below), the usual small downward bias of pooled out-of-fold AUC; the apparent rewiring sat at its nominal 0.05 (`known_answers.rds`). The paired differences the clauses judge are unaffected, since the bias is shared.

## D7. What was not done

- No gene expression from GSE235508 was analysed; only its sample records were read.
- The counts file (25 MB) was not downloaded; it waits for the user's approval.
- One noise distribution (normal), one loading (0.8), four samples per woman.
- No fix for the paired edge test was run (D5).
