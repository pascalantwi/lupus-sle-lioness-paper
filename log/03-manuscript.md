# Drafting log (stage 3)

Status: LIONESS edges, network metrics, and pathway connectivity and crosstalk are, to first order and in practice, fixed quadratic functions of the sample's own expression, so over repeated visits they describe, explain, and predict as the expression scores do; there is no correction (pilot-f) and no fix, only a diagnostic (compare with the expression terms); prediction with pathway metrics and real gene sets has not been piloted.

Workflow 0.2.4-md. Started 2026-11-17.

## Start of the stage

- The user (2026-11-17): "now stage 3 and the graphs incase you can add the lioness diagram for the most if you ploting do not combine any plot. do them individually". Read as the user's word to start stage 3 on pilot-g's `piloted-fail`, accepted with the reading of `pilot-g/decisions.md`, D2 (the one failed clause is a higher-order expression term), and on pilot-d's close call read as a pass (D-S6). The precondition `2, piloted-pass` is not met; the user's decision replaces it, as at D-S6.
- Figure rule from the user, for every deliverable of this project: one plot per figure, no panels combined; LIONESS network diagrams wherever a figure can show one.
- The supervisor's confirmation of the re-scope and of the pathway message is still open; it stays a `[TO CONFIRM]` in the draft.

## Paper test

| Question | Answer | Where |
|---|---|---|
| The procedure | LIONESS single-sample co-expression networks (Pearson), summarized per sample as edges, network metrics, and pathway connectivity (mean edge within a gene set) and crosstalk (mean edge between two sets) | `R/lioness.R`, `R/metrics.R`, `R/pathways.R`; proposal sections 12 to 14 |
| The setting | Repeated blood samples of the same people (GSE65391, 131 children, 897 visits), where the summaries are read as a personal network that tracks disease | `pilot-d/results/data_summary.rds`; proposal section 3 |
| The estimand as a functional | What a user takes a sample's connectivity to estimate: the co-expression of a pathway's genes specific to that sample, a functional of the sample's own joint distribution; what LIONESS returns is a functional of the cohort's correlation matrix and the sample's standardized expression vector | proposal sections 4 and 14.1 |
| The target and the bias term | Connectivity is $\bar r_A + (S_A^2 - Q_A)/(2C_A) - \sum_i R_i^A z_i^2/(2C_A)$ to first order, crosstalk $\bar r_{AB} + S_A S_B/(m_A m_B) - \dots$: the squared pathway score and the product of two scores, so its stability over visits is the expression's | proposal section 14.1; `KnownAnswerChecksG()` |
| The simulation that confirms it | Pilots d to g: edges and metrics equal their first-order form in stability; pathway metrics linear in their terms ($R^2 > 0.99999$); the residual stable only for a very stable, strongly loaded factor, and then a higher-order expression term | `pilot-g/results/`, `log/02-pilot-g-check/` |
| The correction or diagnostic | No correction: within-set adjustment leaves connectivity unchanged (pilot-f, D2), an external score leaves its error (pilot-e). Diagnostic: regress a network summary on its expression terms and report the residual's stability; compare prediction from networks, scores, and the formula | `pilot-f/decisions.md`; proposal section 14 |
| The dataset and what to show | GSE65391: (1) describe how pathway connectivity and crosstalk change over each child's visits; (2) explain them by the expression terms; (3) predict identity and same- and next-visit SLEDAI from networks, scores, and the formula, folds by child and the authors' 91/40 split | proposal section 14; spec section 5 |

## Step 1: the pilot promoted

- `pilot-g/pilot.R` promoted to `analysis/03_simulation.R` (parameter block, a seed per cell, seeded blocks of 20 replicates on three PSOCK workers with the BLAS pinned, `results/run_info.rds`). Named 03, not 01: `analysis/01_factor_shift.R` and `02_factor_shift_markers.R` are the pregnancy scripts, kept as the record. Until the design stop it carries pilot-g's design.
- The script also saves `results/example_networks.rds` (`ExampleNetworks()`): one cohort's pathway metrics and the full LIONESS edges of every visit of two children, for the network diagrams.
- Figure functions in `R/figures.R`, one plot each (the user's rule): `PlotSampleNetwork()` (one sample's LIONESS gene network, genes grouped by pathway), `PlotPathwayNetwork()` (one sample's pathway network: node area connectivity, line width crosstalk), `PlotConnectivityScore()`, `PlotTrajectories()`, `PlotStabilityByPathway()` (one setting per plot). Pathway colours: the first six slots of a categorical palette checked with a colour-vision validator (all pass; three below 3:1 contrast on white, so every pathway is labelled directly).
- `PathwayMetrics()` now returns a one-row matrix for one sample (it returned a vector); found by the figure test.
- Acceptance test `tests/testthat/test-acceptance.R`, gated by `LIONESSPREG_ACCEPTANCE`: with pilot-g's grid and seeds, the package reproduced every replicate (600) and every per-cell summary of `pilot-g/results/` to 1e-12 (`full`, 2026-11-17).
- Tests: `test-pathways.R`, `test-figures.R`; `devtools::check()` 0 errors, 0 warnings, 0 notes.
- Assumptions carried from the pilot's code and decisions: disjoint pathways of equal size; independent factors; person gene levels independent across genes (variance 0.5); normal noise; six equal visits; the first-order terms only (pilot-g D2).

## Step 2: references

Inventory (new for the re-scope; the stage-1 entries are in `references/references.md` and `manuscript/citations.csv`):

| Key | Grade | Section | Why | Origin |
|---|---|---|---|---|
| kuijjer2019lioness | core | intro, methods | LIONESS | stage 1 |
| kuijjer2019lionessr | supporting | intro | Pearson LIONESS in practice | stage 1 |
| yu2017iena | core | intro | single-sample correlation as a product of deviations | stage 1 |
| banchereau2016personalized | supporting | data | GSE65391 | stage 1 (re-scope) |
| petri2019stability | supporting | intro | interferon signature stable within patients | stage 1 (re-scope) |
| finn2015fingerprinting | supporting | methods | fingerprinting accuracy | stage 1 (re-scope) |
| whitney2003individuality, mishra2026temporal | supporting | intro | blood expression is person-specific | stage 1 (re-scope) |
| farahbod2020composition | supporting | intro | cell composition drives co-expression | stage 1 (re-scope) |
| friedman2010glmnet | definitional | methods | elastic net | stage 1 |
| altman2021bloodgen3 | core | application | the gene sets | stage 3 |
| rinchai2021bloodgen3module | definitional | methods | module lists in R | stage 3 |
| liberzon2015hallmark | definitional | methods | fallback gene sets | stage 3 |
| horvath2008geometric | core | intro, discussion | the closest prior result: connectivity and the eigengene in a network built from all samples | stage 3 |
| langfelder2008wgcna | supporting | intro | module eigengenes | stage 3 |

1. Available (16 of the 17 inventory keys): all but devlin1975robust, which is not cited in the re-scope.
2. Core with full text pending: none new. devlin1975robust stays pending from stage 1.
3. Records that could not be opened: none. Two pages could not be reached (www.nature.com and www.gsea-msigdb.org are outside the session's network); the BloodGen3 article was read from its open Europe PMC copy, and the MSigDB data licence is unconfirmed (msigdbr on CRAN is MIT; the data's own terms were not read).
- Gene-set licences: BloodGen3 article CC BY 4.0, module lists in the Bioconductor package BloodGen3Module (GPL-2, version 1.20.0); BloodGen3 is the primary collection. Hallmark is the fallback, its data licence `[TO CONFIRM]`.
- `build-bib.py --check`: 26 entries, one fault (devlin1975robust's full text, not cited here).

## Step 3.1: purpose draft of the Introduction

Written new in `manuscript/manuscript.qmd` (nothing carried from the pregnancy draft): paragraph 1 the procedure and the lionessR example (passage added to `citations.csv`); paragraph 2 the first-order edge, the squared pathway score, and Horvath and Dong's cohort-network result as the closest prior work; paragraph 3 the three questions. Title is the user's; journal BMC Bioinformatics (spec section 6), `biomed-central.csl` copied beside the qmd.

## Step 3.2: design proposal (stop)

The study has two parts. S1 (describe and explain) varies, one at a time around the measured setting, what the pilots showed controls a pathway summary's stability and its remainder beyond the first-order terms. S2 (test) adds a disease-activity outcome and compares prediction from network summaries, pathway scores, and the terms of the closed forms, with type I error and power. Fixed: 131 children; six pathways; the other factor stabilities as measured in GSE65391 (neutrophil 0.41, lymphocyte 0.48, activity 0.39); person gene levels of variance 0.5; normal noise. Starred levels are the base setting.

| Part | Factor | Levels | Reason | Question |
|---|---|---|---|---|
| S1 | Interferon factor stability over visits | 0, 0.25, 0.55\*, 0.75, 0.9 | none to very stable; 0.55 is the IFN5 score's measured value; the residual's stability rose with it (pilot-g) | 1, 2 |
| S1 | Interferon pathway loading | 0.4, 0.6, 0.8\* | a weak to a strong driving factor; the strength sets the higher-order remainder (pilot-g check after the run) | 2 |
| S1 | Interferon pathway size | 5, 15\*, 40 | small to large modules; the real range from the BloodGen3 list `[TO CONFIRM]` | 2 |
| S1 | Personal co-expression spread | 0, 0.1, 0.25\* | none to strong; whether co-expression shows beyond the expression terms | 2 |
| S1 | Visits per child | 3, 6\* | GSE65391 median 6; fewer visits widen the stability estimate | 1 |
| S1 | A null setting | every factor stability 0, no personal co-expression | the known answer: residual stability 0 | 1, 2 |
| S2 | Activity effect on the interferon and plasmablast factors | 0, then four levels calibrated with low-replicate cells from low power to about 80% | type I error and power of the prediction tests; GSE65391 placed by its measured association of IFN5 with SLEDAI | 3 |

- Methods compared. S1: raw and first-order LIONESS; the diagnostic with first-order terms and with fourth-order terms added. S2: A network summaries (6 connectivities and 15 crosstalks), B pathway scores and their squares, C the terms of the closed forms, R0 no expression features (the current SLEDAI for next-visit prediction); Gaussian elastic net, folds by child.
- Measures. S1: $R^2$ on the terms; stability (intraclass correlation over visits) of connectivity, crosstalk, squared score, residual (first and fourth order), and single edges; each with its MCSE. S2: fingerprinting accuracy; concordance for same- and next-visit activity; the differences A minus B and A minus C; rejection rates of two tests used in the application: concordance above 0.5 (one-sided, child-level permutation) and equivalence of A and B within 0.05 (two one-sided tests, child-level bootstrap). Type I error is the first test's rate at effect 0; power its rate at the four levels.
- Cells and replicates. S1: 13 cells (the base, 4 more stabilities, 2 loadings, 2 sizes, 2 spreads, 3 visits, the null); S2: 5 cells; 200 replicates each (MCSE of a stability difference about 0.002; of a rejection rate at most 0.035).
- Runtime (timed 2026-11-17, one replicate on one core): 1.4 s at the typical size (90 genes), 0.8 s at the largest (115 genes, without first-order edges), 0.3 s for nine elastic-net fits. S1 about 13 x 200 x 2 s / 3 cores, about 15 minutes; S2 with five-fold fits, permutation, and bootstrap about 30 minutes. Production about 45 minutes; quick run (20 replicates) about 5 minutes.
- Application plan. GSE65391, the 131 children with two or more visits (897 visits). Gene sets: BloodGen3 modules (download of the module list from Bioconductor `[TO CONFIRM: approval]`), modules with at least 10 genes on the array `[TO CONFIRM]`; crosstalk among pre-specified interferon, neutrophil, lymphocyte, and plasmablast modules `[TO CONFIRM]`. Describe: trajectories over visits, stability per module, gene and pathway network diagrams for two children. Explain: $R^2$ on the terms and residual stability, first and fourth order; each module placed on the stability and loading gradients by its score's measured stability and its first-component share. Test: identity by fingerprinting; same- and next-visit SLEDAI with A, B, C, R0; folds by child and the authors' 91/40 split; the two tests above.

## The design stop answered

- 2026-11-17. After the proposal, the user asked to add degree, betweenness, total edge weight, transitivity, number of edges, and number of nodes, and said the study is LIONESS read through pathways, not a pathway study. The agent restated the design with LIONESS networks at the centre and pathways as the lens (gene choice and the reading of metrics), all metrics in every setting, 100 replicates. The user: "STAGE 3". Read as approval of that design. The BloodGen3 download, the module-size threshold, the crosstalk modules, and the supervisor's confirmation stay open (`[TO CONFIRM]` in `protocol.md`, section 9).
- Refinements made while writing `protocol.md` (before any run): (1) feature set C is the metrics of the first-order network rather than the closed-form terms, because degree, betweenness, transitivity, and the counts have no closed form; the terms remain the diagnostic of aim 2. (2) Graph metrics (degree, betweenness, transitivity, number of nodes) are computed on the thresholded graph, at a fixed threshold ($\lvert e\rvert > 1$, as pilot-e) and a proportional one (the top 5% of pairs), so that no result rests on one threshold; betweenness on the full weighted graph would add about 5 s per replicate. (3) The superiority test uses a child-level bootstrap instead of a permutation, which would refit every model per permutation.
- `protocol.md` committed alone (3007d5b); the code (dd9d60e); `devtools::check()` clean.
- A first full-size S2 replicate during timing (effect 0.3) gave same-visit concordance 0.55 for the pathway scores and 0.48 for the LIONESS metrics: the metrics carry squared expression and lose its sign. A debugging observation, not a result.

## Step 3.3: calibration and quick run

- The first calibration run (2026-11-17) computed every replicate and then stopped with a parse error: the agent edited `analysis/03_simulation.R` (the settings in the run file) while Rscript was still reading it. Nothing was written. Rerun unchanged: 5 effects x 12 replicates in 386 s. Lesson: never edit a script while it runs.
- Calibration result and the chosen effects: `protocol.md`, section 8 (0.15, 0.18, 0.21, 0.24).
- Quick run: `Rscript analysis/03_simulation.R --run quick --reps 10`, 18 settings x 10 replicates in 613 s; `results/run_info.rds` says `run = "quick"`.

## Step 3.4: drafting

- Carried from earlier files: the closed forms (proposal 14.1, pilot-g report), the paper test's framing. Written new: the exact closed form of the LIONESS edge (S1.1, found while writing the Methods and checked numerically to 1e-14, now `LionessEdgesClosedForm()` with a test), its expansion (S1.2), strength and total weight (S1.4, with a test), Methods, Design, Results, Discussion, Abstract, the AI-use statement (moved into the Methods, as BMC's policy asks).
- The manuscript loads the package (`devtools::load_all()`) so its chunks call the figure functions in `R/`; every figure is one plot (the user's rule): two gene-network diagrams, two pathway-network diagrams, two trajectories, metric stability, connectivity against squared score, the stability gradient, power.
- The GSE65391 sizes and stabilities are read from `pilot-d/results/data_summary.rds`; the application script will write them to `results/` (marker in the draft).

## Revision pass

1. Paragraph pass. Results, aim 1: the first sentence described a figure; it now states the answer ("Every LIONESS metric changed from visit to visit, and none was stable"), and the figures follow. Aim 2 and aim 3 opened with a number; each now opens with its answer. Fingerprinting: "identified 0.018 of the people" became "fingerprinting accuracy was ... against ... by chance", since a reader cannot judge the share without chance.
2. Sentence pass. Zero and one rates printed as "at most 0.00" and "at least 1.00"; the sentences now read "in no replicate" and "in every replicate", built from the numbers. "Co-expression is usually summarized" and "from no stable driver to a very stable one" reworded (checker words).
3. Cross-checks. The Abstract's numbers (largest difference, $R^2$, residual stability, power, equivalence) come from the same expressions as the Results. The Introduction's three questions are answered in the Discussion in order (change over visits; expression; prediction). Table 2's labels were garbled (pathway suffixes stripped wrongly); rebuilt from a fixed label map.
4. Second pass not due: no section changed more than three sentences after the first pass.
- Checker: 0 errors; one warning kept: bib entries never cited (pregnancy-scope and stage-1 entries); `build-bib.py --qmd` drops them at the reference review.

## Step 5: render and check

- Citations: 15 distinct keys cited, 13 passage-backed, 2 definitional (`check-manuscript.py`). Deviations from records: none new.
- Markers left for the user (11 in the rendered draft): the authors (YAML), the BloodGen3 download and module size, the application not yet run, the supervisor's confirmation, the repository URL, whether figure legends must name AI-assisted code.
- Literature pass for this stage (records opened, 2026-11-17): Altman 2021 (BloodGen3, open full text: CC BY 4.0, modules in the BloodGen3Module package), Rinchai 2021 (package), Liberzon 2015 (Hallmark), Horvath and Dong 2008 (open full text: approximately factorizable module networks and eigengene-based connectivity, for a network of all samples), Langfelder and Horvath 2008 (WGCNA), Kuijjer 2019 lionessR (open full text: limma on edges between survival groups), Morris 2019 and Marwick 2018 (Crossref records). Two pages were outside the session's network (nature.com, gsea-msigdb.org).
- Time per step (2026-11-17): step 0 and 1 about 1 h with the figure work; references 20 min; design proposal 20 min; protocol, study code, and tests 1.5 h; calibration 2 x 7 min (one lost to the agent's edit of a running script); quick run 10 min; drafting and revision 1.5 h.

## Process feedback

- Editing an R script while Rscript runs it corrupts the run, because Rscript reads the file as it goes; it cost one calibration run.
- The workflow's step 1 names `analysis/01_simulation.R`; a re-scoped project already holding 01 and 02 scripts needs another number. Worth a sentence in `03-manuscript.md`.
- The skeleton's separate "Statement on AI use" section conflicts with BMC's rule that AI use is described in the Methods; the profile's placement was followed, as the stage file says.
- Clone proof (2026-11-17): a fresh clone with `data/raw/` deleted; `renv::restore()` clean; `analysis/00_download.R` refetched every file, checksums equal; tests 185 pass (the acceptance test skipped, as gated); `devtools::check()` 0 errors, 0 warnings, 0 notes; the quick run rerun (623 s) gave per-setting summaries equal to the committed ones to 1e-12; the manuscript rendered. No Python in this project, so no `uv sync`.

## Revision 1 (2026-11-17)

The user answered open points 1 and 2 ("2 and 1"): the BloodGen3 download and the module choice approved.

| Comment | Change | Where |
|---|---|---|
| Download approved | BloodGen3Module 1.20.0 source package (Bioconductor 3.23, GPL-2) added to `analysis/00_download.R`; module lists read from its `R/sysdata.rda` as data | `analysis/00_download.R`, `R/application.R` |
| Module choice approved | Before any edge: the four pathways chosen by BloodGen3 function label (interferon 6 modules, neutrophil 4, lymphocyte 4, plasma cells 5); modules with at least 10 genes on the array | `protocol.md`, section 8 (b713de2) |
| (run) | First run stopped: two genes vary in one sample only, so their LIONESS edges are undefined; dropped (protocol entry). Second run stopped on an unforced argument in a closure sent to workers; fixed with a cluster test. Third run: 382 modules, 535 network genes, 345 s | `analysis/04_application.R`, `results/application_*.rds` |
| (results) | Application section written; Data section rewritten; Abstract and Discussion extended; five figures (two gene networks, two pathway networks of one child, module stabilities) and one table | `manuscript/manuscript.qmd` |

Pre-registered application claims (`protocol.md`, section 7): claim 5 failed narrowly (2 of 382 modules below $R^2 = 0.99$ on the first-order terms, lowest 0.981; all at least 0.9992 with fourth-order terms), stated in the Results with a marker; claim 6 passed (largest stability difference 0.006); claim 7 passed (equivalence of LIONESS and first-order metrics declared for same-visit and next-visit SLEDAI). Not foreseen by the simulation: the residual beyond fourth order, under 0.1 percent of the variance, kept a median stability of 0.08 over modules; marked in the Results and named as a limitation.

Check after the run (not pre-registered), behind one sentence of the Application: in the example child, the gene touched by most of the 150 strongest edges (146 at visit 2, 35 at visit 1) was in both visits the gene with the largest absolute standardized expression (TIPIN at 10.8 at visit 2; CHEK1 at 5.6 at visit 1), as the closed form's product term implies.

Figure functions: `PlotSampleNetwork()` and `PlotPathwayNetwork()` take pathway labels and leave room for side labels; `PlotModuleStability()` added. Checker: 0 errors; the uncited-bib warning kept as before.

## Revision 2 (2026-11-19)

| Comment | Change | Where |
|---|---|---|
| "Where is the prediction results ... to show to Tyler" | The application's identity and SLEDAI prediction results were in the text only, against the rule that every finding has a table or figure; added a table of identity, same-visit and next-visit concordance, the superiority test, and the authors' split, by feature set | `manuscript/manuscript.qmd`, `tbl-app-pred` |
| Model ladder of the first method (R0, R1, E1, N1, F1) | Not done; proposed to the user as a design change (needs a dated protocol entry before any run); waiting for the answer | this log |

## Revision 3 (2026-11-19): the model ladder

| Comment | Change | Where |
|---|---|---|
| "Yes" to the model ladder of the first project | Protocol entry before any ladder fit (aba5f68); `LadderPredictions()`, `BootstrapLadder()` with tests; ladder added to S2 and the application beside A, B, C; fingerprinting of the genes added | `protocol.md` section 8; `R/study.R`; `analysis/03_simulation.R`, `analysis/04_application.R` |
| (runs) | Quick run rerun with the ladder (1,048 s); application rerun (431 s); the S1 numbers are unchanged (the ladder draws its random numbers after them) | `results/` |
| (draft) | Methods: the ladder and the superiority test of N1 over E1; Results: `tbl-ladder`; Application: `tbl-app-ladder` and text | `manuscript/manuscript.qmd` |

Pre-registered ladder claims: (9, GSE65391) passed: N1 not better than E1 (same visit -0.007, next +0.001), N1 and F1 equivalent at both visits. (8, simulation) N1 and F1 equivalent in every replicate (passes); N1 declared better than E1 in 3 of 10 replicates at the largest effect for next-visit activity, above the limit of 1 in 10 (quick run; marked); in all 6 replicates with such a declaration, N1 and E1 were also declared equivalent, and every mean difference was at most 0.002: the bootstrap declares negligible differences between nearly identical nested predictions. The production run decides the claim.

## Revision 4 (2026-11-19), in progress

The user is revising the draft from the Introduction to the end of the Results and will send the changes; the agent waits for them.

| Comment | Change | Where |
|---|---|---|
| "use we in the intro as the writer", then "i mean it the last paragraph" | The last paragraph of the Introduction now uses "we" (we derive, we then ask, in the Methods we give, in a simulation study we vary, finally we follow); a change to the second paragraph's last sentence was reverted | `manuscript/manuscript.qmd`, Introduction |

Proposed rule for this file (not adopted until the user agrees): the paper's own actions are written with "we" throughout, not only in the Introduction.
| "same thing in the method" | The Methods' sentences about what was done now use "we" (graph metrics, the derivation's first sentence, the two checks, Data, software), and so do the included reproducibility paragraph and AI-use statement; definitions and formulas unchanged. The AI statement's "every citation is verified" became "we verify every citation ... before submission", since the reference review has not happened yet | `manuscript/manuscript.qmd` Methods; `manuscript/parts/reproducibility.md`; `manuscript/parts/ai-use-statement.md` |
| Supervisor, via the user: add a section on modelling disease activity and predicting the patients, with the equations stated and every term defined; derivations to the appendix if needed | New Methods subsection "Stability, identification, and prediction of disease activity": the intraclass correlation (equations), fingerprinting (equations), the same-visit and next-visit models and the elastic-net objective (equations), a table of the models (A, B, C, R0, R1, E1, N1, F1), the concordance (equation) and the bootstrap tests; coefficients written $\theta$ so as not to clash with the activity effect $\beta$. Definitions added under every existing Methods equation (the metrics' sums, $r^{(-q)}$, $O(\varepsilon)$, $S_A$, $Q_A$, $R^B_i$, $C_A$, $S$, $R_i$). Supplement S1.6: the concordance from Kendall's tau-b, as the code computes it. Equations for the simulation's data-generating model were drafted and then removed, because the user will revise the simulation section himself later | `manuscript/manuscript.qmd` Methods; `manuscript/_supplement.qmd` |
| "where are the R_0, R_1 ... etc" | Each model now has its own equation in Methods 2.4 (R1, E1, N1, F1, A, B, C for the same visit; R0 and the added current activity for the next visit), with $\mathbf{s}_{pv}$, $\mathbf{s}^2_{pv}$, $\mathbf{g}_{pv}$, $\mathbf{m}_{pv}$, $\mathbf{m}^{(1)}_{pv}$, and the penalized ($\theta$) and unpenalized ($\gamma$) coefficients defined; the table of models stays | `manuscript/manuscript.qmd` Methods 2.4 |
| "add subsection by each models, under the section 2.4" | Section 2.4 split into subsections: stability over visits; identification of a person; the common form of the prediction models (with the elastic-net objective and the models' table); one subsection per model (R0, R1, E1, N1, F1, A, B, C), each with its equation, its terms, and what it adds; scoring and tests. The combined equation block was replaced by the per-model equations | `manuscript/manuscript.qmd` Methods 2.4 |
| R0's row in the models' table: "do I have to add this", then "do it" (option 1, keep and reword) | R0 kept; its base-terms cell reads "today's activity (next-visit prediction only)" instead of "none (current activity only)" | `manuscript/manuscript.qmd`, `tbl-ladder-def` |
| "then remove it" (the models' table, after each model got its own subsection) | Table of the models removed with its cross-reference; every model stays defined in its subsection (2.4.4 to 2.4.11); the other tables renumber | `manuscript/manuscript.qmd` |
| "make the method properly or arrange properly without touching the simulation section" | Methods reordered, no sentence or equation removed (checked line by line against the previous version): an opening overview (two sentences, "we"); LIONESS single-sample networks; network metrics read through pathways; what LIONESS edges and metrics measure; a diagnostic; stability over visits; identification of a person; prediction of disease activity (notation, the common form, the eight models, scoring and tests); data; software and reproducibility; use of AI (its own subsection, as BMC's policy places it in the Methods). Template headings ("The estimand and the procedure", "What the procedure estimates", "The fix or the diagnostic") renamed to say what each part holds; the template's instruction comments in the Methods removed. The two opening sentences of the old 2.4 were merged into the first sentence of "Stability over visits". Simulation section unchanged (checked) | `manuscript/manuscript.qmd` Methods |
| "add write up to each equation section ... make it strong" | A plain-language paragraph after every equation block of the Methods: how to read the LIONESS edge and why everything starts from it (2.1); the two groups of metrics and why the edge's closed form carries over (2.2); why the exact form is the central result, a reading of the three first-order terms, the two consequences of the pathway forms (direction lost; no information beyond the pathway's own expression), and why stability is predictable (2.3); what each check answers (2.4); why the intraclass correlation (2.5); why identification is the stronger test (2.6); the two choices that keep the model comparison fair (2.7.1); what a result means for each model (R0, R1, E1, N1, F1, A, B, C); why concordance, people-level resampling, and equivalence tests (2.7.10). No new citation and no number; checker words reworded ("simply", "can"). Simulation section unchanged (checked) | `manuscript/manuscript.qmd` Methods |
| "define the terms of the equations clearly, number them and use cross reference" | Every display equation numbered (28: the Introduction's first-order edge and 27 in the Methods, each model with its own number), each term defined in the sentence that follows its equation, and the text points to equations by number (`@eq-...`) and to Methods subsections by number (`@sec-...`). Blocks split so that each quantity has its own number: the metrics (connectivity, crosstalk, total weight and strength), the standardization (new, defines $z_{qi}$ and $\varepsilon$), the leave-one-out correlation, the exact edge, the expression terms, the two pathway closed forms, the strength closed form (was inline), the tests (were inline). Symbol clashes removed: the model error is $\xi$ (was $\epsilon$, beside $\varepsilon = 1/(N-1)$); the training sample count is $n_{\mathrm{tr}}$ and its index $t$ (were $n$ and $k$); a person's number of visits is $n_p$ and one feature's value $u_{pv}$ (were $k_p$ and $x_{pv}$); the features under test are $\mathbf{w}_{pv}$ (were $\mathbf{x}_{pv}$, beside expression $x_{qi}$); the bootstrap count is $B^{*}$; one sentence separates the models A, B, C from the pathways $A$, $B$. Cross-references added in the Application and the Discussion to the exact and first-order forms. Simulation section unchanged (checked) | `manuscript/manuscript.qmd` |
| Supervisor, via the user: no text written on the simulation figures; it belongs in the figure names (captions); every figure explains itself with a clear legend; the space around the figures is too large; "remove ... dots genes coloured by pathway, let that be in the caption" | Figure functions (`R/figures.R`): the notes written under the network and pathway diagrams removed; legends added to every figure (edge sign and width; circle size and fill and line colour and width; one sample and the least-squares line; one child's visits; the two Monte Carlo standard error bars; the reference line of a squared factor); no legend titles (their content is in the captions); the network diagrams' frames fitted to the drawing (no outer margins, tighter limits, figures 7 by 5.5 and 7 by 5.9 inches); headroom so legends never cover data; coinciding series in the gradient and power figures drawn side by side (the LIONESS line had been hidden under the first-order line). Captions of the ten simulation figures rewritten as descriptive names built in the setup chunk, carrying everything the removed notes said; the real-data network figures use the same functions and changed the same way. Found while checking the render: the AI-use statement had not been included in any render since 09987c6, because its include line carried a marker on the same line; fixed | `R/figures.R`; `manuscript/manuscript.qmd` |
| "let the color label be described in the caption instead, remove from the figure"; "the trajectory one should add another colored line that indicated if there are no changes" | Colour legends removed from the gene-network, pathway-network, and trajectory figures (an option of the functions, off by default) and the colours described in full in their captions: pathway colours named, red and blue edges and crosstalk, thin and thick lines, filled and hollow circles. The trajectory figures gained an orange dashed no-change line at the mean of all values, the flat path a child with no change between visits follows, described in the captions. The series legends of the stability, gradient, and power figures were kept (they name the series; the user was asked). Network frames shortened to the drawing now that no legend sits below. The real-data figures' captions moved to the setup chunk (they held typed digits). While replacing the caption block, three captions (connectivity against score, gradient, power) were deleted by mistake and restored from the previous commit before any render was committed | `R/figures.R`; `manuscript/manuscript.qmd` |
| "figure 2 and 4 captions are few" | Captions of the visit-2 gene and pathway networks (simulation, and the real child's) written in full like the visit-1 captions, each pointing to its visit-1 figure and stating that the scales are shared | `manuscript/manuscript.qmd` |
| "change the short dash color for fig 6" | Figure 6 (interferon degree over visits): its children's lines are orange, the colour of the no-change line; the no-change line is now dark blue. Its caption had called the children's lines blue; corrected to orange, with the dark blue dashed line described | `manuscript/manuscript.qmd` |
