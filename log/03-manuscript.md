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
