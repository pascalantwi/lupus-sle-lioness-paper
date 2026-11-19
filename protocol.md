# Protocol: LIONESS networks over repeated visits, read through pathways

| | |
|---|---|
| Proposal | `proposal.md`, sections 1 to 14 (the re-scope and revisions 1 to 3) |
| Design chosen | by the user, 2026-11-17: the design proposal of `log/03-manuscript.md`, step 3.2, with LIONESS at the centre and pathways as the lens, the six added metrics (degree, betweenness, total edge weight, transitivity, number of edges, number of nodes), and 100 replicates per setting ("STAGE 3") |
| Pilots | `pilot-d/` (close call read as a pass), `pilot-e/` and `pilot-f/` (failed), `pilot-g/` (failed on one clause, accepted with its reading) |
| Written | 2026-11-17, before any run it governs |
| Code | `analysis/03_simulation.R` (simulation), `R/study.R`, `R/pathways.R`, `R/lioness.R` |

## 1. Aims

1. Describe: how stable over a person's visits is each LIONESS network metric (edges, total edge weight, connectivity, crosstalk, strength, degree, betweenness, transitivity, number of edges, number of nodes), read for the whole network and for the genes of each pathway?
2. Explain: does each metric's stability equal that of the same metric computed from the first-order form of the edges, and are pathway connectivity and crosstalk the expression terms of their closed forms (proposal section 14.1), with a remainder that is itself expression?
3. Test: do LIONESS metrics identify a person and predict same-visit and next-visit disease activity better than the pathway expression scores, and differently from the first-order metrics?

## 2. Data-generating mechanisms

**Fixed in every setting** (`SimulateStudyCohort()`): 131 people (the GSE65391 children with two or more visits); six pathways, of which neutrophil-like, lymphocyte-like, plasmablast-like, personal, and null have 15 genes; the neutrophil-like, lymphocyte-like, and plasmablast-like factors load 0.6 with within-person stabilities 0.41, 0.48, and 0.39 (GSE65391: neutrophil and lymphocyte percentages, SLEDAI); the null pathway's noise correlation 0.3; every gene carries a person level of variance 0.5 and unit normal noise; factors independent; an activity score for every visit, the activity effect $\beta$ times the mean of the interferon and plasmablast factors scaled to unit variance, plus a person level and visit noise with within-person stability 0.39.

**S1, describe and explain** (aims 1 and 2): 13 settings, one factor at a time around the base (starred).

| Factor | Levels | Reason |
|---|---|---|
| Interferon factor stability over visits $\kappa_1$ | 0, 0.25, 0.55\*, 0.75, 0.9 | from no stable driver to a very stable one; 0.55 is the IFN5 score's stability in GSE65391; the remainder's stability rose with it (pilot-g) |
| Interferon loading $\lambda_1$ | 0.4, 0.6, 0.8\* | a weak to a strong driving factor; the remainder grows with it (pilot-g, check after the run) |
| Interferon pathway size $m_1$ | 5, 15\*, 40 | small to large modules |
| Personal co-expression spread $\tau$ | 0, 0.1, 0.25\* | none to strong person-specific co-expression in the personal pathway |
| Visits per person | 3, 6\* | the GSE65391 median is 6 |
| Null | all factor stabilities 0, $\tau = 0$ | the known answer: no stable driver |

**S2, test** (aim 3): the base setting with activity effect $\beta \in \{0, \beta_1, \beta_2, \beta_3, \beta_4\}$. $\beta = 0$ gives the type I error; $\beta_1 < \dots < \beta_4$ are set before the quick run by low-replicate calibration cells (section 6) so that the superiority test for the pathway scores (B) runs from low power to about 80%, and are recorded in section 8 with the date.

**Known answers.** The null setting: residual stabilities near 0. First-order connectivity and crosstalk equal their closed forms exactly; graph metrics of a hand-built graph; the concordance against a literal implementation; the LIONESS edges against leave-one-out LIONESS (`tests/testthat/`).

## 3. Estimands

- $I(M)$: the intraclass correlation over visits (one-way, `EdgeICC()`) of a metric $M$ computed on every sample's LIONESS network; $I_1(M)$ the same on the first-order network.
- $D(M) = I(M) - I_1(M)$: the stability a metric has beyond its first-order form.
- For connectivity and crosstalk: $R^2_1$ and $R^2_4$, the $R^2$ over samples of the metric on its first-order terms and with the fourth-order terms added (`PathwayTerms4()`); $I(\text{res}_1)$ and $I(\text{res}_4)$, the stability of the residuals.
- $I(\bar z_A^2)$, the stability of the squared pathway score; $\bar I(e_A)$, the mean stability of the pathway's single edges.
- Fingerprinting accuracy of each feature set (visit 2 against visit 1).
- $C(S)$: the concordance of out-of-fold predictions of activity from feature set $S$, same visit and next visit; the differences $C(A) - C(B)$, $C(A) - C(C)$, and, for the next visit, $C(A) - C(R_0)$.

## 4. Methods

**LIONESS networks.** Pearson LIONESS edges on all samples of a cohort (`LionessEdges()`, equal to leave-one-out LIONESS), and the first-order edges (`LionessEdgesFirstOrder()`).

**Metrics** (`LionessMetrics()`), per sample. Whole network: total edge weight; with the fixed threshold $\lvert e\rvert > 1$ the number of edges, number of nodes (genes with an edge), and transitivity (global clustering); with the proportional threshold (the 5% of pairs with the largest $\lvert e\rvert$) the number of nodes and transitivity. Per pathway: connectivity (mean edge within), crosstalk (mean edge between each pair), and the mean over the pathway's genes of strength (sum of a gene's edges), degree, and betweenness (shortest paths with distance $1/\lvert e\rvert$ on the thresholded graph), the last two under both thresholds.

**Feature sets** (S2, and fingerprinting in every setting): A, the metrics of the LIONESS network; B, the six pathway scores (mean standardized expression) and their squares; C, the metrics of the first-order network; $R_0$, no expression features. Constant columns dropped; columns standardized.

**Prediction** (`OutOfFold()`): Gaussian elastic net, mixing 0.5, the penalty by inner 5-fold cross-validation with folds by person, `lambda.min`; outer 5-fold cross-validation with folds by person. Next-visit prediction uses the features and activity of visit $v$ for the activity of visit $v + 1$, with the current activity unpenalized; $R_0$ is the linear regression on current activity.

**Tests** (`BootstrapTests()`), the procedure of the application: 200 bootstrap resamples of people; superiority when the 5% quantile of $C(S)$ exceeds 0.5; equivalence of A with B, C, or $R_0$ when the 90% interval of the difference lies within $\pm 0.05$.

## 5. Performance measures

Means over replicates with their Monte Carlo standard errors of every estimand in section 3; rejection rates of the superiority tests (type I error at $\beta = 0$, power at $\beta > 0$) and of the equivalence tests. A difference is called nonzero only when it exceeds three Monte Carlo standard errors.

## 6. Replicates, runs, and seeds

- 100 replicates per setting (production); the MCSE of a stability is then about 0.003 to 0.005, of a rejection rate at most 0.05.
- Quick run: 10 replicates per setting.
- Calibration of $\beta_1, \dots, \beta_4$ before the quick run: the S2 base with $\beta \in \{0.05, 0.1, 0.15, 0.2, 0.25\}$, 12 replicates each; the levels are chosen from that grid, or interpolated, so that the power of B's superiority test for same-visit activity spans about 20% to 80%.
- Seeds: setting $c$ has seed $21461117 + 1000c$; replicate $r$ adds $r$. Three PSOCK workers, the BLAS pinned to one thread. Timed 2026-11-17 on one core: 8 to 11 s per S1 replicate, 19 s per S2 replicate; production about 2 hours on three cores.

## 7. Claims pre-registered before any GSE65391 edge is computed

Simulation (production numbers):

1. For every metric in every S1 setting, $\lvert D(M)\rvert \le \max(0.02, 3\,\text{MCSE})$.
2. $R^2_1 \ge 0.99$ for every connectivity and crosstalk in every S1 setting; $I(\text{res}_4) \le 0.05$ in every S1 setting.
3. In the base setting, no pathway's connectivity is more stable than its squared score by more than $\max(0.05, 3\,\text{MCSE})$.
4. In every S2 setting, equivalence of A and C is declared in at least 80% of replicates for same-visit and next-visit activity.

Application (GSE65391):

5. For every BloodGen3 module analysed, $R^2_1 \ge 0.99$ for its connectivity.
6. For every metric, $\lvert I(M) - I_1(M)\rvert \le 0.02$ over the children.
7. Equivalence of A and C is declared for same-visit and next-visit SLEDAI.

The application (dataset, modules, and procedure) is in section 9.

## 8. Dated entries before production

- 2026-11-17, before the quick run. Calibration (`Rscript analysis/03_simulation.R --calibrate`, `results/calibration.rds`; 12 replicates per effect, seeds offset by 500000 from the study's): power of B's superiority test for same-visit activity 0.00 at 0.05 and 0.10, 0.25 at 0.15, 0.42 at 0.20, 0.92 at 0.25 (MCSE up to 0.15). Interpolated: $\beta_1, \dots, \beta_4 = 0.15, 0.18, 0.21, 0.24$. A and C were declared above one half in no replicate at any effect.
- 2026-11-17, before any GSE65391 LIONESS edge is computed. The user approved the BloodGen3 download and the module choice. The application, made precise (section 9 is replaced by this entry where they differ):
  - Samples: the visits of the children with two or more SLE visits (`DescribeGse65391()`), ordered by cumulative time within child; the LIONESS cohort is these visits alone.
  - Expression: the series matrix values; if their 99th percentile exceeds 50 they are taken as linear and log2-transformed (with an offset making the minimum 1). Probes mapped to gene symbols with the GPL10558 table; for a gene with several probes, the probe with the highest mean expression.
  - Modules: BloodGen3 (`BloodGen3Module` 1.20.0, `Module_listGen3`), each module's genes matched to the array; a module is analysed when at least 10 of its genes match.
  - The network: the genes of the four pathways chosen by their BloodGen3 function label, before any edge was computed: interferon ("Interferon", "Type 1 Interferon"; 6 modules), neutrophil ("Neutrophils"; 4), lymphocyte ("Lymphocytes"; 4), plasmablast ("Plasma cells"; 5); a gene in two pathways is kept in the first. LIONESS and first-order networks over these genes are computed one sample at a time from the exact closed form (`LionessEdgesClosedForm()`'s formula), and every metric of section 4 is taken over the four pathways (connectivity, the six crosstalks, strength, degree, betweenness, the whole-network metrics).
  - Module connectivity: for every analysed module, its LIONESS connectivity, its $R^2$ on the first-order terms and with the fourth-order terms, and the stability of each residual (claim 5).
  - Feature sets: A and C as section 4 over the four-pathway network; B the four pathway scores and their squares. Activity is SLEDAI. Same-visit and next-visit prediction, folds by child, the bootstrap tests of section 4 (claim 7); in addition, fit on the authors' training children and the concordance on their test children.
  - Fingerprinting: each child's first two visits.
  - Gradient placement: for each analysed module, the stability of its score over visits and the share of variance of its first principal component.
  - Figures: the LIONESS gene network and pathway network of two children (the highest and the median mean interferon score), visits 1 and 2.
  - Seed 21461118; script `analysis/04_application.R`, results `results/application_*.rds`.
- 2026-11-17, after a first application run stopped before writing any result. A gene that varies in only one sample has zero variance without it, so its leave-one-out correlation, and its LIONESS edge, is undefined (supplement S1.1: $1 - \varepsilon z_{qi}^2 = 0$). Such genes are dropped, with genes of zero variance. In GSE65391 this removed 2 of 10,086 matched genes, none in the four pathways.
- 2026-11-19, after the quick run and the application, before any ladder model is fitted. The user asked to add the model ladder of the first project (R0, R1, E1, N1, F1) as a stricter baseline: the comparison so far (A, B, C) never gave a model the genes themselves. Added to S2 and to the application, beside A, B, and C, which stay as they are:
  - R0: no expression; for next-visit prediction, linear regression on the current activity; for same-visit prediction it has no predictor and is not fitted.
  - R1: linear regression on the pathway scores and their squares (and the current activity for the next visit).
  - E1: elastic net on all genes of the network, with the R1 terms unpenalized.
  - N1: E1 plus the LIONESS metrics (penalized).
  - F1: E1 plus the first-order metrics (penalized).
  - Elastic net, folds, and outer cross-validation as section 4. Contrasts $\Delta_1 = C(\text{N1}) - C(\text{E1})$ (the network beyond the genes) and $\Delta_2 = C(\text{N1}) - C(\text{F1})$ (the network beyond its first-order form), and $C(\text{E1}) - C(\text{R1})$. Tests by the child-level bootstrap of section 4: superiority of N1 over E1 when the 5% bootstrap quantile of $\Delta_1$ exceeds 0; equivalence when the 90% interval of a contrast lies within $\pm 0.05$. Identity: fingerprinting accuracy of the genes beside A, B, C.
  - Claims, fixed now: (8) in every S2 setting, equivalence of N1 and F1 is declared in at least 80% of replicates, and superiority of N1 over E1 in at most 10%; (9) in GSE65391, equivalence of N1 and F1 is declared, and superiority of N1 over E1 is not, for same-visit and next-visit SLEDAI.

## 9. Application

GSE65391 [@banchereau2016personalized]: the 131 children with two or more visits (897 visits), the series matrix and the GPL10558 table (`analysis/00_download.R`). Gene sets: BloodGen3 modules [@altman2021bloodgen3] with at least 10 genes on the array `[TO CONFIRM]`, from the BloodGen3Module package (download `[TO CONFIRM: approval]`); crosstalk among the interferon, neutrophil, lymphocyte, and plasmablast modules `[TO CONFIRM: the module list]`. LIONESS on all 897 visits; every metric and test of section 4; SLEDAI as the activity; the authors' 91/40 training and test split as a second check. Each module placed on the S1 gradients by its score's stability and its first principal component's share of variance.

## Changes after production

| Date | Change | Reason |
|---|---|---|
