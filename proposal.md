# Proposal: LIONESS edges across pregnancy in SLE and RA

| | |
|---|---|
| From | `idea.md`, 2026-11-13 |
| Written | 2026-11-13, workflow 0.2.4-md |
| Status | final 2026-11-13 (stop 1 waived by the user's request to go straight to the pilot) |

## 1. Summary and verdict

Single-sample networks (LIONESS) are read as each sample's own co-expression network, and a change in edges between time points is read as "network rewiring". To first order a LIONESS edge is a quadratic function of the sample's own standardized expression, so any factor that moves many genes and changes over time enters the edges as its square, with no change in co-expression. Pregnancy is such a setting: cell composition changes, and in systemic lupus erythematosus (SLE) the interferon (IFN) signature rises, while in rheumatoid arthritis (RA) the disease often improves (@lien2024pregnancy). This project asks whether LIONESS edges separate pregnancy from the postpartum period beyond the shared factors, their squares, the gene levels, and their own first-order formula, within SLE and within RA, using one longitudinal study (GSE235508), each disease analysed on its own. The closed forms are the user's earlier result, carried over; what is new is a longitudinal, within-woman setting in which the "rewiring" claim is routinely made. Recommendation: pilot as proposed. The pilot is a simulation that copies the study's structure (women with repeated samples, folds by woman), and the dataset is described, not analysed.

## 2. Question and why it matters

**Question (one sentence):** Do LIONESS edges separate pregnant from non-pregnant samples of the same women beyond the IFN score, cell composition, their squares, the gene levels, and their own first-order features, and does the answer differ between SLE and RA?

**Why it matters:** A user who reads a change in LIONESS edges over pregnancy as rewiring of immune co-expression would, under the mechanism, be reading the squared change of the factors that pregnancy moves. If the pilot's claims hold, such a user would first fit the base model with the squared score and the gene levels and then test whether edges add anything; the first-order formula would reproduce any network result without building a network.

## 3. The paper test

The idea is a **critique**: it asks what a procedure in wide use estimates in a setting it was not built for.

| Step | Answer | What is missing |
|---|---|---|
| 1. A procedure in wide use | LIONESS with Pearson correlation, the default of lionessR (@kuijjer2019lioness; @kuijjer2019lionessr); edges compared between time points or used as classifier features | a survey of longitudinal uses of LIONESS ([TO CONFIRM: a cited example of LIONESS used across time points in human blood]{.tc}) |
| 2. A setting it was not built for | Repeated whole-blood samples of the same women through pregnancy; composition shifts, and in SLE the IFN signature rises (@lien2024pregnancy) | none |
| 3. The estimand, as a functional | $\rho_{ij}(P_t)$, the correlation of genes $i$ and $j$ in the distribution of samples at time $t$; rewiring is $\rho_{ij}(P_{\text{preg}}) - \rho_{ij}(P_{\text{post}})$ | none |
| 4. The target and the bias term | The edge targets $r + \mathrm{IF}(x_q)$; under a shared factor $F$ its conditional mean is linear in $F_q^2$ (section 4), so the time difference in mean edge carries the time difference in $\mathrm{E}[F^2]$ with no change in $\rho_{ij}$ | the remainder under within-woman dependence, sketched only |
| 5. A simulation that confirms it | The pilot of section 7: women with four samples each, a factor that shifts in pregnancy, a module that rewires, folds by woman | to run |
| 6. A correction or a diagnostic | The base with the squared score and the gene levels (E1); the first-order formula as a check on every network result | pilot |
| 7. The application | GSE235508, within SLE and within RA, pregnant against non-pregnant samples, folds by woman | stage 3 |

## 4. Mechanism and theory

**Notation.** Sample $q$ of $N$; genes $i$ and $j$; $r$ the Pearson correlation over the $N$ samples; $z_{qi}$ gene $i$ in sample $q$ standardized by the $N$ samples (divisor $N$); $a = 1/(N - 1)$.

**The exact and first-order edge** (the user's earlier result, checked there against literal leave-one-out to $10^{-12}$). The LIONESS edge is $e_q = N r - (N - 1) r_{(-q)}$, and removing one sample is a rank-one downdate of the moments:

$$
\begin{aligned}
r_{(-q)} &= \frac{r - a\, z_{qi} z_{qj}}{\sqrt{\left(1 - a\, z_{qi}^2\right)\left(1 - a\, z_{qj}^2\right)}}, \\
e_q &= r + z_{qi} z_{qj} - \tfrac{r}{2}\left(z_{qi}^2 + z_{qj}^2\right) + O_p\left(N^{-1}\right).
\end{aligned}
$$

The leading terms are the reference correlation plus the influence function of the correlation at the sample (@devlin1975robust); the edge has the form of a jackknife pseudo-value (@quenouille1956bias; @miller1974jackknife). The product term alone is the iENA edge (@yu2017iena).

**A factor that changes over pregnancy** (sketch). Let $z_{qi} = \lambda_i F_q + \varepsilon_{qi}$ with independent noise of variance $\sigma_i^2$. Then

$$
\mathrm{E}\left[e_q \mid F_q\right] = r + \left[\lambda_i \lambda_j - \tfrac{r}{2}\left(\lambda_i^2 + \lambda_j^2\right)\right] F_q^2 - \tfrac{r}{2}\left(\sigma_i^2 + \sigma_j^2\right) + O\left(N^{-1}\right),
$$

so the difference in mean edge between pregnant and non-pregnant samples is $\left[\lambda_i \lambda_j - \tfrac{r}{2}(\lambda_i^2 + \lambda_j^2)\right]\left(\mathrm{E}_{\text{preg}}[F^2] - \mathrm{E}_{\text{post}}[F^2]\right)$, which is not zero whenever pregnancy moves the mean or the spread of $F$, with $\rho_{ij}$ unchanged within each period. Assumptions: finite fourth moments; $N$ grows with the loadings fixed; the samples of one woman are dependent through her own level of $F$, which changes the variance of the edge's mean but not its expectation. What a full derivation would need: the remainder's bound under within-woman dependence, which the pilot checks numerically rather than derives.

**The fix, read off the formula.** A classifier that holds $F$ and $F^2$ (the score and its square) and the gene levels leaves edges only the co-expression that $F$ does not explain. When the score measures $F$ with error, the gene products of genes that share $F$ estimate $F^2$ more precisely than the squared score, so a residual gain can remain (the E5 model, the genes with their products and squares, measures it).

**Numerical checks.** (i) The exact edge against literal leave-one-out refits; (ii) the first-order remainder shrinking with $N$; (iii) the add-one test edge against direct computation; (iv) a null cell, where every model's AUC is near 0.5; (v) the share of edges called different between periods by a paired test in a cell with no change in co-expression, against its nominal 5% in the null cell. Checks (i) to (iii) are unit tests; (iv) and (v) are pilot measures.

## 5. Prior work and novelty

| Source | What it did | What it did not do | Relation to the idea |
|---|---|---|---|
| @kuijjer2019lioness | LIONESS; suggests sample networks for studying change across time | no closed form for the Pearson edge | the procedure |
| @kuijjer2019lionessr | the R implementation, Pearson default | the same | the default analysed |
| @devlin1975robust | influence function of the correlation | no link to single-sample networks | the first-order term |
| @yu2017iena | iENA, product of standardized deviations as an edge | no LIONESS link, no time | the product term |
| @deschildre2024evaluation | benchmark of single-sample network methods | no closed form, no longitudinal setting | edges near-linear in expression, empirically |
| @kuijjer2026challenges | review of single-sample network modelling; its supplement relates edges to covariance products | no correlation closed form, no squared-factor artifact, no longitudinal test | the closest source on the mechanism |
| @lien2024pregnancy | GSE235508: longitudinal whole blood through pregnancy in SLE, RA, and healthy women; deconvolution and differential expression | no co-expression or single-sample network analysis | the dataset; the IFN rise in SLE pregnancy and the composition change motivate the factor |
| Wright and colleagues, Arthritis Res Ther 2023, pre-pregnancy RA signatures, WGCNA ([TO CONFIRM: record seen only in a search result, PMC10548620; not opened]{.tc}) | population co-expression modules, RA and healthy | no single-sample networks, no SLE | a co-expression use in RA pregnancy |
| Hong and colleagues, J Exp Med 2019, longitudinal blood transcriptome in lupus and healthy pregnancy ([TO CONFIRM: record seen only in a search result; not opened]{.tc}) | SLE and healthy pregnancies, modules over time | no RA, no single-sample networks | a fallback for SLE only |

**Verdict:** the mechanism (the squared shared factor in LIONESS edges) is the user's own earlier result and is not published; no source was found that applies LIONESS or any single-sample network to pregnancy, in SLE, RA, or healthy women. The idea adds a longitudinal, within-woman test of the rewiring claim, with folds by woman, and a disease contrast (SLE, where IFN rises, against RA) inside one study. Searches are logged in `log/01-idea-review.md`.

## 6. Datasets

| Dataset | Accession | What is deposited (samples by group, format, access) | Status | Opened |
|---|---|---|---|---|
| Lien et al. 2024 | GSE235508 | 335 whole-blood RNA-seq samples from 84 pregnancies: SLE 101 samples (25 pregnancies), seropositive RA 93 (23), seronegative RA 48 (15), healthy 93 (21); time points 0 pre-pregnancy, 1 to 3 trimesters, 4 to 6 at 6 weeks, 6 months, and 12 months postpartum; processed counts (`GSE235508_mRNA_counts.txt.gz`, 25 MB) and a series matrix (16 KB) with disease activity, medication, and prednisolone use; raw reads withheld under Norwegian law; open | verified | 2026-11-13 |

**Recommended:** GSE235508, within SLE and within RA separately, pregnant (time points 1 to 3) against non-pregnant (0, 4 to 6) samples of the same women, folds by woman. No fallback holds both diseases in one study; the lupus pregnancy cohort of Hong and colleagues is an SLE-only fallback ([TO CONFIRM: accession not opened]{.tc}). **What the application shows:** whether LIONESS edges add anything over E1 in either disease, whether the edges and their formula agree, and whether the apparent rewiring over pregnancy is larger in SLE, where the IFN signature rises, than in RA.

## 7. Pilot plan and pass criteria

**First simulation.** A cohort of $W$ women, each with four samples: two pregnant and two non-pregnant. Fifty genes: genes 1 to 20 load 0.8 on a shared factor $F$ that changes in pregnancy (composition and IFN), genes 21 to 30 form a module with noise correlation 0.3, genes 31 to 50 load 0.4 on a stable composition factor. Each woman has her own level of both factors (within-woman correlation 0.5 for $F$ and 0.8 for composition) and her own offset on every gene (variance 0.5). Pregnancy moves $F$ by a shift and a scale (moderate: 0.5 and 1.5; large: 1 and 2.25, the levels of the user's earlier pilots), or raises the module's correlation by $\Delta$ (rewiring), or both. Eight cells, with $W$ = 20 (the study's size per disease) or 40. Models: $R_0$ (composition and the IFN score residual on composition, logistic), $R_1$ ($R_0$ with both squares), and elastic nets with the $R_1$ terms unpenalized: E1 (the 45 genes 6 to 50), E5 (E1 with the 1,035 products and squares), F1 (E1 with the 990 first-order edges), N1 (E1 with the 990 LIONESS edges), and N1 with $R_0$. Outer 5-fold cross-validation with each woman's four samples in one fold; out-of-fold AUC. Known-answer cells: the null cell. Replicates: 200 if a timing run projects the run within one hour on three cores. Runtime: about one hour, by analogy with the user's third pilot (nine models, 160 samples, nine cells, 200 replicates in 3,316 seconds).

**The three questions and the pass rule, decided before any run.** $\Delta_1 = A(\text{N1}) - A(\text{E1})$, $\Delta_2 = A(\text{N1}) - A(\text{F1})$, the artifact $G = A(\text{N1 with } R_0) - A(R_0)$; each a paired difference within a replicate, averaged, with its MCSE.

1. Does the effect exist, and do the edges equal their formula over time? (a) $\Delta_1 \ge 0.05$ and more than 3 MCSE in the rewiring cell at $W = 40$; (b) $\lvert\Delta_2\rvert \le \max(0.01, 3\,\text{MCSE})$ in every cell.
2. Is it large enough to matter at the study's size? (a) $G \ge 0.02$ in the large-factor cell at $W = 20$; (b) $\Delta_1 \ge 0.03$ in the rewiring cell ($\Delta = 0.4$) at $W = 20$.
3. Is there a fix? $\lvert\Delta_1\rvert \le \max(0.01, 3\,\text{MCSE})$ in the three factor-only cells.

Pass when every clause passes; fail when 1(a), 1(b), or 3 fails; close call otherwise (a failed 2(a) or 2(b), or a value within 3 MCSE of its threshold), stated both ways.

**Code and data that exist:** the user's first project (`work-flow-for-lioness-paper`) holds tested functions for the exact, first-order, and add-one edges and the elastic-net fits; they are copied into this package. GSE235508 is not yet downloaded; the pilot describes it from its 16 KB series matrix.

## 8. Impact and fit

LIONESS is the most used single-sample network method, and time-course "rewiring" is one of the uses its authors propose (@kuijjer2019lioness). From the mechanism, the artifact is as large as the change in $\mathrm{E}[F^2]$ times the loadings, and pregnancy moves composition in every woman and IFN in SLE (@lien2024pregnancy), so the artifact is expected in exactly the comparison users make. A correction (the base with the squared score and the gene levels) and a diagnostic (the first-order formula) are in sight. Paper type: a critique, with a disease application.

**Target journal:** BMC Bioinformatics, **fallback:** Statistical Applications in Genetics and Molecular Biology, from the journal worksheet in `log/01-idea-review.md`.

## 9. Open questions for the user

- [TO CONFIRM: the comparison in the application: pregnant (trimesters 1 to 3) against non-pregnant (pre-pregnancy and postpartum), or the third trimester against 6 months postpartum only (fewer samples, a sharper contrast)]{.tc}
- [TO CONFIRM: seropositive and seronegative RA analysed together as one RA group (about 38 pregnancies), or seropositive RA alone (23)]{.tc}
- [TO CONFIRM: the healthy pregnancies (no pre-pregnancy or 12-month samples) as a third, descriptive analysis, or left out]{.tc}
- [TO CONFIRM: approval to download `GSE235508_mRNA_counts.txt.gz` (25 MB, GEO) for stage 3; the pilot uses only the 16 KB series matrix]{.tc}

## 10. Reference list

Every source this proposal cites has an entry in `references/references.md`. Core citations whose full text is pending:

| Key | Grade | Identifier | Full text | Needed for |
|---|---|---|---|---|
| devlin1975robust | core | doi 10.1093/biomet/62.3.531 | pending | the influence function of the correlation |

The other keys: kuijjer2019lioness, kuijjer2019lionessr, yu2017iena, deschildre2024evaluation, kuijjer2026challenges, quenouille1956bias, miller1974jackknife, lien2024pregnancy; the methods keys friedman2010glmnet, nadeau2003inference, schuirmann1987tost, mai2021ifn5 are carried for stage 3.

## 11. Revisions

- 2026-11-13. Written in one pass at the user's request to redo the pilot for idea 1 in this repository; stop 1 waived by that request, and the open questions of section 9 go to the pilot's stop.

## 12. Revision 1 (2026-11-13): the user's three aims, and the second pilot

This section governs where it differs from sections 1 to 11, which stand as the record.

### 12.1 Why this revision

The idea the user chose ("idea 1", from the options of 2026-11-13) named three things to predict or compare: **pregnancy stage**, **disease activity at each visit**, and **whether edges change differently in SLE than in RA across pregnancy**. Sections 1 to 11 kept only the first; the second was left out and the third reduced to a phrase, without a reason recorded. Stop 1, where the user would have caught it, was waived. The first pilot (`pilot/`) is a correct test of aim 1 by classification and gave a close call (`pilot/results/verdict.rds`): edges equal their first-order formula and the E1 base removes the artifact, but at 20 women neither the artifact nor a real rewiring moves the AUC much; the paired test of edges between periods called rewiring in about 45% of factor-gene pairs where co-expression did not change (`pilot/results/cell_summaries.rds`). The user chose a re-pilot on the three aims (2026-11-13, "go").

### 12.2 The three aims

1. **Pregnancy stage (within SLE, within RA).** Does the paired test of LIONESS edges between pregnancy and the postpartum period call rewiring that is the squared change of the interferon and composition factors, and do partial-correlation edges (given the IFN score and composition) remove it while keeping a real change? Classification of pregnancy stage is kept as a secondary measure.
2. **Disease activity at each visit.** Do LIONESS edges predict activity beyond the IFN score, composition, their squares, and the gene levels, and beyond their own first-order formula? Main analysis in RA (DAS28: 138 visits, 37 women); in SLE the activity index (LAI-P) is zero at 75 of 97 visits (22 visits from 13 women above zero), so SLE activity is analysed only as any activity against none, as a secondary analysis.
3. **SLE against RA across pregnancy.** Do raw LIONESS edges show a different change over pregnancy in SLE than in RA when the only difference is how the shared factors move in each disease, and does the fix remove that false difference while keeping a real one? Each disease's edges use its own reference; the women of the two diseases come from the one study and are compared, never pooled into one reference.

### 12.3 The real size of the factor change (measured before the pilot)

`analysis/01_factor_shift.R` (output `results/factor_shift.rds`) measured, within each disease and among women with samples in both periods, the change in pregnancy of the IFN5 score (the mean of EPSTI1, IFI44L, LY6E, OAS3, RSAD2; @mai2021ifn5) and of the first principal component of the 2,000 most variable genes (a proxy for the main composition axis), in units of the non-pregnant samples. No edge was computed. Rounded: in SLE the IFN5 score shifts by about 0.2 with a spread ratio of about 1.1 and a within-woman correlation of about 0.9, and the composition proxy barely moves; in RA the IFN5 score falls by about 0.3 with a spread ratio of about 0.6, and the composition proxy shifts by about 0.45 with a spread ratio of about 1.45. These smaller-than-assumed changes set the "real" levels of pilot-b; the first pilot's large level is kept as a stress setting.

### 12.4 Pilot-b design

Fifty genes as before, with two shared factors that may both move in pregnancy: the IFN factor on genes 1 to 20 (score: genes 1 to 5) and the composition factor on genes 31 to 50 (observed without error); module genes 21 to 30. Each factor has its own shift, spread ratio, and within-woman correlation per disease (12.3). Four samples per woman (two in each period) for aims 1 and 3; four visits per woman with an activity score for aim 2, where the activity raises the spread of the IFN factor (artifact cell) or the module's correlation (signal cell).

| Cell | Aim | What changes | Women |
|---|---|---|---|
| 1 | 1 | nothing (null) | 18 |
| 2 | 1 | SLE levels of both factors | 18 |
| 3 | 1 | RA levels of both factors | 27 |
| 4 | 1 | large IFN change (stress: shift 1, spread 2.25) | 18 |
| 5 | 1 | SLE levels and module correlation +0.4 in pregnancy | 18 |
| 6 | 2 | RA levels; activity raises the IFN factor's spread only | 36 |
| 7 | 2 | RA levels; activity raises the module's correlation | 36 |
| 8 | 3 | SLE levels (18 women) against RA levels (27 women); no rewiring | 18 + 27 |
| 9 | 3 | as cell 8, with module +0.4 in SLE pregnancy only | 18 + 27 |
| 10 | 3 | large IFN change in SLE against none in RA (stress) | 18 + 27 |

**Measures.** Aims 1 and 3: the share of gene pairs called changed at the 5% level by a paired test across women (aim 1) or by a two-sample test of the women's changes between diseases (aim 3), for raw LIONESS edges and for partial-correlation LIONESS edges given the IFN score and composition, among the factor-gene pairs (genes 6 to 20) and the module pairs. Aims 1 and 2: the out-of-fold AUC (pregnancy) or the out-of-fold correlation between predicted and observed activity, for $R_0$, $R_1$, E1, F1, N1, and N1 with $R_0$, folds by woman. 200 replicates per cell if the timing run fits one hour on three cores.

### 12.5 Pass rule (to be frozen in `pilot-b/preregistration.md`)

| Q | Clause | Quantity | Threshold | Cells |
|---|---|---|---|---|
| 1 | (a) | raw-edge share called changed, factor pairs | $\ge 0.10$ | 4 |
| 1 | (b) | partial-edge share called changed (aim 1) or differently changed (aim 3), factor pairs | $\le 0.05 + \max(0.01, 3\,\text{MCSE})$ | 2, 3, 4, 8, 10 |
| 1 | (c) | partial-edge share of module pairs called changed, over the raw-edge share | $\ge 0.75$ | 5 |
| 2 | (a) | raw-edge share called changed, factor pairs | $\ge 0.10$ in cell 2 or cell 3 | 2, 3 |
| 2 | (b) | raw-edge share called differently changed, factor pairs | $\ge 0.10$ | 8 |
| 2 | (c) | partial-edge share of module pairs called differently changed | $\ge 0.15$ | 9 |
| 3 | (a) | $\lvert A(\text{N1}) - A(\text{F1}) \rvert$ (AUC or activity correlation) | $\le \max(0.01, 3\,\text{MCSE})$ | 1 to 7 |
| 3 | (b) | activity: $\lvert$N1 $-$ E1$\rvert$ | $\le \max(0.01, 3\,\text{MCSE})$ | 6 |
| 3 | (c) | activity: N1 $-$ E1 | $\ge 0.03$ | 7 |

Pass when every clause passes. Fail when 1(b), 3(a), or 3(b) fails (the fix, the equivalence, and the artifact's removal from activity prediction are the claims the paper rests on). Close call otherwise, stated both ways. A value within 3 MCSE of a one-sided threshold is a close call for that clause.

### 12.6 What the user should check at this stop

- That the three aims read as intended.
- RA as one group (seropositive and seronegative); all pregnancy samples against all non-pregnant ones (defaults, proposal section 9).
- SLE activity as any against none (secondary), RA activity (DAS28) as the main activity analysis.
- [TO CONFIRM: the user's approval of this revision and its pass rule before pilot-b runs]{.tc}

## 13. Revision 2 (2026-11-13): the third pilot

This section governs where it differs from section 12. The user read pilot-b's close call (`pilot-b/results/verdict.rds`: every deciding clause passed; 2(b) and 3(c) failed; 2(a) and 2(c) close) and asked for a third pilot ("go to C"), going straight on without a separate stop. Three design changes, each for a reason pilot-b gave, with **the pass rule of section 12.5 unchanged**:

1. **Every visit.** Five samples per woman, mirroring the visits per woman in GSE235508 (`results/factor_shift_markers.rds`: about 2.7 pregnant and 2.4 non-pregnant in SLE, 2.3 and 2.3 in RA): three pregnant and two non-pregnant in SLE, two and three in RA. Pilot-b used four.
2. **Cell mix from marker genes.** The composition level comes from a neutrophil-minus-lymphocyte marker score (FCGR3B, CXCR2, CSF3R, FPR1 against CD3E, CD3D, CD2, MS4A1, CD79A; [TO CONFIRM: a published source for this marker list]{.tc}), measured by `analysis/02_factor_shift_markers.R`. Rounded: SLE shift 0.6, spread ratio 0.85, within-woman share 0.8; RA shift 1.4, spread ratio 0.85, within-woman share 0.1. Pilot-b's principal-component proxy gave much smaller shifts (about 0 and 0.45).
3. **Activity within woman.** Activity and every feature are centered on the woman's own mean before fitting and scoring (pilot-b, decisions log D3), so the prediction uses each visit's departure from the woman's usual level.

A change considered and dropped: testing a smaller pre-chosen set of edges does not change any pair's test, since each share counts pairs tested one by one.
