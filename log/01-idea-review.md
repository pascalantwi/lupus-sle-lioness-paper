# Log, stage 1: idea review

Workflow 0.2.4-md. Started 2026-11-13.

First line (step 1.4): the question as read is whether LIONESS edges separate pregnant from non-pregnant samples of the same women beyond the IFN score, composition, their squares, the gene levels, and their own first-order features, within SLE and within RA in GSE235508; the paper type is a critique (what LIONESS "rewiring" over time estimates when shared factors change), with a disease application.

## 1. Inventory of `idea.md`

| Element | What `idea.md` says | Stated, implied, or absent |
|---|---|---|
| Question | edges vs shared factors over pregnancy, SLE and RA | stated |
| Procedure | LIONESS, Pearson | stated |
| Setting | repeated whole-blood samples through pregnancy | stated |
| Mechanism | the squared shared factor in the first-order edge | stated |
| Estimand | the correlation at each period; rewiring as its difference | implied |
| Claimed theory | exact and first-order edges (the user's earlier result) | stated |
| Aims | equivalence, artifact and fix, real rewiring | stated |
| Literature | none cited beyond the dataset paper | implied |
| User's prior work | the first LIONESS project | stated |
| Datasets | GSE235508 | stated |
| Journal | none | absent; carried from the first project |
| Constraints | one study, no merging; SLE and RA both | stated by the user in conversation |

`tools/init.R` was run before this stage (package `lionesspreg`).

## 2. Paper test

In `proposal.md`, section 3 (the same table, written once).

## 3. Searches

| # | Date | Source | Query | What it gave | Time |
|---|---|---|---|---|---|
| 1 | 2026-11-13 | GEO | GSE235508 series record and its sample characteristics (series matrix, 16 KB) | 335 samples, groups and time points as in proposal section 6 | 3 min |
| 2 | 2026-11-13 | PubMed | PMID 38049980 (identifier lookup, outside the budget) | Lien 2024 record and abstract | 1 min |
| 3 | 2026-11-13 | web search | LIONESS single-sample network pregnancy longitudinal gene expression rewiring | LIONESS method papers; Wright 2023 healthy-pregnancy RNA-seq (Front Immunol); no LIONESS use in pregnancy | 1 min |
| 4 | 2026-11-13 | web search | GSE235508 pregnancy rheumatoid arthritis lupus co-expression network analysis | Wright 2023 RA pre-pregnancy WGCNA (PMC10548620); no network analysis of GSE235508 | 1 min |
| 5 | 2026-11-13 | web search | pregnancy immune gene co-expression network rewiring blood transcriptome trimester postpartum | Hong 2019 J Exp Med lupus pregnancy; Brummaier 2024 BloodGen3 modules; vitamin D WGCNA; no single-sample networks | 1 min |

The mechanism's prior work (LIONESS, lionessR, Devlin 1975, iENA, the 2024 evaluation, Kuijjer 2026 and its supplement) was searched and read in the first project (its `log/01-idea-review.md`, queries 1 to 25) and is carried, not searched again. Three new queries against a budget of about twenty; the stop rule held at once (no query added a single-sample network use in pregnancy). Wright 2023 and Hong 2019 were seen only in search results; they carry markers in the proposal and no reference entry.

## 4. Theory check

The exact and first-order edges and the conditional mean under a shared factor are the first project's results, checked there by unit tests to $10^{-12}$ and by its pilots. New here: within-woman dependence. It changes the variance of a mean edge, not its expectation, because the conditional-mean step uses only the sample's own $F_q$; the remainder under dependence is not derived and is checked numerically by the pilot (D2).

## 5. Datasets

GSE235508 opened 2026-11-13 (series record and series matrix): see proposal section 6. Women with both pregnant and non-pregnant samples, by group, are counted by `pilot/describe_data.R` at stage 2.

## 6. Journal worksheet

Carried from the first project (same method, same readership): a statistical-methods critique with a genomics application.

- Target: BMC Bioinformatics.
- Fallback: Statistical Applications in Genetics and Molecular Biology.
- Reason: methods readership that uses single-sample networks; open access; no length limit that the simulation would strain.

## 7. Decisions

- D1. Paper type: critique (the procedure's target under time-varying shared factors), with a disease application.
- D2. The pilot copies the study's structure: four samples per woman, woman-level factor levels and gene offsets, folds by woman. Without folds by woman, a woman's other samples in the training fold would leak her offsets into the test score.
- D3. The realistic size is $W = 20$ women per disease (SLE 25 pregnancies; seropositive RA 23), to be confirmed by the data description.
- D4. Models E4 (regularized discriminant analysis) and N2 (partial-correlation edges) of the first project are left out of this pilot: they answered questions of the first project, and the time box goes to the longitudinal structure.
- D5. Stop 1 waived: the user asked for the pilot to be run and reported (2026-11-13); the open questions of proposal section 9 go to the pilot's stop.

## 8. Time

Step 1 to 8: about 40 minutes of agent time.

## 9. Process feedback

- The seed reference list was carried from the first project with each entry's original record and date; the user reviews `references.md` once, at the end.

# Re-scope, 2026-11-16: personal-network stability (GSE65391)

First line (step 1.4): the question as read is whether LIONESS personal networks are personal, that is, how much of their stability over repeated samples, their ability to identify the person, and their prediction of the next disease activity is real co-expression and how much is IFN and cell mix entering the edges squared; a critique with a correction.

## R1. Inventory of the new `idea.md`

| Element | What it says | Stated, implied, or absent |
|---|---|---|
| Question | stability, fingerprinting, next-visit activity of personal LIONESS networks | stated |
| Procedure | LIONESS, Pearson | stated |
| Setting | repeated samples of the same patient | stated |
| Mechanism | the squared shared factor; stability of $F^2$ | stated |
| Estimand | person-level co-expression and its intraclass correlation | implied; stated in the proposal |
| Theory | the closed forms (user's); ICC of $F^2$ (agent's sketch) | stated |
| Data | GSE65391; ILLUMINATE optional | stated |
| Constraints | method and title unchanged (user); SLE against RA dropped (user); simulation mirrors the application (user) | stated |

## R2. Paper test

In `proposal.md`, section 3.

## R3. Searches (2026-11-13 to 2026-11-16; web search unless stated)

| # | Date | Query | What it gave |
|---|---|---|---|
| 1 | 11-14 | GSE88884 ILLUMINATE lupus whole blood SLEDAI interferon | ILLUMINATE design; IFN-positive 75%; Hoffman et al. 2017 |
| 2 | 11-14 | longitudinal autoimmune flare IFN dermatomyositis Sjogren dataset | Stanford DM cohort (no GEO); GSE318067 DM controlled |
| 3 | 11-14 | LIONESS single-sample networks breast cancer TCGA | LIONESS applied to breast cancer subtypes (2024), TNBC response (2025) |
| 4 | 11-14 | longitudinal blood transcriptome pregnancy SLE RA GEO | only GSE235508 holds both |
| 5 | 11-14 | public GEO autoimmune IFN longitudinal myositis Sjogren SSc | MS IFN-beta series GSE41846; a harmonized compendium |
| 6 | 11-15 | tabalumab ILLUMINATE co-expression network GSE88885 | drug effect on B-cell genes; no network re-analysis |
| 7 | 11-15 | cell composition confounds differential co-expression | @farahbod2020composition |
| 8 | 11-15 | single-sample network LIONESS reproducibility test-retest | only reference-set and replicate reproducibility |
| 9 | 11-15 | LIONESS SLE or autoimmune | none found |
| 10 | 11-15 | interferon discordance SLEDAI SACQ | SACQ interferon resembles remission |
| 11 | 11-15 | ILLUMINATE ancestry JCI Insight | 2020 ancestry paper (title only) |
| 12 | 11-15 | sex differences networks LIONESS GTEx | Lopes-Ramos et al. 2020 |
| 13 | 11-15 | lupus subtypes GSE65391 | 3, 5, 7 subtypes reported |
| 14 | 11-15 | longitudinal stability IFN signature lupus | @petri2019stability |
| 15 | 11-15 | single-sample network reliability intraclass correlation | nothing for blood networks |
| 16 | 11-15 | differential co-expression network lupus whole blood | Ntasis et al. 2020 (Sci Rep) |
| 17 | 11-15 | LIONESS reference dependence | acknowledged, not quantified |
| 18 | 11-16 | fingerprinting gene co-expression network blood repeated samples | none for blood networks; @whitney2003individuality; @mishra2026temporal |
| 19 | 11-16 | GSE65391 single-sample network LIONESS SLEDAI prediction | no LIONESS on GSE65391; several ML predictions of activity and flares |
| 20 | 11-16 | Europe PMC records by DOI (identifier lookups, outside the budget) | the six new entries of `references/references.md` |

Budget: 19 queries against about 20. Stop rule: queries 18 and 19 added no study of personal-network stability or blood-network fingerprinting.

## R4. Theory check

The closed forms are the user's (tested in the first project). New, the agent's sketch: for a factor with within-person intraclass correlation $\kappa$, the intraclass correlation of $F^2$ is $\kappa^2$ (from $\mathrm{Cov}(X^2, Y^2) = 2\rho^2$ for a standard bivariate normal pair). Check proposed: a simulation over a grid of $\kappa$, with $\kappa = 0$ as the known-answer case.

## R5. Datasets

GSE65391 and ILLUMINATE opened 2026-11-15; descriptive checks in `log/01-rescope-data-check/` (README there). Neither file is in `data/raw/` yet.

## R6. Decisions

- D-R1. Re-scope approved by the user's requests of 2026-11-14 to 2026-11-16 ("same method", "same title", "the simulation should be same as application", "I do not want to keep SLE vs RA", "do stage 1 ... save it in the new repo").
- D-R2. GSE65391 chosen over ILLUMINATE: SLEDAI and measured blood counts at every visit; ILLUMINATE has clinical values at baseline only and showed little relation of IFN and cell mix to activity.
- D-R3. Ideas set aside, with the data check behind each, in the scratch review (summarized in this section): within-lupus activity, kidney, and flare predictions alone (factors barely change), drug effect (known), ancestry and sex (small differences), lupus against healthy (too easy, the supervisor), subtypes (crowded).
- D-R4. Stop 1 is held this time; the pregnancy proposal's waived stop 1 led to a missed aim.

## R7. Process feedback

- The data check before the proposal (factor stability, group differences) is what the pregnancy pilots lacked; it should come before every proposal.
