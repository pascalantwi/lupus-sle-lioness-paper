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
