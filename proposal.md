# Proposal: LIONESS edges are quadratic functions of expression: closed forms, an interferon artifact, and equivalence with first-order features

| | |
|---|---|
| From | `idea.md`, 2026-11-16 (re-scope) |
| Written | 2026-11-16, workflow 0.2.4-md |
| Status | final 2026-11-16 (the user: "start with the stage 2"; defaults of section 11 taken) |
| Earlier proposal | the pregnancy proposal and its revisions 1 and 2 are the record in git: `git show c677b11:proposal.md` |

## 1. Summary and verdict

LIONESS gives every sample its own gene network, and users read it as the network of that person. If it is, it should be stable over repeated samples of the same person, it should identify her, and it may track her disease. A LIONESS edge is, to first order, a quadratic function of the sample's own standardized expression, so a shared factor such as the interferon (IFN) response or the blood-cell mix enters every edge as its square, and the edge inherits that factor's stability. This project asks how much of a personal network's stability over time, its ability to identify the person ("fingerprinting"), and its ability to predict the next visit's disease activity is real gene coordination, and how much is IFN and cell mix, in children with lupus followed over repeated visits (GSE65391, @banchereau2016personalized). The closed forms are the user's earlier result; the stability of personal LIONESS networks and blood-network fingerprinting were not found in the searches of this stage. A first look at the data shows a moderately stable IFN score and a changing cell mix within each child, so the theory predicts both false stability and visit-to-visit noise in raw edges. Recommendation: pilot as proposed, after the open questions of section 9 are answered.

## 2. Question and why it matters

**Question (one sentence):** Over repeated blood samples of the same patient, how stable are her LIONESS edges, can they identify her, and do they predict her next disease activity, beyond her IFN level, her measured blood-cell counts, her gene levels, and the edges' own first-order formula?

**Why it matters:** A personal network can serve as a biomarker only if it is reproducible and carries something the cheap measures do not. If its stability is mostly the person's stable IFN level, a user would measure IFN directly; if a corrected network keeps a stable personal part, that part is the candidate biomarker. Clinicians need the first property of any blood test, the same answer on repeat testing; researchers using LIONESS need to know whether "personal" is earned.

## 3. The paper test

The idea is a **critique** with a constructive correction.

| Step | Answer | What is missing |
|---|---|---|
| 1. A procedure in wide use | LIONESS with Pearson correlation, the lionessR default (@kuijjer2019lioness; @kuijjer2019lionessr); edges read as each person's network | none |
| 2. A setting it was not built for | Repeated samples of the same person, where shared factors have their own within-person stability: IFN fairly stable, cell mix changing | none |
| 3. The estimand, as a functional | The person-level co-expression $\rho_{ij}(P_w)$ of person $w$; its stability is the share of edge variance between persons, the intraclass correlation of the true personal co-expression | none |
| 4. The target and the bias term | The edge's intraclass correlation equals the share of its variance from stable personal co-expression plus a term from the stability of $F^2$, about $\kappa_F^2$ for a factor with intraclass correlation $\kappa_F$ (section 4) | the remainder under unequal visit numbers, sketched only |
| 5. A simulation that confirms it | Section 7: children with repeated visits, stable personal factors, real or no personal network | to run |
| 6. A correction or a diagnostic | Partial-correlation edges given the IFN score and the measured cell counts; the base with the squared scores; the first-order formula as a check | pilot |
| 7. The application | GSE65391: stability, fingerprinting, and next-visit activity in 131 children with two or more visits | stage 3 |

## 4. Mechanism and theory

**The edge** (the user's earlier result, tested there against literal leave-one-out to $10^{-12}$): with $z_{qi}$ gene $i$ of sample $q$ standardized by the reference (divisor $N$) and $r$ the reference correlation,

$$
\begin{aligned}
e_q &= N r - (N - 1)\, r_{(-q)}, \\
e_q &= r + z_{qi} z_{qj} - \tfrac{r}{2}\left(z_{qi}^2 + z_{qj}^2\right) + O_p\left(N^{-1}\right).
\end{aligned}
$$

The leading terms are the reference correlation plus the influence function of the correlation (@devlin1975robust); the product term is the iENA edge (@yu2017iena).

**Stability under a stable factor** (sketch, labelled as the agent's). Let a visit $v$ of person $w$ have $z_{wvi} = \lambda_i F_{wv} + \varepsilon_{wvi}$, with $F_{wv} = a_w + b_{wv}$, $a_w \sim N(0, \kappa)$ and $b_{wv} \sim N(0, 1 - \kappa)$, so that $F$ has unit variance and intraclass correlation $\kappa$ across visits. Two visits of the same person have $(F_{w1}, F_{w2})$ bivariate normal with correlation $\kappa$, so

$$
\begin{aligned}
\mathrm{Cov}\left(F_{w1}^2, F_{w2}^2\right) &= 2\kappa^2, \qquad \mathrm{Var}\left(F_{wv}^2\right) = 2, \\
\mathrm{ICC}\left(F^2\right) &= \kappa^2 .
\end{aligned}
$$

Line 1 uses $\mathrm{Cov}(X^2, Y^2) = 2\rho^2$ for a standard bivariate normal pair with correlation $\rho$. Since the edge's mean given $F$ is linear in $F^2$ (the user's earlier result), the part of an edge's variance that comes from $F$ is stable with intraclass correlation $\kappa^2$: with no personal co-expression at all, edges of genes that share $F$ look personal to the degree $\kappa^2$ times their share of variance from $F^2$. A factor with low $\kappa$, such as a changing cell mix, adds mostly within-person noise and lowers the edge's stability. Assumptions: normal factor and noise, a fixed reference, equal visits per person. What a full derivation would need: the noise and remainder terms, and unequal visit numbers. **Numerical check:** a simulation over a grid of $\kappa$ comparing the edges' intraclass correlation with $\kappa^2$ times the factor share, and the special case $\kappa = 0$, where the edge's intraclass correlation must be that of the personal co-expression alone.

**The correction.** Partial-correlation edges given the IFN score and the measured cell counts remove $F$ from each gene before the edge is built (the first project's N2), so their stability should be that of the personal co-expression; the score's measurement error leaves a residual, which the pilot measures.

**What the data say about the levels** (agent's data check, `log/01-rescope-data-check/`, not a result of this project): within children with two or more visits in GSE65391, the IFN5 score (@mai2021ifn5) has an intraclass correlation of about 0.55, neutrophil percentage about 0.22, lymphocyte percentage about 0.30, and SLEDAI about 0.39. So $\kappa_F^2 \approx 0.30$ for IFN and about 0.05 to 0.09 for the cell mix. The IFN signature has been reported as the most stable lupus gene signature within patients (@petri2019stability).

## 5. Prior work and novelty

| Source | What it did | What it did not do | Relation to the idea |
|---|---|---|---|
| @kuijjer2019lioness | LIONESS; reproducibility across independent reference sets | no repeated samples of the same person | the procedure; the open question |
| @kuijjer2019lionessr | lionessR, Pearson default | the same | the default analysed |
| @devlin1975robust; @yu2017iena | influence function of the correlation; iENA's product edge | no LIONESS, no repeated samples | the first-order terms |
| @kuijjer2026challenges | review of single-sample networks; covariance products | no correlation closed form, no stability over time | closest on the mechanism |
| @farahbod2020composition | cell composition drives much bulk co-expression in blood; correction gave mixed results | no single-sample edges, no repeated samples | the composition mechanism for population co-expression |
| @whitney2003individuality | donor-specific blood expression over sequential samples, partly traced to cell-subset proportions | no networks | individuality of expression |
| @mishra2026temporal | temporal variation of blood expression in 333 healthy people over six months; immune transcripts stable | no networks, no disease | expression-level stability baseline |
| @petri2019stability | the IFN signature is the most stable lupus signature within patients; its changes do not track activity | no networks | the factor stability the theory uses |
| @finn2015fingerprinting | brain connectivity profiles identify individuals across sessions | not blood, not gene networks | the fingerprinting design |
| @banchereau2016personalized | GSE65391; IFN signature, plasmablast marker of activity, neutrophils toward nephritis | no single-sample networks | the dataset |

**Verdict:** the stability of personal LIONESS networks over repeated samples, and identifying people from blood gene networks, were not found in this stage's searches (logged in `log/01-idea-review.md`, re-scope section). The idea adds a closed-form account of where a personal network's stability comes from, a correction, and a test on 131 children with repeated visits. Predicting lupus activity and flares from GSE65391 has been done with other methods; the contribution there is the comparison of edges with genes and with their formula. [TO CONFIRM: a full search of the fingerprinting and network-reliability literature before "first" is claimed]{.tc}

## 6. Datasets

| Dataset | Accession | What is deposited | Status | Opened |
|---|---|---|---|---|
| Pediatric lupus, longitudinal (@banchereau2016personalized) | GSE65391 | 924 SLE visits from 158 children (131 with two or more visits; median 6, up to 22; median 84 days apart) and 72 samples from 46 healthy children; Illumina HT-12 (GPL10558) normalized matrix (91 MB); SLEDAI and its 24 items at every visit, neutrophil, lymphocyte, and monocyte counts and percentages, anti-dsDNA, C3, C4, kidney class, treatments; the authors' training (91 children) and test (40) split | verified | 2026-11-15 |
| ILLUMINATE-1 and -2, tabalumab trials | GSE88884, GSE88885, GSE88886 | 1,760 adults with SLE and 60 healthy at baseline; 456 patients at baseline, week 16, and week 52 by arm; HTA 2.0 preprocessed matrix (1.3 GB); SLEDAI, anti-dsDNA, C3, C4 at baseline only | verified | 2026-11-15 |

**Recommended:** GSE65391 for every aim. **Fallback and adult check:** ILLUMINATE placebo patients with three visits, for the stability aim only. **What the application shows:** how much of a child's network stability and identifiability is her IFN and cell mix, and whether corrected networks predict her next disease activity.

## 7. Pilot plan and pass criteria

**First simulation (mirrors GSE65391).** $W = 131$ children, six visits each (the median), with a cell of three visits to mirror ILLUMINATE. Fifty genes: genes 1 to 20 load 0.8 on the IFN factor $F$ (score: genes 1 to 5), genes 31 to 50 load 0.6 on the cell-mix factor $C$ (observed, as the measured counts are), genes 21 to 30 form a module. $F$ and $C$ have within-child intraclass correlations $\kappa_F$ and $\kappa_C$; each child's module correlation is $0.3 + \tau u_w$ with $u_w$ uniform on $(-1, 1)$, the same at every visit (a real personal network when $\tau > 0$). Activity has within-child intraclass correlation 0.39 and is either unrelated to expression or raises the module correlation at that visit.

| Cell | $\kappa_F$ | $\kappa_C$ | $\tau$ | Visits | Activity | Purpose |
|---|---|---|---|---|---|---|
| 1 | 0 | 0 | 0 | 6 | unrelated | null: nothing personal |
| 2 | 0.55 | 0.22 | 0 | 6 | unrelated | real levels, no personal network |
| 3 | 0.9 | 0.22 | 0 | 6 | unrelated | stress: very stable IFN |
| 4 | 0.55 | 0.22 | 0.25 | 6 | unrelated | real personal network |
| 5 | 0.55 | 0.22 | 0.25 | 3 | unrelated | ILLUMINATE-like visits |
| 6 | 0.55 | 0.22 | 0 | 6 | raises the module correlation | prediction signal |

**Measures.** The mean intraclass correlation of edges over IFN pairs, module pairs, and cell-mix pairs, for raw LIONESS edges, first-order edges (F1), and partial-correlation edges; fingerprinting accuracy, the share of children whose visit-2 edge vector is closest (Pearson correlation) to their own visit-1 vector among all 131; next-visit activity: the concordance index of $R_0$, $R_1$, E1, F1, N1, N2 (folds by child). 200 replicates per cell if the time box allows; MCSE with every value.

**The three questions and the pass rule** (draft numbers, fixed in the pre-registration):

1. Does the effect exist? (a) In cell 3, the mean raw intraclass correlation over IFN pairs $\ge 0.2$ (the sketch predicts about $0.81$ times the IFN share); (b) in every cell and pair group, raw and first-order edges' intraclass correlations differ by at most 0.02.
2. Is it large enough at realistic sizes? In cell 2, (a) raw fingerprinting accuracy $\ge 0.05$ (chance is $1/131 \approx 0.008$); (b) mean raw intraclass correlation over IFN pairs $\ge 0.10$.
3. Is there a fix? (a) In cells 2 and 3, partial-edge intraclass correlation over IFN pairs $\le 0.05$ and partial-edge fingerprinting $\le 0.02$; (b) in cell 4, the partial module intraclass correlation at least 0.75 times the raw, and partial fingerprinting on module pairs $\ge 0.05$; (c) in cell 6, $C(\text{N1}) - C(\text{E1}) \ge 0.03$, and in cells 1 to 5 $\lvert C(\text{N1}) - C(\text{E1}) \rvert \le \max(0.01, 3\,\text{MCSE})$.

Pass when every clause passes; fail when 1(b), 3(a), or 3(b) fails; close call otherwise, stated both ways; a value within 3 MCSE of a one-sided threshold is a close call for its clause. Reasons: 0.05 fingerprinting is six times chance; 0.02 for the corrected fingerprinting allows about two times chance for the score's measurement error; 0.75 keeps three quarters of a real signal, as in the earlier pilots.

**Code and data that exist.** In this repository: the exact, first-order, add-one, and partial-correlation edges with tests; simulators with person-level factors and repeated samples (`SimulateCohortB()`), folds by person, the concordance index, and elastic-net fits. The data-check scripts are in `log/01-rescope-data-check/`. GSE65391 is not yet in `data/raw/`.

## 8. Impact and fit

LIONESS is the most widely used single-sample network method, and "a network for each person" is its promise; no test of that promise over repeated samples was found. The sketch predicts that at the IFN stability measured in GSE65391 a sizeable part of edge stability is $F^2$, and that the changing cell mix adds noise, so both a false personal signal and a hidden one are expected. A correction and a diagnostic (the first-order formula) are in sight. Paper type: a critique with a correction.

**Target journal:** BMC Bioinformatics, **fallback:** Statistical Applications in Genetics and Molecular Biology (carried from the first project's worksheet; same readership).

## 9. Open questions for the user

- [TO CONFIRM: the supervisor agrees to drop the SLE against RA application and to redesign the simulation to mirror the repeated-visit data, with the method and title unchanged]{.tc}
- [TO CONFIRM: aim 5 (next-visit disease activity) kept, with SLEDAI at the next visit as the outcome rather than a flare]{.tc}
- [TO CONFIRM: ILLUMINATE placebo patients as an adult check of the stability aim (two more downloads, about 300 MB each), or GSE65391 alone]{.tc}
- [TO CONFIRM: approval to download GSE65391 (91 MB) and the GPL10558 table into `data/raw/` at stage 2, not committed]{.tc}

## 10. Reference list

Every source this proposal cites has an entry in `references/references.md`. Core citation with full text pending: devlin1975robust (doi 10.1093/biomet/62.3.531). The other keys: kuijjer2019lioness, kuijjer2019lionessr, yu2017iena, kuijjer2026challenges, mai2021ifn5, banchereau2016personalized, petri2019stability, farahbod2020composition, whitney2003individuality, finn2015fingerprinting, mishra2026temporal.

## 11. Revisions

- 2026-11-16. Re-scope from the pregnancy project after its third pilot's close call and a review of candidate ideas (`log/01-idea-review.md`, re-scope section). Stop 1 is held for the user and the supervisor.
- 2026-11-16, stop 1 answered by the user's instruction to start stage 2, without answers to section 9; the agent took these defaults, recorded as D-R5 in `log/01-idea-review.md`: (1) the user confirms the re-scope with the supervisor; (2) aim 5 kept, with the next visit's SLEDAI as the outcome; (3) GSE65391 alone, ILLUMINATE not used; (4) the GSE65391 download into `data/raw/` approved, not committed. The reference check (step 9) found every cited record opened and one core full text pending (devlin1975robust).
- 2026-11-16, at stage 2. Section 4's neutrophil and lymphocyte stabilities (0.22, 0.30) were wrong: the scratch script misaligned values when a count was missing. The tested package code gives 0.41 and 0.48 (`pilot-d/results/data_summary.rds`); the IFN (0.55) and SLEDAI (0.39) values stand. So $\kappa_C^2$ is about 0.17 to 0.23, not 0.05 to 0.09, and the cell mix adds both stable and changing parts to edges.

## 12. Revision 1 (2026-11-16): network metrics, after the supervisor's suggestion

This section governs where it differs from sections 1 to 11. Pilot-d (`pilot-d/results/verdict.rds`, a close call the user read as a pass for the method claims) found single-sample edges unreliable over visits (intraclass correlation 0.07 to 0.22). The supervisor suggested that edges are too noisy and that metrics built from each person's LIONESS network (connectivity, betweenness, distance, and others) may be stable enough to study over time. The aims become: (1) what these metrics measure, through the closed form where one exists (module connectivity is close to the square of the module's mean standardized expression; a gene's strength is close to its level times the person's summed expression); (2) their stability over visits, for raw LIONESS, first-order, and partial-correlation edges; (3) whether that stability is the expression score's; (4) whether the correction keeps a real personal module. Metrics: module connectivity and between-module connectivity; node strength, betweenness, eigenvector centrality, and local clustering; average shortest path length, global efficiency, and density; hub overlap between visits (reported). Distances are the inverse absolute edge; thresholded metrics keep edges with absolute value above 1. Prior work: LIONESS degree has been used as a survival biomarker (PMC12064794, seen in a search result, not opened), and the test-retest reliability of graph metrics is a large literature for brain networks; neither was found for blood LIONESS networks over repeated samples. The pilot is `pilot-e/`. Prediction of the next visit's activity is left to stage 3.

## 13. Revision 2 (2026-11-16): the own-module correction

This section governs where it differs from section 12. Pilot-e failed its deciding clause 3(a) (`pilot-e/results/verdict.rds`): conditioning on a five-gene interferon score did not correct a module's connectivity, because the score measures the factor with error and carries the person's own gene levels (`pilot-e/decisions.md`, D2). The user chose a re-pilot with a correction aimed at that cause ("option 1"): each gene is regressed on the mean standardized expression of the other genes of its own set (a leave-one-gene-out set score) and on the measured cell mix, and LIONESS edges are built from the residuals. The question is whether this removes a shared factor's stability from a set's connectivity while keeping a real personal module, which the same regression could also remove. The pilot is `pilot-f/`; the pass rule is pilot-e's with this correction in place of the partial-correlation one.

## 14. Revision 3 (2026-11-16): a pathway approach

This section governs where it differs from sections 12 and 13. The supervisor asked for LIONESS in a pathway approach; the user: "do it". Pilots e and f (`pilot-e/`, `pilot-f/`) showed that a gene set's connectivity is its squared expression score and that no linear adjustment separates a network part from it. The pathway approach summarizes each person's LIONESS network by pathway: the connectivity within each pathway (mean edge) and the crosstalk between each pair of pathways (mean edge between their genes). The first-order form gives both in closed form (agent's derivation; section 14.1), and the pilot (`pilot-g/`) tests it with six simulated pathways that mirror blood biology; real gene-set collections (blood transcriptional modules, MSigDB Hallmark, or Reactome, after a licence check) are for stage 3.

### 14.1 The closed forms (sketch, the agent's)

For pathway $A$ with $m_A$ genes, $S_A = \sum_{i \in A} z_i$, $Q_A = \sum_{i \in A} z_i^2$, $R_i^{B} = \sum_{j \in B,\, j \ne i} r_{ij}$, and $C_A = m_A (m_A - 1) / 2$, averaging the first-order edge over pairs gives

$$
\begin{aligned}
\text{within}_A &= \bar r_A + \frac{S_A^2 - Q_A}{2 C_A} - \frac{1}{2 C_A} \sum_{i \in A} R_i^{A} z_i^2, \\
\text{crosstalk}_{AB} &= \bar r_{AB} + \frac{S_A S_B}{m_A m_B} - \frac{1}{2 m_A m_B} \left( \sum_{i \in A} R_i^{B} z_i^2 + \sum_{j \in B} R_j^{A} z_j^2 \right).
\end{aligned}
$$

Line 1 uses $\sum_{i < j} z_i z_j = (S_A^2 - Q_A)/2$ and collects the $-\tfrac{r_{ij}}{2}(z_i^2 + z_j^2)$ terms by gene; line 2 averages over the $m_A m_B$ cross pairs. The reference correlations are the same for every sample, so a sample's pathway connectivity is, up to $O(N^{-1})$, a fixed linear function of three expression terms ($S_A^2$, $Q_A$, and $\sum_i R_i^A z_i^2$), and its crosstalk of three ($S_A S_B$ and the two weighted sums of squares). Prediction: LIONESS pathway connectivity and crosstalk carry nothing beyond these expression terms, even when a pathway has real personal co-expression; that co-expression shows through the terms themselves. Numerical check: the regression of LIONESS connectivity on its three terms should have $R^2$ near 1, and its residual no stability over visits.
