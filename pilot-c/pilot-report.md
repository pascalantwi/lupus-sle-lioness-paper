# Third pilot report: every visit, measured cell mix, and activity within woman

2026-11-13

# Summary

**Scope.** This document reports the third pilot[^1] of the project, for
revision 2 of the proposal (`proposal.md`, section 13), on the user’s
three aims: pregnancy stage, disease activity at each visit, and whether
edges change differently in SLE than in RA across pregnancy. It changes
three things from the second pilot (`pilot-b/`): every visit (five per
woman), the cell mix measured from marker genes, and activity analysed
within each woman; the pass rule is unchanged. It covers the idea, how
the pilot was run, what it found, the project spec, and the decisions
made along the way. It was built from `project-spec.md`, `decisions.md`,
the saved results in `pilot-c/results/`, the factor shifts measured on
GSE235508 (`results/factor_shift.rds`,
`results/factor_shift_markers.rds`), and, read-only, the earlier pilots’
results in `pilot/results/` and `pilot-b/results/`; every number is
computed from those files when the document renders.

**Status.** In `STATUS.md` the third pilot is marked
`piloted-close-call`[^2]. Every clause the rule names as deciding passed
again, and the detection of a real SLE-only change now passed; the two
clauses on the size of the artifact at the study’s levels and the clause
on activity-linked co-expression failed.

| Pilot question | Answer | Key evidence |
|----|----|----|
| 1\. Does the fix work? | pass, pass, pass | Raw edges called 0.179 (MCSE 0.0032) of factor pairs changed at the stress level; partial edges gave 0.0466 to 0.0496 across the five factor-only cells (nominal 0.05), and kept 0.277 (MCSE 0.0081) of a real module change against 0.281 (MCSE 0.0082) for raw edges. |
| 2\. Does it matter at the study’s levels? | fail, fail, pass | False rewiring: 0.046 (MCSE 0.0008) at the SLE level, 0.059 (MCSE 0.0013) at the RA level; false SLE-against-RA difference 0.059 (MCSE 0.0011); a real SLE-only change found in 0.180 (MCSE 0.0056) of module pairs by partial edges. |
| 3\. Prediction | pass, pass, fail | N1 against F1 ranged over -0.0001 to 0.0003 in cells 1 to 7; for activity, N1 against E1 was 0.0013 (MCSE 0.0009) when only the IFN factor’s spread tracked activity, and 0.0027 (MCSE 0.0012) when the module’s correlation did. |

**What needs your decision** (details in [Action items and
limitations](#action-items-and-limitations)): the status call; how the
pregnancy work enters the one paper, given that the artifact is small at
this dataset’s levels; whether aim 2 stays; and the open `[TO CONFIRM]`
items in the spec (highlighted in yellow).

# The idea

## In plain language

LIONESS[^3] gives every blood sample a value for every pair of genes,
called an edge. When the edges of the same women differ between
pregnancy and the months after birth, the difference is read as the
immune network “rewiring”. But an edge is, almost exactly, a simple
product of the sample’s own two gene values, so anything that moves many
genes together, such as the interferon[^4] response or the mix of white
cells, changes the edges by its square. Suppose two genes that follow
the white-cell mix are each 1 unit above average after birth and 1.5
units above in pregnancy: their product rises from 1 to 2.25, and the
edge “rewires” with no change in how the genes work together.

The fix tested here removes the interferon score and the cell mix from
each gene before the edge is built (a partial correlation[^5]). The
pilot simulates women followed through pregnancy, with interferon and
cell-mix changes as large as those measured in the real study, and asks
three things: whether the usual before-and-after test of edges calls
false rewiring and whether the fix stops it; whether edges predict
pregnancy stage or disease activity beyond simpler models; and whether a
false difference between lupus and arthritis appears and disappears with
the fix. No edge was computed on patient data.

## Potential impact

LIONESS builds one network per sample and has an R implementation (M. L.
Kuijjer et al. 2019; Marieke L. Kuijjer et al. 2019); its authors
propose sample networks for studying change across time. The first-order
term of the edge is the influence function of the correlation (Devlin et
al. 1975), the product term is the iENA edge (Yu et al. 2017), and a
recent review relates sample edges to covariance products (Kuijjer et
al. 2026). The study behind the data described pregnancy in SLE and RA
by differential expression and deconvolution (Lien et al. 2024); no
single-sample network analysis of pregnancy was found at stage 1.

The problem is the reading of a change in edges between time points, or
between diseases, as rewiring, when the factors that move many genes
change too. With twice the replicates of the second pilot, the fix held
in every cell: at a large factor change raw edges called 0.179 (MCSE
0.0032) of the pairs whose co-expression never changed “rewired”, and
the partial-correlation edges 0.048 (MCSE 0.0006), the test’s own level;
across the five factor-only cells the fixed share stayed within 0.0466
to 0.0496 (<a href="#tbl-shares" class="quarto-xref">Table 5</a>). The
fix kept real change: 0.277 (MCSE 0.0081) of module pairs within a
disease, and 0.180 (MCSE 0.0056) for a real SLE-only change between
diseases.

The pilot also sharpens what drives the artifact. Measured with marker
genes, pregnancy shifts the cell mix by about 1.4 standard deviations in
RA with a slightly narrower spread
(<a href="#tbl-markers" class="quarto-xref">Table 3</a>), and raw edges
then called only 0.059 (MCSE 0.0013) of factor pairs changed, against
0.107 in the second pilot, whose proxy had a wider spread. A shift of
the mean largely cancels in a before-and-after comparison; what remains
follows the first project’s unequal-groups law and the change in spread
(decisions log, D2). For users, the message is that rewiring over time
is inflated by changes in the spread of a shared factor and by
unbalanced visits, and that the partial-correlation edges remove it in
either case.

## Application: the study behind the data

### The study

The data come from a Norwegian longitudinal study of pregnancy in
rheumatic disease (Lien et al. 2024),
https://doi.org/10.1136/ard-2023-224644, GEO series GSE235508[^6]. Women
with systemic lupus erythematosus (SLE)[^7], rheumatoid arthritis
(RA)[^8], and healthy women gave whole-blood samples for RNA sequencing
before pregnancy, in each trimester, and at 6 weeks, 6 months, and 12
months after birth, with disease activity recorded at each visit
(DAS28[^9] in RA, LAI-P[^10] in SLE). The paper reports that the
interferon signature was augmented in SLE pregnancy and that the altered
cell type composition strongly influences the profile.
<a href="#tbl-data" class="quarto-xref">Table 1</a> counts the women
with samples in both periods, and
<a href="#tbl-shift" class="quarto-xref">Table 2</a> and
<a href="#tbl-markers" class="quarto-xref">Table 3</a> give the size of
the factor changes measured in this project from the counts: the
principal-component proxy of the second pilot, and the marker-gene
scores used here.

<div id="tbl-data">

Table 1: GSE235508 by group, from `../pilot/results/data_summary.rds`
(read-only). Usable pregnancies have at least one pregnant and one
non-pregnant sample.

<div class="cell-output-display">

| Group           | Samples | Pregnancies | Usable pregnancies |
|:----------------|--------:|------------:|-------------------:|
| Healthy         |      93 |          21 |                 20 |
| SLE             |     101 |          25 |                 18 |
| Seronegative RA |      48 |          15 |                  9 |
| Seropositive RA |      93 |          23 |                 18 |
| RA, both types  |     141 |          38 |                 27 |

</div>

</div>

<div id="tbl-shift">

Table 2: The change in pregnancy of the IFN5 interferon score and of a
composition proxy (the first principal component of the 2,000 most
variable genes), in units of each disease’s non-pregnant samples, from
`../results/factor_shift.rds`. These set the SLE and RA levels of the
pilot.

<div class="cell-output-display">

| Disease | Score | Shift | Spread ratio | Within-woman share of variance | Women |
|:--------|:------|:------|:-------------|:-------------------------------|------:|
| SLE     | IFN5  | 0.19  | 1.08         | 0.93                           |    18 |
| SLE     | PC1   | -0.11 | 1.00         | 0.56                           |    18 |
| RA      | IFN5  | -0.29 | 0.63         | 0.38                           |    27 |
| RA      | PC1   | 0.45  | 1.45         | 0.02                           |    27 |

</div>

</div>

<div id="tbl-markers">

Table 3: The change in pregnancy of the marker-gene scores, in units of
each disease’s non-pregnant samples, from
`../results/factor_shift_markers.rds`; the composition axis (neutrophil
minus lymphocyte score) sets the composition levels of this pilot.

<div class="cell-output-display">

| Disease | Score | Shift | Spread ratio | Within-woman share of variance | Women |
|:---|:---|:---|:---|:---|---:|
| SLE | neutrophil | 0.43 | 0.92 | 0.76 | 18 |
| SLE | lymphocyte | -0.54 | 0.92 | 0.85 | 18 |
| SLE | composition axis | 0.58 | 0.84 | 0.80 | 18 |
| RA | neutrophil | 0.65 | 0.98 | 0.43 | 27 |
| RA | lymphocyte | -1.40 | 1.03 | 0.29 | 27 |
| RA | composition axis | 1.37 | 0.86 | 0.08 | 27 |

</div>

</div>

### Why the share of pairs called changed

The quantity that decides the first and third aims is the share of gene
pairs that a test calls changed: across pregnancy within the same women
(aim 1), or differently in SLE than in RA (aim 3). With no change in
co-expression the share should be 5%, the test’s level; a higher share
is false rewiring. The second aim is judged by the concordance
index[^11] of predicted against observed activity, the continuous
counterpart of the AUC[^12].

### Did the original study compute it?

No. The original study used differential expression and cell-type
deconvolution, and reported no expression pattern attributable to
disease activity; it computed no co-expression or single-sample network
and made no test of rewiring.

## In symbols

The standardized expression of gene $i$ in sample $q$ ($z_{qi}$), the
reference correlation ($r$), the shared factors ($F_q$, $C_q$, with
loadings $\lambda$), and the conditioning set ($S_q$, the IFN score and
composition) determine the LIONESS edge ($e_q$), its first-order form,
its mean given the factor, and the partial-correlation edge
($e^{\text{p}}_q$):

$$
\begin{aligned}
e_q &= N r - (N-1)\, r_{(-q)} = r + z_{qi} z_{qj} - \frac{r}{2}\left(z_{qi}^2 + z_{qj}^2\right) + O_p\left(N^{-1}\right), \\
\mathrm{E}\left[e_q \mid F_q\right] &= \text{const} + \left[\lambda_i \lambda_j - \frac{r}{2}\left(\lambda_i^2 + \lambda_j^2\right)\right] F_q^2 + O\left(N^{-1}\right), \\
e^{\text{p}}_q &= N \rho - (N-1)\, \rho_{(-q)}, \quad \rho = \mathrm{cor}\left(x_i - \hat x_i(S), \, x_j - \hat x_j(S)\right).
\end{aligned}
$$

| Symbol | Meaning |
|----|----|
| $z_{qi}$ | expression of gene $i$ in sample $q$, standardized by the reference |
| $r$, $r_{(-q)}$ | correlation of genes $i$ and $j$ in the reference, and without sample $q$ |
| $F_q$, $\lambda_i$ | shared factor in sample $q$ and gene $i$’s loading on it |
| $S_q$ | the conditioning set: the IFN score and composition |
| $\rho$, $e^{\text{p}}_q$ | partial correlation given $S$, and its LIONESS edge |

The target is the change in correlation between periods; the raw edge’s
change adds the bracket times the change in $\mathrm{E}[F^2]$. The
partial edge removes $F$ from each gene first, so its change carries no
$F^2$ term, up to the error with which the score measures $F$.

<div id="fig-idea">

![](pilot-report_files/figure-commonmark/fig-idea-1.png)

Figure 1: The idea, an illustration on made-up numbers (not pilot
output): 40 simulated women at the stress level, no change in
co-expression. Each line is one woman’s mean edge for two IFN genes
outside and in pregnancy; left, the raw LIONESS edge; right, the
partial-correlation edge.

</div>

<div id="fig-data">

![](pilot-report_files/figure-commonmark/fig-data-1.png)

Figure 2: The real data (GSE235508, from `../results/factor_shift.rds`
and `../results/factor_shift_markers.rds`): the change in pregnancy of
the IFN5 score and of the marker-gene composition axis, within SLE and
within RA. Bars are the shift in units of the non-pregnant standard
deviation; the numbers above are the spread ratios.

</div>

## Critiquing the model

Any model is a simplification, and the useful question is not whether it
is true but where it bends and whether the bend matters for the question
being asked. The model above was chosen so that the problem could be
written as a formula and checked by simulation, and each simplification
that made that possible is a place where real data may behave
differently. The table lists the main ones. Work through it one row at a
time. The first column names the feature and the symbols in the
equations that carry it, so find them there first. Then decide in your
own words how the feature could fail for real data before reading the
second column, and finally read the third column for what should be done
to check it. Several of these checks are natural first tasks for the
project. Which of them the pilot examined, and which it left open, is
for you to find when you read the Methods and Results that follow.

| Feature of the model | Why it may not be realistic | What should be done to check it |
|----|----|----|
| Two shared factors with equal loadings (0.8 IFN, 0.6 composition) | Blood has several cell types and IFN subtypes, with different loadings per gene | Estimate the factors and loadings from GSE235508 and plant changes in real data (a plasmode) |
| Composition from one neutrophil-minus-lymphocyte marker score | Monocytes, eosinophils, and platelets also change in pregnancy | Add a score per cell type, or a published deconvolution |
| Composition observed without error | Composition is estimated, so the partial edge removes it only partly | Add measurement error to $C$ |
| Five samples per woman (three and two, or two and three) | The study’s visits per woman vary from two to seven | Mirror each woman’s real time points |
| Activity linked to the module through a smooth function of the visit’s activity | Real activity may change co-expression in bursts (flares) | Simulate flares as rare visits with a large change |
| Normal noise | Log counts have heavier tails | Resample residuals from the real data |

# Methods

## Design

The pilot simulates women followed through pregnancy in 10 cells, three
aims. Each woman gives five samples (three pregnant and two not at the
SLE, null, and stress levels; two and three at the RA levels). Fifty
genes: genes 1 to 20 follow the IFN factor $F$ (the IFN score is the
mean of genes 1 to 5), genes 21 to 30 form a co-expressed module, and
genes 31 to 50 follow the composition factor $C$. Both factors may shift
and spread in pregnancy, by the amounts measured in SLE or in RA
(<a href="#tbl-shift" class="quarto-xref">Table 2</a>), or by a large
stress amount, with the composition levels from the marker-gene axis
(<a href="#tbl-markers" class="quarto-xref">Table 3</a>). In aim 1 and
aim 3 cells, LIONESS edges (raw) and partial-correlation edges (the fix)
are computed on each disease’s samples, and each woman’s change in edge
between periods is tested. In aim 1 and aim 2 cells, six models predict
pregnancy stage or activity with 5-fold cross-validation[^13] by woman;
for activity, the activity and every feature are first centered on each
woman’s own mean, so the models predict each visit’s departure from her
usual level.

<div id="fig-design">

![](pilot-report_files/figure-commonmark/fig-design-1.png)

Figure 3: The design at a glance. Each row is one part of the pilot,
read from left to right: what is simulated, which factors are varied,
what is estimated, what is measured, and how the result is shown. The
shares are in <a href="#tbl-shares" class="quarto-xref">Table 5</a> and
<a href="#fig-shares" class="quarto-xref">Figure 4</a>, the prediction
results in <a href="#tbl-models" class="quarto-xref">Table 6</a> and
<a href="#fig-models" class="quarto-xref">Figure 5</a>, the verdict in
<a href="#tbl-verdict" class="quarto-xref">Table 7</a>, and the
known-answer checks in
<a href="#tbl-known-answer" class="quarto-xref">Table 4</a>.

</div>

Each figure and table in the Results traces back to one row of this
diagram.

| Part | Cells | Main factors | Replicates |
|----|----|----|----|
| A | 1 to 5 | levels none, SLE, RA, stress; module +0.4 | 400 per cell |
| B | 6 and 7 | activity linked to the IFN factor’s spread or to the module | 400 per cell |
| C | 8 to 10 | SLE against RA; module +0.4 in SLE; stress | 400 per cell |

The six models are $R_0$ (composition and the residual IFN score), $R_1$
($R_0$ with both squares), and elastic nets[^14] with the $R_1$ terms
unpenalized: E1 (the 45 non-score genes), F1 (E1 with the 990
first-order edges), N1 (E1 with the 990 LIONESS edges), and N1 with
$R_0$ (the edges on $R_0$, to show the artifact). Pregnancy stage is
scored by the AUC, activity by the concordance index of the within-woman
activity, with the correlation, out-of-fold $R^2$, and root mean squared
error reported.

## Known-answer checks

Known-answer checks come before any result, because every later number
is only as good as the simulation that produced it. Each check sets up a
case where the correct answer is known from theory, runs the simulation
code on it, and compares the two. Agreement says nothing about the
scientific question or about real data. It says that the code does what
the model says. Disagreement means stopping to look for a bug before
reading anything else, and the same checks apply to any code you write
for the project.

Read each row of the table in four steps: the setting that was
simulated, the known answer for that setting, what the simulation gave,
and whether the gap between the two is small compared with the Monte
Carlo standard error (MCSE). The MCSE is how much a simulated result
moves from run to run because only a finite number of replicates was
used. A gap within about two MCSE is noise, and a gap of several MCSE is
a problem.

<div id="tbl-known-answer">

Table 4: Known-answer checks, from `results/known_answers.rds`.

<div class="cell-output-display">

| Check | Setting | Known answer | Simulation gave | How to read it |
|----|----|----|----|----|
| Exact LIONESS Pearson edge equals leave-one-out LIONESS | N = 12, 6 genes, one outlier; largest absolute difference | 0.000 | 8.88e-15 | 8.882e-15 from the known value, against a tolerance of 1e-08, so the check passes. |
| Exact add-one edge equals (N + 1) r(+q) - N r | reference N = 12, 5 new samples; largest absolute difference | 0.000 | 5.33e-15 | 5.329e-15 from the known value, against a tolerance of 1e-08, so the check passes. |
| Exact LIONESS Pearson edge equals leave-one-out LIONESS | N = 60, 6 genes, one outlier; largest absolute difference | 0.000 | 1.42e-14 | 1.421e-14 from the known value, against a tolerance of 1e-08, so the check passes. |
| Exact add-one edge equals (N + 1) r(+q) - N r | reference N = 60, 5 new samples; largest absolute difference | 0.000 | 1.07e-14 | 1.066e-14 from the known value, against a tolerance of 1e-08, so the check passes. |
| First-order edge approaches the exact edge as N grows | N = 4000, r = 0.5; largest absolute difference | 0.000 | 0.0123 | 0.01235 from the known value, against a tolerance of 0.05, so the check passes. |
| AUC of perfectly separated scores | 5 against 5 | 1.000 | 1 | 0 from the known value, against a tolerance of 1e-08, so the check passes. |
| AUC of tied scores | all scores equal | 0.500 | 0.5 | 0 from the known value, against a tolerance of 1e-08, so the check passes. |
| Folds keep each woman’s samples together | 20 women, 4 samples each; women split across folds | 0.000 | 0 | 0 from the known value, against a tolerance of 1e-08, so the check passes. |
| Partial correlation equals the inverse-covariance formula | 30 samples, 4 genes, 2 conditioning variables; largest difference | 0.000 | 2.07e-16 | 2.073e-16 from the known value, against a tolerance of 1e-08, so the check passes. |
| Partial-correlation edges with no conditioning equal Pearson edges | 30 samples, 4 genes; largest difference | 0.000 | 5.33e-15 | 5.329e-15 from the known value, against a tolerance of 1e-08, so the check passes. |
| Concordance of a perfect prediction | 10 distinct outcomes | 1.000 | 1 | 0 from the known value, against a tolerance of 1e-08, so the check passes. |
| Concordance of a constant prediction | 10 distinct outcomes | 0.500 | 0.5 | 0 from the known value, against a tolerance of 1e-08, so the check passes. |
| Concordance equals the AUC for a binary outcome | 5 against 5 | 0.000 | 0 | 0 from the known value, against a tolerance of 1e-08, so the check passes. |
| Share called changed (raw factor) in the null cell | cell 1, 400 replicates | 0.050 | 0.0471 (MCSE 0.0008) | 3.8 MCSE from the known value, so the gap is too large: look for a bug. |
| Share called changed (raw module) in the null cell | cell 1, 400 replicates | 0.050 | 0.0496 (MCSE 0.0019) | 0.2 MCSE from the known value, so the gap is noise. |
| Share called changed (par factor) in the null cell | cell 1, 400 replicates | 0.050 | 0.0479 (MCSE 0.0006) | 3.3 MCSE from the known value, so the gap is too large: look for a bug. |
| Share called changed (par module) in the null cell | cell 1, 400 replicates | 0.050 | 0.0494 (MCSE 0.0019) | 0.3 MCSE from the known value, so the gap is noise. |
| AUC of R0 in the null cell | cell 1, 400 replicates | 0.500 | 0.4896 (MCSE 0.0033) | 3.2 MCSE from the known value, so the gap is too large: look for a bug. |
| AUC of R1 in the null cell | cell 1, 400 replicates | 0.500 | 0.4892 (MCSE 0.0036) | 3.0 MCSE from the known value, so the gap is too large: look for a bug. |
| AUC of E1 in the null cell | cell 1, 400 replicates | 0.500 | 0.4907 (MCSE 0.0036) | 2.6 MCSE from the known value, so the gap is too large: look for a bug. |
| AUC of F1 in the null cell | cell 1, 400 replicates | 0.500 | 0.4899 (MCSE 0.0035) | 2.8 MCSE from the known value, so the gap is too large: look for a bug. |
| AUC of N1 in the null cell | cell 1, 400 replicates | 0.500 | 0.4902 (MCSE 0.0036) | 2.8 MCSE from the known value, so the gap is too large: look for a bug. |
| AUC of N1R0 in the null cell | cell 1, 400 replicates | 0.500 | 0.4888 (MCSE 0.0033) | 3.4 MCSE from the known value, so the gap is too large: look for a bug. |

</div>

</div>

0 of 10 rows with a standard error are more than 2 MCSE from the known
value.

# Results

## How to read the outputs

Each output traces back to a row of the design diagram (the captions
name the part where one applies). Read each in the same four steps.
First, find its row in the diagram, which says what was simulated and
what was varied. Second, work out what one point or one table entry is,
usually one combination of the varied factors, averaged over replicates.
Third, identify the reference value that would hold if nothing were
wrong, shown as a dashed or grey line wherever possible. Fourth, look
for where the results leave that reference, and for which settings they
do. The outputs report three kinds of quantity.

- **Bias** is the average of an estimate minus its true value. When the
  truth is no difference, any average difference from 0 is a spurious
  difference created by the problem under study.
- **Level** is how often a test rejects when there is truly nothing to
  find. A correct test at 0.05 rejects 5% of the time, so a rate above
  0.05 means too many false alarms.
- **Power** is how often a test rejects when there is a real effect.
  More is better, but only when the level is not inflated, because a
  test that rejects too often also looks powerful.

All rates and averages come from simulation, so each carries Monte Carlo
error. That is why the figures show bars or a grey band, and why small
differences should not be over-read.

This pilot reports shares of gene pairs called changed (a rate whose
correct value with no change in co-expression is 0.05), the AUC, the
concordance index, and paired differences between models. A paired
difference is computed within each replicate and then averaged; its MCSE
is the replicate standard deviation of the difference divided by the
square root of the number of replicates.

## False rewiring and the fix

**Figure guide (parts A and C).** Each row is one cell; the left panel
shows factor pairs (genes whose co-expression never changes), the right
panel module pairs (genes whose co-expression changes in cells 5 and 9).
Red circles are raw LIONESS edges, blue triangles partial-correlation
edges; the dashed line is the 5% a correct test should give with no
change. Look for: whether the red circles sit right of the line in the
factor panel; whether the blue triangles sit on it; whether the blue
triangles stay close to the red circles in the module panel where the
change is real.

<div id="fig-shares">

![](results/shares.png)

Figure 4: Third pilot result (parts A and C): share of gene pairs called
changed, raw against partial-correlation edges, means over replicates,
from `results/cell_summaries.rds`.

</div>

## Prediction

**Figure guide (parts A and B).** Each row is one prediction cell; the
points are paired differences in AUC (cells 1 to 5) or within-woman
concordance (cells 6 and 7) averaged over 400 replicates, with bars of
two MCSE: N1 with $R_0$ against $R_0$ (the artifact), N1 against E1, N1
against F1. Look for: whether the purple points sit on the line in every
row; where the red points sit away from it; whether the blue points move
away from the line in cell 7.

<div id="fig-models">

![](results/model_gains.png)

Figure 5: Third pilot result (parts A and B): paired model differences
by cell, from `results/cell_summaries.rds`.

</div>

## Numbered findings

Each finding is a claim followed by the output behind it. Read the
claim, open the output, and check that you can see the claim there
before accepting it.

1.  **Raw edges call false rewiring when the factor’s spread changes,
    and little when only its mean moves.** With no change in
    co-expression, the share of factor pairs called changed was 0.046
    (MCSE 0.0008) at the SLE level, 0.059 (MCSE 0.0013) at the RA level,
    and 0.179 (MCSE 0.0032) at the stress level, where the IFN factor’s
    spread more than doubles
    (<a href="#tbl-shares" class="quarto-xref">Table 5</a>, red circles
    in <a href="#fig-shares" class="quarto-xref">Figure 4</a>).
2.  **The fix holds the test’s level everywhere.** Partial-correlation
    edges gave 0.0466 to 0.0496 in every factor-only cell, within a
    disease and between diseases
    (<a href="#tbl-verdict" class="quarto-xref">Table 7</a>, clause
    1(b)).
3.  **The fix keeps a real change.** Partial edges called 0.277 (MCSE
    0.0081) of module pairs changed within a disease, raw edges 0.281
    (MCSE 0.0082); between diseases, a real SLE-only change was found in
    0.180 (MCSE 0.0056) of module pairs, now above the 0.15 threshold
    (clause 2(c)).
4.  **At the measured levels a false SLE-against-RA difference is
    small.** Raw edges called 0.059 (MCSE 0.0011) of factor pairs
    differently changed, below the 0.10 threshold (clause 2(b) fail),
    and partial edges 0.049 (MCSE 0.0006).
5.  **Edges equal their formula in prediction.** N1 and F1 differed by
    -0.0001 to 0.0003 in AUC or concordance across cells 1 to 7 (clause
    3(a)).
6.  **Within woman, the squared-score base removes the activity
    artifact, and edges find no activity-linked co-expression.** When
    activity raised only the IFN factor’s spread, N1 with $R_0$ gained
    0.026 (MCSE 0.0019) over $R_0$ and N1 added 0.0013 (MCSE 0.0009)
    over E1; when activity raised the module’s correlation, N1 added
    0.0027 (MCSE 0.0012) (clause 3(c) fail), with every model’s
    concordance near 0.5
    (<a href="#tbl-models" class="quarto-xref">Table 6</a>).
7.  **The known-answer checks hold**
    (<a href="#tbl-known-answer" class="quarto-xref">Table 4</a>).

## Tables

**Table guide (shares).** One row per testing cell: the mean share of
pairs called changed, with its MCSE, for raw and partial edges, factor
and module pairs. The reference value is 0.05 wherever co-expression
does not change.

<div id="tbl-shares">

Table 5: Share of gene pairs called changed (MCSE), from
`results/cell_summaries.rds`. Cells 1 to 5 use the paired test within a
disease, cells 8 to 10 the test between diseases.

<div class="cell-output-display">

| Cell | Raw, factor | Fix, factor | Raw, module | Fix, module |
|:---|:---|:---|:---|:---|
| cell 1 (aim 1: none levels, 18 women) | 0.047 (MCSE 0.0008) | 0.048 (MCSE 0.0006) | 0.050 (MCSE 0.0019) | 0.049 (MCSE 0.0019) |
| cell 2 (aim 1: SLE levels, 18 women) | 0.046 (MCSE 0.0008) | 0.047 (MCSE 0.0006) | 0.049 (MCSE 0.0020) | 0.047 (MCSE 0.0019) |
| cell 3 (aim 1: RA levels, 27 women) | 0.059 (MCSE 0.0013) | 0.050 (MCSE 0.0006) | 0.049 (MCSE 0.0020) | 0.049 (MCSE 0.0020) |
| cell 4 (aim 1: stress levels, 18 women) | 0.179 (MCSE 0.0032) | 0.048 (MCSE 0.0006) | 0.050 (MCSE 0.0019) | 0.050 (MCSE 0.0019) |
| cell 5 (aim 1: SLE levels, 18 women; module +0.4 in pregnancy) | 0.047 (MCSE 0.0007) | 0.047 (MCSE 0.0006) | 0.281 (MCSE 0.0082) | 0.277 (MCSE 0.0081) |
| cell 8 (aim 3: SLE levels (18 women) against RA levels (27 women)) | 0.059 (MCSE 0.0011) | 0.049 (MCSE 0.0006) | 0.048 (MCSE 0.0018) | 0.047 (MCSE 0.0018) |
| cell 9 (aim 3: SLE levels (18 women) against RA levels (27 women); module +0.4 in SLE pregnancy) | 0.062 (MCSE 0.0011) | 0.049 (MCSE 0.0006) | 0.181 (MCSE 0.0057) | 0.180 (MCSE 0.0056) |
| cell 10 (aim 3: stress levels (18 women) against none levels (27 women)) | 0.102 (MCSE 0.0020) | 0.048 (MCSE 0.0006) | 0.051 (MCSE 0.0019) | 0.050 (MCSE 0.0019) |

</div>

</div>

**Table guide (prediction).** One row per prediction cell; the columns
are the mean out-of-fold AUC (cells 1 to 5) or concordance index (cells
6 and 7) of each model. The reference value is 0.5 for no prediction.

<div id="tbl-models">

Table 6: Mean out-of-fold AUC (cells 1 to 5) or concordance index (cells
6 and 7) by model, from `results/cell_summaries.rds`. The largest MCSE
of any entry is given below the table.

<div class="cell-output-display">

| Cell | R0 | R1 | E1 (genes) | F1 (E1 + first-order edges) | N1 (E1 + edges) | N1 with R0 |
|:---|:---|:---|:---|:---|:---|:---|
| cell 1 (aim 1: none levels, 18 women) | 0.490 | 0.489 | 0.491 | 0.490 | 0.490 | 0.489 |
| cell 2 (aim 1: SLE levels, 18 women) | 0.669 | 0.652 | 0.651 | 0.650 | 0.650 | 0.666 |
| cell 3 (aim 1: RA levels, 27 women) | 0.852 | 0.850 | 0.848 | 0.848 | 0.848 | 0.850 |
| cell 4 (aim 1: stress levels, 18 women) | 0.624 | 0.692 | 0.691 | 0.691 | 0.691 | 0.632 |
| cell 5 (aim 1: SLE levels, 18 women; module +0.4 in pregnancy) | 0.667 | 0.653 | 0.652 | 0.657 | 0.657 | 0.668 |
| cell 6 (aim 2: RA levels, 36 women; activity raises the IFN factor’s spread) | 0.496 | 0.559 | 0.556 | 0.557 | 0.557 | 0.522 |
| cell 7 (aim 2: RA levels, 36 women; activity raises the module’s correlation) | 0.497 | 0.498 | 0.497 | 0.500 | 0.500 | 0.499 |

</div>

</div>

The largest MCSE of an entry in the table is 0.0036. In the activity
cells, the out-of-fold $R^2$ of N1 for the within-woman activity was
-0.018 and -0.055, and its root mean squared error 0.636 and 0.648 (the
within-woman activity has variance about 0.4).

**Table guide (the verdict).** One row per clause of the pre-registered
rule, read by `pilot-b/verdict.R`. Check that each outcome follows from
its value and threshold, and that the combination gives the status.

<div id="tbl-verdict">

Table 7: The pass rule applied clause by clause, with the second pilot’s
value and outcome beside each, from `results/verdict.rds`.

<div class="cell-output-display">

| Q | Clause | Cell | Quantity | Value | MCSE | Threshold | Outcome | Pilot-b value | Pilot-b outcome |
|---:|:---|---:|:---|:---|:---|:---|:---|:---|:---|
| 1 | \(a\) | 4 | raw share, factor pairs (paired) | 0.1786 | 0.0032 | \>= 0.10 | pass | 0.1697 | pass |
| 1 | \(b\) | 2 | partial share, factor pairs (paired) | 0.0466 | 0.0006 | \<= 0.0600 | pass | 0.0476 | pass |
| 1 | \(b\) | 3 | partial share, factor pairs (paired) | 0.0496 | 0.0006 | \<= 0.0600 | pass | 0.0485 | pass |
| 1 | \(b\) | 4 | partial share, factor pairs (paired) | 0.0479 | 0.0006 | \<= 0.0600 | pass | 0.0480 | pass |
| 1 | \(b\) | 8 | partial share, factor pairs (between) | 0.0486 | 0.0006 | \<= 0.0600 | pass | 0.0485 | pass |
| 1 | \(b\) | 10 | partial share, factor pairs (between) | 0.0484 | 0.0006 | \<= 0.0600 | pass | 0.0483 | pass |
| 1 | \(c\) | 5 | partial minus 0.75 raw share, module pairs (paired) | 0.0654 | 0.0030 | \>= 0 | pass | 0.0492 | pass |
| 2 | \(a\) | 3 | raw share, factor pairs (paired), larger of cells 2, 3 | 0.0592 | 0.0013 | \>= 0.10 | fail | 0.1075 | close |
| 2 | \(b\) | 8 | raw share, factor pairs (between) | 0.0589 | 0.0011 | \>= 0.10 | fail | 0.0789 | fail |
| 2 | \(c\) | 9 | partial share, module pairs (between) | 0.1801 | 0.0056 | \>= 0.15 | pass | 0.1573 | close |
| 3 | \(a\) | 1 | AUC(N1) - AUC(F1) | 0.0003 | 0.0002 | \|D\| \<= 0.0100 | pass | -0.0003 | pass |
| 3 | \(a\) | 2 | AUC(N1) - AUC(F1) | 0.0002 | 0.0001 | \|D\| \<= 0.0100 | pass | 0.0000 | pass |
| 3 | \(a\) | 3 | AUC(N1) - AUC(F1) | 0.0001 | 0.0000 | \|D\| \<= 0.0100 | pass | 0.0004 | pass |
| 3 | \(a\) | 4 | AUC(N1) - AUC(F1) | 0.0001 | 0.0002 | \|D\| \<= 0.0100 | pass | 0.0004 | pass |
| 3 | \(a\) | 5 | AUC(N1) - AUC(F1) | -0.0001 | 0.0001 | \|D\| \<= 0.0100 | pass | -0.0007 | pass |
| 3 | \(a\) | 6 | C(N1) - C(F1) | 0.0001 | 0.0001 | \|D\| \<= 0.0100 | pass | -0.0002 | pass |
| 3 | \(a\) | 7 | C(N1) - C(F1) | -0.0001 | 0.0001 | \|D\| \<= 0.0100 | pass | -0.0004 | pass |
| 3 | \(b\) | 6 | C(N1) - C(E1) | 0.0013 | 0.0009 | \|E\| \<= 0.0100 | pass | 0.0020 | pass |
| 3 | \(c\) | 7 | C(N1) - C(E1) | 0.0027 | 0.0012 | \>= 0.03 | fail | 0.0034 | fail |

</div>

</div>

**Table guide (numbers quoted in the spec).** Each row recomputes a
number that section 3 of the spec quotes for the third pilot.

<div id="tbl-spec-numbers">

Table 8: Numbers quoted in the spec’s section 3 for the third pilot,
recomputed from `results/cell_summaries.rds` and `results/verdict.rds`.

<div class="cell-output-display">

| Quantity | Value |
|:---|:---|
| 1(a): raw share, factor pairs, cell 4 | 0.179 (MCSE 0.0032) |
| 1(b): range of partial shares, cells 2, 3, 4, 8, 10 | 0.0466 to 0.0496 |
| 1(c): partial and raw module shares, cell 5 | 0.277 (MCSE 0.0081) and 0.281 (MCSE 0.0082) |
| 2(a): raw share, factor pairs, cells 2 and 3 | 0.046 (MCSE 0.0008) and 0.059 (MCSE 0.0013) |
| 2(b): raw between-disease share, factor pairs, cell 8 | 0.059 (MCSE 0.0011) |
| 2(c): partial between-disease share, module pairs, cell 9 | 0.180 (MCSE 0.0056) |
| 3(a): range of N1 - F1, cells 1 to 7 | -0.0001 to 0.0003 |
| 3(b): C(N1) - C(E1), cell 6 | 0.0013 (MCSE 0.0009) |
| 3(c): C(N1) - C(E1), cell 7 | 0.0027 (MCSE 0.0012) |
| Status | piloted-close-call |

</div>

</div>

# The draft spec

This is `project-spec.md`, read from the file when the document renders.
Anything the agent[^15] left for the user to decide is highlighted in
yellow and starts with “To confirm”.

|  |  |
|----|----|
| Proposal | `proposal.md`, final 2026-11-13 |
| Authors | as in `DESCRIPTION`, set by `tools/authors.R` before the final render |
| Written | 2026-11-13 |
| Status | draft; section 3 completed 2026-11-13 after the verdict |

## 1. Question and why it matters

**Question (one sentence):** Do LIONESS edges separate pregnant from
non-pregnant samples of the same women beyond the IFN score, cell
composition, their squares, the gene levels, and their own first-order
features, and does the answer differ between SLE and RA?

**Why it matters:** A change in LIONESS edges over pregnancy is read as
rewiring of immune co-expression. If the mechanism holds, much of that
change is the squared change of the factors pregnancy moves (composition
in every woman, IFN in SLE), and a user would fit the base with the
squared score and the gene levels before claiming that edges add
anything, and would check every network result against the first-order
formula.

## 2. Mechanism we suspect

A LIONESS edge equals, to first order,
$r + z_{qi} z_{qj} - \tfrac{r}{2}(z_{qi}^2 + z_{qj}^2)$. Under a shared
factor $F$ with loadings $\lambda$, its conditional mean is linear in
$F_q^2$ with coefficient
$\lambda_i \lambda_j - \tfrac{r}{2}(\lambda_i^2 + \lambda_j^2)$, so the
difference in mean edge between pregnancy and the postpartum period is
that coefficient times the difference in $\mathrm{E}[F^2]$, with the
period-specific correlation unchanged. The gap is zero when pregnancy
moves neither the mean nor the spread of $F$ (the null cell), and grows
with the shift and the scale of $F$ (the factor cells). If the mechanism
is wrong, edges would gain over E1 in factor-only cells, or edges and
the first-order formula would classify differently.

## 3. Pilot

**Code:** `pilot/pilot.R` (`Rscript pilot/pilot.R --reps <n>`), then
`Rscript pilot/verdict.R`; data description
`Rscript pilot/describe_data.R`.

**Pass rule:** `pilot/preregistration.md`, section 2, copied from
proposal section 7.

**Verdict file:** `pilot/results/verdict.rds`, after the run.

**Plot:** `pilot/results/auc_gains.png`; pilot-b:
`pilot-b/results/shares.png`, `pilot-b/results/model_gains.png`;
pilot-c: `pilot-c/results/shares.png`,
`pilot-c/results/model_gains.png`.

**Second pilot (`pilot-b/`, proposal section 12, the user’s three
aims):** `piloted-close-call` (2026-11-13), 200 replicates per cell, 10
cells, three cores, 965 seconds; verdict file
`pilot-b/results/verdict.rds`. Numbers rounded from
`pilot-b/results/cell_summaries.rds`.

1.  Does the fix work? Yes: raw edges called 0.17 of unchanged factor
    pairs changed at the stress level; partial-correlation edges gave
    0.044 to 0.049 in every factor-only cell, within and between
    diseases, and kept a real module change (0.22 against 0.23).
2.  Does it matter at the study’s levels? Modestly: no false rewiring at
    the SLE level, 0.107 at the RA level (close call), a false
    SLE-against-RA difference of 0.079 (below 0.10, fail), and a real
    SLE-only change found in 0.157 of module pairs (close call).
3.  Prediction: edges equalled their formula (N1 against F1 within
    0.0007); the base with the squared score removed the activity
    artifact; edges did not detect activity-linked co-expression across
    women (fail).

**Third pilot (`pilot-c/`, proposal section 13: every visit, marker-gene
cell mix, activity within woman; pass rule unchanged):**
`piloted-close-call` (2026-11-13), 400 replicates per cell, 10 cells,
three cores, 2,233 seconds; verdict file `pilot-c/results/verdict.rds`.
Numbers rounded from `pilot-c/results/cell_summaries.rds`.

1.  Does the fix work? Yes, again: partial-correlation edges gave 0.047
    to 0.050 in every factor-only cell (raw 0.18 at the stress level),
    and kept real change within a disease (0.28 against 0.28) and
    between diseases.
2.  Does it matter at the study’s levels? Little: with the cell mix
    measured from marker genes, raw edges called 0.046 (SLE level) and
    0.059 (RA level) of factor pairs changed, and 0.059 differently
    changed between diseases (2(a), 2(b) fail); a real SLE-only change
    was found in 0.18 of module pairs (2(c) pass).
3.  Prediction: edges equalled their formula (within 0.0003); the
    squared-score base removed the activity artifact within woman; edges
    found no activity-linked co-expression within woman (3(c) fail).

**Status:** `piloted-close-call` (2026-11-13), 200 replicates per cell,
8 cells, three cores, 1,622 seconds.

**The three pilot answers** (`pilot/results/verdict.rds`; numbers
rounded from `pilot/results/cell_summaries.rds`):

1.  Does the effect exist? Edges equal their formula: N1 and F1 differed
    by at most 0.0011 in AUC in any cell (1(b) pass). With a module
    change at 40 women, N1 gained 0.043 (MCSE 0.0039) over E1, within 3
    MCSE of 0.05 (1(a) close call).
2.  Is it large enough at realistic sizes? Not by classification at 20
    women: the artifact was 0.005 (MCSE 0.0035), below 0.02 (2(a) fail),
    and real rewiring added 0.009 (MCSE 0.0034) over E1, below 0.03
    (2(b) fail). By the paired edge test, yes: with no change in
    co-expression, 0.45 of factor-gene pairs were called rewired at 20
    women under a large factor change, against 0.047 in the null cell.
3.  Is there a fix? For classification, yes: with the E1 base edges
    added between -0.0004 and 0.0015 in the factor-only cells (3(a)
    pass). For the paired edge test, none was tested.

## 4. ADEMP skeleton

**Aims.** (1) Do LIONESS edges and their first-order formula classify
pregnancy stage alike? (2) How much do edges gain from a factor that
changes over pregnancy, and does the E1 base remove it? (3) When
co-expression truly changes, do edges detect it beyond E1 at the study’s
size?

**Data-generating mechanisms.**

- Factors varied, with levels: what pregnancy changes (nothing; the
  factor $F$, moderate or large; the module’s correlation by 0.2 or 0.4;
  both); the number of women (20, 40).
- Factors fixed, and why: four samples per woman (two in each period,
  the median the study provides per disease); 50 genes (as in the first
  project, so results are comparable); loadings 0.8 and 0.4; woman-level
  gene offsets of variance 0.5 (the between-person part of whole-blood
  expression); normal noise (heavy tails were tested in the first
  project).
- A cell where the answer is known: the null cell, every AUC 0.5.

**Estimands.** The out-of-sample AUC of each model for the period of a
new woman’s sample; the paired AUC differences $\Delta_1$, $\Delta_2$,
and the artifact $G$.

**Methods compared.** $R_0$, $R_1$ (logistic); E1, E5, F1, N1, N1 with
$R_0$ (elastic net, `glmnet` `cv.glmnet`, mixing 0.5, 5 inner folds
grouped by woman, binomial deviance, `lambda.min`, base terms
unpenalized). Outer 5-fold cross-validation by woman.

**Performance measures.** Mean out-of-fold AUC and paired differences
with MCSE; the share of gene pairs called different between periods by a
paired test (apparent rewiring).

**Replicates and Monte Carlo SE target.** 200 per cell in the pilot if
the time box allows, an MCSE of about 0.003 for a difference with
standard deviation 0.04.

## 5. Application

**Dataset:** GSE235508 (Lien et al. 2024), GEO; the series matrix in
`data/raw/` for the pilot’s description, the counts at stage 3 after
approval.

**Comparison:** within SLE and within RA separately, pregnant
(trimesters 1 to 3) against non-pregnant (pre-pregnancy and postpartum)
samples of the same women, folds by woman; counts per group from
`pilot/results/data_summary.rds`.

**What to show:** whether the apparent rewiring over pregnancy, and any
gain of edges over the base, disappears once the squared IFN and
composition scores and the gene levels are in the base, and whether that
differs between SLE and RA.

## 6. Target journal

**Target journal:** BMC Bioinformatics

**Fallback:** Statistical Applications in Genetics and Molecular Biology

**The target’s word limit and AI policy:** <span class="tc">**To
confirm:** copied from the journal profile at stage 3</span>

## 7. Milestones

| Stage | Ends with | Due |
|----|----|----|
| 2\. Pilot | `pilot/pilot-report.html`, the verdict file, `pilot/decisions.md` | 2026-11-13 |
| 3\. Manuscript, design stop | `protocol.md` from the approved design | <span class="tc">**To confirm:** date</span> |
| 3\. Manuscript, first draft | `manuscript/manuscript.html` in draft mode | <span class="tc">**To confirm:** date</span> |
| 3\. Manuscript, settled | the draft after the revision rounds | <span class="tc">**To confirm:** date</span> |
| 4\. Background report | `background/background.html` | <span class="tc">**To confirm:** date</span> |
| 5\. Mock review | the last `review/round-N/` | <span class="tc">**To confirm:** date</span> |
| Submission | the final render and the checker in final mode | <span class="tc">**To confirm:** date</span> |

## 8. Out of scope

- A second dataset merged with GSE235508 (the user’s rule: no studies
  combined).
- Pregnancy outcomes or flare prediction (a predictive project of its
  own).
- Other single-sample network methods (SSN, iENA, SWEET, BONOBO):
  covered in the first project.
- Partial-correlation edges and regularized discriminant analysis (the
  first project’s N2 and E4).
- Mixed models of the edges over time: a different estimand, parked as
  its own idea.

## 9. Reading

1.  1.  After reading, a reader can answer: how is a LIONESS sample
        network defined, and what use across time do the authors
        propose?
2.  2.  After reading, a reader can answer: what happens to the IFN
        signature and to cell composition over pregnancy in SLE and in
        RA?
3.  3.  After reading, a reader can answer: what is the influence
        function of the correlation coefficient?

# Decisions log

This is the agent’s `decisions.md` for the third pilot, read from the
file when the document renders.

## D0. The literature pass

As for pilot-b (`pilot-b/decisions.md`, D0). New since: the marker-gene
composition score (`analysis/02_factor_shift_markers.R`); its marker
list carries a `[TO CONFIRM]` marker in proposal section 13 until a
published source is opened.

## D1. The status call: `piloted-close-call`, with both readings

Every deciding clause passed again (1(b), 3(a), 3(b)), and with 400
replicates clause 2(c) moved from a close call to a pass. Clauses 2(a),
2(b), and 3(c) failed. Values, with pilot-b’s beside them, in
`pilot-c/results/verdict.rds`.

- **The reading that passes.** The paper’s claims held in every cell,
  with twice the replicates of pilot-b: partial-correlation edges kept
  the share of factor pairs called changed at the test’s level in every
  factor-only cell, within and between diseases; they kept a real change
  in co-expression within a disease and a real SLE-only change between
  diseases (2(c) now passes); edges equalled their first-order formula
  in every prediction cell; and the base with the squared score removed
  the artifact from activity prediction. The failures of 2(a) and 2(b)
  say that, at the factor changes measured in GSE235508 with marker
  genes, raw edges call little false rewiring: a property of this
  dataset, which the theory explains (D2), not a failure of the method.
- **The reading that fails.** The real-data application of aims 1 and 3
  will show little or no false rewiring to remove in GSE235508, and aim
  2 cannot be answered: edges did not detect an activity-linked change
  in co-expression, across women (pilot-b) or within women (pilot-c).
  The pregnancy data then illustrate the fix weakly.

The agent’s reading: a pass for the method claims and for detecting real
SLE-only change; a fail for a strong real-data illustration of the
artifact in this dataset and for activity prediction. The user decides.

## D2. Why the false rewiring is small at the measured levels: the unequal-groups law

With the marker score, pregnancy shifts the composition axis by about
0.6 (SLE) and 1.4 (RA) standard deviations with a spread ratio of about
0.85 (`results/factor_shift_markers.rds`). To first order, a woman’s
change in edge for two genes loading on a factor $C$ has mean
proportional to $(1 - 2\pi)\Delta^2 + (\sigma_1^2 - \sigma_0^2)$, where
$\pi$ is the share of pregnant samples, $\Delta$ the shift, and
$\sigma_1, \sigma_0$ the spreads: the first project’s unequal-groups law
(its supplement, S1.5) plus the spread term. A pure mean shift cancels
in a balanced design; at the RA level ($\pi = 0.4$) the shift term
($0.2 \times 1.96$) is partly offset by the narrower spread
($0.85^2 - 1$), so the false rewiring is small. Pilot-b’s
principal-component proxy gave a spread ratio of 1.45 at the RA level,
which is where its larger false rewiring came from. This is a prediction
the real data can check at stage 3: the artifact should follow the
spread change and the imbalance, not the shift alone. The derivation
here is a sketch; its numerical check is the comparison of cells 2, 3, 4
across pilots b and c, and a dedicated cell at stage 3.

## D3. Activity within woman

Centering activity and features on each woman removed the downward bias
of pilot-b’s across-woman analysis (concordances now near 0.5 rather
than below it in cell 7), but edges still found no activity-linked
change in co-expression (`cell_summaries.rds`, cell 7). A sample’s edge
is a single noisy product; a change in correlation of at most 0.4 spread
over a woman’s visits is too small to see in single samples at 36 women.
Aim 2 by single-sample edges is not supported at this size.

## D4. Seeds and the run

400 replicates per cell, 2,233 seconds on three workers, one pass
(`pilot-c/results/run_info.rds`, `pilot-c/run.log`). The
pre-registration and the amendment were committed and pushed before the
run.

## D5. What was not done

- No edge was computed on GSE235508.
- The marker list has no opened published source yet.
- Composition observed without error; two factors; normal noise.
- No dedicated cell for the spread term of D2.

# Action items and limitations

## Action items

1.  **Decide the status.** The rule gives `piloted-close-call` again;
    both readings are in the decisions log (D1). The method claims
    passed with 400 replicates; the artifact is small at this dataset’s
    measured levels, and edges do not detect activity-linked
    co-expression.
2.  **Decide how the pregnancy work enters the one paper.** One reading:
    the pregnancy application shows the fix holds in a longitudinal
    design, that the artifact follows the spread change and the
    imbalance of visits (decisions log, D2), and that a real SLE-only
    change remains detectable after the fix.
3.  **Decide whether aim 2 (activity) stays.** Neither across-woman nor
    within-woman prediction found activity-linked co-expression at 36
    women.
4.  **Resolve the yellow `[TO CONFIRM]` items** in the spec, the marker
    list’s source (proposal section 13), and verify the citations.

## Limitations

- Everything here is simulation; no edge was computed on GSE235508.
- The marker list has no opened published source; composition is
  observed without error in the simulation.
- Two factors with equal loadings, normal noise, five samples per woman.
- The spread-and-imbalance account of D2 is a sketch, checked only by
  comparing cells across pilots.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-devlin1975robust" class="csl-entry">

Devlin, Susan J., R. Gnanadesikan, and J. R. Kettenring. 1975. “Robust
Estimation and Outlier Detection with Correlation Coefficients.”
*Biometrika* 62 (3): 531-45. <https://doi.org/10.1093/biomet/62.3.531>.

</div>

<div id="ref-kuijjer2019lioness" class="csl-entry">

Kuijjer, M. L., M. G. Tung, G. Yuan, J. Quackenbush, and K. Glass. 2019.
“Estimating Sample-Specific Regulatory Networks.” *iScience* 14: 226-40.
<https://doi.org/10.1016/j.isci.2019.03.021>.

</div>

<div id="ref-kuijjer2019lionessr" class="csl-entry">

Kuijjer, Marieke L., Ping-Han Hsieh, John Quackenbush, and Kimberly
Glass. 2019. “<span class="nocase">lionessR</span>: Single Sample
Network Inference in R.” *BMC Cancer* 19 (1): 1003.
<https://doi.org/10.1186/s12885-019-6235-7>.

</div>

<div id="ref-kuijjer2026challenges" class="csl-entry">

Kuijjer, Marieke, Margherita De Marzio, and Kimberly Glass. 2026.
*Challenges and Opportunities in Single-Sample Network Modeling*.
<https://doi.org/10.64898/2026.02.27.708608>.

</div>

<div id="ref-lien2024pregnancy" class="csl-entry">

Lien, Hilde Johanne Tokle, Tina T. Pedersen, Bente Jakobsen, et al.
2024. “Single-Cell Resolution of Longitudinal Blood Transcriptome
Profiles in Rheumatoid Arthritis, Systemic Lupus Erythematosus and
Healthy Control Pregnancies.” *Annals of the Rheumatic Diseases* 83 (3):
300-311. <https://doi.org/10.1136/ard-2023-224644>.

</div>

<div id="ref-yu2017iena" class="csl-entry">

Yu, Xiangtian, Jingsong Zhang, Shaoyan Sun, Xin Zhou, Tao Zeng, and
Luonan Chen. 2017. “Individual-Specific Edge-Network Analysis for
Disease Prediction.” *Nucleic Acids Research* 45 (20): e170.
<https://doi.org/10.1093/nar/gkx787>.

</div>

</div>

[^1]: A pilot is a quick, small first run of the proposal’s simulation,
    done before the project goes further, to check that the idea is
    worth the time. It asks three questions: does the effect exist, is
    it large enough at realistic sample sizes, and is there an obvious
    fix or diagnostic? The proposal (`proposal.md`) holds the pilot plan
    and its pass criteria.

[^2]: `piloted-pass` means all three pilot questions got a clear yes.
    `piloted-fail` means at least one clearly did not, so the project
    should not go on as proposed. `piloted-close-call` means at least
    one answer could be read either way, so the call is stated both
    ways: the reading that passes and the reading that fails.

[^3]: LIONESS (Linear Interpolation to Obtain Network Estimates for
    Single Samples) builds a network from all samples and again without
    one sample, and scales the difference to give that sample’s network.

[^4]: Interferons are signalling proteins of the immune system; many
    genes respond to them together, and their summed expression, the
    interferon score, is high in many patients with lupus.

[^5]: A partial correlation is the correlation of two genes after the
    part of each that follows other variables (here the interferon score
    and the cell mix) has been removed by regression.

[^6]: The Gene Expression Omnibus (GEO) is the public archive of the US
    National Center for Biotechnology Information for gene expression
    data; a series (GSE) holds one study’s samples.

[^7]: Systemic lupus erythematosus is an autoimmune disease that can
    affect skin, joints, kidneys, and blood, and often flares in
    pregnancy.

[^8]: Rheumatoid arthritis is an autoimmune disease of the joints;
    seropositive patients carry rheumatoid factor or anti-citrullinated
    protein antibodies, seronegative patients do not. It often improves
    in pregnancy.

[^9]: DAS28 (Disease Activity Score in 28 joints) combines tender and
    swollen joint counts, a blood inflammation marker, and the patient’s
    global assessment; above 5.1 is high activity, below 2.6 remission.

[^10]: LAI-P is a lupus activity index adapted for pregnancy, which
    discounts features that pregnancy itself causes; 0 means no
    activity.

[^11]: The concordance index is the share of pairs of visits with
    different activity in which the visit with higher activity gets the
    higher prediction; 0.5 is chance and 1 is perfect, and for a
    yes-or-no outcome it equals the AUC.

[^12]: The area under the ROC curve is the probability that a randomly
    chosen pregnant sample gets a higher score than a randomly chosen
    non-pregnant one; 0.5 is chance and 1 is perfect separation.

[^13]: In cross-validation the data are split into folds, each model is
    fitted without one fold and scored on it, and the scores are pooled;
    here every sample of a woman is in the same fold, so a woman’s own
    samples never help predict each other.

[^14]: The elastic net is logistic regression with a penalty that
    shrinks coefficients and sets many to zero, so it can use more
    features than samples; glmnet is the R package that fits it.

[^15]: The AI agent that ran the pilot, wrote the files this document is
    built from, and kept a log of its decisions.
