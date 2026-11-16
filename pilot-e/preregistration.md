# Pre-registration: `pilot-e/`

| | |
|---|---|
| Proposal | `proposal.md`, section 12 (revision 1, network metrics), approved by the user ("pilot e", 2026-11-16) |
| Written | 2026-11-16 |
| Committed | in the commit "Stage 2, step 2 (pilot-e): pre-registration", before any pilot-e code or run |

## 1. The three questions

1. Do the closed forms hold? Is a LIONESS module's connectivity the square of its expression score, and do raw LIONESS and first-order edges give the same stability for every metric?
2. Are metrics more stable than single edges, and is their stability the expression score's? At the interferon stability measured in GSE65391, is the interferon module's connectivity more stable than its single edges, and as stable as its squared expression score?
3. Does the correction separate real structure? Do partial-correlation edges remove the interferon module's excess stability and keep a real personal module's?

## 2. The pass rule

For each sample, with edges over genes 6 to 50: the connectivity of a gene set is its mean edge (IFN set: genes 6 to 20; module: 21 to 30; cell mix: 31 to 50; cross: one gene of 6 to 20 and one of 21 to 30). The IFN score is the mean of genes 6 to 20 standardized over the samples. $I(m)$ is the intraclass correlation of metric $m$ over visits (`EdgeICC()`); for a node metric, the mean over the genes of the set of the per-gene intraclass correlations. Means over replicates, MCSE as before.

| Q | Clause | Quantity | Threshold | Cells | Within its MCSE |
|---|---|---|---|---|---|
| 1 | (a) | correlation over samples of raw IFN connectivity with the squared IFN score | $\ge 0.9$ | 1 to 4, each | $\lvert v - 0.9\rvert \le 3$ MCSE: close call |
| 1 | (b) | $\lvert I_{\text{raw}}(m) - I_{\text{F1}}(m)\rvert$ for every metric of section 3 | $\le \max(0.02, 3\,\text{MCSE})$ | 1 to 4 | $\max(0.02, 3\,\text{MCSE}) < v \le 0.02 + 3\,\text{MCSE}$: close call |
| 2 | (a) | $I_{\text{raw}}(\text{IFN connectivity}) - \bar I_{\text{raw}}(\text{IFN edges})$ | $\ge 0.05$ | 2 | $\lvert v - 0.05\rvert \le 3$ MCSE: close call |
| 2 | (b) | $I_{\text{raw}}(\text{IFN connectivity}) - I(\text{IFN score}^2)$ | $\lvert v\rvert \le \max(0.05, 3\,\text{MCSE})$ | 2 | $\max(0.05, 3\,\text{MCSE}) < \lvert v\rvert \le 0.05 + 3\,\text{MCSE}$: close call |
| 3 | (a) | $I_{\text{partial}}(\text{IFN connectivity}) - I_{\text{partial}}(\text{cross connectivity})$ | $\lvert v\rvert \le \max(0.05, 3\,\text{MCSE})$ | 2, 3 | as 2(b) |
| 3 | (b) | $I_{\text{partial}}(\text{module connectivity}) - I_{\text{partial}}(\text{cross connectivity})$ | $\ge 0.05$ | 4 | as 2(a) |

**How the answers combine.** Pass when every clause passes. Fail when 1(a), 1(b), or 3(a) fails. Close call otherwise.

## 3. Design

- **Mechanism.** `SimulateRepeatedCohort()` as in pilot-d: 131 people, six visits, cell-mix stability 0.45, activity unrelated ($\gamma = 0$).
- **Cells (4).** 1: IFN stability 0, no personal module; 2: 0.55, none; 3: 0.9, none; 4: 0.55, personal module spread 0.25.
- **Edge types.** Raw LIONESS, first-order, partial-correlation given the IFN score of genes 1 to 5 and the observed cell mix, on all samples.
- **Metrics, per sample and edge type** (`NetworkMetrics()`): connectivity of the IFN, module, cell-mix, and cross sets; node strength (signed sum of a gene's edges); weighted betweenness and average shortest path length with distance $1 / \lvert e\rvert$; global efficiency (mean inverse distance); eigenvector centrality with weight $\lvert e\rvert$; on the graph of edges with $\lvert e\rvert > 1$, density and local clustering (0 for a gene with fewer than two neighbours). Node metrics are summarized by set (IFN, module, cell-mix genes). Reported, not judged: hub overlap (Jaccard of the five strongest genes between visits 1 and 2, within a person minus between people), and the stability of the IFN score itself.
- **Not tested here.** Prediction of the next visit's activity (stage 3); path metrics have no closed form, so they are judged only in clause 1(b).
- **Replicates.** From a timing run of ten replicates per cell with the workers loaded: 200 if the projected run fits one hour on three cores, otherwise the largest multiple of 20 that fits, not below 60.
- **Seeds.** Cell $c$ has seed $21161116 + 1000c$; replicate $r$ the cell seed plus $r$.
- **Known answers.** Cell 1 for the stability excesses; tests: metric functions against hand-computed values on small graphs (path length, efficiency, density, clustering), and module connectivity against its closed form on a small example.
- **Time box.** One hour on three cores.

## 4. What is reported whatever the outcome

Every metric's intraclass correlation by cell and edge type, the closed-form correlations, hub overlap, the verdict table, the known-answer table, and one figure of stability by metric, in `pilot-e/results/`.

## 5. Readings of open cases

- 2026-11-16, before the run. Thresholds were set from the sketch in proposal section 4 with the noise share computed (pilot-d's lesson): the IFN score's square has an intraclass correlation of about 0.27 at the measured level, against about 0.14 for single IFN edges.

## 6. Amendments before the run

- 2026-11-16, after the timing run and before the run. Replicates per cell: 120. The timing run (10 replicates per cell, three workers loaded, the run's own scheduling; `pilot-e/results/timing.rds`) took 284 seconds, which projects 200 replicates to about 5,700 seconds, over the one-hour box; 120 replicates project to about 3,400 seconds, the largest multiple of 20 that fits.
