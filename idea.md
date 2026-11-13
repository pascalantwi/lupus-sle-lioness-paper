# Idea: Does LIONESS find network rewiring over pregnancy in SLE and RA, or only interferon and cell composition?

## The idea in a paragraph

Pregnancy changes the immune system, and the two closest rheumatic diseases respond to it differently: rheumatoid arthritis (RA) often improves, while systemic lupus erythematosus (SLE) is prone to flare. Single-sample network methods such as LIONESS are used to claim that gene networks "rewire" over time. A LIONESS edge is, to first order, a quadratic function of the patient's own standardized expression, so any factor that moves many genes and changes over time (cell composition in pregnancy, and interferon in SLE) enters the edges as its square. This project asks whether LIONESS edges track pregnancy beyond the interferon score, cell composition, their squares, the gene levels, and their own first-order formula, in SLE and in RA, using one longitudinal study of whole blood through pregnancy (GSE235508), each disease analysed on its own and no studies merged.

## Question

Do LIONESS edges distinguish late pregnancy from the postpartum period, within the same women, beyond the IFN score, cell composition, their squares, the gene levels, and their own first-order features, and does the answer differ between SLE and RA?

## Mechanism

The LIONESS edge $e_q = N r - (N - 1) r_{(-q)}$ equals, to first order, $r + z_{qi} z_{qj} - \tfrac{r}{2}(z_{qi}^2 + z_{qj}^2)$. A factor $P$ that moves genes $i$ and $j$ gives $z_{qi} z_{qj} \approx \lambda_i \lambda_j P_q^2$, so a factor that changes over pregnancy changes the edges through its square, with no change in co-expression.

## Aims

1. Equivalence over time: do LIONESS edges and their first-order formula classify pregnancy stage alike?
2. Artifact and fix: how much do edges gain from a shared factor that changes over pregnancy, and does a base holding its square remove the gain?
3. Real rewiring: when co-expression truly changes over pregnancy, do edges detect it beyond the gene levels, and is the formula equally good?

## Data

GSE235508 (Lien et al. 2024): 335 whole-blood RNA-seq samples from 84 pregnancies (SLE, RA, healthy), before pregnancy, at each trimester, and 6 weeks, 6 months, and 12 months postpartum; one study.

## Your own prior work

The first LIONESS project (`pascalantwi/work-flow-for-lioness-paper`): the exact and first-order forms of the LIONESS edge, the squared-factor artifact and its fix, three pilots, and a first draft on SLE against RA.

## Notes

Supervised by Dr Tyler Grimes. Chosen by the user on 2026-11-13 as "idea 1" of the options for a more interesting application with the same method; the user asked for the pilot to be run in this repository.
