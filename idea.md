# Idea: LIONESS edges are quadratic functions of expression: closed forms, an interferon artifact, and equivalence with first-order features

## The idea in a paragraph

LIONESS promises a gene network for each person. If that network is a property of the person, it should stay much the same when her blood is drawn again, it should identify her among other people, and it may track how her disease changes. A LIONESS edge is, to first order, a quadratic function of the person's own standardized expression, so a factor that moves many genes, such as the interferon response or the blood-cell mix, enters every edge as its square. A person's network can then look stable and personal only because her interferon level is stable, and look noisy because her cell mix changes from visit to visit. This project asks how much of a personal network's stability, its ability to identify the person, and its ability to predict her next disease activity is real gene coordination, and how much is interferon and cell mix, in children with lupus followed over repeated visits (GSE65391).

## Question

Are LIONESS personal networks personal? Over repeated blood samples of the same patient, how stable are her edges, can they identify her, and do they predict her next disease activity, beyond her interferon level, her measured blood-cell counts, her gene levels, and the edges' own first-order formula?

## Mechanism

The edge $e_q = N r - (N-1) r_{(-q)}$ equals $r + z_{qi} z_{qj} - \tfrac{r}{2}(z_{qi}^2 + z_{qj}^2)$ to first order. Under a shared factor $F$ with loadings $\lambda$, its mean given $F_q$ is linear in $F_q^2$, so the edge inherits the stability of $F^2$ over visits, not only of co-expression.

## Aims

1. Theory: the edge's stability under a stable shared factor, in closed form.
2. Simulation matching the data: patients with repeated visits, stable personal factors, changing cell mix, and a real personal network or none.
3. Stability of personal networks over time, raw against corrected.
4. Fingerprinting: identifying a patient from her network at another visit.
5. Prediction of the next visit's disease activity.
6. Equivalence of edges and their first-order formula throughout.

## Data

GSE65391: 158 children with lupus, 924 visits, SLEDAI and measured blood counts at every visit, kidney class, treatments; 46 healthy children. ILLUMINATE (GSE88884 to GSE88886) as an optional adult check.

## Your own prior work

The first LIONESS project (`work-flow-for-lioness-paper`): the exact and first-order forms of the edge, the squared-factor artifact and its fix, the unequal-groups law, three pilots, and a draft whose title this project keeps. In this repository, three pilots on pregnancy in SLE and RA (GSE235508), all close calls, which showed that within-disease clinical contrasts move interferon and cell mix little.

## Notes

The user's re-scope of 2026-11-16, after a review of candidate ideas: the method and the title stay; the SLE against RA application is not kept; the simulation must mirror the application.
