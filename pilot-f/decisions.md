# Decisions log: pilot-f (the own-module correction)

## D0. The literature pass

As pilot-e. Nothing read since.

## D1. The status call: `piloted-fail`

Clause 3(a) failed in both cells and clause 3(b) failed (`pilot-f/results/verdict.rds`): with the own-module correction, the interferon set's connectivity stayed exactly as stable as with raw edges, and the personal module's connectivity exceeded the cross connectivity by less than 0.05. Clauses 1(a), 1(b), 2(a), and 2(b) passed again, with pilot-e's values: a module's connectivity correlated about 0.97 with its squared expression score; raw and first-order edges gave the same stability for all 19 metrics; connectivity was steadier than single edges and as steady as the squared score. There is no reading under which the rule passes.

## D2. Why: a set's connectivity cannot be corrected from within the set

The own-module correction regresses each gene on the mean of the other genes of its set. The mean of those residuals over the set is then a fixed multiple of the set's original mean standardized expression (each gene's standardized value enters $m - 1$ of the leave-one-out means, so their average is the set mean), and LIONESS standardizes the residuals again. Connectivity is, to first order, the squared set mean of standardized expression plus small terms (pilot-e report, "In symbols"), so it is unchanged up to scale, and an intraclass correlation does not depend on scale. The interferon set's connectivity stability was about 0.24 raw and 0.25 corrected at the measured level, and 0.62 for both when interferon was very stable (`cell_summaries.rds`). The correction did steady the metrics that do not reduce to one set's mean (density, path length, strength, and the cell-mix set, whose correction includes the measured cell mix).

This is a result, not a defect of the code: within-set adjustment cannot separate a set's connectivity from its expression score, and an external score (pilot-e) leaves the score's measurement error. Together, pilots e and f say that a gene set's LIONESS connectivity is its squared expression score, and that no linear adjustment recovers a separate "network" part from it. The derivation above is the agent's sketch; its numerical check is the equality of raw and corrected stabilities across the four cells.

## D3. The personal module

A real personal module raised module connectivity's stability from about 0.03 to 0.095 for raw and corrected edges alike; the excess over cross pairs was about 0.036, below 0.05 (3(b) fail). The module's connectivity is again its squared module score, whose stability rises when the module's correlation is personal; the metric sees the personal module only through the module's expression.

## D4. The run

120 replicates per cell (amendment before the run), 3,136 seconds on three workers, one pass.

## D5. What was not done

As pilot-e; and no correction from outside both the set and the score (for example a held-out score from other samples of the same person), which would need its own design.
