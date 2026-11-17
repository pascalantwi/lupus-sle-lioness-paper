# Scratch check after the pilot-g run (not pre-registered): is the stable
# residual specific to the interferon pathway, or to a very stable factor?
# Cell 2 as pre-registered, but the neutrophil-like factor made very stable
# (0.9) while interferon stays at the measured 0.55. Ten replicates.
devtools::load_all(here::here(), quiet = TRUE)
row <- PilotGDesign()[2, ]
row$kappa_neut <- 0.9
sets <- PathwaySetsG()
blocks <- PathwayBlocks(sets, 90)
out <- do.call(rbind, lapply(1:10, function(rep) {
  set.seed(row$seed + 500 + rep)
  coh <- SimulatePathwayCohort(row$people, row$visits, row$kappa_ifn,
                               row$kappa_neut, row$kappa_lymph, row$kappa_pb,
                               row$tau)
  met <- PathwayMetrics(LionessEdges(coh$x), blocks)
  tr <- PathwayTerms(coh$x, sets)
  sapply(c(ifn = "within_ifn", neut = "within_neut", lymph = "within_lymph"),
         function(b) EdgeICC(R2OnTerms(met[, b], tr[[b]])$resid, coh$person))
}))
cat("Residual stability (ICC) after the first-order terms, mean of 10:\n")
print(round(colMeans(out), 3))
cat("MCSE:\n")
print(round(apply(out, 2, sd) / sqrt(nrow(out)), 4))
