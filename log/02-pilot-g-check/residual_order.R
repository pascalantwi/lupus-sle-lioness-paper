# Scratch check after the pilot-g run (not pre-registered): is the stable
# residual of IFN connectivity in cell 3 a higher-order expression term?
# Five replicates of cell 3: the residual's share of variance, and its
# stability before and after adding fourth-order terms of the pathway's
# standardized expression (S^4, S^2 Q, Q^2, and the sum of z^4).
devtools::load_all(here::here(), quiet = TRUE)
row <- PilotGDesign()[3, ]
sets <- PathwaySetsG()
blocks <- PathwayBlocks(sets, 90)
out <- do.call(rbind, lapply(1:5, function(rep) {
  set.seed(row$seed + rep)
  coh <- SimulatePathwayCohort(row$people, row$visits, row$kappa_ifn,
                               row$kappa_neut, row$kappa_lymph, row$kappa_pb,
                               row$tau)
  met <- PathwayMetrics(LionessEdges(coh$x), blocks)[, "within_ifn"]
  tr <- PathwayTerms(coh$x, sets)$within_ifn
  z <- StandardizeN(coh$x)$z[, sets$ifn]
  f1 <- R2OnTerms(met, tr)
  hi <- cbind(tr, tr[, "s2"]^2, tr[, "s2"] * tr[, "q"], tr[, "q"]^2,
              rowSums(z^4))
  f2 <- R2OnTerms(met, hi)
  data.frame(rep = rep, one_minus_r2 = 1 - f1$r2,
             icc_resid = EdgeICC(f1$resid, coh$person),
             one_minus_r2_4th = 1 - f2$r2,
             icc_resid_4th = EdgeICC(f2$resid, coh$person))
}))
print(out, digits = 3)
print(colMeans(out[, -1]), digits = 3)
