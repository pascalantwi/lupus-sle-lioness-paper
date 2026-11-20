# Figure functions. Each draws exactly one plot on the current device (the
# project's rule: no figure combines panels), so the manuscript, the
# background report, the slides, and the poster draw the same figures.

#' Colors of the six pathways
#'
#' The first six slots of a categorical palette checked for colour-vision
#' deficiency; three fall below 3:1 contrast on white, so every figure labels
#' pathways directly as well.
#'
#' @return A named character vector of hex colors, names as [PathwaySetsG()].
#' @export
PathwayColors <- function() {
  c(ifn = "#2a78d6", neut = "#eb6834", lymph = "#1baf7a", pb = "#eda100",
    personal = "#e87ba4", null = "#008300")
}

#' Labels of the six pathways
#'
#' @return A named character vector, names as [PathwaySetsG()].
#' @export
PathwayLabels <- function() {
  c(ifn = "Interferon", neut = "Neutrophil-like", lymph = "Lymphocyte-like",
    pb = "Plasmablast-like", personal = "Personal", null = "Null")
}

#' Node positions for a gene network grouped by pathway
#'
#' Each pathway's genes sit on a small circle, and the pathways sit on a large
#' circle, so the within-pathway edges are short and the crosstalk edges cross
#' the middle.
#'
#' @param sets A named list of disjoint gene index sets.
#' @param inner Radius of each pathway's circle relative to the large one.
#' @return A matrix, one row per gene (row `i` is gene `i`), columns `x`, `y`;
#'   attribute `centers` holds the pathways' centres.
#' @export
PathwayLayout <- function(sets, inner = 0.28) {
  k <- length(sets)
  p <- max(unlist(sets))
  ang <- pi / 2 - 2 * pi * (seq_len(k) - 1) / k
  centers <- cbind(x = cos(ang), y = sin(ang))
  rownames(centers) <- names(sets)
  xy <- matrix(NA_real_, p, 2, dimnames = list(NULL, c("x", "y")))
  for (j in seq_len(k)) {
    g <- sets[[j]]
    a <- 2 * pi * (seq_along(g) - 1) / length(g)
    xy[g, ] <- cbind(centers[j, 1] + inner * cos(a),
                     centers[j, 2] + inner * sin(a))
  }
  attr(xy, "centers") <- centers
  xy
}

#' Write labels just outside points on a circle
#'
#' Right-hand labels start at their point, left-hand ones end there, and the
#' top and bottom ones are centred, so no label covers its node.
#'
#' @param xy Points, one row each, around the origin.
#' @param labels Text.
#' @param gap Distance from the point, outward.
#' @return `NULL`, invisibly.
#' @keywords internal
LabelOutside <- function(xy, labels, gap) {
  r <- sqrt(rowSums(xy^2))
  ux <- xy[, 1] / r
  uy <- xy[, 2] / r
  hadj <- ifelse(ux > 0.3, 0, ifelse(ux < -0.3, 1, 0.5))
  for (i in seq_along(labels)) {
    graphics::text(xy[i, 1] + gap * ux[i], xy[i, 2] + gap * uy[i],
                   labels[i], adj = c(hadj[i], 0.5), cex = 0.8,
                   col = "grey20", font = 2)
  }
  invisible(NULL)
}

#' Edge color by sign
#'
#' @param e Edge weights.
#' @param alpha Opacity.
#' @return Colors: red for positive, blue for negative edges.
#' @keywords internal
EdgeColors <- function(e, alpha = 0.55) {
  base <- ifelse(e > 0, "#e34948", "#1c5cab")
  alpha <- rep_len(alpha, length(e))
  # adjustcolor() takes one opacity per call
  vapply(seq_along(e), function(i) grDevices::adjustcolor(base[i], alpha[i]),
         character(1))
}

#' Draw one sample's LIONESS gene network
#'
#' The strongest edges of one sample's LIONESS network, with genes grouped by
#' pathway ([PathwayLayout()]). Edge width is proportional to the absolute
#' edge, on a fixed scale so that networks of different visits compare; red
#' edges are positive, blue negative.
#'
#' @param e One sample's edges over all gene pairs ([PairIndex()] order).
#' @param sets A named list of gene index sets.
#' @param top Number of strongest edges drawn (absolute value).
#' @param scale_max Absolute edge drawn at the widest line; defaults to the
#'   largest drawn edge.
#' @param main Plot title.
#' @param node_cex Size of the gene points.
#' @param labels Pathway labels, named as `sets`.
#' @param legend Draw the legend of edge sign and width.
#' @return The drawn edges' indices, invisibly.
#' @export
PlotSampleNetwork <- function(e, sets, top = 120L, scale_max = NULL,
                              main = "", node_cex = 1.3,
                              labels = PathwayLabels(), legend = TRUE) {
  p <- max(unlist(sets))
  pr <- PairIndex(p)
  keep <- order(-abs(e))[seq_len(min(top, length(e)))]
  keep <- keep[order(abs(e[keep]))]  # the strongest drawn last, on top
  if (is.null(scale_max)) scale_max <- max(abs(e[keep]))
  xy <- PathwayLayout(sets)
  memb <- rep(NA_character_, p)
  for (k in names(sets)) memb[sets[[k]]] <- k
  graphics::par(mar = c(0, 0, if (nzchar(main)) 2 else 0, 0))
  graphics::plot(NA, xlim = c(-2.35, 2.35), ylim = c(-2.1, 1.6), asp = 1,
                 xaxs = "i", yaxs = "i",
                 axes = FALSE, xlab = "", ylab = "", main = main)
  w <- pmin(abs(e[keep]) / scale_max, 1)
  graphics::segments(xy[pr[keep, 1], 1], xy[pr[keep, 1], 2],
                     xy[pr[keep, 2], 1], xy[pr[keep, 2], 2],
                     col = EdgeColors(e[keep], 0.25 + 0.5 * w),
                     lwd = 0.3 + 3.2 * w)
  cols <- PathwayColors()[memb]
  cols[is.na(cols)] <- "grey50"
  graphics::points(xy[, 1], xy[, 2], pch = 21, bg = cols, col = "white",
                   cex = node_cex, lwd = if (node_cex < 1) 0.3 else 1)
  cen <- attr(xy, "centers")
  lab <- labels[names(sets)]
  lab[is.na(lab)] <- names(sets)[is.na(lab)]
  LabelOutside(cen, lab, gap = 0.4)
  if (legend) {
    graphics::legend("bottom", inset = 0.01, ncol = 2, cex = 0.75,
                     legend = c("Positive edge", "Negative edge",
                                "Weak edge (thin)", "Strong edge (thick)"),
                     col = c("#e34948", "#1c5cab", "grey45", "grey45"),
                     lwd = c(2.5, 2.5, 0.6, 3.5), bty = "n")
  }
  invisible(keep)
}

#' Draw one sample's pathway network
#'
#' A LIONESS network summarized by pathway: one node per pathway, its area
#' proportional to the pathway's connectivity (the mean edge within it), and
#' one line per pair of pathways, its width proportional to the absolute
#' crosstalk (the mean edge between them); red is positive, blue negative.
#' Nodes with negative connectivity are drawn hollow.
#'
#' @param conn Named connectivities, names `within_<set>` or `conn_<set>`.
#' @param cross Named crosstalks, names `cross_<A>_<B>`.
#' @param scale_conn,scale_cross Values drawn at the largest node and the
#'   widest line; default to the largest absolute values given.
#' @param main Plot title.
#' @param labels Pathway labels, named by pathway.
#' @param legend Draw the legend of circle size and fill and line colour
#'   and width.
#' @return `NULL`, invisibly.
#' @export
PlotPathwayNetwork <- function(conn, cross, scale_conn = NULL,
                               scale_cross = NULL, main = "",
                               labels = PathwayLabels(), legend = TRUE) {
  nm <- sub("^(within|conn)_", "", names(conn))
  k <- length(nm)
  ang <- pi / 2 - 2 * pi * (seq_len(k) - 1) / k
  xy <- cbind(cos(ang), sin(ang))
  rownames(xy) <- nm
  if (is.null(scale_conn)) scale_conn <- max(abs(conn))
  if (is.null(scale_cross)) scale_cross <- max(abs(cross))
  graphics::par(mar = c(0, 0, if (nzchar(main)) 2 else 0, 0))
  graphics::plot(NA, xlim = c(-2.3, 2.3), ylim = c(-2.3, 1.62), asp = 1,
                 xaxs = "i", yaxs = "i",
                 axes = FALSE, xlab = "", ylab = "", main = main)
  pairs <- strsplit(sub("^cross_", "", names(cross)), "_")
  for (i in seq_along(cross)) {
    ab <- pairs[[i]]
    graphics::segments(xy[ab[1], 1], xy[ab[1], 2], xy[ab[2], 1], xy[ab[2], 2],
                       col = EdgeColors(cross[i], 0.7),
                       lwd = 0.5 + 9 * pmin(abs(cross[i]) / scale_cross, 1))
  }
  cols <- PathwayColors()[nm]
  cols[is.na(cols)] <- "grey50"
  size <- 1.5 + 6 * sqrt(pmin(abs(conn) / scale_conn, 1))
  graphics::points(xy[, 1], xy[, 2], pch = 21, cex = size,
                   bg = ifelse(conn >= 0, cols, "white"), col = cols, lwd = 2)
  lab <- labels[nm]
  lab[is.na(lab)] <- nm[is.na(lab)]
  LabelOutside(xy, lab, gap = 0.42)
  if (legend) {
    graphics::legend("bottom", inset = 0.01, ncol = 2, cex = 0.75, bty = "n",
                     legend = c("High connectivity (large circle)",
                                "Low connectivity (small circle)",
                                "Negative connectivity (hollow)",
                                "Positive crosstalk", "Negative crosstalk",
                                "Thicker line: stronger crosstalk"),
                     pch = c(21, 21, 21, NA, NA, NA),
                     pt.cex = c(2.6, 1.2, 1.8, NA, NA, NA),
                     pt.bg = c("grey45", "grey45", "white", NA, NA, NA),
                     col = c("grey45", "grey45", "grey45", "#e34948",
                             "#1c5cab", "grey45"),
                     lwd = c(NA, NA, NA, 3, 3, 6))
  }
  invisible(NULL)
}

#' Draw connectivity against the squared pathway score
#'
#' @param conn Connectivity, one value per sample.
#' @param score2 Squared pathway score, one value per sample.
#' @param col Point color.
#' @param xlab,ylab,main Labels.
#' @return The correlation, invisibly.
#' @export
PlotConnectivityScore <- function(conn, score2, col = "#2a78d6",
                                  xlab = "Squared pathway score",
                                  ylab = "Pathway connectivity", main = "") {
  graphics::par(mar = c(4.2, 4.4, if (nzchar(main)) 2.5 else 1, 1))
  graphics::plot(score2, conn, pch = 19, cex = 0.55,
                 col = grDevices::adjustcolor(col, 0.5), xlab = xlab,
                 ylab = ylab, main = main, las = 1, bty = "l")
  fit <- stats::lm(conn ~ score2)
  graphics::abline(fit, col = "grey30", lwd = 1.5)
  graphics::legend("topleft", c("One sample (one child at one visit)",
                                "Least-squares line"),
                   pch = c(19, NA), lwd = c(NA, 1.5),
                   col = c(grDevices::adjustcolor(col, 0.7), "grey30"),
                   cex = 0.75, bg = "white", box.col = "grey80")
  invisible(stats::cor(conn, score2))
}

#' Draw a quantity over each child's visits
#'
#' @param value One value per sample.
#' @param person,visit Person and visit of each sample.
#' @param ids People drawn; defaults to all.
#' @param col Line color.
#' @param ylab,main Labels.
#' @return `NULL`, invisibly.
#' @export
PlotTrajectories <- function(value, person, visit, ids = unique(person),
                             col = "#2a78d6", ylab = "", main = "") {
  use <- person %in% ids
  graphics::par(mar = c(4.2, 4.4, if (nzchar(main)) 2.5 else 1, 1))
  yr <- range(value[use])
  # headroom above the data for the legend
  graphics::plot(NA, xlim = range(visit[use]),
                 ylim = c(yr[1], yr[2] + 0.22 * diff(yr)),
                 xlab = "Visit", ylab = ylab, main = main, las = 1, bty = "l",
                 xaxt = "n")
  graphics::axis(1, at = sort(unique(visit[use])))
  for (id in ids) {
    w <- which(person == id)
    w <- w[order(visit[w])]
    graphics::lines(visit[w], value[w], col = grDevices::adjustcolor(col, 0.6),
                    lwd = 1.5)
    graphics::points(visit[w], value[w], pch = 19, cex = 0.6, col = col)
  }
  graphics::legend("topright", c("One child: the value at each visit",
                                 "joined in visit order"),
                   pch = c(19, NA), lwd = c(NA, 1.5),
                   col = c(col, grDevices::adjustcolor(col, 0.6)),
                   cex = 0.75, bg = "white", box.col = "grey80")
  invisible(NULL)
}

#' Draw stability by pathway for one simulated setting
#'
#' One row per pathway; one point per quantity (connectivity, squared score,
#' single edges, residual), the stability over visits with a bar of two Monte
#' Carlo standard errors.
#'
#' @param summaries Output of [SummarizeCellsG()].
#' @param cell The setting drawn.
#' @param kinds Named quantities to draw (prefix of `icc_<kind>_<pathway>`)
#'   with their legend labels.
#' @param main Plot title.
#' @return The plotted means, invisibly.
#' @export
PlotStabilityByPathway <- function(summaries, cell,
                                   kinds = c(conn = "Connectivity",
                                             score2 = "Squared score",
                                             edges = "Single edges (mean)",
                                             resid = "Residual on the terms"),
                                   main = "") {
  pw <- names(PathwayLabels())
  cols <- c("#e34948", "grey35", "#7d4fc4", "#2a78d6")[seq_along(kinds)]
  pch <- c(19, 17, 15, 18)[seq_along(kinds)]
  get <- function(q, f) {
    v <- summaries[summaries$cell == cell & summaries$quantity == q, f]
    if (length(v)) v else NA_real_
  }
  m <- sapply(names(kinds), function(k) {
    sapply(pw, function(p) get(sprintf("icc_%s_%s", k, p), "mean"))
  })
  se <- sapply(names(kinds), function(k) {
    sapply(pw, function(p) get(sprintf("icc_%s_%s", k, p), "mcse"))
  })
  yy <- rev(seq_along(pw))
  off <- seq(0.24, -0.24, length.out = length(kinds))
  graphics::par(mar = c(4.2, 8.5, if (nzchar(main)) 2.5 else 1, 1))
  graphics::plot(NA, xlim = range(c(0, m), na.rm = TRUE) + c(-0.02, 0.02),
                 ylim = c(0.5, length(pw) + 0.5), yaxt = "n", ylab = "",
                 xlab = "Stability over visits (intraclass correlation)",
                 main = main, bty = "l", las = 1)
  graphics::axis(2, at = yy, labels = PathwayLabels(), las = 1, cex.axis = 0.8)
  graphics::abline(v = 0, col = "grey80")
  for (k in seq_along(kinds)) {
    graphics::segments(m[, k] - 2 * se[, k], yy + off[k],
                       m[, k] + 2 * se[, k], yy + off[k], col = cols[k])
    graphics::points(m[, k], yy + off[k], pch = pch[k], col = cols[k],
                     cex = 0.9)
  }
  graphics::legend("bottomright", kinds, col = cols, pch = pch, cex = 0.7,
                   bg = "white", box.col = "grey80")
  invisible(m)
}

#' Labels of the LIONESS metrics
#'
#' @param m Metric names as in [LionessMetrics()].
#' @return Readable labels.
#' @export
MetricLabels <- function(m) {
  pw <- PathwayLabels()
  lab <- m
  base <- c(total_weight = "Total edge weight", n_edges = "Number of edges",
            n_nodes = "Number of nodes", transitivity = "Transitivity",
            n_nodes_p = "Number of nodes (top 5%)",
            transitivity_p = "Transitivity (top 5%)")
  hit <- m %in% names(base)
  lab[hit] <- base[m[hit]]
  kinds <- c(conn = "connectivity", strength = "strength",
             degree_p = "degree (top 5%)", degree = "degree",
             betweenness_p = "betweenness (top 5%)",
             betweenness = "betweenness")
  for (k in names(kinds)) {
    for (p in names(pw)) {
      lab[m == paste0(k, "_", p)] <- paste(pw[[p]], kinds[[k]])
    }
  }
  cr <- grepl("^cross_", m)
  lab[cr] <- vapply(strsplit(sub("^cross_", "", m[cr]), "_"), function(ab) {
    paste(pw[[ab[1]]], "-", pw[[ab[2]]], "crosstalk")
  }, character(1))
  lab
}

#' Draw the stability of LIONESS metrics against their first-order versions
#'
#' One row per metric; a filled point for the LIONESS network and an open
#' point for the first-order network, each with a bar of two Monte Carlo
#' standard errors.
#'
#' @param summaries Output of [SummarizeStudy()].
#' @param cell The setting drawn.
#' @param metrics Metric names, as in [LionessMetrics()].
#' @param main Plot title.
#' @return The plotted means, invisibly.
#' @export
PlotMetricStability <- function(summaries, cell, metrics, main = "") {
  get <- function(q, f) {
    v <- summaries[summaries$cell == cell & summaries$quantity == q, f]
    if (length(v)) v else NA_real_
  }
  raw <- vapply(metrics, function(m) get(paste0("icc_raw_", m), "mean"), 0)
  fo <- vapply(metrics, function(m) get(paste0("icc_fo_", m), "mean"), 0)
  se <- vapply(metrics, function(m) get(paste0("icc_raw_", m), "mcse"), 0)
  yy <- rev(seq_along(metrics))
  graphics::par(mar = c(4.2, 15, if (nzchar(main)) 2.5 else 1, 1))
  graphics::plot(NA, xlim = range(c(0, raw, fo), na.rm = TRUE) + c(-0.02, 0.02),
                 ylim = c(0.5, length(metrics) + 0.5), yaxt = "n", ylab = "",
                 xlab = "Stability over visits (intraclass correlation)",
                 main = main, bty = "l", las = 1)
  graphics::axis(2, at = yy, labels = MetricLabels(metrics), las = 1,
                 cex.axis = 0.7)
  graphics::abline(v = 0, col = "grey80")
  graphics::segments(raw - 2 * se, yy, raw + 2 * se, yy, col = "#e34948")
  graphics::points(raw, yy, pch = 19, col = "#e34948", cex = 0.9)
  graphics::points(fo, yy, pch = 1, col = "grey20", cex = 1.3, lwd = 1.2)
  graphics::legend("bottomright", c("LIONESS network", "First-order network",
                                    "Two Monte Carlo standard errors"),
                   pch = c(19, 1, NA), lwd = c(NA, NA, 1.5),
                   col = c("#e34948", "grey20", "#e34948"), cex = 0.7,
                   bg = "white", box.col = "grey80")
  invisible(cbind(raw = raw, fo = fo))
}

#' Draw quantities along one factor of the design
#'
#' One line per quantity, with bars of two Monte Carlo standard errors.
#'
#' @param summaries Output of [SummarizeStudy()].
#' @param design Output of [StudyDesign()].
#' @param factor Design column varied (for example `kappa_ifn`).
#' @param cells The settings on the gradient.
#' @param quantities Named quantities: names are summary quantities, values
#'   legend labels.
#' @param xlab,ylab,main Labels.
#' @param ref Optional reference line `y = f(x)`, a function.
#' @param ref_label Legend label of the reference line.
#' @return The plotted means, invisibly.
#' @export
PlotGradient <- function(summaries, design, factor, cells, quantities,
                         xlab = factor, ylab = "", main = "", ref = NULL,
                         ref_label = "Reference line") {
  x <- design[[factor]][match(cells, design$cell)]
  o <- order(x)
  x <- x[o]
  cells <- cells[o]
  get <- function(q, f) {
    vapply(cells, function(c) {
      v <- summaries[summaries$cell == c & summaries$quantity == q, f]
      if (length(v)) v else NA_real_
    }, 0)
  }
  m <- sapply(names(quantities), get, f = "mean")
  se <- sapply(names(quantities), get, f = "mcse")
  m <- matrix(m, length(cells))
  se <- matrix(se, length(cells))
  cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#7d4fc4", "#e34948",
            "grey35")[seq_along(quantities)]
  pch <- c(19, 17, 15, 18, 1, 2)[seq_along(quantities)]
  graphics::par(mar = c(4.2, 4.6, if (nzchar(main)) 2.5 else 1, 1))
  yl <- range(c(m - 2 * se, m + 2 * se, 0), na.rm = TRUE)
  # headroom above the data for the legend
  yl[2] <- yl[2] + 0.3 * diff(yl)
  # series that coincide are drawn side by side, so that none is hidden
  dodge <- (seq_along(quantities) - (length(quantities) + 1) / 2) *
    0.012 * diff(range(x))
  graphics::plot(NA, xlim = range(x) + c(-1, 1) * 2 * max(abs(dodge)),
                 ylim = yl, xlab = xlab, ylab = ylab,
                 main = main, bty = "l", las = 1, xaxt = "n")
  graphics::axis(1, at = x)
  if (!is.null(ref)) {
    xs <- seq(min(x), max(x), length.out = 100)
    graphics::lines(xs, ref(xs), col = "grey60", lty = 2)
  }
  graphics::abline(h = 0, col = "grey85")
  for (k in seq_along(quantities)) {
    xk <- x + dodge[k]
    graphics::lines(xk, m[, k], col = cols[k], lwd = 2)
    graphics::segments(xk, m[, k] - 2 * se[, k], xk, m[, k] + 2 * se[, k],
                       col = cols[k])
    graphics::points(xk, m[, k], pch = pch[k], col = cols[k])
  }
  leg <- c(unname(quantities), "Bars: two Monte Carlo standard errors")
  leg_col <- c(cols, "grey40")
  leg_pch <- c(pch, NA)
  leg_lty <- c(rep(1, length(quantities)), 1)
  leg_lwd <- c(rep(2, length(quantities)), 1)
  if (!is.null(ref)) {
    leg <- c(leg, ref_label)
    leg_col <- c(leg_col, "grey60")
    leg_pch <- c(leg_pch, NA)
    leg_lty <- c(leg_lty, 2)
    leg_lwd <- c(leg_lwd, 1)
  }
  graphics::legend("topleft", leg, col = leg_col, pch = leg_pch,
                   lty = leg_lty, lwd = leg_lwd, cex = 0.7, bg = "white",
                   box.col = "grey80")
  invisible(m)
}

#' Draw one stability against another, one point per module
#'
#' @param x,y Stabilities, one per module.
#' @param highlight Logical, modules drawn in the accent color.
#' @param xlab,ylab,main Labels.
#' @return The correlation, invisibly.
#' @export
PlotModuleStability <- function(x, y, highlight = rep(FALSE, length(x)),
                                xlab = "", ylab = "", main = "") {
  graphics::par(mar = c(4.2, 4.6, if (nzchar(main)) 2.5 else 1, 1))
  lim <- range(c(x, y), na.rm = TRUE)
  graphics::plot(x, y, xlim = lim, ylim = lim, pch = 19, cex = 0.6,
                 col = ifelse(highlight, "#2a78d6",
                              grDevices::adjustcolor("grey40", 0.6)),
                 xlab = xlab, ylab = ylab, main = main, las = 1, bty = "l")
  graphics::abline(0, 1, col = "grey60", lty = 2)
  if (any(highlight)) {
    graphics::legend("topleft", c("Interferon modules", "Other modules"),
                     pch = 19, col = c("#2a78d6", "grey40"), cex = 0.7,
                     bg = "white", box.col = "grey80")
  }
  invisible(stats::cor(x, y))
}
