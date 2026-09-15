# ------------------------------------------------------------------
# 23: Figure 1 -- conceptual diagram. Two panels share the same observed
#     inconsistent triad (A-B positive, A-C negative, B-C positive).
#
#     Panel A (structural balance alone): uniform nodes, uniform grey
#     edges, three equal arrows to three candidate outcomes. Balance
#     theory supplies the choice set and the pressure to adjust, not a
#     reason to prefer one candidate over another.
#
#     Panel B (asymmetric dependence added) visualizes the MECHANISM,
#     not just the prediction -- the earlier all-arrows version showed
#     H1's conclusion (thicker arrow toward the more asymmetric edge)
#     without showing why asymmetry should matter, which is exactly the
#     causal step Section 3 argues for. Here:
#       - node size = each state's illustrative total trade (A large, C
#         mid-sized, B small): the generative source of the asymmetry.
#       - edge width = schematic RESISTANCE TO REORIENTATION, the
#         inverse of that edge's dependence asymmetry -- a thin edge is
#         *not* a weak or unstable relationship (H2 shows asymmetry can
#         instead protect a relationship from outright termination);
#         it is a relationship with less room, on its more-exposed side,
#         to absorb pressure without its political sign giving way.
#         A-B, spanning large A and small B, is the most asymmetric tie
#         (.57) and is drawn thinnest; A-C, between two closer-sized
#         economies (.15), is drawn thickest.
#       - a single callout arrow marks H1's prediction on the thin,
#         low-resistance A-B edge, rather than three graduated arrows
#         implying a precise probability mapping H1 does not claim.
#     A-B is deliberately a positive (alliance) edge, not the
#     pre-existing negative A-C tie, so the figure does not visually
#     conflate "negative" with "changes."
#
#     Grayscale throughout; all values are illustrative round numbers,
#     not data. Base R graphics only.
# ------------------------------------------------------------------


A <- c(0, 1); B <- c(-0.87, -0.5); C <- c(0.87, -0.5)
nodes <- list(A = A, B = B, C = C)
centroid <- (A + B + C) / 3
sign <- c(AB = 1, AC = -1, BC = 1)
ends <- list(AB = list(A, B), AC = list(A, C), BC = list(B, C))
edge_label <- c(AB = "A-B changes", AC = "A-C changes", BC = "B-C changes")
asym <- c(AB = 0.57, AC = 0.15, BC = 0.29)

plot_base <- function(title, subtitle, legend_lines = NULL) {
  plot(NA, xlim = c(-2.5, 2.5), ylim = c(-1.85, 1.9), axes = FALSE, xlab = "", ylab = "",
       asp = 1, main = "", cex.main = 1.8, font.main = 2)
  mtext(title, side = 3, line = 4.4, cex = 1.55, font = 2)
  mtext(subtitle, side = 3, line = 2.6, cex = 1.0, font = 3)
  if (!is.null(legend_lines)) {
    for (i in seq_along(legend_lines)) {
      mtext(legend_lines[i], side = 3, line = 2.6 - 1.35 * i, cex = 0.8, font = 3, col = "grey30")
    }
  }
}

draw_nodes <- function(node_cex) {
  for (nm in names(nodes)) {
    p <- nodes[[nm]]
    points(p[1], p[2], pch = 21, bg = "white", col = "black", cex = node_cex[[nm]], lwd = 1.5)
    text(p[1], p[2], nm, cex = 1.2, font = 2)
  }
}

# ---- Panel A: structural balance alone ----
draw_panel_a <- function() {
  plot_base("A. Structural balance alone",
            "Creates pressure to adjust, not where adjustment lands")
  for (tn in names(sign)) {
    p1 <- ends[[tn]][[1]]; p2 <- ends[[tn]][[2]]
    lty <- if (sign[[tn]] > 0) 1 else 5
    lines(c(p1[1], p2[1]), c(p1[2], p2[2]), col = "grey55", lty = lty, lwd = 2.2, lend = 1)
  }
  draw_nodes(setNames(rep(3.0, 3), names(nodes)))
  for (tn in names(sign)) {
    p1 <- ends[[tn]][[1]]; p2 <- ends[[tn]][[2]]
    mid <- (p1 + p2) / 2
    dir <- (mid - centroid); dir <- dir / sqrt(sum(dir^2))
    arrows(mid[1] + dir[1]*0.34, mid[2] + dir[2]*0.34, mid[1] + dir[1]*0.85, mid[2] + dir[2]*0.85,
           lwd = 3.2, length = 0.09, angle = 20, col = "grey40", lend = 1)
    lab_pos <- mid + dir * 1.15
    text(lab_pos[1], lab_pos[2], edge_label[[tn]], cex = 1.14)
  }
}

# ---- Panel B: asymmetric dependence -- node size generates edge resistance,
#      which in turn structures arrow weight (relative likelihood) ----
draw_panel_b <- function() {
  plot_base("B. Asymmetric dependence added",
            "More asymmetric relationships are more likely to absorb the adjustment",
            legend_lines = c("Thinner ties are more asymmetric, and more likely to absorb the reorientation"))

  node_total <- c(A = 100, B = 8, C = 35)
  node_cex <- setNames(1.8 + 6.2 * sqrt(node_total / max(node_total)), names(node_total))

  # R's built-in dashed lty (5) uses a short dash/gap pattern that scales
  # with lwd, so a thick dashed line looks nearly solid. A custom, much
  # longer dash pattern ("8686": long dash, long gap) stays clearly dashed
  # at any line weight, so one line per edge can carry both sign and
  # resistance -- no separate halo layer needed.
  max_asym <- max(asym)
  edge_lwd <- c()
  for (tn in names(sign)) {
    p1 <- ends[[tn]][[1]]; p2 <- ends[[tn]][[2]]
    lty <- if (sign[[tn]] > 0) "solid" else "22"
    resistance <- (max_asym - asym[[tn]]) / max_asym   # 0 = A-B (thinnest), 1 = most resistant
    lwd <- 1.2 + 16 * resistance
    edge_lwd[tn] <- lwd
    lines(c(p1[1], p2[1]), c(p1[2], p2[2]), col = "black", lty = lty, lwd = lwd, lend = 1)
  }
  draw_nodes(node_cex)

  # arrows: dark/thick toward the low-resistance A-B edge (H1's prediction),
  # faint/thin toward the other two -- a two-tier contrast, not a precise
  # probability mapping across all three, since H1 claims a disproportionate
  # concentration, not a specific probability ratio.
  is_favored <- names(sign) == "AB"
  for (tn in names(sign)) {
    p1 <- ends[[tn]][[1]]; p2 <- ends[[tn]][[2]]
    mid <- (p1 + p2) / 2
    dir <- (mid - centroid); dir <- dir / sqrt(sum(dir^2))
    favored <- tn == "AB"
    # value label offset is the same fixed distance for all three edges
    # (just enough to clear the thickest edge, A-C), so the numbers sit
    # close to the edge and at a consistent distance across all three
    val_pos <- mid + dir * 0.22
    text(val_pos[1], val_pos[2], sprintf("%.2f", asym[[tn]]), cex = 0.98)
    # arrow starts close to the edge itself (just past the value label)
    # rather than floating well out toward the label
    start <- mid + dir * 0.40
    if (favored) {
      # arrow travels straight along the edge's own radial direction; the
      # "H1: most likely" / "A-B changes" block sits stacked and centered
      # directly above the arrow's tip, so the arrow points at the middle
      # of the block from below
      tip <- mid + dir * 1.0
      tag_pos <- tip + c(0, 0.16)
      lab_pos <- tag_pos + c(0, 0.24)
      arrows(start[1], start[2], tip[1], tip[2], lwd = 6.2, length = 0.14, angle = 20, col = "black", lend = 1)
      text(lab_pos[1], lab_pos[2], edge_label[[tn]], cex = 1.4, font = 2, col = "black")
      text(tag_pos[1], tag_pos[2], "H1: most likely", cex = 0.88, font = 4, col = "black")
    } else {
      lab_pos <- mid + dir * 1.05
      tip <- mid + dir * 0.75
      arrows(start[1], start[2], tip[1], tip[2], lwd = 1.3, length = 0.06, angle = 20, col = "grey72", lend = 1)
      text(lab_pos[1], lab_pos[2], edge_label[[tn]], cex = 1.02, font = 1, col = "grey40")
    }
  }
}

render_panels <- function() {
  layout(matrix(1:2, nrow = 1))
  par(mar = c(1, 1, 8.5, 1), oma = c(2, 0, 3.5, 0))
  draw_panel_a()
  draw_panel_b()
  mtext("Where does adjustment land?", side = 3, outer = TRUE, cex = 1.75, font = 2, line = 1.7)
  mtext("Structural imbalance creates the choice; dependence asymmetry allocates it.",
        side = 1, outer = TRUE, cex = 1.15, line = 0.6)
}

png("results/fig1_triad_diagram.png", width = 3000, height = 1700, res = 220)
render_panels()
dev.off()

pdf("results/fig1_triad_diagram.pdf", width = 13.6, height = 7.7)
render_panels()
dev.off()

cat("23_figure1_triad_diagram.R complete. Saved results/fig1_triad_diagram.png and .pdf\n")
