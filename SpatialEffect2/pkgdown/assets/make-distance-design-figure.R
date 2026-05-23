assets_dir <- "/Users/zheqiao/Documents/cursor_projects/SpatialEffectPackage/SpatialEffect2/pkgdown/assets"
dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)

base_col <- "#2F6EA9"
axis_col <- "#E7EBF0"
text_col <- "#334155"
title_col <- "#0F172A"
fill_col <- adjustcolor(base_col, alpha.f = 0.20)

circle_poly <- function(r, n = 500, cx = 0, cy = 0) {
  theta <- seq(0, 2 * pi, length.out = n)
  list(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

draw_base <- function(panel_title, panel_desc) {
  plot(
    NA, NA,
    xlim = c(-1.18, 1.18), ylim = c(-1.20, 1.18),
    xlab = "", ylab = "", axes = FALSE, asp = 1
  )
  abline(h = 0, v = 0, col = axis_col, lwd = 2)
  title(main = panel_title, col.main = title_col, font.main = 2, cex.main = 2.0, line = 0.9)
  # Keep blue subtitle clearly below the main title.
  mtext(panel_desc, side = 3, line = -0.95, col = base_col, cex = 1.28, font = 2)
  points(0, 0, pch = 21, bg = "#0F172A", col = "#0F172A", cex = 2.0)
  text(0, -1.11, "intervention node i", cex = 1.25, col = text_col)
}

draw_edge <- function(path) {
  png(path, width = 1250, height = 900, res = 180)
  par(mar = c(2.8, 2.2, 3.6, 1.2), bg = "white")
  r_outer <- 0.82
  draw_base('A. cType = "edge"', "boundary at distance d")
  edge <- circle_poly(r_outer)
  lines(edge$x, edge$y, col = base_col, lwd = 7)
  text(r_outer + 0.08, 0.10, "d", cex = 1.35, col = base_col, font = 2)
  dev.off()
}

draw_disk <- function(path) {
  png(path, width = 1250, height = 900, res = 180)
  par(mar = c(2.8, 2.2, 3.6, 1.2), bg = "white")
  r_outer <- 0.82
  draw_base('B. cType = "disk"', "cumulative region [0, d]")
  disk <- circle_poly(r_outer)
  polygon(disk$x, disk$y, col = fill_col, border = base_col, lwd = 7)
  points(0, 0, pch = 21, bg = "#0F172A", col = "#0F172A", cex = 2.0)
  text(r_outer + 0.08, 0.10, "d", cex = 1.35, col = base_col, font = 2)
  dev.off()
}

draw_donut <- function(path) {
  png(path, width = 1250, height = 900, res = 180)
  par(mar = c(2.8, 2.2, 3.6, 1.2), bg = "white")
  r_inner <- 0.42
  r_outer <- 0.82
  draw_base('C. cType = "donut"', "adjacent-ring difference [d-kappa, d]")
  outer <- circle_poly(r_outer)
  inner <- circle_poly(r_inner)
  polygon(outer$x, outer$y, col = fill_col, border = base_col, lwd = 7)
  polygon(inner$x, inner$y, col = "white", border = "white")
  lines(inner$x, inner$y, col = base_col, lwd = 5, lty = 2)
  points(0, 0, pch = 21, bg = "#0F172A", col = "#0F172A", cex = 2.0)
  text(r_outer + 0.08, 0.10, "d", cex = 1.35, col = base_col, font = 2)
  text(-(r_inner + 0.10), 0.24, "d-kappa", cex = 1.15, col = base_col, font = 2)
  dev.off()
}

edge_path <- file.path(assets_dir, "distance-design-edge.png")
disk_path <- file.path(assets_dir, "distance-design-disk.png")
donut_path <- file.path(assets_dir, "distance-design-donut.png")

draw_edge(edge_path)
draw_disk(disk_path)
draw_donut(donut_path)

cat("Wrote:\n", edge_path, "\n", disk_path, "\n", donut_path, "\n")
