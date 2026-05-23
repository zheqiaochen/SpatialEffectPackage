assets_dir <- "/Users/zheqiao/Documents/cursor_projects/SpatialEffectPackage/SpatialEffect2/pkgdown/assets"
dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(terra))

set.seed(1)
ras <- rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40)
values(ras) <- rnorm(ncell(ras))

nz <- 60
Zdata <- data.frame(
  x = runif(nz, 2, 38),
  y = runif(nz, 2, 38),
  treat = rbinom(nz, 1, 0.5)
)

cols <- hcl.colors(60, "YlGnBu", rev = TRUE)

png(file.path(assets_dir, "outcome-surface-ras.png"), width = 1200, height = 820, res = 170)
plot(ras, main = "Example Outcome Surface (ras)", col = cols)
dev.off()

png(file.path(assets_dir, "outcome-surface-with-nodes.png"), width = 1200, height = 860, res = 170)
plot(ras, main = "Outcome Surface with Intervention Nodes", col = cols)
points(
  Zdata$x, Zdata$y,
  pch = 21,
  bg = ifelse(Zdata$treat == 1, "#D94841", "#2C7FB8"),
  col = "white",
  cex = 0.95
)
legend(
  "topright",
  legend = c("treated", "control"),
  pch = 21,
  pt.bg = c("#D94841", "#2C7FB8"),
  col = "white",
  bty = "n"
)
dev.off()

cat("Wrote:\n")
cat(file.path(assets_dir, "outcome-surface-ras.png"), "\n")
cat(file.path(assets_dir, "outcome-surface-with-nodes.png"), "\n")
