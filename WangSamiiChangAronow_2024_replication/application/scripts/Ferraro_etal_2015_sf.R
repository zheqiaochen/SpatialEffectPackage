## Ferraro et al. (2015) replication using SpatialEffect2 (sf/terra based)
## Adapted from the original raster/sp version

rm(list = ls())

library(SpatialEffect2)
library(sf)
library(terra)
library(fields)

set.seed(2024)

## Detect script directory robustly
.get_script_dir <- function() {
  # Rscript via commandArgs
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  # source()
  for (i in sys.nframe():1) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  # fallback
  return(getwd())
}
script_dir <- .get_script_dir()
data_path <- file.path(script_dir, "..", "data")
graph_path <- file.path(script_dir, "..", "graphs")

# ---- 1. Load spatial data ----
protected_areas <- st_read(file.path(data_path, "Costa Rica/PAs_Reproj/PAs_before80_1.shp"), quiet = TRUE)
forest <- st_read(file.path(data_path, "Costa Rica/Deforestation and Carbon/Spatial/final_pixels_w_covs.shp"), quiet = TRUE)
forest <- st_transform(forest, st_crs(protected_areas))

# Fix invalid geometries
invalid_pa <- st_is_valid(protected_areas, NA_on_exception = TRUE)
cat("Number of invalid PA geometries:", sum(!invalid_pa, na.rm = TRUE), "\n")
protected_areas <- st_make_valid(protected_areas)

crs_proj4 <- st_crs(protected_areas)$proj4string

# ---- 2. Create buffers and rasterize (using terra) ----
buffer_distance <- 5000
pa_buffers <- st_buffer(protected_areas, dist = buffer_distance)

full_extent_sf <- st_union(st_union(protected_areas, pa_buffers))
bbox <- st_bbox(full_extent_sf)
full_ext <- ext(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"])

res_coarse <- 5000
blank_rast <- rast(full_ext, res = res_coarse, crs = crs_proj4)

# Align CRS to the raster's CRS
st_crs(protected_areas) <- crs(blank_rast)
st_crs(pa_buffers) <- crs(blank_rast)

pa_rast <- rasterize(vect(protected_areas), blank_rast, field = 1, fun = "mean", background = 0)
buf_rast <- rasterize(vect(pa_buffers), blank_rast, field = 1, fun = "mean", background = 0)
pa_rast[pa_rast == 0] <- NA
buf_rast[buf_rast == 0] <- NA
buf_rast[buf_rast == 1] <- 0

combined_rast <- cover(pa_rast, buf_rast)
cat("Combined raster values:\n")
print(table(values(combined_rast)[, 1], useNA = "ifany"))

# ---- 3. Boundary detection via focal ----
boundary_function <- function(x) {
  if (length(unique(x[!is.na(x)])) > 1) 1 else NA
}

window_mat <- matrix(1, 3, 3)
boundary_rast <- focal(combined_rast, w = window_mat, fun = boundary_function, na.policy = "omit")
boundary_rast[is.na(combined_rast)] <- NA
boundary_rast[combined_rast == 0] <- 0

cat("Boundary raster values:\n")
print(table(values(boundary_rast)[, 1], useNA = "ifany"))

# ---- 4. Outcome raster ----
mean(forest$defor97[forest$prot_b80 == 1]) - mean(forest$defor97[forest$prot_b80 == 0])

res_fine <- 1000
rast_template <- rast(ext(vect(forest)), res = res_fine, crs = crs_proj4)

ras <- rasterize(vect(forest), rast_template, field = "defor97", fun = mean, na.rm = TRUE)
ras[is.na(ras)] <- 0

# ---- 5. Covariate rasters ----
forest$luc_high <- as.numeric(forest$luc1 == 1 | forest$luc2 == 1 | forest$luc3 == 1)
ras_cov1 <- rasterize(vect(forest), rast_template, field = "luc_high", fun = mean, na.rm = TRUE)
ras_cov1[is.na(ras_cov1)] <- 0
ras_cov2 <- rasterize(vect(forest), rast_template, field = "dist_rd_69", fun = mean, na.rm = TRUE)
ras_cov2[is.na(ras_cov2)] <- 0
ras_cov3 <- rasterize(vect(forest), rast_template, field = "dist_mcity", fun = mean, na.rm = TRUE)
ras_cov3[is.na(ras_cov3)] <- 0

values_cov1 <- values(ras_cov1)[, 1]
values_cov2 <- values(ras_cov2)[, 1]
values_cov3 <- values(ras_cov3)[, 1]

# ---- 6. Build Zdata from boundary polygons ----
ras_Z_poly <- as.polygons(boundary_rast)
ras_Z_sf <- st_as_sf(ras_Z_poly)

# Get cell indices of boundary tiles
com_index <- which(!is.na(values(boundary_rast)[, 1]))

cov1 <- cov2 <- cov3 <- c()
num_tiles <- c()
boundary_vals <- values(boundary_rast)[, 1]

for (i in seq_len(nrow(ras_Z_sf))) {
  one_poly <- vect(ras_Z_sf[i, ])
  one_rasterized <- rasterize(one_poly, ras)
  one_vals <- values(one_rasterized)[, 1]
  valid_cells <- which(!is.na(one_vals) & !is.nan(one_vals) & one_vals != 0)

  if (length(valid_cells) == 0) {
    # Mark as NA in boundary raster for this tile
    # (skip this polygon)
    next
  }
  coords <- xyFromCell(ras, valid_cells)
  num_tiles <- c(num_tiles, nrow(coords))
  grid_idx1 <- na.omit(cellFromXY(ras_cov1, coords))
  grid_idx2 <- na.omit(cellFromXY(ras_cov2, coords))
  grid_idx3 <- na.omit(cellFromXY(ras_cov3, coords))
  cov1 <- c(cov1, mean(values_cov1[grid_idx1], na.rm = TRUE))
  cov2 <- c(cov2, mean(values_cov2[grid_idx2], na.rm = TRUE))
  cov3 <- c(cov3, mean(values_cov3[grid_idx3], na.rm = TRUE))
}

# ---- 7. Build Zdata and propensity scores ----
Zdata <- data.frame(treatment = values(boundary_rast)[, 1])
Zdata <- Zdata[!is.na(Zdata$treatment), , drop = FALSE]
names(Zdata)[1] <- "treatment"
Zdata$cov1 <- cov1
Zdata$cov2 <- cov2
Zdata$cov3 <- cov3

pscore_fit <- glm(treatment ~ cov1 + cov2 + cov3, data = Zdata, family = "binomial")
Zdata$prob_treatment <- predict(pscore_fit, type = "response")

# Propensity score density plot
pdf(file = file.path(graph_path, "F2015-forest-pscore-sf.pdf"), height = 8, width = 8)
par(mfrow = c(1, 1), mar = c(5, 5, 4, 3))
plot(density(Zdata$prob_treatment),
     cex.lab = 2, cex.main = 2, cex.axis = 2,
     xlab = "Estimated propensity score",
     main = "Density of propensity score estimates")
dev.off()

# Refresh ras_Z_sf from boundary_rast
ras_Z_sf <- st_as_sf(as.polygons(boundary_rast))

# ---- 8. Set up SpatialEffect2 parameters ----
dVec <- seq(from = 0, to = 20000, by = 500)

# ---- 9. Run SpatialEffect2 analysis (polygon intervention, donut) ----
cat("\n=== Running SpatialEffect2 (polygon, donut) ===\n")
result.list <- SpatialEffect(
  ras = ras,
  Zdata = Zdata,
  ras_Z = ras_Z_sf,
  treatment = "treatment",
  prob_treatment = "prob_treatment",
  dVec = dVec,
  numpts = 1000,
  only.unique = 0,
  smooth = 1,
  per.se = 1,
  conley.se = 1,
  cutoff = 15000,
  alpha = 0.05,
  edf = FALSE,
  bw = NULL,
  bw_debias = NULL,
  bias_correction = TRUE,
  nPerms = 500,
  smooth.conley.se = 1,
  dist.metric = "Euclidean",
  cType = "donut",
  kernel = "tri",
  n_threads = 4L
)

summary(result.list, dVec.range = c(1000, 5000))

# ---- 10. Extract results ----
AMR_est <- result.list[["AMR_est"]]
Per.CI <- result.list[["Per.CI"]]
Conley.SE <- result.list[["Conley.SE"]]
Conley.CI <- result.list[["Conley.CI"]]
AMR_est_smoothed <- result.list[["AMR_est_smoothed"]]
smoothed.Conley.SE <- result.list[["smoothed.Conley.SE"]]
smoothed.Conley.CI <- result.list[["smoothed.Conley.CI"]]
smoothed.Conley.CB <- result.list[["smoothed.Conley.CB"]]

save(result.list, file = file.path(dirname(data_path), "data", "Ferraro_etal_2015_sf.RData"))

# ---- 11. Hypothesis test ----
test.result <- SpatialEffectTest(result.list, test.range = c(1000, 5000), smooth = 0)
cat("\nTest statistic:", test.result$test.stat, "\n")
cat("Test CI:", test.result$test.CI, "\n")

# ---- 12. Plots ----
dVec_real <- dVec / 1000

# Plot using built-in method
pdf(file = file.path(graph_path, "F2015_AMR_sf_plot.pdf"), height = 8, width = 8)
plot(result.list, ci.type = "both")
dev.off()

# Unsmoothed AMR plot
pdf(file = file.path(graph_path, "Ferraro_etal_AMR_polygon_unsmoothed_sf.pdf"), height = 8, width = 8)
par(mar = c(5, 5, 5, 5))
plot(AMR_est[, 2] ~ dVec_real,
     ylab = "Estimates", xlab = "Distance (kilometers)",
     type = "l", ylim = c(-0.04, 0.04),
     main = "Hajek Estimator with Its CIs",
     cex.lab = 2, cex.main = 2, cex.axis = 2, lwd = 3)
points(Conley.CI[, 1] ~ dVec_real, type = "l", col = "black", lty = "dotted", lwd = 3)
points(Conley.CI[, 2] ~ dVec_real, type = "l", col = "black", lty = "dotted", lwd = 3)
points(Per.CI[, 1] ~ dVec_real, type = "l", col = "blue", lty = "dotted", lwd = 3)
points(Per.CI[, 2] ~ dVec_real, type = "l", col = "blue", lty = "dotted", lwd = 3)
abline(h = 0, lty = 2, lwd = 3, col = "gray")
legend("topright", col = c("black", "black", "blue"), lty = c(1, 3, 3), lwd = 3,
       legend = c("AME estimates", "Spatial HAC 95% CI", "Quantiles under sharp null"),
       cex = 2, bty = "n")
dev.off()

# Smoothed AMR plot
if (!is.null(AMR_est_smoothed)) {
  pdf(file = file.path(graph_path, "Ferraro_etal_AMR_polygon_smoothed_sf.pdf"), height = 8, width = 8)
  par(mar = c(5, 5, 5, 5))
  plot(AMR_est_smoothed[, 2] ~ dVec_real,
       ylab = "Estimates", xlab = "Distance (kilometers)",
       type = "l", ylim = c(-0.04, 0.04),
       main = "Kernel Regression Estimator with Its CIs",
       cex.lab = 2, cex.main = 2, cex.axis = 2, lwd = 3)
  if (!is.null(smoothed.Conley.CI)) {
    points(smoothed.Conley.CI[, 1] ~ dVec_real, type = "l", col = "black", lty = "dotted", lwd = 3)
    points(smoothed.Conley.CI[, 2] ~ dVec_real, type = "l", col = "black", lty = "dotted", lwd = 3)
  }
  abline(h = 0, lty = 2, lwd = 3, col = "gray")
  legend("topright", col = c("black", "black"), lty = c(1, 3), lwd = 3,
         legend = c("AME estimates", "Spatial HAC 95% CI"), cex = 2, bty = "n")
  dev.off()
}

cat("\n=== Done ===\n")
