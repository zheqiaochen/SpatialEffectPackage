## Jayachandran et al. (2017) replication using SpatialEffect2

rm(list = ls())

library(SpatialEffect2)
library(sf)
library(sp)
library(raster)
library(terra)
library(fields)
library(ggplot2)

set.seed(2024)

# ---- 0. Paths ----
.get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  for (i in sys.nframe():1) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(dirname(normalizePath(ofile)))
  }
  getwd()
}

script_dir <- .get_script_dir()
data_path <- file.path(script_dir, "..", "data")
graph_path <- file.path(script_dir, "..", "graphs")
graph_subdir <- file.path(graph_path, "J2017_graphs")
dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
dir.create(graph_subdir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Import data ----
# Outcome raster (2013 deforestation layer), coarsened as in original script.
ras <- raster(file.path(data_path, "rct_uganda_gfc_updated.tif"), band = 3)
ras <- aggregate(ras, 20, fun = mean)

# Polygon interventions with treatment label z \in {0,1}.
map <- st_read(file.path(data_path, "digitized_rct", "rct_studyarea_vectorized.shp"), quiet = TRUE)
map <- st_make_valid(map)

# ---- 2. Maps in original geographic CRS ----
ras_Z_geo <- st_as_sf(as(map, "Spatial"))
tr_colors <- ifelse(map$z == 1, "goldenrod1", "darkturquoise")

pdf(file = file.path(graph_path, "J2017-forest-map-polygon-circles-treated.pdf"), height = 8, width = 8)
par(mfrow = c(1, 1))
plot(ras, bty = "n", box = FALSE, col = c("white", "darkgoldenrod4"), legend = FALSE,
     cex.lab = 2, cex.main = 2, cex.axis = 2)
lines(ras_Z_geo[map$z == 0, ], col = "darkturquoise", lwd = 2)
lines(ras_Z_geo[map$z == 1, ], col = "goldenrod1", lwd = 2)
for (i in which(map$z == 1)) {
  buffer_treated <- buffer(as(ras_Z_geo[i, ], "Spatial"), width = 0.015 * 111000)
  lines(buffer_treated, col = "goldenrod1", lwd = 2, lty = 2)
}
legend("bottomright", col = c("darkturquoise", "goldenrod1", "darkgoldenrod4"),
       legend = c("Control", "Treated", "Deforestation"),
       lty = c(1, 1, NA), pch = c(NA, NA, 19), cex = 1.5, pt.cex = c(.8, .8, .4), bty = "n")
dev.off()

pdf(file = file.path(graph_path, "J2017-forest-map-polygon-circles-control.pdf"), height = 8, width = 8)
par(mfrow = c(1, 1))
plot(ras, bty = "n", box = FALSE, col = c("white", "darkgoldenrod4"), legend = FALSE,
     cex.lab = 2, cex.main = 2, cex.axis = 2)
lines(ras_Z_geo[map$z == 0, ], col = "darkturquoise", lwd = 2)
lines(ras_Z_geo[map$z == 1, ], col = "goldenrod1", lwd = 2)
for (i in which(map$z == 0)) {
  buffer_control <- buffer(as(ras_Z_geo[i, ], "Spatial"), width = 0.015 * 111000)
  lines(buffer_control, col = "darkturquoise", lwd = 2, lty = 2)
}
legend("bottomright", col = c("darkturquoise", "goldenrod1", "darkgoldenrod4"),
       legend = c("Control", "Treated", "Deforestation"),
       lty = c(1, 1, NA), pch = c(NA, NA, 19), cex = 1.5, pt.cex = c(.8, .8, .4), bty = "n")
dev.off()

# ---- 3. Reproject to kilometer units for Euclidean donut distances ----
sr <- "+proj=utm +zone=15 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
ras_proj <- projectRaster(ras, crs = sr)
map_sp <- as(map, "Spatial")
map_proj_sp <- spTransform(map_sp, crs(ras_proj))
ras_Z <- st_as_sf(map_proj_sp)

centr <- st_coordinates(st_centroid(ras_Z))
Zdata <- data.frame(
  x_coord = as.numeric(centr[, 1]),
  y_coord = as.numeric(centr[, 2]),
  treatment = as.numeric(ras_Z$z)
)

cat("Intervention nodes:", nrow(Zdata), "\n")
cat("Treated:", sum(Zdata$treatment == 1), "Control:", sum(Zdata$treatment == 0), "\n")

# ---- 4. SpatialEffect2 analysis ----
dVec <- seq(from = 0, to = 30, by = 1) # kilometers

result.list <- SpatialEffect(
  ras = ras_proj,
  Zdata = Zdata,
  ras_Z = ras_Z,
  treatment = "treatment",
  dVec = dVec,
  numpts = 1000,
  only.unique = 0,
  smooth = 1,
  per.se = 1,
  conley.se = 1,
  cutoff = 6,
  alpha = 0.1,
  edf = FALSE,
  bw = 5,
  bw_debias = 5,
  bias_correction = TRUE,
  nPerms = 500,
  smooth.conley.se = 1,
  dist.metric = "Euclidean",
  cType = "donut",
  kernel = "uni",
  n_threads = 4L
)

save(result.list, file = file.path(data_path, "forestResult_reproj.RData"))

summary(result.list, dVec.range = c(1, 15))

AMR_est <- result.list[["AMR_est"]]
Per.CI <- result.list[["Per.CI"]]
Conley.SE <- result.list[["Conley.SE"]]
Conley.CI <- result.list[["Conley.CI"]]
AMR_est_smoothed <- result.list[["AMR_est_smoothed"]]
smoothed.Conley.SE <- result.list[["smoothed.Conley.SE"]]
smoothed.Conley.CI <- result.list[["smoothed.Conley.CI"]]

test.result <- SpatialEffectTest(result.list, test.range = c(1, 5), smooth = 0, alpha = 0.1)
cat("Test statistic:", test.result$test.stat, "\n")
cat("Test CI:", test.result$test.CI, "\n")

# ---- 5. Main plots (match original naming) ----
show_n <- min(17, nrow(AMR_est))
dVec_real <- dVec[seq_len(show_n)]

pdf(file = file.path(graph_path, "J2017_AME_polygon_unsmoothed_reproj_90CI.pdf"), height = 8, width = 8)
par(mar = c(5, 5, 5, 5))
plot(AMR_est[seq_len(show_n), 2] ~ dVec_real,
     ylab = "Estimates", xlab = "Distance (kilometers)",
     type = "l", ylim = c(-0.04, 0.04),
     main = "Hajek Estimator with Its CIs",
     cex.lab = 2, cex.main = 2, cex.axis = 2, lwd = 3)
points(Conley.CI[seq_len(show_n), 1] ~ dVec_real, type = "l", col = "black", lty = "dotted", lwd = 3)
points(Conley.CI[seq_len(show_n), 2] ~ dVec_real, type = "l", col = "black", lty = "dotted", lwd = 3)
points(Per.CI[seq_len(show_n), 1] ~ dVec_real, type = "l", col = "blue", lty = "dotted", lwd = 3)
points(Per.CI[seq_len(show_n), 2] ~ dVec_real, type = "l", col = "blue", lty = "dotted", lwd = 3)
abline(h = 0, lty = 2, lwd = 3, col = "gray")
legend("topright", col = c("black", "black", "blue"), lty = c(1, 3, 3), lwd = c(3, 3, 3),
       legend = c("AME estimates", "Spatial HAC 90% CI", "Quantiles under sharp null"), cex = 2, bty = "n")
dev.off()

pdf(file = file.path(graph_path, "J2017_AME_polygon_smoothed_reproj_90CI.pdf"), height = 8, width = 8)
par(mar = c(5, 5, 5, 5))
plot(AMR_est_smoothed[seq_len(show_n), 2] ~ dVec_real,
     ylab = "Estimates", xlab = "Distance (kilometers)",
     type = "l", ylim = c(-0.04, 0.04),
     main = "Kernel Regression Estimator with Its CIs",
     cex.lab = 2, cex.main = 2, cex.axis = 2, lwd = 3)
points(smoothed.Conley.CI[seq_len(show_n), 1] ~ dVec_real, type = "l", col = "black", lty = "dotted", lwd = 3)
points(smoothed.Conley.CI[seq_len(show_n), 2] ~ dVec_real, type = "l", col = "black", lty = "dotted", lwd = 3)
abline(h = 0, lty = 2, lwd = 3, col = "gray")
legend("topright", col = c("black", "black"), lty = c(1, 3), lwd = c(3, 3),
       legend = c("AME estimates", "Spatial HAC 90% CI"), cex = 2, bty = "n")
dev.off()

# ---- 6. Coefficient-style figures ----
ha_ests <- AMR_est[seq_len(show_n), 2]
ha_CI95_upper <- Conley.CI[seq_len(show_n), 2]
ha_CI95_lower <- Conley.CI[seq_len(show_n), 1]
ha_CI90_upper <- AMR_est[seq_len(show_n), 2] + qnorm(0.95) * Conley.SE[seq_len(show_n)]
ha_CI90_lower <- AMR_est[seq_len(show_n), 2] - qnorm(0.95) * Conley.SE[seq_len(show_n)]

ha.data <- data.frame(
  estimates = ha_ests * 100,
  CI95_u = ha_CI95_upper * 100,
  CI95_l = ha_CI95_lower * 100,
  CI90_u = ha_CI90_upper * 100,
  CI90_l = ha_CI90_lower * 100,
  names = dVec_real
)
ha.data <- ha.data[
  is.finite(ha.data$estimates) &
    is.finite(ha.data$CI95_u) &
    is.finite(ha.data$CI95_l) &
    is.finite(ha.data$CI90_u) &
    is.finite(ha.data$CI90_l),
  , drop = FALSE
]
ha.data$names <- factor(ha.data$names, levels = rev(rev(ha.data$names)))
colors <- rep("red", nrow(ha.data))

coefs_p <- ggplot(ha.data) +
  geom_pointrange(aes(x = names, y = estimates, ymin = CI95_l, ymax = CI95_u), size = 0.3, colour = colors, lty = "dashed") +
  geom_pointrange(aes(x = names, y = estimates, ymin = CI90_l, ymax = CI90_u), size = 0.6, colour = colors, lty = "solid") +
  theme_bw() + geom_hline(aes(yintercept = 0), colour = "black", lty = 2) + xlab("Distance (kilometers)") +
  scale_linetype_manual(labels = c("90% CI", "95% CI"), values = c(1, 2)) + ylim(-3, 1) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 25),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        legend.position = "bottom", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + ylab("Estimates")
ggsave(file.path(graph_subdir, "ha_coefs_reproj_90CI.pdf"), coefs_p, device = "pdf", height = 8, width = 8)

smooth_ests <- AMR_est_smoothed[seq_len(show_n), 2]
smooth_CI95_upper <- smoothed.Conley.CI[seq_len(show_n), 2]
smooth_CI95_lower <- smoothed.Conley.CI[seq_len(show_n), 1]
smooth_CI90_upper <- AMR_est_smoothed[seq_len(show_n), 2] + qnorm(0.95) * smoothed.Conley.SE[seq_len(show_n)]
smooth_CI90_lower <- AMR_est_smoothed[seq_len(show_n), 2] - qnorm(0.95) * smoothed.Conley.SE[seq_len(show_n)]

smooth.data <- data.frame(
  estimates = smooth_ests * 100,
  CI95_u = smooth_CI95_upper * 100,
  CI95_l = smooth_CI95_lower * 100,
  CI90_u = smooth_CI90_upper * 100,
  CI90_l = smooth_CI90_lower * 100,
  names = dVec_real
)
smooth.data <- smooth.data[
  is.finite(smooth.data$estimates) &
    is.finite(smooth.data$CI95_u) &
    is.finite(smooth.data$CI95_l) &
    is.finite(smooth.data$CI90_u) &
    is.finite(smooth.data$CI90_l),
  , drop = FALSE
]
smooth.data$names <- factor(smooth.data$names, levels = rev(rev(smooth.data$names)))

smooth_p <- ggplot(smooth.data) +
  geom_pointrange(aes(x = names, y = estimates, ymin = CI95_l, ymax = CI95_u), size = 0.3, colour = "red", lty = "dashed") +
  geom_pointrange(aes(x = names, y = estimates, ymin = CI90_l, ymax = CI90_u), size = 0.6, colour = "red", lty = "solid") +
  theme_bw() + geom_hline(aes(yintercept = 0), colour = "black", lty = 2) + xlab("Distance (kilometers)") +
  scale_linetype_manual(labels = c("90% CI", "95% CI"), values = c(1, 2)) + ylim(-3, 1) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 25),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        legend.position = "bottom", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + ylab("Estimates")
ggsave(file.path(graph_subdir, "smoothed_coefs_reproj_90CI.pdf"), smooth_p, device = "pdf", height = 8, width = 8)

cat("\n=== Done (SpatialEffect2) ===\n")
