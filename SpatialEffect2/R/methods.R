#' Print method for SpatialEffect objects
#' @param x A SpatialEffect object
#' @param dVec.range Optional 2-element vector to filter displayed distances
#' @param ... Ignored
#' @export
print.SpatialEffect <- function(x, dVec.range = NULL, ...) {
  cat("Call: SpatialEffect\n\n")
  .display_results(x, dVec.range)
}

#' Summary method for SpatialEffect objects
#' @param object A SpatialEffect object
#' @param dVec.range Optional 2-element vector to filter displayed distances
#' @param ... Ignored
#' @export
summary.SpatialEffect <- function(object, dVec.range = NULL, ...) {
  cat("Call: SpatialEffect\n\n")
  .display_results(object, dVec.range)
}

#' Plot method for SpatialEffect objects
#' @param x A SpatialEffect object
#' @param smooth Plot smoothed estimates if available (default TRUE)
#' @param ci.type Which CI to plot: "conley", "permutation", or "both"
#' @param ... Additional arguments passed to plot()
#' @export
plot.SpatialEffect <- function(x, smooth = TRUE, ci.type = "conley", ...) {
  params <- x[["Parameters"]]
  dVec <- params$dVec
  est <- x[["AMR_est"]]$taud_est

  # Determine y-axis range
  ylim <- range(est, na.rm = TRUE)
  if (!is.null(x[["Conley.CI"]])) ylim <- range(ylim, x[["Conley.CI"]], na.rm = TRUE)
  if (!is.null(x[["Per.CI"]])) ylim <- range(ylim, x[["Per.CI"]], na.rm = TRUE)
  ylim <- ylim * 1.1

  plot(dVec, est, type = "b", pch = 16, col = "black",
       xlab = "Distance", ylab = "AME Estimate", ylim = ylim,
       main = "Average Marginalized Effect (AME)", ...)
  abline(h = 0, lty = 2, col = "gray50")

  if (ci.type %in% c("conley", "both") && !is.null(x[["Conley.CI"]])) {
    lines(dVec, x[["Conley.CI"]][, 1], lty = 2, col = "blue")
    lines(dVec, x[["Conley.CI"]][, 2], lty = 2, col = "blue")
  }
  if (ci.type %in% c("permutation", "both") && !is.null(x[["Per.CI"]])) {
    lines(dVec, x[["Per.CI"]][, 1], lty = 3, col = "red")
    lines(dVec, x[["Per.CI"]][, 2], lty = 3, col = "red")
  }

  if (smooth && !is.null(x[["AMR_est_smoothed"]])) {
    lines(x[["AMR_est_smoothed"]][, 1], x[["AMR_est_smoothed"]][, 2],
          col = "darkgreen", lwd = 2)
    if (!is.null(x[["smoothed.Conley.CI"]])) {
      lines(x[["AMR_est_smoothed"]][, 1], x[["smoothed.Conley.CI"]][, 1],
            lty = 2, col = "darkgreen")
      lines(x[["AMR_est_smoothed"]][, 1], x[["smoothed.Conley.CI"]][, 2],
            lty = 2, col = "darkgreen")
    }
    if (!is.null(x[["smoothed.Conley.CB"]])) {
      lines(x[["AMR_est_smoothed"]][, 1], x[["smoothed.Conley.CB"]][, 1],
            lty = 4, col = "purple")
      lines(x[["AMR_est_smoothed"]][, 1], x[["smoothed.Conley.CB"]][, 2],
            lty = 4, col = "purple")
    }
  }

  # Legend
  legend_labels <- "AME"
  legend_cols <- "black"
  legend_lty <- 1
  if (ci.type %in% c("conley", "both") && !is.null(x[["Conley.CI"]])) {
    legend_labels <- c(legend_labels, "Conley CI")
    legend_cols <- c(legend_cols, "blue")
    legend_lty <- c(legend_lty, 2)
  }
  if (ci.type %in% c("permutation", "both") && !is.null(x[["Per.CI"]])) {
    legend_labels <- c(legend_labels, "Permutation CI")
    legend_cols <- c(legend_cols, "red")
    legend_lty <- c(legend_lty, 3)
  }
  if (smooth && !is.null(x[["AMR_est_smoothed"]])) {
    legend_labels <- c(legend_labels, "Smoothed")
    legend_cols <- c(legend_cols, "darkgreen")
    legend_lty <- c(legend_lty, 1)
  }
  legend("topright", legend = legend_labels, col = legend_cols,
         lty = legend_lty, bty = "n", cex = 0.8)
}

#' Internal function to format and display results
#' @keywords internal
.display_results <- function(x, dVec.range = NULL) {
  output <- x[["AMR_est"]]
  cnames <- c("dVec", "AMR_est")

  if (x[["Parameters"]]$conley.se == 1 && !is.null(x[["Conley.CI"]])) {
    output <- cbind(output, x[["Conley.CI"]])
    cnames <- c(cnames, "Conley.CI.l", "Conley.CI.u")
  }
  if (x[["Parameters"]]$per.se == 1 && !is.null(x[["Per.CI"]])) {
    output <- cbind(output, x[["Per.CI"]])
    cnames <- c(cnames, "Per.CI.l", "Per.CI.u")
  }
  if (x[["Parameters"]]$smooth == 1 && !is.null(x[["AMR_est_smoothed"]])) {
    output <- cbind(output, x[["AMR_est_smoothed"]][, 2])
    cnames <- c(cnames, "AMR_smoothed")
  }
  if (x[["Parameters"]]$smooth.conley.se == 1 && !is.null(x[["smoothed.Conley.CI"]])) {
    output <- cbind(output, x[["smoothed.Conley.CI"]])
    cnames <- c(cnames, "sm.Conley.CI.l", "sm.Conley.CI.u")
  }

  names(output) <- cnames
  output <- apply(as.matrix(output), 2, function(col) round(col, 3))

  if (!is.null(dVec.range)) {
    if (length(dVec.range) != 2) stop("dVec.range must have exactly 2 values")
    if (dVec.range[1] > dVec.range[2]) dVec.range <- rev(dVec.range)
    keep <- output[, 1] >= dVec.range[1] & output[, 1] <= dVec.range[2]
    output <- output[keep, , drop = FALSE]
  }
  print(output)
}
