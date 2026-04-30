# Internal helper: prefer AME names, fallback to legacy AMR names
.get_est_component <- function(x, base_name) {
  if (!is.null(x[[base_name]])) return(x[[base_name]])
  legacy <- sub("^AME", "AMR", base_name)
  x[[legacy]]
}

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
#'
#' The default rendering (\code{style = "lines"})
#' is fully black-and-white-print safe: confidence intervals are drawn as
#' paired lines with distinct dash patterns, so they remain legible in
#' grayscale printing. Set \code{style = "shade"} for translucent bands.
#'
#' @param x A SpatialEffect object
#' @param smooth Plot smoothed estimates if available (default TRUE)
#' @param ci.type Which CI to plot: "conley", "permutation", or "both"
#' @param style Rendering style for confidence intervals: \code{"lines"}
#'   (default, B&W-print safe dashed lines) or \code{"shade"} (translucent
#'   polygons).
#' @param main Optional figure title. Defaults to no title (AER convention;
#'   titles typically appear as figure captions).
#' @param xlab,ylab Axis labels.
#' @param ... Additional arguments passed to plot()
#' @export
plot.SpatialEffect <- function(x,
                               smooth = TRUE,
                               ci.type = "conley",
                               style = c("lines", "shade"),
                               main = NULL,
                               xlab = "Distance",
                               ylab = "Average marginalized effect",
                               ...) {
  style <- match.arg(style)
  params <- x[["Parameters"]]
  dVec <- params$dVec
  ame_est <- .get_est_component(x, "AME_est")
  est <- ame_est$taud_est

  finite_est <- is.finite(est)
  if (!any(finite_est)) {
    warning("No finite AME estimates available for plotting.", call. = FALSE)
    return(invisible(NULL))
  }

  has_conley <- ci.type %in% c("conley", "both") && !is.null(x[["Conley.CI"]])
  has_perm   <- ci.type %in% c("permutation", "both") && !is.null(x[["Per.CI"]])
  ame_est_smoothed <- .get_est_component(x, "AME_est_smoothed")
  has_smooth <- isTRUE(smooth) && !is.null(ame_est_smoothed)
  has_sm_ci  <- has_smooth && !is.null(x[["smoothed.Conley.CI"]])
  has_sm_cb  <- has_smooth && !is.null(x[["smoothed.Conley.CB"]])

  ylim <- range(est[finite_est], na.rm = TRUE)
  if (has_conley) ylim <- range(ylim, x[["Conley.CI"]], na.rm = TRUE)
  if (has_perm)   ylim <- range(ylim, x[["Per.CI"]], na.rm = TRUE)
  if (has_sm_ci)  ylim <- range(ylim, x[["smoothed.Conley.CI"]], na.rm = TRUE)
  if (has_sm_cb)  ylim <- range(ylim, x[["smoothed.Conley.CB"]], na.rm = TRUE)
  pad <- diff(ylim) * 0.08
  ylim <- c(ylim[1] - pad, ylim[2] + pad)

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(family = "serif",
      mar = c(4.2, 4.6, if (is.null(main)) 1.2 else 2.8, 12.5),
      mgp = c(2.7, 0.6, 0),
      tcl = -0.3,
      cex.axis = 0.9,
      cex.lab = 1.0,
      las = 1)

  plot(dVec, est, type = "n",
       xlab = xlab, ylab = ylab, ylim = ylim,
       main = main, bty = "l", xaxs = "i", yaxs = "i", ...)

  # Line-type codes chosen so every CI is distinguishable in B&W print
  lty_conley <- 2   # dashed
  lty_perm   <- 3   # dotted
  lty_sm_ci  <- 5   # long dash
  lty_sm_cb  <- 4   # dot-dash

  abline(h = 0, col = "gray40", lwd = 0.8)

  if (style == "shade") {
    shade_band <- function(xs, lo, up, col) {
      ok <- is.finite(lo) & is.finite(up)
      if (!any(ok)) return(invisible(NULL))
      polygon(c(xs[ok], rev(xs[ok])), c(lo[ok], rev(up[ok])),
              col = col, border = NA)
    }
    if (has_perm)   shade_band(dVec, x[["Per.CI"]][, 1],    x[["Per.CI"]][, 2],
                               grDevices::adjustcolor("black", alpha.f = 0.08))
    if (has_conley) shade_band(dVec, x[["Conley.CI"]][, 1], x[["Conley.CI"]][, 2],
                               grDevices::adjustcolor("black", alpha.f = 0.18))
    if (has_sm_cb)  shade_band(ame_est_smoothed[, 1],
                               x[["smoothed.Conley.CB"]][, 1],
                               x[["smoothed.Conley.CB"]][, 2],
                               grDevices::adjustcolor("black", alpha.f = 0.10))
    if (has_sm_ci)  shade_band(ame_est_smoothed[, 1],
                               x[["smoothed.Conley.CI"]][, 1],
                               x[["smoothed.Conley.CI"]][, 2],
                               grDevices::adjustcolor("black", alpha.f = 0.22))
  } else {
    draw_ci_lines <- function(xs, lo, up, lty, lwd = 1.1, col = "black") {
      lines(xs, lo, lty = lty, lwd = lwd, col = col)
      lines(xs, up, lty = lty, lwd = lwd, col = col)
    }
    if (has_conley) draw_ci_lines(dVec, x[["Conley.CI"]][, 1], x[["Conley.CI"]][, 2],
                                  lty = lty_conley)
    if (has_perm)   draw_ci_lines(dVec, x[["Per.CI"]][, 1],    x[["Per.CI"]][, 2],
                                  lty = lty_perm)
    if (has_sm_ci)  draw_ci_lines(ame_est_smoothed[, 1],
                                  x[["smoothed.Conley.CI"]][, 1],
                                  x[["smoothed.Conley.CI"]][, 2],
                                  lty = lty_sm_ci, lwd = 1.2, col = "gray25")
    if (has_sm_cb)  draw_ci_lines(ame_est_smoothed[, 1],
                                  x[["smoothed.Conley.CB"]][, 1],
                                  x[["smoothed.Conley.CB"]][, 2],
                                  lty = lty_sm_cb, lwd = 1.2, col = "gray25")
  }

  lines(dVec, est, col = "black", lwd = 1.2)
  points(dVec, est, pch = 21, bg = "white", col = "black", cex = 0.7, lwd = 0.8)

  if (has_smooth) {
    lines(ame_est_smoothed[, 1], ame_est_smoothed[, 2],
          col = "black", lwd = 2.0)
  }

  box(bty = "l", lwd = 0.8)

  # Build legend
  lg_labels <- "Point estimate"
  lg_lty    <- 1
  lg_lwd    <- 1.2
  lg_pch    <- 21
  lg_col    <- "black"
  lg_ptbg   <- "white"
  add_entry <- function(label, lty = NA, lwd = NA, pch = NA,
                        col = "black", ptbg = NA) {
    lg_labels <<- c(lg_labels, label)
    lg_lty    <<- c(lg_lty, lty)
    lg_lwd    <<- c(lg_lwd, lwd)
    lg_pch    <<- c(lg_pch, pch)
    lg_col    <<- c(lg_col, col)
    lg_ptbg   <<- c(lg_ptbg, ptbg)
  }

  if (style == "shade") {
    if (has_conley) add_entry("Conley CI",       pch = 22, ptbg = grDevices::adjustcolor("black", 0.18))
    if (has_perm)   add_entry("Permutation CI",  pch = 22, ptbg = grDevices::adjustcolor("black", 0.08))
    if (has_smooth) add_entry("Smoothed estimate", lty = 1, lwd = 2.0)
    if (has_sm_ci)  add_entry("Smoothed CI",     pch = 22, ptbg = grDevices::adjustcolor("black", 0.22))
    if (has_sm_cb)  add_entry("Smoothed CB",     pch = 22, ptbg = grDevices::adjustcolor("black", 0.10))
  } else {
    if (has_conley) add_entry("Conley CI",       lty = lty_conley, lwd = 1.1)
    if (has_perm)   add_entry("Permutation CI",  lty = lty_perm,   lwd = 1.1)
    if (has_smooth) add_entry("Smoothed estimate", lty = 1, lwd = 2.0)
    if (has_sm_ci)  add_entry("Smoothed CI",     lty = lty_sm_ci, lwd = 1.2, col = "gray25")
    if (has_sm_cb)  add_entry("Smoothed CB",     lty = lty_sm_cb, lwd = 1.2, col = "gray25")
  }

  usr <- par("usr")
  op_xpd <- par(xpd = NA)
  on.exit(par(op_xpd), add = TRUE)
  legend(x         = usr[2] + (usr[2] - usr[1]) * 0.02,
         y         = usr[4] - (usr[4] - usr[3]) * 0.02,
         legend    = lg_labels,
         lty       = lg_lty,
         lwd       = lg_lwd,
         pch       = lg_pch,
         col       = lg_col,
         pt.bg     = lg_ptbg,
         pt.cex    = 1.0,
         bty       = "n",
         cex       = 0.8,
         x.intersp = 0.6,
         y.intersp = 1.15,
         seg.len   = 1.6)

  invisible(NULL)
}

#' Internal function to format and display results
#' @keywords internal
.display_results <- function(x, dVec.range = NULL) {
  ame_est <- .get_est_component(x, "AME_est")
  ame_est_smoothed <- .get_est_component(x, "AME_est_smoothed")
  output <- ame_est
  cnames <- c("dVec", "AME_est")

  if (x[["Parameters"]]$conley.se == 1 && !is.null(x[["Conley.CI"]])) {
    output <- cbind(output, x[["Conley.CI"]])
    cnames <- c(cnames, "Conley.CI.l", "Conley.CI.u")
  }
  if (x[["Parameters"]]$per.se == 1 && !is.null(x[["Per.CI"]])) {
    output <- cbind(output, x[["Per.CI"]])
    cnames <- c(cnames, "Per.CI.l", "Per.CI.u")
  }
  if (x[["Parameters"]]$smooth == 1 && !is.null(ame_est_smoothed)) {
    output <- cbind(output, ame_est_smoothed[, 2])
    cnames <- c(cnames, "AME_smoothed")
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
