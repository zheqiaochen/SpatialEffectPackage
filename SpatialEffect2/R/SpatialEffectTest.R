#' Sharp null randomization test for spatial treatment effects
#'
#' Tests whether the AME is significantly different from zero over a specified
#' range of distances, using Fisher's sharp null hypothesis.
#'
#' @param result.list A SpatialEffect object from \code{\link{SpatialEffect}}.
#' @param test.range Numeric vector of length 2 specifying the distance range to test.
#' @param smooth Whether to use smoothed estimates (0 or 1).
#' @param alpha Significance level (default 0.05).
#'
#' @return A list with:
#' \describe{
#'   \item{test.stat}{The observed test statistic (sum of AME over test range)}
#'   \item{test.CI}{Permutation confidence interval under the null}
#' }
#'
#' @export
SpatialEffectTest <- function(result.list, test.range, smooth = 0, alpha = 0.05) {

  params <- result.list[["Parameters"]]
  dVec <- params$dVec
  Zup <- params$Zup
  Ybards <- params$Ybards
  Sdata <- params$Sdata
  blockvar <- params$blockvar
  clustvar <- params$clustvar
  nPerms <- params$nPerms
  n_threads <- if (!is.null(params$n_threads)) params$n_threads else 1L
  AMR_est <- result.list[["AMR_est"]]
  nz <- length(Zup)

  if (is.null(test.range)) {
    warning("test.range is not specified, using full dVec range")
    test.range <- c(min(dVec), max(dVec))
  }
  if (length(test.range) != 2) {
    warning("test.range must have exactly 2 values, using range bounds")
    test.range <- range(test.range)
  }
  if (test.range[1] > test.range[2])
    test.range <- rev(test.range)
  test.range[1] <- max(test.range[1], min(dVec))
  test.range[2] <- min(test.range[2], max(dVec))

  in_range <- AMR_est$d >= test.range[1] & AMR_est$d <= test.range[2]

  if (smooth == 1) {
    AMR_est_smoothed <- result.list[["AMR_est_smoothed"]]
    test.stat <- sum(AMR_est_smoothed[AMR_est_smoothed[, 1] >= test.range[1] &
                                       AMR_est_smoothed[, 1] <= test.range[2], 2])
    Sdata$zIndex <- rep(seq_len(nz), length(dVec))
  } else {
    test.stat <- sum(AMR_est$taud_est[in_range])
  }

  permMat <- .gen_perms(Zup, blockvar = blockvar, clustvar = clustvar,
                         maxiter = nPerms)
  test.per <- rep(NA_real_, nPerms)

  for (i in seq_len(ncol(permMat))) {
    if (smooth == 0) {
      z_perm <- permMat[, i]
      AMR_per <- as.numeric(z_perm %*% Ybards / sum(z_perm) -
                             (1 - z_perm) %*% Ybards / sum(1 - z_perm))
      test.per[i] <- sum(AMR_per[in_range])
    } else {
      treatment_per <- data.frame(zIndex = seq_len(nz), treatment.per = permMat[, i])
      Sdata_per <- merge(Sdata, treatment_per, by = "zIndex")
      Sdata_per$treatment <- Sdata_per$treatment.per
      bw <- params$bw
      bw_debias <- params$bw_debias
      kernel <- "tri"  # default for smoothed test
      AMR_per <- LocalReg(dVec, Sdata_per, bw, bw_debias, Zup, xevals = NULL,
                           smooth.conley.se = 0, kernel = kernel, cutoff = 0,
                           dist.metric = "Euclidean",
                           bias_correction = !is.null(bw_debias),
                           n_threads = n_threads)$coefs[, 2]
      test.per[i] <- sum(AMR_per[AMR_est_smoothed[, 1] >= test.range[1] &
                                   AMR_est_smoothed[, 1] <= test.range[2]])
    }
  }

  test.CI <- quantile(test.per, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  list(test.stat = test.stat, test.CI = test.CI)
}
