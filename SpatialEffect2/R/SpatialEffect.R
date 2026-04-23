#' Estimate the Average Marginalized Effect (AME) Curve
#'
#' Estimates how treatment effects vary with distance in space for spatial
#' experiments under unknown interference, following Wang, Samii, Chang and
#' Aronow (2024).
#'
#' @param ras A \code{\link[terra]{SpatRaster}} of the outcome surface, or NULL
#'   to use kriging interpolation.
#' @param Ydata Data frame of outcome observations (used with kriging when
#'   outcome locations differ from intervention nodes).
#' @param outcome Character name of the outcome variable column.
#' @param x_coord_Y,y_coord_Y Character names of coordinate columns in Ydata.
#' @param ras_Z An \code{\link[sf]{sf}} object for polygon interventions, or
#'   NULL for point interventions. Can also be a SpatRaster or SpatVector which
#'   will be converted to sf.
#' @param Zdata Data frame of intervention node attributes.
#' @param x_coord_Z,y_coord_Z Character names of coordinate columns in Zdata
#'   (required for point interventions).
#' @param treatment Character name of the binary treatment column in Zdata.
#' @param covs Optional matrix of covariates (nz rows or nz*length(dVec) rows).
#' @param prob_treatment Character name of propensity score column in Zdata for
#'   IPW (observational studies).
#' @param dVec Numeric vector of distance values at which to evaluate the AME.
#' @param dist.metric "Euclidean" (default) or "Geodesic".
#' @param cType Circle type: "edge" (default), "disk", or "donut".
#' @param numpts Number of sampling points on each circle (auto-selected if NULL).
#' @param evalpts Multiplier for auto-selected numpts (default 1).
#' @param only.unique Deduplicate raster cells per node (0 or 1).
#' @param smooth Whether to apply local polynomial smoothing (0 or 1).
#' @param bw Bandwidth for smoothing (auto-selected via CV if NULL).
#' @param bw_debias Bandwidth for bias correction (auto-selected if NULL).
#' @param bias_correction Use bias-corrected smoothing (default TRUE).
#' @param smooth.conley.se Compute Conley SE for smoothed estimates (0 or 1).
#' @param conf.band Compute uniform confidence bands for smoothed AME (0 or 1).
#' @param per.se Compute permutation-based CIs (0 or 1, default 1).
#' @param blockvar Block variable for stratified permutation.
#' @param clustvar Cluster variable for clustered permutation.
#' @param conley.se Compute Conley spatial HAC SEs (0 or 1, default 1).
#' @param kernel Kernel for HAC: "uni"/"uniform", "tri"/"triangular",
#'   "epa"/"epanechnikov".
#' @param cutoff Spatial bandwidth for Conley SE (0 = EHW).
#' @param alpha Significance level (default 0.05).
#' @param edf Use effective degrees of freedom adjustment (default FALSE).
#' @param m Polynomial order for kriging (default 2).
#' @param lambda Smoothing parameter for kriging (default 0.02).
#' @param nPerms Number of permutations (default 1000).
#' @param n_threads Number of threads for parallel C++ computation (default 1).
#'
#' @return An S3 object of class "SpatialEffect" containing:
#' \describe{
#'   \item{AMR_est}{Data frame of distance and AME estimates}
#'   \item{Per.CI}{Permutation confidence intervals (if per.se=1)}
#'   \item{Conley.SE}{Conley standard errors (if conley.se=1)}
#'   \item{Conley.CI}{Conley confidence intervals (if conley.se=1)}
#'   \item{AMR_est_smoothed}{Smoothed estimates (if smooth=1)}
#'   \item{smoothed.Conley.SE}{Smoothed Conley SEs (if smooth=1 and smooth.conley.se=1)}
#'   \item{smoothed.Conley.CI}{Smoothed Conley CIs}
#'   \item{smoothed.Conley.CB}{Uniform confidence bands (if conf.band=1)}
#'   \item{Parameters}{List of stored parameters for downstream use}
#' }
#'
#' @export
SpatialEffect <- function(ras = NULL, Ydata = NULL, outcome = NULL,
                           x_coord_Y = NULL, y_coord_Y = NULL,
                           ras_Z = NULL, Zdata = NULL,
                           x_coord_Z = NULL, y_coord_Z = NULL,
                           treatment, covs = NULL, prob_treatment = NULL,
                           dVec, dist.metric = "Euclidean",
                           cType = "edge", numpts = NULL, evalpts = 1,
                           only.unique = 0,
                           smooth = 0, bw = NULL, bw_debias = NULL,
                           bias_correction = TRUE, smooth.conley.se = 1,
                           conf.band = 0,
                           per.se = 1, blockvar = NULL, clustvar = NULL,
                           conley.se = 1, kernel = "uni", cutoff = 0,
                           alpha = 0.05, edf = FALSE,
                           m = 2, lambda = 0.02, nPerms = 1000,
                           n_threads = 1L) {

  # ---- Input validation ----
  if (is.null(Zdata) || is.null(treatment))
    stop("No treatment specified")

  nz <- nrow(Zdata)
  Zup <- Zdata[[treatment]]

  if (!dist.metric %in% c("Euclidean", "Geodesic"))
    stop("Unrecognized distance metric: use 'Euclidean' or 'Geodesic'")
  if (!cType %in% c("edge", "disk", "donut"))
    stop("Unrecognized circle type: use 'edge', 'disk', or 'donut'")
  if (is.null(dVec) || length(dVec) == 0)
    stop("No distance values provided")

  # ---- Handle polygon interventions (sf-based) ----
  if (!is.null(ras_Z)) {
    message("Polygon intervention")
    if (inherits(ras_Z, "SpatRaster")) {
      ras_Z <- sf::st_as_sf(terra::as.polygons(ras_Z))
    } else if (inherits(ras_Z, "SpatVector")) {
      ras_Z <- sf::st_as_sf(ras_Z)
    } else if (inherits(ras_Z, "SpatialPolygonsDataFrame")) {
      ras_Z <- sf::st_as_sf(ras_Z)
    } else if (!inherits(ras_Z, "sf")) {
      stop("ras_Z must be an sf, SpatRaster, SpatVector, or SpatialPolygonsDataFrame")
    }
    centroids <- sf::st_coordinates(sf::st_centroid(ras_Z))
    Zdata[["x_coord_Z"]] <- centroids[, 1]
    Zdata[["y_coord_Z"]] <- centroids[, 2]
    x_coord_Z <- "x_coord_Z"
    y_coord_Z <- "y_coord_Z"
  } else {
    message("Point intervention")
    if (is.null(x_coord_Z) || is.null(y_coord_Z))
      stop("Coordinates of intervention nodes are not detected")
  }

  # ---- Propensity score weights ----
  if (is.null(prob_treatment)) {
    pweight <- rep(1, nz)
  } else {
    p <- Zdata[[prob_treatment]]
    pweight <- Zup / p + (1 - Zup) / (1 - p)
  }

  # ---- Covariates ----
  if (!is.null(covs)) {
    covs <- as.matrix(covs)
    n_covs <- ncol(covs)
    if (nrow(covs) == nz) {
      covs <- matrix(rep(t(covs), length(dVec)), ncol = n_covs, byrow = TRUE)
    }
    if (!is.null(blockvar)) {
      blocks <- model.matrix(~ blockvar - 1)
      blocks <- matrix(rep(t(blocks), length(dVec)), ncol = ncol(blocks), byrow = TRUE)
      covs <- cbind(covs, blocks)
      n_covs <- ncol(covs)
    }
    colnames(covs) <- paste0("X_", seq_len(n_covs))
  }

  # ---- Distance unit for disk/donut ----
  dist_unit <- min(abs(diff(dVec)), na.rm = TRUE)

  # ---- Outcome surface ----
  if (is.null(ras)) {
    if (is.null(outcome)) stop("Outcome variable is not specified")
    if (is.null(Ydata) && !is.null(Zdata[[outcome]])) {
      ras <- fields::Krig(cbind(Zdata[[x_coord_Z]], Zdata[[y_coord_Z]]),
                          Zdata[[outcome]], m = m, lambda = lambda)
    } else if (!is.null(Ydata) && !is.null(Ydata[[outcome]])) {
      ras <- fields::Krig(cbind(Ydata[[x_coord_Y]], Ydata[[y_coord_Y]]),
                          Ydata[[outcome]], m = m, lambda = lambda)
    } else {
      stop("Outcome variable is not specified")
    }
    Yobs <- NULL
    gridRes <- 1
    dtype <- "krig"
  } else {
    # Convert raster package objects to terra
    if (inherits(ras, "RasterLayer")) {
      ras <- terra::rast(ras)
    }
    gridRes <- terra::res(ras)[1]
    Yobs <- as.numeric(terra::values(ras)[, 1])
    dtype <- "raster"
  }

  # ---- Compute circle averages for each distance ----
  Sdata.list <- list()
  cov.list <- list()
  AMR_est <- data.frame(d = dVec, taud_est = NA_real_)
  Ybards <- Ybard_sums <- Ybard_lens <- matrix(NA_real_, nz, length(dVec))

  message("Computing circle averages for each distance value...")
  for (d in seq_along(dVec)) {
    dUp <- dVec[d]

    rim <- RimAvg(ras, Yobs, ras_Z, nz, Zdata, x_coord_Z, y_coord_Z,
                  dUp, numpts, gridRes, evalpts, only.unique, dtype,
                  dist.metric, cType, dist_unit, n_threads)

    Ybard <- rim$Ybard
    Ybard_sum <- rim$Ybard_sum
    Ybard_len <- rim$Ybard_len

    if (cType == "donut" && d > 1) {
      Ybard_sums[, d] <- Ybard_sum
      Ybard_lens[, d] <- Ybard_len
      Ybards[, d] <- (Ybard_sums[, d] - Ybard_sums[, d - 1]) /
                      (Ybard_lens[, d] - Ybard_lens[, d - 1])
    } else {
      Ybards[, d] <- Ybard
      Ybard_sums[, d] <- Ybard_sum
      Ybard_lens[, d] <- Ybard_len
    }

    # OLS regression of circle average on treatment
    if (is.null(covs)) {
      one_data <- cbind(outcome = Ybards[, d], treatment = Zup,
                        dVec = rep(dUp, nz), pweight = pweight)
      Sdata.list[[d]] <- one_data
      ols_fit <- lm(outcome ~ treatment, data = as.data.frame(one_data), weights = pweight)
      AMR_est$taud_est[d] <- coef(ols_fit)[2]
    } else {
      one_covs <- as.matrix(covs[((d - 1) * nz + 1):(d * nz), , drop = FALSE])
      for (cc in seq_len(ncol(one_covs))) {
        one_covs[, cc] <- one_covs[, cc] - mean(one_covs[, cc])
      }
      one_data <- cbind(outcome = Ybards[, d], treatment = Zup,
                        dVec = rep(dUp, nz), one_covs, pweight = pweight)
      Sdata.list[[d]] <- one_data
      cov_names <- colnames(one_covs)
      ols_formula <- as.formula(
        paste0("outcome ~ ", paste0("treatment * ", cov_names, collapse = " + "))
      )
      ols_fit <- lm(ols_formula, data = as.data.frame(one_data), weights = pweight)
      AMR_est$taud_est[d] <- coef(ols_fit)[2]
      cov.list[[d]] <- model.matrix(ols_fit)
    }
    message(sprintf("  d = %g (%d/%d)", dUp, d, length(dVec)))
  }

  Sdata <- as.data.frame(do.call(rbind, Sdata.list))
  names(Sdata)[1:4] <- c("outcome", "treatment", "dVec", "pweight")

  # ---- Results ----
  result.list <- list()
  result.list[["AMR_est"]] <- AMR_est
  c_n <- NULL

  # ---- Permutation inference ----
  if (per.se == 1) {
    permMat <- .gen_perms(Zup, blockvar = blockvar, clustvar = clustvar,
                          maxiter = nPerms)
    VCE_per <- matrix(NA_real_, nrow = ncol(permMat), ncol = length(dVec))
    for (i in seq_len(ncol(permMat))) {
      z_perm <- permMat[, i]
      VCE_per[i, ] <- as.numeric(z_perm %*% Ybards / sum(z_perm) -
                                  (1 - z_perm) %*% Ybards / sum(1 - z_perm))
    }
    Per.CI <- apply(VCE_per, 2, function(x) quantile(x, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE))
    result.list[["Per.CI"]] <- t(Per.CI)
  }

  # ---- Conley SE ----
  if (conley.se == 1) {
    metric_int <- if (dist.metric == "Euclidean") 1L else 2L
    x_coord <- as.numeric(Zdata[[x_coord_Z]])
    y_coord <- as.numeric(Zdata[[y_coord_Z]])

    if (metric_int == 1L) {
      dist_mat <- DistMatEuclidean(x_coord, y_coord, n_threads)
    } else {
      # For geodesic, use sf::st_distance
      pts <- sf::st_as_sf(data.frame(x = x_coord, y = y_coord),
                          coords = c("x", "y"), crs = 4326)
      dist_mat <- as.matrix(sf::st_distance(pts))
    }

    k <- .parse_kernel(kernel)
    Conley.SE <- rep(NA_real_, length(dVec))
    Conley.CI <- matrix(NA_real_, length(dVec), 2)

    X_mat <- diag(sqrt(pweight)) %*% cbind(rep(1, nz), Zup)
    W_meat <- X_mat
    XX_mat_inv <- solve(crossprod(X_mat))
    z_vec <- c(0, 1) %*% XX_mat_inv %*% t(X_mat)
    M_mat <- diag(1, nz) - X_mat %*% XX_mat_inv %*% t(X_mat)
    dofs <- rep(NA_real_, length(dVec))

    for (d in seq_along(dVec)) {
      mu_d <- diag(sqrt(pweight)) %*% Ybards[, d]
      mu_d[is.nan(mu_d)] <- 0

      if (!is.null(covs)) {
        X_mat <- cov.list[[d]]
        W_meat <- X_mat
        XX_mat_inv <- solve(crossprod(X_mat))
        z_coef <- rep(0, ncol(X_mat)); z_coef[2] <- 1
        z_vec <- z_coef %*% XX_mat_inv %*% t(X_mat)
        M_mat <- diag(1, nz) - X_mat %*% XX_mat_inv %*% t(X_mat)
      }

      beta_d <- solve(crossprod(X_mat)) %*% (t(X_mat) %*% mu_d)
      res_d <- mu_d - X_mat %*% beta_d

      c_val <- if (cutoff > 0) cutoff + d else cutoff
      conley_result <- ConleySE(as.numeric(res_d), W_meat, dist_mat,
                                 XX_mat_inv, c_val, k, 0L, 0L, n_threads)
      VCE_d <- t(XX_mat_inv) %*% conley_result$VCE_meat %*% XX_mat_inv
      dist_kernel <- conley_result$Dist_kernel

      if (edf) {
        dof <- nz / ((nz - 2) * XX_mat_inv[2, 2]) *
          sum(diag(M_mat %*% ((t(z_vec) %*% z_vec) * dist_kernel) %*% M_mat))
        Conley.SE[d] <- sqrt(VCE_d[2, 2] / dof)
        dofs[d] <- dof
      } else {
        Conley.SE[d] <- sqrt(VCE_d[2, 2])
      }
    }

    Conley.CI <- cbind(AMR_est$taud_est - qnorm(1 - alpha / 2) * Conley.SE,
                       AMR_est$taud_est + qnorm(1 - alpha / 2) * Conley.SE)
    result.list[["Conley.SE"]] <- Conley.SE
    result.list[["Conley.CI"]] <- Conley.CI
  }

  # ---- Smoothing ----
  if (smooth == 1) {
    needs_cv <- is.null(bw) || (bias_correction && is.null(bw_debias))
    if (needs_cv) {
      bw_result <- CrossValidation(Sdata, outcome = "outcome", treatment = "treatment",
                                    dVec, grid = NULL, nfold = 5, block_cv = TRUE,
                                    parallel = TRUE, metric = "MSPE", kernel = kernel,
                                    bias_correction = bias_correction,
                                    Zdata = Zdata, x_coord_Z = x_coord_Z,
                                    y_coord_Z = y_coord_Z, n_threads = n_threads)
      if (is.null(bw)) bw <- bw_result$opt.bw
      if (bias_correction && is.null(bw_debias)) bw_debias <- bw_result$opt.bw.debias
    }
    wls_results <- LocalReg(dVec, Sdata, bw, bw_debias, Zup, xevals = NULL,
                             smooth.conley.se = smooth.conley.se, kernel = kernel,
                             cutoff = cutoff, dist = if (exists("dist_mat")) dist_mat else NULL,
                             dist.metric = dist.metric, bias_correction = bias_correction,
                             Zdata = Zdata, x_coord_Z = x_coord_Z, y_coord_Z = y_coord_Z,
                             n_threads = n_threads)
    result.list[["AMR_est_smoothed"]] <- wls_results$coefs
    if (smooth.conley.se == 1) {
      if (edf) {
        result.list[["smoothed.Conley.SE"]] <- wls_results$ses / sqrt(dofs)
      } else {
        result.list[["smoothed.Conley.SE"]] <- wls_results$ses
      }
      result.list[["smoothed.Conley.CI"]] <- cbind(
        wls_results$coefs[, 2] - qnorm(1 - alpha / 2) * result.list[["smoothed.Conley.SE"]],
        wls_results$coefs[, 2] + qnorm(1 - alpha / 2) * result.list[["smoothed.Conley.SE"]]
      )
      if (conf.band == 1) {
        dVec_range <- max(dVec) - min(dVec)
        a_n <- 2 * log(dVec_range / bw) + 2 * log(0.5^0.5 / (2 * pi))
        c_n <- sqrt(a_n - 2 * log(log((1 - alpha)^(-0.5))))
        result.list[["smoothed.Conley.CB"]] <- cbind(
          wls_results$coefs[, 2] - c_n * result.list[["smoothed.Conley.SE"]],
          wls_results$coefs[, 2] + c_n * result.list[["smoothed.Conley.SE"]]
        )
      }
    }
  }

  # ---- Store parameters for downstream use ----
  result.list[["Parameters"]] <- list(
    dVec = dVec, Zup = Zup, Ybards = Ybards, Sdata = Sdata,
    blockvar = blockvar, clustvar = clustvar, nPerms = nPerms,
    conley.se = conley.se, per.se = per.se, smooth = smooth,
    smooth.conley.se = smooth.conley.se, bw = bw, bw_debias = bw_debias,
    c_n = c_n, n_threads = n_threads
  )
  result.list$call <- match.call()
  class(result.list) <- "SpatialEffect"
  return(result.list)
}


# ---- Internal helpers ----

#' Parse kernel string to integer code
#' @keywords internal
.parse_kernel <- function(kernel) {
  kernel <- tolower(kernel)
  if (kernel %in% c("uni", "uniform")) return(1L)
  if (kernel %in% c("tri", "triangular")) return(2L)
  if (kernel %in% c("epa", "epanechnikov")) return(3L)
  stop("Unrecognized kernel: use 'uni', 'tri', or 'epa'")
}

#' Generate permutation matrix
#'
#' Simple Bernoulli-style permutation matrix generator.
#' Falls back to the ri package if available and blockvar/clustvar are specified.
#' @keywords internal
.gen_perms <- function(Zup, blockvar = NULL, clustvar = NULL, maxiter = 1000) {
  nz <- length(Zup)
  n_treat <- sum(Zup)

  if (!is.null(blockvar) || !is.null(clustvar)) {
    # Try to use ri package for blocked/clustered permutations
    if (requireNamespace("ri", quietly = TRUE)) {
      return(ri::genperms(Zup, blockvar = blockvar, clustvar = clustvar,
                          maxiter = maxiter))
    }
    warning("Package 'ri' not available; ignoring blockvar/clustvar for permutations")
  }

  # Simple permutation: reshuffle treatment labels
  perm_mat <- matrix(NA_real_, nrow = nz, ncol = maxiter)
  for (i in seq_len(maxiter)) {
    perm_mat[, i] <- sample(Zup, nz, replace = FALSE)
  }
  return(perm_mat)
}
