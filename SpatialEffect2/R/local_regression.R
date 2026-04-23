#' Local polynomial regression for smoothing AME estimates
#'
#' Implements the local linear regression smoother with optional bias correction
#' via higher-order polynomial terms, following Hainmueller et al. (2018).
#'
#' @param dVec Distance values
#' @param Sdata Stacked data frame with outcome, treatment, dVec, pweight columns
#' @param bw Bandwidth
#' @param bw_debias Bandwidth for bias correction polynomial
#' @param Zup Treatment vector
#' @param xevals Evaluation points (defaults to dVec)
#' @param smooth.conley.se Compute Conley SE for smoothed estimates
#' @param kernel Kernel type string
#' @param cutoff Conley SE bandwidth
#' @param dist Pre-computed distance matrix
#' @param dist.metric Distance metric
#' @param bias_correction Use bias correction
#' @param Zdata Intervention data frame
#' @param x_coord_Z,y_coord_Z Coordinate column names
#' @param n_threads Number of threads
#' @return List with coefs, coefs_all, ses, coefs_debias
#' @keywords internal
LocalReg <- function(dVec, Sdata, bw, bw_debias = NULL, Zup = NULL,
                     xevals = NULL, smooth.conley.se = 1,
                     kernel = "uni", cutoff = 0, dist = NULL,
                     dist.metric = "Euclidean", bias_correction = TRUE,
                     Zdata = NULL, x_coord_Z = NULL, y_coord_Z = NULL,
                     n_threads = 1L) {

  if (is.null(xevals)) xevals <- dVec

  all_names <- names(Sdata)
  wls_formula <- as.formula("outcome ~ treatment + dVec_demeaned + interaction")
  wls_formula_debias <- as.formula(
    "outcome ~ treatment + dVec_demeaned + interaction + dVec_demeaned_2 + interaction_2"
  )
  num_coefs <- 4
  num_debias_coefs <- 6

  if (length(all_names) > 4) {
    cov_names <- all_names[5:length(all_names)]
    cov_formula <- paste0("treatment * ", cov_names, collapse = " + ")
    wls_formula <- as.formula(
      paste0("outcome ~ treatment + dVec_demeaned + interaction + ", cov_formula)
    )
    wls_formula_debias <- as.formula(
      paste0("outcome ~ treatment + dVec_demeaned + interaction + dVec_demeaned_2 + interaction_2 + ", cov_formula)
    )
    num_coefs <- num_coefs + 2 * length(cov_names)
    num_debias_coefs <- num_debias_coefs + 2 * length(cov_names)
  }

  wls_coefs <- wls_ses <- rep(NA_real_, length(xevals))
  wls_coefs_all <- matrix(NA_real_, length(xevals), num_coefs)
  wls_coefs_debias_all <- matrix(NA_real_, length(xevals), num_debias_coefs)

  # Pre-compute spatial distance for Conley SE of smoothed estimates
  if (smooth.conley.se == 1) {
    metric_int <- if (dist.metric == "Euclidean") 1L else 2L
    if (is.null(dist)) {
      if (is.null(Zdata) || is.null(x_coord_Z) || is.null(y_coord_Z))
        stop("Zdata and coordinates required when smooth.conley.se = 1")
      x_c <- as.numeric(Zdata[[x_coord_Z]])
      y_c <- as.numeric(Zdata[[y_coord_Z]])
      if (metric_int == 1L) {
        dist <- DistMatEuclidean(x_c, y_c, n_threads)
      } else {
        pts <- sf::st_as_sf(data.frame(x = x_c, y = y_c),
                            coords = c("x", "y"), crs = 4326)
        dist <- as.matrix(sf::st_distance(pts))
      }
    }
    # Instead of building the full (nz*nD)^2 replicated distance matrix,
    # we build the subset dist matrix on-the-fly for each evaluation point.
    # This saves O((nz*nD)^2) memory.
    nz_local <- nrow(dist)
  }

  k <- .parse_kernel(kernel)

  for (i in seq_along(xevals)) {
    xeval <- xevals[i]
    dat1 <- Sdata

    # Remove NaN outcomes upfront to keep all indices aligned
    valid_rows <- !is.nan(dat1$outcome) & !is.na(dat1$outcome)
    dat1 <- dat1[valid_rows, ]

    dat1$dVec_demeaned <- dat1$dVec - xeval
    dat1$interaction <- dat1$treatment * dat1$dVec_demeaned

    if (bw == 0) {
      dat1$w <- as.numeric(dat1$dVec_demeaned == 0) * dat1$pweight
    } else {
      dat1$w <- (1 - abs(dat1$dVec_demeaned / bw)) *
                (abs(dat1$dVec_demeaned / bw) <= 1) * dat1$pweight
    }
    w_index <- (dat1$w > 0)

    if (bias_correction && !is.null(bw_debias)) {
      dat1$w_debias <- (1 - abs(dat1$dVec_demeaned / bw_debias)) *
                       (abs(dat1$dVec_demeaned / bw_debias) <= 1) * dat1$pweight
      w_index <- (dat1$w_debias > 0)
      dat1$dVec_demeaned_2 <- dat1$dVec_demeaned^2
      dat1$interaction_2 <- dat1$interaction^2
      wls_reg_debias <- lm(wls_formula_debias, data = dat1, weights = w_debias)
      wls_coefs_debias_all[i, ] <- coef(wls_reg_debias)
    }

    # Skip if insufficient positive-weight observations or treatment has no variation
    n_pos_w <- sum(w_index)
    if (n_pos_w < 3 || length(unique(dat1$treatment[w_index])) < 2) {
      wls_coefs[i] <- NA_real_
      next
    }

    result_i <- tryCatch({
      wls_reg <- lm(wls_formula, data = dat1, weights = w)
      wls_coef <- coef(wls_reg)
      wls_coef[is.na(wls_coef)] <- 0

      # Use subset of rows with positive weight
      dat1_w <- dat1[w_index, ]
      W_K <- diag(sqrt(dat1_w$w))
      res_vec <- as.numeric(W_K %*% wls_reg$residuals[w_index])
      X_mat <- W_K %*% as.matrix(model.matrix(wls_reg))[w_index, ]
      W_meat_local <- X_mat
      XX_mat_inv_local <- solve(crossprod(X_mat))

      if (bias_correction && !is.null(bw_debias)) {
        W_K_debias <- diag(sqrt(dat1_w$w_debias))
        res_debias <- as.numeric(W_K_debias %*% wls_reg_debias$residuals[w_index])
        X_mat_debias_all <- W_K_debias %*% as.matrix(model.matrix(wls_reg_debias))[w_index, ]
        X_mat_debias <- X_mat_debias_all[, 5:6, drop = FALSE]
        nu <- matrix(0, 2, num_debias_coefs)
        nu[1, 5] <- 1; nu[2, 6] <- 1
        W_meat_local <- t(t(X_mat) - 0.5 * t(X_mat) %*% X_mat_debias %*%
                     (nu %*% solve(crossprod(X_mat_debias_all)) %*% t(X_mat_debias_all)))
        coef_debias <- 0.5 * XX_mat_inv_local %*%
          (t(X_mat) %*% X_mat_debias) %*% wls_coefs_debias_all[i, 5:6]
        wls_coef[2] <- wls_coef[2] - coef_debias[2]
        res_vec <- res_debias
      }

      se_val <- NA_real_
      if (smooth.conley.se == 1) {
        c_val <- if (cutoff > 0) cutoff + xeval else cutoff
        orig_idx <- as.integer(rownames(dat1_w))
        node_idx <- ((orig_idx - 1) %% nz_local) + 1
        dist_sub <- dist[node_idx, node_idx, drop = FALSE]
        d_vals <- dat1_w$dVec
        d_sep <- abs(outer(d_vals, d_vals, "-"))
        dist_sub <- dist_sub + d_sep
        conley_result <- ConleySE(res_vec, W_meat_local, dist_sub,
                                   XX_mat_inv_local, c_val, k, 0L, 0L, n_threads)
        VCE <- t(XX_mat_inv_local) %*% conley_result$VCE_meat %*% XX_mat_inv_local
        se_val <- sqrt(max(0, VCE[2, 2]))
      }

      list(coef = wls_coef, coef2 = wls_coef[2], se = se_val)
    }, error = function(e) NULL)

    if (is.null(result_i)) {
      wls_coefs[i] <- NA_real_
      next
    }
    wls_coefs[i] <- result_i$coef2
    wls_coefs_all[i, ] <- result_i$coef
    wls_ses[i] <- result_i$se
  }

  wls_results <- list(
    coefs = cbind(xevals, wls_coefs),
    coefs_all = wls_coefs_all,
    ses = wls_ses
  )
  if (bias_correction) wls_results$coefs_debias <- wls_coefs_debias_all
  return(wls_results)
}

#' Cross-validation for bandwidth selection
#'
#' Two-stage grid search with K-fold (optionally spatial block) cross-validation.
#'
#' @param Sdata Stacked data frame
#' @param outcome,treatment Column name strings
#' @param dVec Distance vector
#' @param grid Bandwidth grid (auto-generated if NULL)
#' @param nfold Number of folds
#' @param block_cv Use spatial block CV
#' @param parallel Use future.apply parallelism
#' @param metric "MSPE" or "MAPE"
#' @param kernel Kernel type
#' @param bias_correction Use bias correction
#' @param Zdata,x_coord_Z,y_coord_Z Intervention data (required for block_cv)
#' @param n_threads C++ threads
#' @return List with CV.out, opt.bw, opt.bw.debias
#' @keywords internal
CrossValidation <- function(Sdata, outcome, treatment, dVec, grid = NULL,
                             nfold = 5, block_cv = TRUE, parallel = TRUE,
                             metric = "MSPE", kernel = "uni",
                             bias_correction = FALSE,
                             Zdata = NULL, x_coord_Z = NULL, y_coord_Z = NULL,
                             n_threads = 1L) {

  getError <- function(train, test, bw) {
    one_reg <- tryCatch(
      LocalReg(dVec, train, bw = bw, bw_debias = bw, xevals = NULL,
               smooth.conley.se = 0, kernel = kernel, cutoff = 0,
               dist.metric = "Euclidean", bias_correction = bias_correction,
               n_threads = n_threads),
      error = function(e) NULL
    )
    if (is.null(one_reg)) return(c(Inf, Inf, Inf, Inf))

    predict_from_coefs <- function(coef_mat, test_data) {
      esCoef <- function(x) {
        Xnew <- abs(dVec - x)
        d1 <- min(Xnew)
        label1 <- which.min(Xnew)
        Xnew[label1] <- Inf
        d2 <- min(Xnew)
        label2 <- which.min(Xnew)
        if (d1 == 0) coef_mat[label1, ]
        else if (d2 == 0) coef_mat[label2, ]
        else (coef_mat[label1, ] * d2 + coef_mat[label2, ] * d1) / (d1 + d2)
      }
      Knn <- t(sapply(test_data$dVec, esCoef))
      test_Y <- test_data[[outcome]]
      test_X <- cbind(1, test_data[[treatment]])
      test_Y - rowSums(test_X * Knn)
    }

    coef <- one_reg$coefs_all[, 1:2]
    coef[is.na(coef)] <- 0
    error <- predict_from_coefs(coef, test)

    if (bias_correction) {
      coef_debias <- one_reg$coefs_debias[, 1:2]
      coef_debias[is.na(coef_debias)] <- 0
      error_debias <- predict_from_coefs(coef_debias, test)
      c(mean(abs(error)), mean(error^2), mean(abs(error_debias)), mean(error_debias^2))
    } else {
      c(mean(abs(error)), mean(error^2), 0, 0)
    }
  }

  n <- nrow(Sdata)
  nz <- n / length(dVec)
  range_dVec <- max(dVec) - min(dVec)
  user_grid <- !is.null(grid)
  if (!user_grid)
    grid <- exp(seq(log(range_dVec / 10), log(range_dVec), length.out = 15))

  # Generate folds
  if (block_cv) {
    if (is.null(Zdata) || is.null(x_coord_Z) || is.null(y_coord_Z))
      stop("Zdata and coordinates required for block CV")
    fold <- rep(0L, nz)
    x_vec <- as.numeric(Zdata[[x_coord_Z]])
    y_vec <- as.numeric(Zdata[[y_coord_Z]])
    nfold_sqr <- round(sqrt(nfold))
    x_q <- quantile(x_vec, seq(0, 1, 1 / nfold_sqr))
    y_q <- quantile(y_vec, seq(0, 1, 1 / nfold_sqr))
    for (pp in 2:length(x_q)) {
      for (qq in 2:length(y_q)) {
        block_ind <- (pp - 2) * (length(x_q) - 1) + qq - 1
        cond <- x_vec >= x_q[pp - 1] & x_vec <= x_q[pp] &
                y_vec >= y_q[qq - 1] & y_vec <= y_q[qq]
        fold[cond] <- block_ind
      }
    }
    fold <- rep(sample(fold, nz, replace = FALSE), length(dVec))
    nfold <- nfold_sqr^2
  } else {
    fold <- rep(sample((0:(nz - 1)) %% nfold + 1, nz, replace = FALSE), length(dVec))
  }

  cv_one <- function(bw_val) {
    mse <- matrix(NA_real_, nfold, 4)
    for (j in seq_len(nfold)) {
      testid <- which(fold == j)
      mse[j, ] <- getError(Sdata[-testid, ], Sdata[testid, ], bw = bw_val)
    }
    c(bw_val, colMeans(mse))
  }

  runGrid <- function(g) {
    if (parallel) {
      do.call(rbind, future.apply::future_lapply(g, cv_one, future.seed = TRUE))
    } else {
      do.call(rbind, lapply(g, cv_one))
    }
  }

  message("Cross-validating bandwidth...")
  Error <- runGrid(grid)

  if (!user_grid) {
    # Stage 2: refine around optimum
    col_idx <- if (metric == "MAPE") 2 else 3
    best_idx <- which.min(Error[, col_idx])
    lo <- max(1, best_idx - 1)
    hi <- min(length(grid), best_idx + 1)
    grid2 <- exp(seq(log(grid[lo]), log(grid[hi]), length.out = 15))
    grid2 <- grid2[!grid2 %in% grid]
    if (length(grid2) > 0) {
      Error2 <- runGrid(grid2)
      Error <- rbind(Error, Error2)
      grid <- c(grid, grid2)
      ord <- order(grid)
      grid <- grid[ord]
      Error <- Error[ord, , drop = FALSE]
    }
  }

  colnames(Error) <- c("bandwidth", "MAPE", "MSPE", "MAPE_debias", "MSPE_debias")
  col_idx <- if (metric == "MAPE") 2 else 3

  best_idx <- which.min(Error[, col_idx])
  opt.bw <- if (length(best_idx) > 0) grid[best_idx] else grid[ceiling(length(grid) / 2)]

  if (bias_correction) {
    best_debias_idx <- which.min(Error[, col_idx + 2])
    opt.bw.debias <- if (length(best_debias_idx) > 0) grid[best_debias_idx] else opt.bw
  } else {
    opt.bw.debias <- NULL
  }

  message(sprintf("Bandwidth = %.3f", opt.bw))
  if (bias_correction && !is.null(opt.bw.debias))
    message(sprintf("Bandwidth for bias correction = %.3f", opt.bw.debias))

  list(CV.out = round(Error, 3), opt.bw = opt.bw, opt.bw.debias = opt.bw.debias)
}
