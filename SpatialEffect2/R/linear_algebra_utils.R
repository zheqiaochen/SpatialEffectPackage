#' Safe inverse of cross-product matrix
#'
#' Uses a tiny ridge penalty fallback when \code{crossprod(X)} is singular or
#' nearly singular. If ridge escalation still fails, falls back to an SVD-based
#' pseudo-inverse.
#'
#' @param X Numeric matrix.
#' @param ridge_init Initial ridge value.
#' @param max_tries Maximum number of ridge escalation attempts.
#' @return A list with \code{inv} and \code{ridge}.
#' @keywords internal
.safe_solve_crossprod <- function(X, ridge_init = 1e-8, max_tries = 8L) {
  XtX <- crossprod(X)
  p <- ncol(XtX)
  ridge <- 0

  for (i in seq_len(max_tries + 1L)) {
    mat <- if (ridge > 0) XtX + diag(ridge, p) else XtX
    inv <- tryCatch(solve(mat), error = function(e) NULL)
    if (!is.null(inv) && all(is.finite(inv))) {
      return(list(inv = inv, ridge = ridge))
    }
    ridge <- if (ridge == 0) ridge_init else ridge * 10
  }

  sv <- svd(XtX)
  tol <- max(dim(XtX)) * max(sv$d) * .Machine$double.eps
  d_inv <- ifelse(sv$d > tol, 1 / sv$d, 0)
  inv <- sv$v %*% (d_inv * t(sv$u))
  list(inv = inv, ridge = NA_real_)
}
