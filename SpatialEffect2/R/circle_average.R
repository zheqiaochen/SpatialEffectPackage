#' Generate sampling points on a circle (edge), disk, or donut
#'
#' @param center Numeric vector c(x, y)
#' @param radius Distance radius
#' @param numpts Number of points around the circle
#' @param dist.metric "Euclidean" or "Geodesic"
#' @param cType "edge", "disk", or "donut"
#' @param dist_unit Step size for disk/donut radii
#' @return Matrix of (x, y) sampling coordinates
#' @keywords internal
GenSamplingPoints <- function(center, radius, numpts, dist.metric = "Euclidean",
                              cType = "edge", dist_unit = NULL) {
  if (cType == "edge") {
    if (dist.metric == "Euclidean") {
      angles <- ((1:numpts) / numpts) * 2 * pi
      return(cbind(center[1] + cos(angles) * radius,
                   center[2] + sin(angles) * radius))
    } else if (dist.metric == "Geodesic") {
      if (!requireNamespace("geosphere", quietly = TRUE))
        stop("Package 'geosphere' is required for geodesic distance")
      return(geosphere::destPoint(center, b = ((1:numpts) / numpts) * 360, d = radius))
    }
  } else if (cType %in% c("disk", "donut")) {
    if (is.null(dist_unit)) dist_unit <- radius
    radii <- seq(0, radius, by = dist_unit)
    if (dist.metric == "Euclidean") {
      angles <- ((1:numpts) / numpts) * 2 * pi
      all_points <- lapply(radii, function(r) {
        cbind(center[1] + cos(angles) * r,
              center[2] + sin(angles) * r)
      })
      return(do.call(rbind, all_points))
    } else if (dist.metric == "Geodesic") {
      if (!requireNamespace("geosphere", quietly = TRUE))
        stop("Package 'geosphere' is required for geodesic distance")
      all_points <- lapply(radii, function(r) {
        geosphere::destPoint(center, b = ((1:numpts) / numpts) * 360, d = r)
      })
      return(do.call(rbind, all_points))
    }
  }
}

#' Generate sampling points by buffering a polygon (sf-based)
#'
#' @param geom An sf geometry (single polygon)
#' @param ras A SpatRaster object
#' @param radius Buffer distance
#' @param cType "edge", "disk", or "donut"
#' @return Matrix of (x, y) sampling coordinates on valid raster cells
#' @keywords internal
GenSamplingPointsBuffer <- function(geom, ras, radius, cType) {
  buffered <- sf::st_buffer(geom, dist = radius)
  if (cType == "edge") {
    # Extract boundary of the buffer
    boundary <- sf::st_cast(buffered, "MULTILINESTRING")
    buffer_rast <- terra::rasterize(terra::vect(boundary), ras, touches = TRUE)
  } else {
    buffer_rast <- terra::rasterize(terra::vect(buffered), ras)
  }
  vals <- terra::values(buffer_rast)[, 1]
  all_coords <- terra::xyFromCell(ras, 1:terra::ncell(ras))
  valid <- !is.na(vals)
  return(all_coords[valid, , drop = FALSE])
}

#' Compute circle averages for all intervention nodes at a given distance
#'
#' @param ras SpatRaster or Krig object
#' @param Yobs Raster values vector (NULL for krig)
#' @param ras_Z sf object for polygon interventions (NULL for point)
#' @param nz Number of intervention nodes
#' @param Zdata Data frame of intervention data
#' @param x_coord_Z Name of x-coordinate column
#' @param y_coord_Z Name of y-coordinate column
#' @param dUp Distance value
#' @param numpts Number of sampling points
#' @param gridRes Raster resolution
#' @param evalpts Multiplier for numpts
#' @param only.unique Deduplicate grid cells?
#' @param dtype "raster" or "krig"
#' @param dist.metric "Euclidean" or "Geodesic"
#' @param cType "edge", "disk", or "donut"
#' @param dist_unit Step size for disk/donut
#' @param n_threads Number of threads for C++ computation
#' @return List(Ybard, Ybard_sum, Ybard_len)
#' @keywords internal
RimAvg <- function(ras, Yobs = NULL, ras_Z = NULL, nz, Zdata, x_coord_Z, y_coord_Z,
                   dUp, numpts = NULL, gridRes = NULL, evalpts = 10,
                   only.unique = 0, dtype = "raster", dist.metric = "Euclidean",
                   cType = "edge", dist_unit = NULL, n_threads = 1L) {

  if (is.null(numpts)) {
    if (is.null(gridRes)) gridRes <- 1
    numpts <- ceiling(ceiling(dUp / (gridRes / sqrt(2)) + 1) * pi) * evalpts
  }

  # --- Kriging path ---
  if (dtype == "krig") {
    Ybard <- Ybard_sum <- Ybard_len <- rep(NA_real_, nz)
    # Generate all sampling coords at once for batch prediction
    centers <- cbind(Zdata[[x_coord_Z]], Zdata[[y_coord_Z]])
    offset_template <- GenSamplingPoints(c(0, 0), dUp, numpts, dist.metric, cType, dist_unit)
    n_offsets <- nrow(offset_template)

    all_coords <- cbind(
      rep(centers[, 1], each = n_offsets) + rep(offset_template[, 1], times = nz),
      rep(centers[, 2], each = n_offsets) + rep(offset_template[, 2], times = nz)
    )
    group_ids <- rep(seq_len(nz), each = n_offsets) - 1L  # 0-based for C++

    # Batch kriging prediction
    preds <- predict(ras, all_coords)
    result <- CircleAvgFromPredictions(as.numeric(preds), as.integer(group_ids), as.integer(nz))
    return(result)
  }

  # --- Raster path ---
  # Fast C++ path: Euclidean + point interventions
  if (is.null(ras_Z) && dist.metric == "Euclidean") {
    centers_x <- as.numeric(Zdata[[x_coord_Z]])
    centers_y <- as.numeric(Zdata[[y_coord_Z]])

    offset_template <- GenSamplingPoints(c(0, 0), dUp, numpts, dist.metric, cType, dist_unit)

    # Get raster properties from terra
    raster_vals <- as.numeric(terra::values(ras)[, 1])
    ext <- terra::ext(ras)
    xmin <- ext[1]
    ymax <- ext[4]
    xres <- terra::res(ras)[1]
    yres <- terra::res(ras)[2]
    nc <- terra::ncol(ras)
    nr <- terra::nrow(ras)

    result <- CircleAvgRaster(
      raster_vals, nc, nr, xmin, ymax, xres, yres,
      centers_x, centers_y,
      offset_template[, 1], offset_template[, 2],
      as.integer(only.unique), as.integer(n_threads)
    )
    return(result)
  }

  # --- Polygon intervention path (sf-based, parallelized in R) ---
  if (!is.null(ras_Z)) {
    one_rim_stat <- function(i) {
      geom_i <- sf::st_geometry(ras_Z)[i]
      coords <- GenSamplingPointsBuffer(geom_i, ras, dUp, cType)
      cell_idx <- terra::cellFromXY(ras, coords)
      if (only.unique == 1) cell_idx <- unique(na.omit(cell_idx))
      else cell_idx <- na.omit(cell_idx)
      grids <- Yobs[cell_idx]
      grids <- grids[!is.na(grids)]
      c(mean(grids), sum(grids), length(grids))
    }
    rim_stats <- future.apply::future_lapply(seq_len(nz), one_rim_stat,
                                              future.seed = TRUE)
    rim_mat <- do.call(rbind, rim_stats)
    return(list(
      Ybard = rim_mat[, 1],
      Ybard_sum = rim_mat[, 2],
      Ybard_len = rim_mat[, 3]
    ))
  }

  # --- Geodesic point intervention fallback (loop) ---
  Ybard <- Ybard_sum <- Ybard_len <- rep(NA_real_, nz)
  for (i in seq_len(nz)) {
    coords <- GenSamplingPoints(
      center = c(Zdata[[x_coord_Z]][i], Zdata[[y_coord_Z]][i]),
      radius = dUp, numpts = numpts, dist.metric = dist.metric,
      cType = cType, dist_unit = dist_unit
    )
    cell_idx <- terra::cellFromXY(ras, coords)
    if (only.unique == 1) cell_idx <- unique(na.omit(cell_idx))
    else cell_idx <- na.omit(cell_idx)
    grids <- Yobs[cell_idx]
    grids <- grids[!is.na(grids)]
    Ybard[i] <- if (length(grids) > 0) mean(grids) else NaN
    Ybard_sum[i] <- sum(grids)
    Ybard_len[i] <- length(grids)
  }
  list(Ybard = Ybard, Ybard_sum = Ybard_sum, Ybard_len = Ybard_len)
}
