\name{SpatialEffect2-workflow}
\alias{SpatialEffect2-workflow}
\title{SpatialEffect2 Workflow Guide: Detailed Inputs, Parameters, and Practice}
\description{
Practical operations manual for \code{SpatialEffect2}. This page documents
input templates, argument meanings, expected output formats, and step-by-step
analysis practice.
}

\section{Core functions}{
\itemize{
  \item \code{\link{SpatialEffect}}: estimate distance-indexed AME curves.
  \item \code{\link{SpatialEffectTest}}: randomization test on a distance range.
  \item S3 methods for class \code{"SpatialEffect"}: \code{print()}, \code{summary()}, \code{plot()}.
}
}

\section{Input templates (recommended formats)}{
\strong{Template A: point interventions}
\preformatted{
# one row per intervention node
Zdata <- data.frame(
  id = 1:nz,
  x = ...,   # numeric x coordinate
  y = ...,   # numeric y coordinate
  treat = ... # 0/1 treatment
)

# required call pieces
x_coord_Z = "x"
y_coord_Z = "y"
treatment = "treat"
}

\strong{Template B: polygon interventions}
\preformatted{
# ras_Z must be sf / SpatRaster / SpatVector
ras_Z <- intervention_sf
Zdata <- data.frame(treat = c(1,0,1,...))
# node coordinates are derived internally from centroids
}

\strong{Template C: kriging mode (no raster)}
\preformatted{
ras = NULL
Ydata <- data.frame(x = ..., y = ..., outcome = ...)
outcome = "outcome"
x_coord_Y = "x"
y_coord_Y = "y"
}
}

\section{Argument reference for \code{SpatialEffect}}{
\strong{Always required}
\describe{
  \item{\code{Zdata}}{Data frame of intervention nodes; one row per node.}
  \item{\code{treatment}}{Character name of binary treatment column in \code{Zdata}.}
  \item{\code{dVec}}{Numeric distance grid where AME is evaluated.}
}

\strong{Outcome and geometry inputs}
\describe{
  \item{\code{ras}}{\code{terra::SpatRaster} outcome surface. If \code{NULL}, kriging interpolation is used.}
  \item{\code{Ydata}}{Outcome point data frame used in kriging mode.}
  \item{\code{outcome}}{Outcome column name in \code{Ydata} (or in \code{Zdata} if \code{Ydata = NULL}).}
  \item{\code{x_coord_Y}, \code{y_coord_Y}}{Coordinate column names in \code{Ydata}.}
  \item{\code{ras_Z}}{Polygon interventions as \code{sf}, \code{SpatRaster}, or \code{SpatVector}; \code{NULL} means point intervention mode.}
  \item{\code{x_coord_Z}, \code{y_coord_Z}}{Coordinate column names in \code{Zdata}; required in point mode.}
}

\strong{Distance definition and sampling}
\describe{
  \item{\code{cType}}{Distance set type: \code{"edge"}, \code{"disk"}, or \code{"donut"}.}
  \item{\code{dist.metric}}{\code{"Euclidean"} or \code{"Geodesic"}.}
  \item{\code{numpts}}{Points sampled per circle boundary; auto-selected if \code{NULL}.}
  \item{\code{evalpts}}{Multiplier for auto-selected \code{numpts}.}
  \item{\code{only.unique}}{If \code{1}, deduplicate sampled raster cells per node.}
}

\strong{Design adjustment}
\describe{
  \item{\code{covs}}{Optional covariate matrix. Accepts either \code{nz x p} or \code{(nz * length(dVec)) x p}.}
  \item{\code{prob_treatment}}{Propensity-score column name for IPW adjustment.}
  \item{\code{blockvar}}{Optional block variable for permutation design.}
  \item{\code{clustvar}}{Optional cluster variable for permutation design.}
}

\strong{Inference controls}
\describe{
  \item{\code{per.se}}{If \code{1}, compute permutation interval \code{Per.CI}.}
  \item{\code{conley.se}}{If \code{1}, compute Conley HAC \code{Conley.SE} and \code{Conley.CI}.}
  \item{\code{kernel}}{Conley kernel: \code{"uni"}/\code{"uniform"}, \code{"tri"}/\code{"triangular"}, \code{"epa"}/\code{"epanechnikov"}.}
  \item{\code{cutoff}}{Spatial kernel bandwidth for Conley HAC.}
  \item{\code{alpha}}{Significance level used for intervals.}
  \item{\code{edf}}{Use effective-DOF adjustment for Conley SE.}
  \item{\code{nPerms}}{Number of permutations.}
}

\strong{Smoothing controls}
\describe{
  \item{\code{smooth}}{If \code{1}, run local-polynomial smoothing on AME curve.}
  \item{\code{bw}}{Smoothing bandwidth; if \code{NULL}, selected by CV.}
  \item{\code{bw_debias}}{Bias-correction bandwidth; if \code{NULL}, selected by CV when needed.}
  \item{\code{bias_correction}}{If \code{TRUE}, compute debiased smoother.}
  \item{\code{smooth.conley.se}}{If \code{1}, compute uncertainty for smoothed curve.}
  \item{\code{conf.band}}{If \code{1}, compute uniform confidence band for smoothed curve.}
}

\strong{Kriging and compute controls}
\describe{
  \item{\code{m}, \code{lambda}}{Kriging controls used when \code{ras = NULL}.}
  \item{\code{n_threads}}{Thread count for C++ distance and HAC kernels.}
}
}

\section{Argument reference for \code{SpatialEffectTest}}{
\describe{
  \item{\code{result.list}}{A \code{"SpatialEffect"} object returned by \code{SpatialEffect}.}
  \item{\code{test.range}}{Length-2 numeric vector \code{c(d_min, d_max)} for the tested distance interval.}
  \item{\code{smooth}}{If \code{1}, test uses smoothed curve; otherwise unsmoothed AME.}
  \item{\code{alpha}}{Significance level for permutation interval.}
}
}

\section{Output format (example structure)}{
\strong{\code{SpatialEffect(...)} returns class \code{"SpatialEffect"}.}
\preformatted{
$AME_est               # data.frame with columns d, taud_est
$Per.CI                # optional matrix [length(dVec), 2]
$Conley.SE             # optional numeric vector [length(dVec)]
$Conley.CI             # optional matrix [length(dVec), 2]
$AME_est_smoothed      # optional matrix [length(dVec), 2]
$smoothed.Conley.SE    # optional numeric vector
$smoothed.Conley.CI    # optional matrix [length(dVec), 2]
$smoothed.Conley.CB    # optional matrix [length(dVec), 2]
$Parameters            # list used by methods/tests
$call                  # stored call
}

\strong{\code{SpatialEffectTest(...)} returns:}
\preformatted{
$test.stat             # scalar statistic on selected range
$test.CI               # length-2 permutation interval
}
}

\section{Practical workflow checklist}{
\enumerate{
  \item Check CRS alignment and spatial overlap first.
  \item Start with conservative baseline: \code{smooth = 0}, \code{per.se = 1}, \code{conley.se = 1}.
  \item Use \code{dVec} step size similar to raster resolution for \code{"donut"} to reduce empty-ring \code{NA}.
  \item Add smoothing only after unsmoothed curve is interpretable.
  \item Test policy-relevant ranges with \code{SpatialEffectTest}.
  \item Run sensitivity checks over \code{cType}, \code{cutoff}, \code{kernel}, and \code{dVec} granularity.
}
}

\section{Common issues and fixes}{
\strong{Issue: many \code{NA} values in AME}
\itemize{
  \item CRS mismatch between outcome and intervention objects.
  \item No or weak spatial overlap.
  \item Misalignment between \code{Zdata} rows and intervention geometries.
  \item \code{"donut"} rings too narrow relative to raster resolution.
}

\strong{Issue: unstable smoothing / singular design}
\itemize{
  \item Increase \code{bw} and \code{bw_debias}.
  \item Reduce collinearity in \code{covs}.
  \item Use less dense \code{dVec}.
}

\strong{Issue: runtime too slow}
\itemize{
  \item Increase \code{n_threads}.
  \item Decrease \code{nPerms}.
  \item Reduce sampling load via \code{numpts} and \code{evalpts}.
}
}

\examples{
\dontrun{
library(SpatialEffect2)
library(terra)

set.seed(1)

# 1) outcome raster
ras <- rast(nrows = 60, ncols = 60, xmin = 0, xmax = 60, ymin = 0, ymax = 60)
values(ras) <- rnorm(ncell(ras))

# 2) intervention nodes
nz <- 100
Zdata <- data.frame(
  x = runif(nz, 2, 58),
  y = runif(nz, 2, 58),
  treat = rbinom(nz, 1, 0.5)
)

# 3) estimate AME
res <- SpatialEffect(
  ras = ras,
  Zdata = Zdata,
  x_coord_Z = "x",
  y_coord_Z = "y",
  treatment = "treat",
  dVec = seq(0, 10, by = 1),
  cType = "donut",
  dist.metric = "Euclidean",
  smooth = 1,
  per.se = 1,
  conley.se = 1,
  nPerms = 500,
  n_threads = 2L
)

summary(res, dVec.range = c(1, 5))
plot(res, ci.type = "both")
SpatialEffectTest(res, test.range = c(1, 5), smooth = 0)
}
}

\author{
Package authors: Ye Wang, Cyrus Samii, Haoge Chang, and P. M. Aronow.

This workflow page is maintained as a practical user reference.
}
