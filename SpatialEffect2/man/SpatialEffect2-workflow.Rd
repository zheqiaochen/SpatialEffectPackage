\name{SpatialEffect2-workflow}
\alias{SpatialEffect2-workflow}
\title{SpatialEffect2 Analysis Workflow}
\description{
A decision guide for planning a SpatialEffect2 analysis. The tutorial shows how
to run a complete example; this page helps decide what to pass, which estimand to
use, and which checks to run before interpreting results.
}

\section{How to use this page}{
Use this page before or during an empirical analysis. It is not a second
tutorial. It is a checklist for choices that affect the estimand, computation,
and uncertainty.
}

\section{Decision 1: what is the outcome surface?}{
\describe{
  \item{Raster outcome}{Use \code{ras} when the outcome is already a gridded
    surface, such as forest loss, pollution, temperature, or population density.
    This is the most direct workflow.}
  \item{Point outcome}{Use \code{Ydata}, \code{outcome}, \code{x_coord_Y}, and
    \code{y_coord_Y} when outcomes are point measurements. Set \code{ras = NULL}
    so SpatialEffect2 interpolates an outcome surface.}
}

\strong{Check:} the outcome and intervention objects should use compatible
coordinate reference systems. For \code{dist.metric = "Euclidean"}, distances
are interpreted in the units of the projected CRS.
}

\section{Decision 2: what is the intervention geometry?}{
\describe{
  \item{Point interventions}{Use \code{Zdata} with one row per node plus
    coordinate columns supplied through \code{x_coord_Z} and \code{y_coord_Z}.}
  \item{Polygon interventions}{Pass polygons through \code{ras_Z}. Keep one row
    in \code{Zdata} per polygon, in the same order as \code{ras_Z}. Coordinates
    are taken from polygon centroids internally.}
}

\strong{Check:} \code{nrow(Zdata)} must match the number of intervention nodes.
For polygons, row order matters.
}

\section{Decision 3: which distance design answers the question?}{
\describe{
  \item{\code{cType = "edge"}}{Use the outcome boundary at exactly distance
    \code{d}. This is closest to asking, "What happens right at this distance?"}
  \item{\code{cType = "disk"}}{Use the cumulative region from \code{0} through
    \code{d}. This is useful for total neighborhood effects up to a radius.}
  \item{\code{cType = "donut"}}{Use adjacent rings. This is often easiest to
    interpret when asking how effects change as we move outward.}
}

\strong{Check:} for \code{"donut"}, choose a \code{dVec} step that is not much
smaller than the raster resolution. Very thin rings can create empty-ring
\code{NA} values.
}

\section{Decision 4: experimental or observational study?}{
\describe{
  \item{Experimental study}{If treatment probabilities are known by design and
    no weighting is needed, omit \code{prob_treatment}. Use \code{blockvar} or
    \code{clustvar} if randomization was blocked or clustered.}
  \item{Observational study}{Estimate or supply a propensity score in
    \code{Zdata}, then pass its column name through \code{prob_treatment}.
    Inspect overlap before interpreting the weighted estimates.}
}

\strong{Check:} propensity scores near 0 or 1 can dominate weighted estimates.
Plot or summarize the propensity-score distribution.
}

\section{Decision 5: how should uncertainty be shown?}{
\describe{
  \item{\code{conley.se = 1}}{Adds spatial HAC standard errors and
    \code{Conley.CI}. Use \code{cutoff} and \code{kernel} for sensitivity
    checks.}
  \item{\code{per.se = 1}}{Adds permutation intervals \code{Per.CI}. Increase
    \code{nPerms} for final results.}
  \item{\code{smooth = 1}}{Adds a smoothed AME curve. Use it after the
    unsmoothed curve is understandable.}
  \item{\code{smooth.conley.se = 1}}{Adds uncertainty for the smoothed curve.
    This is more computationally demanding.}
}

\strong{Practical default:} start with \code{smooth = 0}, inspect the raw AME
curve, then add smoothing and sensitivity checks.
}

\section{Templates}{
\strong{Point interventions}
\preformatted{
res <- SpatialEffect(
  ras = ras,
  Zdata = Zdata,
  x_coord_Z = "x",
  y_coord_Z = "y",
  treatment = "treat",
  dVec = dVec,
  cType = "donut"
)
}

\strong{Polygon interventions}
\preformatted{
res <- SpatialEffect(
  ras = ras,
  ras_Z = ras_Z,
  Zdata = Zdata,
  treatment = "treat",
  dVec = dVec,
  cType = "donut"
)
}

\strong{Observational study with propensity weights}
\preformatted{
res <- SpatialEffect(
  ras = ras,
  ras_Z = ras_Z,
  Zdata = Zdata,
  treatment = "treat",
  prob_treatment = "p_treat",
  dVec = dVec,
  cType = "donut"
)
}
}

\section{After estimation}{
\preformatted{
summary(res, dVec.range = c(d_min, d_max))
plot(res, ci.type = "conley")
SpatialEffectTest(res, test.range = c(d_min, d_max), smooth = 0)
}

\describe{
  \item{\code{summary()}}{Use for a compact numeric table at a relevant distance
    range.}
  \item{\code{plot()}}{Use for the AME curve and available intervals.}
  \item{\code{SpatialEffectTest()}}{Use for a policy-relevant distance band, such
    as 0--5 km.}
}
}

\section{Troubleshooting checklist}{
\describe{
  \item{Many \code{NA} estimates}{Check CRS alignment, spatial overlap, row
    order between \code{ras_Z} and \code{Zdata}, and whether donut rings are too
    thin relative to raster resolution.}
  \item{Runtime is slow}{Reduce \code{numpts}, reduce \code{nPerms}, start with
    \code{smooth = 0}, or increase \code{n_threads}.}
  \item{Confidence intervals are very wide}{Check treatment balance,
    propensity-score overlap, distance bands with few outcome cells, and
    sensitivity to \code{cutoff}.}
  \item{Results change with \code{dVec}}{Use a distance grid aligned with the
    raster resolution and report sensitivity to reasonable grids.}
}
}

