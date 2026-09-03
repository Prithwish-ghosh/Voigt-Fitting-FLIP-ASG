## univariate_slice_and_effective_support.R
## ---------------------------------------------------------------------------
## NOTE ON PROVENANCE: flipping.R (the file you uploaded) does
##   source("univariate slice and effective support.R")
## but that file itself was not part of your upload, so it could not be
## ported -- only reconstructed. What follows is a from-scratch, documented
## implementation of the two functions flipping.R actually calls from it
## (kernel_from_nlk() and effective.support()), built purely from how they
## are used in flipping.R. If you still have your original file, use that
## instead -- it will be functionally interchangeable as far as the sampler
## is concerned, but may differ in numerical details (step-out schedule,
## refinement tolerance, etc.).
##
## Importantly: for the rbvfit application in this delivery, bounds are
## supplied directly and explicitly to the sampler (physically motivated
## lb/ub from set_bounds(), matching the Python fit), so
## multivariate_gibbs_sample_ASG_log_freshu_flip() never actually calls
## effective.support() in that workflow -- it only needs this file to exist
## so that its `source()` call at the top does not error out. The
## reconstruction below is provided so the whole pipeline is self-contained
## and also usable in the no-bounds-supplied mode.
## ---------------------------------------------------------------------------

kernel_from_nlk <- function(nlk) {
  ## The "kernel" is just the unnormalised log-density implied by nlk
  ## (negative log kernel): g(theta) = -nlk(theta). This mirrors exactly
  ## how `g` is defined inline in multivariate_gibbs_sample_ASG_log_freshu_flip().
  function(theta) -nlk(theta)
}

effective.support <- function(g, tol = 0.01, scale0 = 1, x0 = 0, max_expand = 100L) {
  ## Finds an interval [lower, upper] outside of which the 1-D log-kernel
  ## g(x) has dropped to a fraction `tol` of its peak value (i.e. where
  ## exp(g(x) - g(mode)) < tol), for use as slice-sampler bounds when the
  ## caller does not supply `bounds` directly.
  ##
  ## Algorithm: (1) crude coordinate hill-climb from x0 to locate a local
  ## mode of g with geometrically shrinking step size (robust to g being
  ## non-smooth, at the cost of only finding *a* local mode, not
  ## necessarily global); (2) march outward from the mode in growing
  ## (x1.6) steps until g drops below the tol threshold or a candidate is
  ## non-finite; (3) bisect between the last accepted point and the first
  ## rejected point to refine the crossing. Repeated for both directions.

  ## ---- 1. locate a mode near x0 ----
  gx0 <- g(x0)
  if (!is.finite(gx0)) {
    found <- FALSE
    for (step in scale0 * c(1, -1, 2, -2, 5, -5, 10, -10, 20, -20, 50, -50)) {
      cand <- x0 + step
      if (is.finite(g(cand))) { x0 <- cand; gx0 <- g(cand); found <- TRUE; break }
    }
    if (!found) stop("effective.support: could not find a finite starting point near x0")
  }

  best_x <- x0; best_g <- gx0; step <- max(scale0, 1e-8)
  for (rep in seq_len(80L)) {
    improved <- FALSE
    for (cand in c(best_x - step, best_x + step)) {
      gc <- g(cand)
      if (is.finite(gc) && gc > best_g) { best_g <- gc; best_x <- cand; improved <- TRUE }
    }
    if (!improved) step <- step / 2
    if (step < scale0 * 1e-6) break
  }
  mode_x <- best_x; mode_g <- best_g

  ## ---- 2 & 3. expand outward then bisect to the tol-crossing, each side ----
  log_tol <- log(tol)  # negative; crossing is where mode_g - g(x) == -log_tol
  expand_and_refine <- function(sign) {
    inside <- mode_x; outside <- NA_real_; w <- scale0
    for (i in seq_len(max_expand)) {
      cand <- inside + sign * w
      gc <- g(cand)
      if (!is.finite(gc) || (mode_g - gc) > -log_tol) { outside <- cand; break }
      inside <- cand
      w <- w * 1.6
    }
    if (is.na(outside)) outside <- inside + sign * w  # gave up expanding; use last step
    for (i in seq_len(50L)) {
      mid <- (inside + outside) / 2
      gm <- g(mid)
      if (is.finite(gm) && (mode_g - gm) <= -log_tol) inside <- mid else outside <- mid
    }
    outside
  }

  upper <- expand_and_refine(1)
  lower <- expand_and_refine(-1)
  list(lower = lower, upper = upper, mode = mode_x)
}
