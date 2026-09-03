## rbvfit_likelihood.R
## ---------------------------------------------------------------------------
## R port of the statistical model in vfit_mcmc.py :: vfit (lnprior, lnlike,
## lnprob). Multi-instrument aware: `instrument_data` is a named list, each
## element list(line_arrays=, wave=, flux=, error=, kernel=, voigt_method=),
## exactly mirroring the Python `instrument_data` dict of
## {'model':, 'wave':, 'flux':, 'error':}. Requires voigt_model.R to be
## sourced first (for evaluate_voigt_model()).
## ---------------------------------------------------------------------------

make_lnprior <- function(lb, ub) {
  ## Uniform prior within [lb, ub], -Inf outside -- identical to
  ## vfit_mcmc.py :: vfit.lnprior().
  function(theta) {
    if (any(theta < lb) || any(theta > ub)) return(-Inf)
    0.0
  }
}

make_lnlike <- function(instrument_data) {
  ## Gaussian log-likelihood summed over instruments -- identical in form
  ## to vfit_mcmc.py :: vfit.lnlike() (the "fast" variant that is actually
  ## used for MCMC in the Python code; the commented-out "more accurate but
  ## slower" alternative differs only by an additive constant and was not
  ## used by the original for sampling).
  function(theta) {
    lnlike_total <- 0.0
    for (nm in names(instrument_data)) {
      d <- instrument_data[[nm]]
      model_flux <- tryCatch(
        evaluate_voigt_model(theta, d$wave, d$line_arrays, kernel = d$kernel,
                              voigt_method = d$voigt_method %||% "wofz"),
        error = function(e) NULL
      )
      if (is.null(model_flux) || any(!is.finite(model_flux))) return(-Inf)
      inv_sigma2 <- d$inv_sigma2
      log_inv_sigma2 <- d$log_inv_sigma2
      lnlike_instrument <- -0.5 * sum((d$flux - model_flux)^2 * inv_sigma2 - log_inv_sigma2)
      lnlike_total <- lnlike_total + lnlike_instrument
    }
    lnlike_total
  }
}
`%||%` <- function(a, b) if (is.null(a)) b else a

prepare_instrument_data <- function(raw_instrument_data) {
  ## raw_instrument_data: named list of list(line_arrays=, wave=, flux=,
  ## error=, kernel=NULL, voigt_method="wofz"). Precomputes inv_sigma2 /
  ## log_inv_sigma2 once, mirroring vfit.py's _compile_models().
  out <- list()
  for (nm in names(raw_instrument_data)) {
    d <- raw_instrument_data[[nm]]
    out[[nm]] <- list(
      line_arrays = d$line_arrays,
      wave = d$wave, flux = d$flux, error = d$error,
      kernel = d$kernel, voigt_method = d$voigt_method %||% "wofz",
      inv_sigma2 = 1.0 / d$error^2,
      log_inv_sigma2 = log(1.0 / d$error^2)
    )
  }
  out
}

make_lnprob <- function(instrument_data, lb, ub) {
  ## Log posterior = lnprior + lnlike, mirroring vfit_mcmc.py :: vfit.lnprob().
  lnprior <- make_lnprior(lb, ub)
  lnlike  <- make_lnlike(instrument_data)
  function(theta) {
    lp <- lnprior(theta)
    if (!is.finite(lp)) return(-Inf)
    lp + lnlike(theta)
  }
}

make_nlk <- function(instrument_data, lb, ub) {
  ## Negative log kernel for the flip-move slice sampler in flipping_sampler.R
  ## (which expects g(theta) = -nlk(theta) to be an unnormalised log-density).
  lnprob <- make_lnprob(instrument_data, lb, ub)
  function(theta) {
    lp <- lnprob(theta)
    if (!is.finite(lp)) return(Inf)
    -lp
  }
}
