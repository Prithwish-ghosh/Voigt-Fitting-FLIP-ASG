## run_rbvfit_flip_mcmc.R
## ---------------------------------------------------------------------------
## Driver script: rbvfit's Voigt-profile absorption-line model (ported from
## the Python rbvfit-master package you uploaded) fitted with YOUR flip-move
## slice sampler (flipping.R) instead of rbvfit's own emcee/zeus backend.
##
## Run this file directly (Rscript run_rbvfit_flip_mcmc.R, or source it in
## an R session) to see:
##   PART A: a self-contained synthetic-data recovery test. This is exactly
##           the check used while building this port -- it simulates a
##           two-component SiII absorption profile, adds noise, and confirms
##           the sampler recovers the injected N/b/v values. Numbers will
##           differ slightly run to run (new random noise draw + random
##           seed you can change) but should always recover the truth to
##           within a couple of posterior standard deviations.
##   PART B: a template showing exactly where to plug in your own spectrum
##           (wavelength / flux / error vectors) and absorber configuration,
##           mirroring the rbvfit-master/src/rbvfit/examples/*.py tutorials.
##
## Required files (same folder): voigt_lines.R, atom_full.dat,
## rbvfit_config.R, voigt_model.R, rbvfit_likelihood.R,
## univariate_slice_and_effective_support.R, flipping_sampler.R.
## ---------------------------------------------------------------------------

source("voigt_lines.R")
source("rbvfit_config.R")
source("voigt_model.R")
source("rbvfit_likelihood.R")
source("univariate_slice_and_effective_support.R")
source("flipping_sampler.R")

linelist <- load_atom_linelist("atom_full.dat")

## =============================================================================
## PART A -- synthetic recovery test
## =============================================================================

cat("\n================ PART A: synthetic recovery test ================\n\n")

set.seed(20260902)

## ---- 1. Absorber configuration: one system, one ion, two velocity clouds ----
## Mirrors:  config = FitConfiguration(); config.add_system(z=0, ion='SiII',
##           transitions=[1190.5, 1193.5], components=2)
config <- new_fit_configuration()
config <- add_system(config, z = 0.0, ion = "SiII",
                      transitions = c(1190.5, 1193.5), components = 2)
line_arrays <- build_line_arrays(config, linelist)

## theta layout is c(N_1, N_2, b_1, b_2, v_1, v_2), same order rbvfit uses
## (all N's, then all b's, then all v's, in ion/component order).
true_theta <- c(13.6, 13.3, 18, 12, -25, 40)

## ---- 2. Simulate an observed spectrum ----
wave <- c(seq(1189.5, 1191.5, by = 0.02), seq(1192.3, 1194.5, by = 0.02))
FWHM_pixels <- 2.7                                  # instrumental resolution, in pixels
kernel <- gaussian_kernel1d(FWHM_pixels / 2.355)

true_flux <- evaluate_voigt_model(true_theta, wave, line_arrays,
                                   kernel = kernel, voigt_method = "wofz")
sigma <- 0.02
noisy_flux <- true_flux + rnorm(length(wave), sd = sigma)
error <- rep(sigma, length(wave))

## ---- 3. Bounds (ion-aware, matching vfit_mcmc.py :: set_bounds) ----
bnds <- set_bounds(nguess = c(13.6, 13.3), bguess = c(18, 12), vguess = c(-25, 40),
                    ions = c("SiII", "SiII"))
lb <- bnds$lb; ub <- bnds$ub
bounds_list <- lapply(seq_along(lb), function(i) c(a = lb[i], b = ub[i]))

## ---- 4. Likelihood / posterior (mirrors vfit_mcmc.py :: vfit) ----
instrument_data <- prepare_instrument_data(list(
  SIM = list(line_arrays = line_arrays, wave = wave, flux = noisy_flux, error = error,
             kernel = kernel, voigt_method = "wofz")
))
nlk <- make_nlk(instrument_data, lb, ub)  # nlk = -log posterior; g = -nlk inside the sampler

## ---- 5. Starting point (perturbed away from truth, like a real fit would be) ----
theta_init <- true_theta + c(0.3, -0.3, 3, -2, 5, -5)
theta_init <- pmin(pmax(theta_init, lb + 1e-6), ub - 1e-6)

## ---- 6. Run YOUR flip-move slice sampler ----
## Bounds are supplied explicitly (bounds_list, from the physically-motivated
## ion table), so effective.support() is never invoked here -- the sampler
## goes straight to slice sampling within [lb, ub].
res <- multivariate_gibbs_sample_ASG_log_freshu_flip(
  nlk = nlk, n_samples = 3000, burn_in = 1000, thin = 1,
  theta_init = theta_init, bounds = bounds_list,
  param_names = c("N1", "N2", "b1", "b2", "v1", "v2"),
  make_plots = FALSE,   # set TRUE if you have ggplot2 installed
  diagnose = TRUE, flip_prob = 0.3
)

post_mean <- colMeans(res$samples)
post_sd   <- apply(res$samples, 2, sd)

cat("\nTrue theta:      ", round(true_theta, 3), "\n")
cat("Posterior mean:  ", round(post_mean, 3), "\n")
cat("Posterior sd:    ", round(post_sd, 3), "\n")
cat("(z-scores)       ", round((post_mean - true_theta) / post_sd, 2), "\n\n")

## NOTE on the flip move for this application: the flip proposes
## theta -> -theta for the WHOLE vector at once. Column densities N
## (log10 cm^-2, ~11-22) and Doppler parameters b (km/s, ~5-150) are
## always positive and far from zero, so reflecting them through the
## origin almost always lands outside [lb, ub] and gets rejected --
## exactly the "harmless, no correctness cost" case the sampler's own
## header comment anticipates (see flip_accept_rate in the output above,
## typically ~0%). The flip move only pays off for parameters with a
## genuine +/- mirror-image degeneracy passing through zero (velocity v
## alone might qualify in some setups, but here it's flipped jointly with
## N and b, which blocks it). If you want the flip move to do real work
## for this kind of fit, it would need to be restricted to the velocity
## block only (theta[v_idx] -> -theta[v_idx] with N, b left alone) --
## that is a change to flipping_sampler.R's proposal, not something this
## script does for you, since it changes the sampler's semantics.

## ---- 7. Simple diagnostics (base R, no extra packages required) ----
if (interactive() || Sys.getenv("RBVFIT_SAVE_PLOTS") == "1") {
  png_path <- "rbvfit_flip_diagnostics.png"
  png(png_path, width = 1000, height = 700)
  par(mfrow = c(2, 3))
  for (i in seq_along(res$param_names)) {
    plot(res$chain[, i], type = "l", main = res$param_names[i],
         xlab = "iteration", ylab = "value")
    abline(v = res$burn_in, col = "red", lty = 2)
    abline(h = true_theta[i], col = "blue", lty = 3)
  }
  dev.off()
  cat("Saved trace plots to", png_path, "\n")
}

## =============================================================================
## PART B -- template for your own data
## =============================================================================
## Uncomment and edit. Everything below mirrors
## rbvfit-master/src/rbvfit/examples/example_voigt_fitter.py line for line.
##
## ---- 1. Load your spectrum -------------------------------------------------
## # dat   <- read.csv("my_spectrum.csv")   # needs columns wave, flux, error
## # wave  <- dat$wave
## # flux  <- dat$flux
## # error <- dat$error
## # (If your data is FITS, the `FITSio` package can read it: FITSio::readFITS())
##
## ---- 2. Define the absorber(s) ---------------------------------------------
## # config <- new_fit_configuration()
## # config <- add_system(config, z = 0.0,       ion = "SiII", transitions = c(1190.5, 1193.5), components = 1)
## # config <- add_system(config, z = 0.162005,  ion = "HI",   transitions = c(1025.7),          components = 1)
## # line_arrays <- build_line_arrays(config, linelist)
##
## ---- 3. Initial guesses + bounds --------------------------------------------
## # nguess <- c(14.2, 14.5); bguess <- c(40, 30); vguess <- c(0, 0)
## # theta_init <- c(nguess, bguess, vguess)
## # bnds <- set_bounds(nguess, bguess, vguess,
## #                     Nlow = c(12,12), Nhi = c(17,17),
## #                     blow = c(5,5),   bhi = c(100,100),
## #                     vlow = c(-300,-300), vhi = c(300,300))
## # lb <- bnds$lb; ub <- bnds$ub
## # bounds_list <- lapply(seq_along(lb), function(i) c(a = lb[i], b = ub[i]))
##
## ---- 4. Instrument resolution + likelihood ----------------------------------
## # FWHM_vel_kms <- 18.0
## # FWHM_pixels  <- mean_fwhm_pixels(FWHM_vel_kms, wave)
## # kernel <- gaussian_kernel1d(FWHM_pixels / 2.355)
## # instrument_data <- prepare_instrument_data(list(
## #   COS = list(line_arrays = line_arrays, wave = wave, flux = flux, error = error,
## #              kernel = kernel, voigt_method = "wofz")
## # ))
## # nlk <- make_nlk(instrument_data, lb, ub)
##
## ---- 5. Run the sampler ------------------------------------------------------
## # res <- multivariate_gibbs_sample_ASG_log_freshu_flip(
## #   nlk = nlk, n_samples = 5000, burn_in = 2000, thin = 2,
## #   theta_init = theta_init, bounds = bounds_list,
## #   param_names = c("N_SiII","N_HI","b_SiII","b_HI","v_SiII","v_HI"),
## #   make_plots = TRUE, diagnose = TRUE, flip_prob = 0.3
## # )
## # colMeans(res$samples); apply(res$samples, 2, sd)
## =============================================================================

setwd(old_wd)

source("rbvfit_plots.R")
best_theta <- colMeans(res$samples)

rbvfit_corner_plot(res$samples, param_names = res$param_names, truths = true_theta,
                   save_path = "corner.png")
rbvfit_correlation_plot(res$samples, param_names = res$param_names,
                        save_path = "correlation.png")
rbvfit_velocity_plot(best_theta, config, line_arrays, wave, noisy_flux, error,
                     kernel = kernel, voigt_method = "wofz",
                     velocity_range = c(-150, 150), show_components = TRUE,
                     save_path = "velocity.png")
