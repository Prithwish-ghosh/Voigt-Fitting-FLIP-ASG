## voigt_model.R
## ---------------------------------------------------------------------------
## R port of core/voigt_model.py: builds the flat atomic-parameter arrays and
## theta index maps for a config (mirrors VoigtModel._cache_atomic_parameters
## and ._setup_fast_mapping), computes the vectorized multi-line Voigt optical
## depth, and evaluates model flux with an optional Gaussian instrumental
## convolution. Requires voigt_lines.R (rb_setline / voigt_hjerting) to be
## sourced first.
## ---------------------------------------------------------------------------

mean_fwhm_pixels <- function(FWHM_vel_kms, wave_obs_grid) {
  ## Direct port of voigt_model.py :: mean_fwhm_pixels().
  if (any(wave_obs_grid <= 0)) stop("Wavelength grid must be strictly positive.")
  if (length(wave_obs_grid) < 2) stop("Wavelength grid must have at least two points.")
  c_kms <- 299792.458
  delta_lambda <- c(wave_obs_grid[2] - wave_obs_grid[1],
                     (wave_obs_grid[3:length(wave_obs_grid)] - wave_obs_grid[1:(length(wave_obs_grid) - 2)]) / 2,
                     wave_obs_grid[length(wave_obs_grid)] - wave_obs_grid[length(wave_obs_grid) - 1])
  fwhm_lambda <- wave_obs_grid * FWHM_vel_kms / c_kms
  fwhm_pixels <- fwhm_lambda / delta_lambda
  mean(fwhm_pixels)
}

gaussian_kernel1d <- function(stddev) {
  ## Matches astropy.convolution.Gaussian1DKernel's default array (odd
  ## length ~ floor(8*stddev+1), normalised to sum 1). Verified against
  ## astropy numerically to machine precision for stddev in [0.5, 10].
  size <- floor(8 * stddev + 1)
  if (size %% 2 == 0) size <- size + 1
  half <- (size - 1) / 2
  offsets <- (-half):half
  k <- exp(-offsets^2 / (2 * stddev^2))
  k / sum(k)
}

convolve1d_nearest <- function(x, kernel) {
  ## Matches scipy.ndimage.convolve1d(x, kernel, mode='nearest') for a
  ## symmetric kernel: same-length output, edges replicated for padding.
  ## Verified against scipy+astropy numerically to machine precision.
  n <- length(x); ksize <- length(kernel); half <- (ksize - 1L) %/% 2L
  padded <- c(rep(x[1], half), x, rep(x[n], half))
  out <- numeric(n)
  for (i in seq_len(n)) out[i] <- sum(padded[i:(i + ksize - 1L)] * kernel)
  out
}

build_line_arrays <- function(config, linelist) {
  ## Mirrors VoigtModel._cache_atomic_parameters() + ._setup_fast_mapping().
  ## Returns a list with the flat per-line atomic arrays and the N/b/v
  ## index vectors (0-based conceptually, but stored as 1-based R indices
  ## into a theta vector laid out as c(N_1..N_p, b_1..b_p, v_1..v_p)).
  lambda0 <- c(); gamma <- c(); f_osc <- c(); z_factor <- c()
  N_idx <- integer(0)
  global_comp_idx <- 0L

  for (system in config$systems) {
    for (ion_group in system$ion_groups) {
      for (wl in ion_group$transitions) {
        line_info <- rb_setline(wl, "closest", linelist)
        for (comp in seq_len(ion_group$components)) {
          lambda0   <- c(lambda0, line_info$wave)
          gamma     <- c(gamma, line_info$gamma)
          f_osc     <- c(f_osc, line_info$fval)
          z_factor  <- c(z_factor, 1.0 + system$z)
          N_idx     <- c(N_idx, global_comp_idx + comp)   # 1-based within this ion group's block
        }
      }
      global_comp_idx <- global_comp_idx + ion_group$components
    }
  }

  total_components <- global_comp_idx
  list(
    lambda0 = lambda0, gamma = gamma, f_osc = f_osc, z_factor = z_factor,
    N_indices = N_idx,                          # 1-based index into theta[1:total_components]
    b_indices = N_idx + total_components,       # 1-based index into theta[(tc+1):(2tc)]
    v_indices = N_idx + 2L * total_components,  # 1-based index into theta[(2tc+1):(3tc)]
    n_lines = length(lambda0),
    total_components = total_components
  )
}

vectorized_voigt_tau <- function(lambda0, gamma, f_osc, N_linear, b_values, wave_rest, voigt_method = "wofz") {
  ## Direct port of voigt_model.py :: _vectorized_voigt_tau().
  ## wave_rest is an (n_lines x n_wavelengths) matrix (one row per line).
  c_freq <- 2.99792458e18
  atomic_constant <- 4.48898479507e3

  n_lines <- length(lambda0)
  b_f   <- b_values / lambda0 * 1e13
  freq0 <- c_freq / lambda0
  constant <- atomic_constant / (freq0 * b_values)

  a <- gamma / (4 * pi * b_f)

  freq <- c_freq / wave_rest    # matrix, n_lines x n_wave
  x <- sweep(freq, 1, freq0, "-")
  x <- sweep(x, 1, b_f, "/")

  H <- matrix(0, nrow = nrow(wave_rest), ncol = ncol(wave_rest))
  for (i in seq_len(n_lines)) {
    H[i, ] <- voigt_hjerting(x[i, ], a[i], method = voigt_method)
  }

  tau_all <- sweep(H, 1, N_linear * f_osc * constant, "*")
  tau_all
}

evaluate_voigt_model <- function(theta, wavelength, line_arrays, kernel = NULL, voigt_method = "wofz") {
  ## Direct port of voigt_model.py :: _evaluate_compiled_model(), returning
  ## flux only (return_components = FALSE path). `kernel` is either NULL
  ## (no instrumental convolution) or a numeric vector from
  ## gaussian_kernel1d().
  N_linear <- 10^theta[line_arrays$N_indices]
  b_values <- theta[line_arrays$b_indices]
  v_values <- theta[line_arrays$v_indices]

  c_kms <- 299792.458
  z_total <- line_arrays$z_factor * (1 + v_values / c_kms) - 1

  n_lines <- line_arrays$n_lines
  n_wave  <- length(wavelength)
  wave_rest <- matrix(rep(wavelength, each = n_lines), nrow = n_lines) /
    (1 + matrix(rep(z_total, n_wave), nrow = n_lines))

  tau_all <- vectorized_voigt_tau(line_arrays$lambda0, line_arrays$gamma, line_arrays$f_osc,
                                   N_linear, b_values, wave_rest, voigt_method)
  tau_total <- colSums(tau_all)
  flux <- exp(-tau_total)

  if (!is.null(kernel)) flux <- convolve1d_nearest(flux, kernel)
  flux
}
