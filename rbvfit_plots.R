## rbvfit_plots.R
## Base-R equivalents of rbvfit's core/results_plot.py:
##   corner_plot() / correlation_plot() / velocity_plot()
## No extra packages required. Source voigt_lines.R, rbvfit_config.R,
## voigt_model.R first.

## ---- 1. corner_plot(): posterior pairs plot with 1-D marginals ----
rbvfit_corner_plot <- function(samples, param_names = colnames(samples), truths = NULL,
                               save_path = NULL) {
  samples <- as.matrix(samples)
  p <- ncol(samples)
  if (!is.null(save_path)) png(save_path, width = 180 * p, height = 180 * p, res = 150)
  op <- par(mfrow = c(p, p), mar = c(2, 2, 0.5, 0.5), oma = c(3, 3, 3, 1))
  on.exit({ par(op); if (!is.null(save_path)) dev.off() }, add = TRUE)
  
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (i == j) {
        d <- density(samples[, i])
        plot(d, main = "", xlab = "", ylab = "", axes = FALSE)
        box()
        if (!is.null(truths)) abline(v = truths[i], col = "red", lwd = 2)
        qs <- quantile(samples[, i], c(0.16, 0.5, 0.84))
        title(main = sprintf("%s = %.3f (+%.3f/-%.3f)", param_names[i], qs[2],
                             qs[3] - qs[2], qs[2] - qs[1]), cex.main = 0.8)
      } else if (i > j) {
        plot(samples[, j], samples[, i], pch = 16, cex = 0.3,
             col = adjustcolor("steelblue", alpha.f = 0.3), xlab = "", ylab = "")
        if (!is.null(truths)) { abline(v = truths[j], col = "red", lty = 2); abline(h = truths[i], col = "red", lty = 2) }
      } else {
        plot.new()  # blank upper triangle, like corner.corner's default
      }
      if (i == p) mtext(param_names[j], side = 1, line = 2.2, cex = 0.7)
      if (j == 1) mtext(param_names[i], side = 2, line = 2.2, cex = 0.7)
    }
  }
  mtext("Parameter Posterior Distributions", outer = TRUE, cex = 1.1, line = 1)
}

## ---- 2. correlation_plot(): correlation matrix heatmap ----
rbvfit_correlation_plot <- function(samples, param_names = colnames(samples), save_path = NULL) {
  cm <- cor(samples)
  p <- ncol(cm)
  if (!is.null(save_path)) png(save_path, width = 900, height = 800, res = 150)
  op <- par(mar = c(6, 6, 3, 5))
  on.exit({ par(op); if (!is.null(save_path)) dev.off() }, add = TRUE)
  
  image(1:p, 1:p, t(cm)[, p:1], col = colorRampPalette(c("blue", "white", "red"))(100),
        zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "", main = "Parameter Correlation Matrix")
  axis(1, at = 1:p, labels = param_names, las = 2)
  axis(2, at = 1:p, labels = rev(param_names), las = 2)
  for (i in 1:p) for (j in 1:p) {
    text(i, p - j + 1, sprintf("%.2f", cm[j, i]),
         col = if (abs(cm[j, i]) < 0.5) "black" else "white", cex = 0.8)
  }
  box()
}

## ---- 3. velocity_plot(): data + best-fit model in velocity space, per transition ----
list_transitions <- function(config) {
  rows <- list()
  for (system in config$systems) {
    for (ion_group in system$ion_groups) {
      for (wl in ion_group$transitions) {
        rows[[length(rows) + 1]] <- data.frame(
          system_z = system$z, ion_name = ion_group$ion_name,
          rest_wavelength = wl, obs_wavelength = wl * (1 + system$z)
        )
      }
    }
  }
  do.call(rbind, rows)
}

evaluate_voigt_model_components <- function(theta, wavelength, line_arrays, kernel = NULL, voigt_method = "wofz") {
  N_linear <- 10^theta[line_arrays$N_indices]
  b_values <- theta[line_arrays$b_indices]
  v_values <- theta[line_arrays$v_indices]
  c_kms <- 299792.458
  z_total <- line_arrays$z_factor * (1 + v_values / c_kms) - 1
  n_lines <- line_arrays$n_lines; n_wave <- length(wavelength)
  wave_rest <- matrix(rep(wavelength, each = n_lines), nrow = n_lines) /
    (1 + matrix(rep(z_total, n_wave), nrow = n_lines))
  tau_all <- vectorized_voigt_tau(line_arrays$lambda0, line_arrays$gamma, line_arrays$f_osc,
                                  N_linear, b_values, wave_rest, voigt_method)
  flux_total <- exp(-colSums(tau_all))
  if (!is.null(kernel)) flux_total <- convolve1d_nearest(flux_total, kernel)
  components <- lapply(seq_len(n_lines), function(i) {
    fc <- exp(-tau_all[i, ])
    if (!is.null(kernel)) fc <- convolve1d_nearest(fc, kernel)
    fc
  })
  list(flux = flux_total, components = components,
       lambda0 = line_arrays$lambda0, v_value = v_values)
}

rbvfit_velocity_plot <- function(best_theta, config, line_arrays, wave, flux, error,
                                 kernel = NULL, voigt_method = "wofz",
                                 velocity_range = c(-500, 500), y_range = c(0, 1.2),
                                 show_components = TRUE, save_path = NULL) {
  transitions <- list_transitions(config)
  n_trans <- nrow(transitions)
  c_kms <- 299792.458
  
  model_out <- evaluate_voigt_model_components(best_theta, wave, line_arrays, kernel, voigt_method)
  
  if (!is.null(save_path)) png(save_path, width = 900, height = 300 * n_trans, res = 150)
  op <- par(mfrow = c(n_trans, 1), mar = c(4, 4, 2, 1))
  on.exit({ par(op); if (!is.null(save_path)) dev.off() }, add = TRUE)
  
  for (t in seq_len(n_trans)) {
    obs_wl <- transitions$obs_wavelength[t]
    vel <- c_kms * (wave / obs_wl - 1)
    mask <- vel >= velocity_range[1] & vel <= velocity_range[2]
    
    plot(vel[mask], flux[mask], type = "s", col = "black", lwd = 1,
         xlim = velocity_range, ylim = y_range, xlab = "Velocity (km/s)", ylab = "Normalized flux",
         main = sprintf("%s %.1f (z=%.5f)", transitions$ion_name[t], transitions$rest_wavelength[t], transitions$system_z[t]))
    polygon(c(vel[mask], rev(vel[mask])),
            c(flux[mask] - error[mask], rev(flux[mask] + error[mask])),
            col = adjustcolor("gray", alpha.f = 0.3), border = NA)
    lines(vel[mask], model_out$flux[mask], col = "red", lwd = 2)
    abline(h = 1, lty = 3, col = "gray50")
    abline(v = 0, lty = 3, col = "gray50")
    
    if (show_components) {
      comp_colors <- rainbow(length(model_out$lambda0))
      for (i in seq_along(model_out$lambda0)) {
        if (abs(model_out$lambda0[i] - transitions$rest_wavelength[t]) < 0.1) {
          lines(vel[mask], model_out$components[[i]][mask], col = comp_colors[i], lty = 2, lwd = 1)
        }
      }
    }
    legend("bottomright", legend = c("Data", "Best-fit model"), col = c("black", "red"),
           lty = 1, lwd = c(1, 2), bty = "n", cex = 0.7)
  }
}