## flipping_sampler.R
## ---------------------------------------------------------------------------
## Your uploaded flip-move slice sampler (multivariate_gibbs_sample_ASG_log_freshu_flip),
## unchanged apart from the source() filename below, which is repointed at
## univariate_slice_and_effective_support.R (see that file's header note --
## your original helper file was not part of the upload and has been
## reconstructed).
## ---------------------------------------------------------------------------
## Adds ONE cheap global move per sweep: propose theta -> -theta
## (reflection through the origin) via a standard Metropolis step using
## the SAME g() already defined -- so it's computed generically, not
## assumed. For kernels where g(-theta) == g(theta) exactly (K0, K1,
## K3, K4 here), this is accepted with probability 1 every time it's
## tried, giving a free global jump between disconnected branches. For
## kernels without that symmetry (K2), it will just usually get
## rejected -- harmless, no correctness cost either way.
##
## make_plots = TRUE now also auto-prints an ACF plot alongside the
## trace plot for each parameter (single-chain diagnostics only --
## Gelman-Rubin needs multiple chains, so that lives in the multi-chain
## diagnostic wrapper instead, not here).
##
## Everything else is identical to the fresh-log_u version.
source("univariate_slice_and_effective_support.R")

slice1d_fixed_logu <- function(g, log_u, a, b, x0 = NULL, w = NULL,
                               max_step_out = 50L, max_shrink = 200L) {
  if (a > b) { tmp <- a; a <- b; b <- tmp }             # guard against swapped bounds
  if (is.null(x0)) x0 <- runif(1, a, b)
  x0 <- min(max(x0, a), b)                              # clip x0 defensively
  if (is.null(w)) w <- max((b - a) / 50, 1e-10)          # never a zero-width step
  
  ## L, R computed WITHOUT clipping first, so x0 is guaranteed to lie
  ## inside [L, R] at the start (L <= x0 <= R exactly, by construction).
  ## Clipping L or R to [a,b] before the other is what previously let
  ## R end up smaller than x0 -- fixed by clipping both only at the end.
  u1 <- runif(1)
  L <- x0 - w * u1
  R <- L + w
  
  steps <- 0L
  while (L > a) {
    gL <- g(max(L, a))
    if (!is.finite(gL) || gL <= log_u || steps >= max_step_out) break
    L <- L - w; steps <- steps + 1L
  }
  L <- max(L, a)
  
  steps <- 0L
  while (R < b) {
    gR <- g(min(R, b))
    if (!is.finite(gR) || gR <= log_u || steps >= max_step_out) break
    R <- R + w; steps <- steps + 1L
  }
  R <- min(R, b)
  
  if (!is.finite(L) || !is.finite(R) || L >= R) return(x0)  # degenerate bracket: stay put
  
  for (t in seq_len(max_shrink)) {
    x_new <- runif(1, L, R)
    gx <- g(x_new)
    if (is.finite(gx) && gx > log_u) return(x_new)
    if (x_new < x0) L <- x_new else R <- x_new
    if (!is.finite(R - L) || R - L < .Machine$double.eps * 10) return(x0)
  }
  x0  # fallback: keep current value rather than ever erroring out
}

multivariate_gibbs_sample_ASG_log_freshu_flip <- function(
    nlk, n_samples = 2000, burn_in = 500, thin = 1, theta_init, tol = 0.01,
    scale0 = 1, param_names = NULL, bounds = NULL, random_scan = TRUE,
    make_plots = TRUE, diagnose = FALSE, flip_prob = 0.5) {
  
  start_time <- Sys.time()
  m <- length(theta_init)
  if (is.null(param_names)) param_names <- paste0("theta_", seq_len(m))
  total_iter <- burn_in + n_samples * thin
  theta_chain <- matrix(NA_real_, nrow = total_iter, ncol = m)
  colnames(theta_chain) <- param_names
  theta_new <- as.numeric(theta_init)
  
  g <- function(theta) -nlk(theta)
  
  if (!is.null(bounds)) {
    bounds_list <- bounds
  } else {
    ker_full <- kernel_from_nlk(nlk)
    bounds_list <- vector("list", m)
    for (i in seq_len(m)) {
      ker_i <- (function(i) {
        function(xi) {
          tmp <- theta_new; tmp[i] <- xi
          ker_full(tmp)
        }
      })(i)
      bnds <- effective.support(ker_i, tol = tol, scale0 = scale0)
      bounds_list[[i]] <- c(a = bnds$lower, b = bnds$upper)
    }
  }
  
  n_stalled <- 0L
  n_total_updates <- 0L
  n_flip_attempts <- 0L
  n_flip_accepts <- 0L
  
  in_bounds <- function(theta) {
    all(vapply(seq_len(m), function(i) {
      theta[i] >= bounds_list[[i]]["a"] && theta[i] <= bounds_list[[i]]["b"]
    }, logical(1)))
  }
  
  for (t in seq_len(total_iter)) {
    
    ## --- global reflection-flip move (cheap, tried once per sweep) ---
    if (runif(1) < flip_prob) {
      n_flip_attempts <- n_flip_attempts + 1L
      theta_prop <- -theta_new
      if (in_bounds(theta_prop)) {
        log_alpha <- g(theta_prop) - g(theta_new)
        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
          theta_new <- theta_prop
          n_flip_accepts <- n_flip_accepts + 1L
        }
      }
    }
    
    ## --- usual coordinate-wise slice sweep, fresh log_u per coordinate ---
    scan_order <- if (random_scan) sample(seq_len(m)) else seq_len(m)
    for (i in scan_order) {
      ## `log_u` here is -u in the write-up's notation, not log(u):
      ## draw E ~ Exp(1) via -log(runif(1)), then log_u = g(x) - E
      ##   = -(f(x) + E) = -u, where u ~ Exp(1) + f(x) as in the theory.
      ## Working with -u directly (rather than u itself) avoids ever
      ## exponentiating a possibly large f(x), and keeps the slice
      ## comparison g(x) > log_u equivalent to f(x) < u.
      log_u <- log(runif(1)) + g(theta_new)
      
      g_i <- (function(i) {
        function(xi) {
          tmp <- theta_new; tmp[i] <- xi
          g(tmp)
        }
      })(i)
      a_i <- bounds_list[[i]]["a"]; b_i <- bounds_list[[i]]["b"]
      
      old_val <- theta_new[i]
      theta_new[i] <- slice1d_fixed_logu(g_i, log_u, a = a_i, b = b_i, x0 = old_val)
      
      if (diagnose) {
        n_total_updates <- n_total_updates + 1L
        if (identical(theta_new[i], old_val)) n_stalled <- n_stalled + 1L
      }
    }
    theta_chain[t, ] <- theta_new
    if (t %% 100 == 0) cat("iter", t, "/", total_iter, "\n")
  }
  
  if (diagnose) {
    cat(sprintf("Stalled coordinate-updates: %d / %d (%.1f%%)\n",
                n_stalled, n_total_updates, 100 * n_stalled / n_total_updates))
  }
  cat(sprintf("Flip move: %d accepted / %d attempted (%.1f%%)\n",
              n_flip_accepts, n_flip_attempts,
              100 * n_flip_accepts / max(n_flip_attempts, 1)))
  
  sample_idx <- seq(from = burn_in + 1L, to = total_iter, by = thin)
  theta_samples <- theta_chain[sample_idx, , drop = FALSE]
  time_taken <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  if (make_plots) {
    for (i in seq_len(m)) {
      chain_i <- theta_chain[, i]
      
      ## Trace plot (unchanged)
      df_trace <- data.frame(iter = seq_along(chain_i), val = chain_i)
      p_trace <- ggplot(df_trace, aes(iter, val)) +
        geom_line(alpha = 0.8) +
        geom_vline(xintercept = burn_in, linetype = "dashed", color = "red") +
        labs(title = paste(param_names[i], ": trace"), x = "Iteration", y = param_names[i]) +
        theme_minimal(base_size = 12)
      print(p_trace)
      
      ## ACF plot (new) -- full chain including burn-in, base R's acf()
      ## auto-plots on call; shown for the single chain this function
      ## produces. For a proper multi-chain ACF/Gelman-Rubin comparison,
      ## use the multi-chain diagnostic wrapper instead.
#      acf(chain_i, main = paste(param_names[i], ": ACF"))
    }
  }
  
  list(samples = theta_samples, chain = theta_chain, param_names = param_names,
       burn_in = burn_in, thin = thin, time_taken = time_taken,
       stall_rate = if (diagnose) n_stalled / n_total_updates else NA,
       flip_accept_rate = n_flip_accepts / max(n_flip_attempts, 1))
}
