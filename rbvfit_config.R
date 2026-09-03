## rbvfit_config.R
## ---------------------------------------------------------------------------
## R port of the *parts of* core/fit_configuration.py and vfit_mcmc.py's
## set_bounds() that are actually needed to build a theta vector, its
## N/b/v index maps, and parameter bounds. This intentionally reproduces
## just the data model (systems -> ion groups -> transitions/components),
## not the full validation / GUI-support machinery of the original
## ParameterManager and FitConfiguration classes (mirroring, plotting
## helpers, .rbv save files, etc.), since none of that changes the fitted
## numbers.
##
## Config representation
## ----------------------
## A config is a plain R list of "systems". Each system is
##   list(z = <redshift>, ion_groups = list(ion_group, ion_group, ...))
## and each ion_group is
##   list(ion_name = <string>, transitions = <numeric vector, rest-frame A>,
##        components = <integer>)
##
## Build one with new_fit_configuration() + add_system(), exactly mirroring
## the Python call pattern:
##   config <- new_fit_configuration()
##   config <- add_system(config, z = 0.0,       ion = "SiII", transitions = c(1190.5, 1193.5), components = 1)
##   config <- add_system(config, z = 0.162005,  ion = "HI",   transitions = c(1025.7),          components = 1)
## ---------------------------------------------------------------------------

new_fit_configuration <- function() {
  list(systems = list())
}

add_system <- function(config, z, ion, transitions, components = 1L, merge = FALSE) {
  ## Mirrors FitConfiguration.add_system() -> AbsorptionSystem.add_ion().
  ## Systems are grouped by redshift (tolerance 1e-6); within a system,
  ## ion groups are grouped by ion name. `ion` must be given explicitly
  ## here (the R port does not replicate the Python 'auto' ion-detection
  ## convenience -- pass the ion name you want, e.g. via rb_setline()'s
  ## $name field if you need to look it up from a wavelength).
  z_tol <- 1e-6
  sys_idx <- NULL
  for (i in seq_along(config$systems)) {
    if (abs(config$systems[[i]]$z - z) < z_tol) { sys_idx <- i; break }
  }
  if (is.null(sys_idx)) {
    config$systems[[length(config$systems) + 1L]] <- list(z = z, ion_groups = list())
    sys_idx <- length(config$systems)
  }

  ion_groups <- config$systems[[sys_idx]]$ion_groups
  grp_idx <- NULL
  for (j in seq_along(ion_groups)) {
    if (identical(ion_groups[[j]]$ion_name, ion)) { grp_idx <- j; break }
  }

  if (!is.null(grp_idx)) {
    if (!merge) {
      stop(sprintf("Ion %s already exists in system at z=%g. Use merge=TRUE to add transitions.", ion, z))
    }
    if (ion_groups[[grp_idx]]$components != components) {
      stop(sprintf("Cannot merge ion %s: component mismatch (existing: %d, new: %d)",
                    ion, ion_groups[[grp_idx]]$components, components))
    }
    existing <- ion_groups[[grp_idx]]$transitions
    ion_groups[[grp_idx]]$transitions <- unique(c(existing, transitions))
  } else {
    ion_groups[[length(ion_groups) + 1L]] <- list(
      ion_name = ion, transitions = as.numeric(transitions), components = as.integer(components)
    )
  }
  config$systems[[sys_idx]]$ion_groups <- ion_groups
  config
}

config_parameter_count <- function(config) {
  sum(vapply(config$systems, function(s) {
    sum(vapply(s$ion_groups, function(g) g$components * 3L, integer(1)))
  }, integer(1)))
}

## ---- Ion-specific bounds table (direct port of vfit_mcmc.py) --------------

ION_BOUNDS_TABLE <- list(
  HI    = list(N = c(12.0, 22.0), b = c(5.0, 100.0), v = c(-500.0, 500.0)),
  CIV   = list(N = c(12.0, 16.0), b = c(5.0, 80.0),  v = c(-200.0, 200.0)),
  OVI   = list(N = c(13.0, 16.0), b = c(10.0, 100.0), v = c(-300.0, 300.0)),
  SiIV  = list(N = c(11.0, 15.0), b = c(5.0, 60.0),  v = c(-150.0, 150.0)),
  MgII  = list(N = c(11.0, 16.0), b = c(5.0, 80.0),  v = c(-100.0, 100.0)),
  FeII  = list(N = c(11.0, 16.0), b = c(5.0, 60.0),  v = c(-100.0, 100.0)),
  AlIII = list(N = c(11.0, 15.0), b = c(5.0, 60.0),  v = c(-100.0, 100.0)),
  NV    = list(N = c(12.0, 15.0), b = c(10.0, 80.0), v = c(-200.0, 200.0)),
  OI    = list(N = c(13.0, 16.0), b = c(5.0, 50.0),  v = c(-100.0, 100.0)),
  SiII  = list(N = c(11.0, 16.0), b = c(5.0, 60.0),  v = c(-100.0, 100.0)),
  AlII  = list(N = c(11.0, 15.0), b = c(5.0, 60.0),  v = c(-100.0, 100.0)),
  CII   = list(N = c(13.0, 17.0), b = c(5.0, 50.0),  v = c(-100.0, 100.0)),
  NII   = list(N = c(13.0, 16.0), b = c(5.0, 60.0),  v = c(-100.0, 100.0)),
  SiIII = list(N = c(11.0, 15.0), b = c(5.0, 60.0),  v = c(-100.0, 100.0)),
  CIII  = list(N = c(13.0, 16.0), b = c(5.0, 80.0),  v = c(-150.0, 150.0)),
  NiII  = list(N = c(11.0, 15.0), b = c(5.0, 60.0),  v = c(-100.0, 100.0)),
  MnII  = list(N = c(11.0, 15.0), b = c(5.0, 60.0),  v = c(-100.0, 100.0)),
  CrII  = list(N = c(11.0, 15.0), b = c(5.0, 60.0),  v = c(-100.0, 100.0)),
  TiII  = list(N = c(11.0, 15.0), b = c(5.0, 60.0),  v = c(-100.0, 100.0)),
  ZnII  = list(N = c(11.0, 15.0), b = c(5.0, 60.0),  v = c(-100.0, 100.0))
)

set_bounds <- function(nguess, bguess, vguess, ions = NULL, ion_bounds = list(),
                        Nlow = NULL, blow = NULL, vlow = NULL,
                        Nhi = NULL, bhi = NULL, vhi = NULL) {
  ## R port of vfit_mcmc.py :: set_bounds(). nguess/bguess/vguess are the
  ## per-component initial guesses (length = total_components, in the same
  ## order as the N (or b, or v) block of theta). Returns list(lb, ub)
  ## suitable for concatenation as c(lb) / c(ub) matching
  ## theta = c(nguess, bguess, vguess).
  nguess <- as.numeric(nguess); bguess <- as.numeric(bguess); vguess <- as.numeric(vguess)
  m <- length(nguess)

  if (!is.null(ions)) {
    if (length(ions) != m) stop("length(ions) must match number of components")
    Nlow_v <- numeric(m); NHI_v <- numeric(m)
    blow_v <- numeric(m); bHI_v <- numeric(m)
    vlow_v <- numeric(m); vHI_v <- numeric(m)
    for (i in seq_len(m)) {
      ion <- ions[i]
      ion_data <- if (!is.null(ion_bounds[[ion]])) ion_bounds[[ion]] else ION_BOUNDS_TABLE[[ion]]
      if (is.null(ion_data)) {
        warning(sprintf("Ion '%s' not found in bounds table, using defaults", ion))
        Nlow_v[i] <- nguess[i] - 2.0;              NHI_v[i] <- nguess[i] + 2.0
        blow_v[i] <- max(2.0, bguess[i] - 40.0);   bHI_v[i] <- min(150.0, bguess[i] + 40.0)
        vlow_v[i] <- vguess[i] - 50.0;             vHI_v[i] <- vguess[i] + 50.0
      } else {
        Nlow_v[i] <- ion_data$N[1]; NHI_v[i] <- ion_data$N[2]
        blow_v[i] <- ion_data$b[1]; bHI_v[i] <- ion_data$b[2]
        vlow_v[i] <- ion_data$v[1]; vHI_v[i] <- ion_data$v[2]
      }
    }
  } else {
    Nlow_v <- nguess - 2.0;                    NHI_v <- nguess + 2.0
    blow_v <- pmax(bguess - 40.0, 2.0);         bHI_v <- pmin(bguess + 40.0, 150.0)
    vlow_v <- vguess - 50.0;                    vHI_v <- vguess + 50.0
  }

  if (!is.null(Nlow)) Nlow_v <- as.numeric(Nlow)
  if (!is.null(blow)) blow_v <- as.numeric(blow)
  if (!is.null(vlow)) vlow_v <- as.numeric(vlow)
  if (!is.null(Nhi))  NHI_v  <- as.numeric(Nhi)
  if (!is.null(bhi))  bHI_v  <- as.numeric(bhi)
  if (!is.null(vhi))  vHI_v  <- as.numeric(vhi)

  lb <- c(Nlow_v, blow_v, vlow_v)
  ub <- c(NHI_v, bHI_v, vHI_v)
  list(lb = lb, ub = ub)
}
