# rbvfit-in-R: rbvfit's Voigt-profile model fitted with your flip-move slice sampler

This is an R port of the core scientific pieces of the `rbvfit-master` Python
package you uploaded (Voigt-profile absorption-line fitting for
quasar/ISM/IGM spectra), rewired to use **your** sampler
(`multivariate_gibbs_sample_ASG_log_freshu_flip`, from `flipping.R`) instead
of rbvfit's own `emcee`/`zeus` backend.

Run `run_rbvfit_flip_mcmc.R` first — it simulates data, fits it, and prints
recovered vs. true parameters, so you can confirm everything works before
pointing it at real data.

## What was ported, and how it was checked

Only the pieces that actually determine the fitted numbers were ported —
not the GUI, results-plotting, `.rbv`/HDF5 save formats, or the
`ParameterManager`/`FitConfiguration` validation machinery, none of which
change what gets fit.

| Python source (rbvfit-master) | R file | What it does |
|---|---|---|
| `rb_setline.py`, `lines/atom_full.dat` | `voigt_lines.R`, `atom_full.dat` | Atomic line lookup (wavelength, oscillator strength, damping constant) |
| `core/voigt_approx.py` | `voigt_lines.R` (`H_tepper_garcia`) | Fast Tepper-García (2006) Voigt-Hjerting approximation |
| `scipy.special.wofz` (used inside `core/voigt_model.py`) | `voigt_lines.R` (`humlicek_w4`) | Exact Voigt-Hjerting function. R has no built-in Faddeeva function, so this implements the classic Humlicek (1982) 4-region rational approximation from scratch. |
| `core/fit_configuration.py` (systems/ion groups) | `rbvfit_config.R` | Build up the absorber configuration (`add_system`) |
| `vfit_mcmc.py :: set_bounds`, `ION_BOUNDS_TABLE` | `rbvfit_config.R` | Ion-aware parameter bounds |
| `core/voigt_model.py` (`_cache_atomic_parameters`, `_setup_fast_mapping`, `_vectorized_voigt_tau`, `_evaluate_compiled_model`, kernel convolution) | `voigt_model.R` | Building the flat parameter arrays + evaluating model flux, including the instrumental Gaussian convolution |
| `vfit_mcmc.py :: vfit.lnprior/lnlike/lnprob` | `rbvfit_likelihood.R` | The statistical model (uniform prior in bounds, Gaussian likelihood) |
| your `flipping.R` | `flipping_sampler.R` | The MCMC engine itself — **unchanged** except the `source()` line now points at the reconstructed helper file (see below) |

### Numerical validation performed in this session

I installed the actual Python `rbvfit` package (plus `numpy`/`scipy`/`astropy`)
and R itself in the sandbox, and cross-checked the R port directly against
the real Python code:

- **Line lookup**: `rb_setline()` in R returns identical wavelength/f-value/
  damping-constant to the Python version for every line tested.
- **Voigt-Hjerting function**: `humlicek_w4()` checked against
  `scipy.special.wofz` across weak/intermediate/strong-damping test points —
  agreement better than 1e-4 relative error everywhere, consistent with the
  algorithm's known accuracy.
- **Gaussian instrumental kernel + convolution**: `gaussian_kernel1d()` +
  `convolve1d_nearest()` checked against `astropy.convolution.Gaussian1DKernel`
  + `scipy.ndimage.convolve1d(..., mode='nearest')` — agreement to machine
  precision (~1e-15).
- **Full model flux** (multi-line, multi-system, with convolution): checked
  against `VoigtModel.evaluate()` from the real Python package for a
  2-system (SiII + HI) configuration — max relative error ~7e-5 across the
  spectrum, i.e. essentially just the Voigt-Hjerting approximation error
  above.
- **End-to-end pipeline**: simulated a 2-component SiII spectrum, added
  noise, ran the flip-move sampler through the R-ported likelihood, and
  recovered all 6 injected parameters (N₁, N₂, b₁, b₂, v₁, v₂) to within
  ~1.5 posterior standard deviations. This is `run_rbvfit_flip_mcmc.R`
  PART A, so you get this same check every time you run it.

## One thing that was *not* in your upload

`flipping.R` does `source("univariate slice and effective support.R")`,
but that file wasn't part of what you uploaded, so `kernel_from_nlk()` and
`effective.support()` could not be ported — only reconstructed from how
`flipping.R` calls them. `univariate_slice_and_effective_support.R` in this
delivery is that reconstruction, clearly marked as such in its header.

In practice this matters less than it sounds: in both the demo and the
real-data template, `bounds` is passed to the sampler explicitly (the
ion-aware `lb`/`ub` from `set_bounds()`), so the sampler's own code never
actually calls `effective.support()` — it only needs the file to exist so
the `source()` line doesn't error. If you still have your original helper
file, drop it in instead (same filename) and everything else is unaffected.

## A note on the flip move for this application

The flip move proposes `theta -> -theta` for the whole parameter vector at
once. Column densities `N` (log₁₀ cm⁻², ~11–22) and Doppler parameters `b`
(km/s, ~5–150) are always positive and bounded well away from zero, so
reflecting them through the origin almost always lands outside `[lb, ub]`
and gets rejected — exactly the "harmless, no correctness cost" case your
own header comment anticipates (you'll see `flip_accept_rate` come out near
0% in the output). It only pays off for parameters with a genuine ± mirror
degeneracy through zero; velocity `v` alone might qualify in some setups,
but here it's flipped jointly with `N` and `b`, which blocks it. If you want
the flip move to do real work here, it would need to act on the velocity
block only — that's a change to `flipping_sampler.R`'s proposal step, which
this delivery leaves untouched since it changes the sampler's semantics
rather than just wiring it up to a new model.

## Files in this delivery

- `flipping_sampler.R` — your sampler (only the `source()` path was fixed)
- `univariate_slice_and_effective_support.R` — reconstructed helper (see above)
- `voigt_lines.R` — atomic line lookup + Voigt-Hjerting function (exact + fast)
- `atom_full.dat` — the atomic line list data (copied from the repo)
- `rbvfit_config.R` — absorber configuration builder + ion-aware bounds
- `voigt_model.R` — model flux evaluation (Voigt optical depth + convolution)
- `rbvfit_likelihood.R` — prior / likelihood / posterior / `nlk` for the sampler
- `run_rbvfit_flip_mcmc.R` — driver script: synthetic self-test (Part A) +
  template for your own data (Part B)

## Requirements

Base R only for everything except plotting: `make_plots = TRUE` in the
sampler (as in your original file) needs `ggplot2`. `run_rbvfit_flip_mcmc.R`
defaults to `make_plots = FALSE` and instead writes a base-R trace-plot PNG
when run non-interactively (`RBVFIT_SAVE_PLOTS=1 Rscript run_rbvfit_flip_mcmc.R`),
so nothing extra is required to reproduce the validation run.

## Quick start

```r
source("run_rbvfit_flip_mcmc.R")   # runs Part A automatically
```

Then edit Part B of that script (clearly marked, line-for-line mirroring
`examples/example_voigt_fitter.py` from the original package) to point at
your own wavelength/flux/error vectors and absorber configuration.
