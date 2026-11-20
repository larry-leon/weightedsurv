# weightedsurv 0.2.0

This is an update of weightedsurv 0.1.0 (on CRAN).

## Test environments

* macOS 27.2.1 (aarch64-apple-darwin20), R 4.5.2: `R CMD check --as-cran` on the tarball built by
  `R CMD build` from a clean `git clone`; also run with the suggested packages not installed
  (`_R_CHECK_FORCE_SUGGESTS_=false`): 0 errors, 0 warnings, 1 note (HTML Tidy, below).
* GitHub Actions (r-lib/actions): windows-latest, macos-latest, ubuntu-latest, R release:
  0 errors, 0 warnings.

## R CMD check results

0 errors | 0 warnings | 2 notes (local macOS check)

Tests: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1691 ]`.

* NOTE, "checking CRAN incoming feasibility": one possibly invalid URL,

  ```
  URL: https://larry-leon.github.io/weightedsurv/articles/articles/estimated_weights_methodology.html
    From: README.md
    Status: 404
  ```

  This is a new article of the package website, added in 0.2.0. The local check ran before the
  website for 0.2.0 was deployed; the site has since been deployed and the URL now resolves.

* NOTE, "checking HTML version of manual": "Skipping checking HTML validation: 'tidy' doesn't look
  like recent enough HTML Tidy." This comes from the HTML Tidy installed on the test machine, not
  from the package.

## Reverse dependencies

`tools::package_dependencies("weightedsurv", reverse = TRUE, which = "all")` gives one package on
CRAN: forestsearch (it imports weightedsurv).

Reverse dependency checked: forestsearch, CRAN version 0.5.0 (source tarball from CRAN). `R CMD
check` was run on it twice, in separate temporary libraries:

* with weightedsurv 0.1.0 (CRAN): Status OK; tests `[ FAIL 0 | WARN 0 | SKIP 163 | PASS 42 ]`;
* with weightedsurv 0.2.0 (this submission): Status OK; tests `[ FAIL 0 | WARN 0 | SKIP 163 | PASS 42 ]`.

No error, warning or note in either run; examples and the vignette rebuild are OK in both.
forestsearch calls `df_counting()`, `plot_weighted_km()` and `plotKM.band_subgroups()` without
subject weights.

## Changes since 0.1.0

The full list is in `NEWS.md`. In summary:

* **Weighted standard errors and tests change.** In 0.1.0, subject weights (for example
  inverse-probability weights) were treated in the variance formulas as if they were counts, so
  standard errors and test statistics depended on how the weights were scaled. In 0.2.0, when the
  weights are not all equal, the standard error of the weighted Kaplan-Meier curve is the robust
  (infinitesimal-jackknife) standard error, equal to that of
  `survival::survfit(weights = w, robust = TRUE)`, and the variance of the weighted log-rank score
  and of the weighted Cox estimate is built from per-subject score contributions, equal for
  Fleming-Harrington (0, 0) weights to the robust results of
  `survival::coxph(weights = w, robust = TRUE, ties = "breslow")`. The resampled curves, the
  simultaneous bands and the RMST bands of weighted analyses use the same per-subject
  contributions. Results no longer depend on the scale of the weights.
* Estimates (survival curves, differences, the RMST point estimate, the hazard ratio) do not
  change. Unweighted results, and results with all weights equal, do not change in these
  variances and tests.
* New: `ps_weights()` and the arguments `weight_variance` / `weight_model`, which carry the
  estimation of logistic propensity-score weights into standard errors, bands and tests (the
  default, `weight_variance = "known"`, treats the weights as fixed); `df_counting()$ess`
  (effective sample size by arm); `cumulative_rmst_bands(mean_rule = )`.
* Changed defaults: `plotKM.band_subgroups()` `draws` 20 to 0 and `draws.band` 20 to 1000;
  `KM_diff()` `show_resamples` `TRUE` to `FALSE`.
* Input checks now give errors that name the column or arm; Kaplan-Meier quantiles equal
  `quantile(survfit())`; `KM_diff()` returns the Kaplan-Meier value at every requested time
  point; plots draw each curve from time zero; several bug fixes (see `NEWS.md`).
* Documentation: the vignette and the website articles were revised, and a new website article
  gives the statistical details of estimated weights.
