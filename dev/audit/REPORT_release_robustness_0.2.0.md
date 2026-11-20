# REPORT: release robustness and documentation corrections — weightedsurv 0.2.0

Task: `dev/tasks/TASK_release_robustness_0.2.0.md`. Branch `release-0.2.0-robust`, created from `master` at
`672829126e59f5b608b3e4cb860f3d9af3b7f271` (`6728291`). Version unchanged (0.2.0). Nothing pushed, merged,
tagged, submitted or deployed. No simulation, bootstrap or rerun of a stored simulation. Outputs of the
gates are in `dev/audit/release_robust/`.

Environment: R 4.5.2 (2025-10-31), aarch64-apple-darwin20, macOS 27.2.1; pandoc from quarto
(`/Applications/quarto/bin/tools` on the `PATH` and as `RSTUDIO_PANDOC`); `DISPLAY` unset.

## 1. Gates

| # | Result | Evidence |
|---|---|---|
| G1 | PASS | `g1_paths.txt`. Changed: `.github/workflows/pkgdown.yaml` (new), `DESCRIPTION`, `NEWS.md`, `README.md`, `dev/` (task copy, `dev/scripts/check_doc_deps.R`, `dev/audit/release_robust/`, this report), six files under `vignettes/` and the CSV of section 5. `git diff 6728291..HEAD -- R/ man/ NAMESPACE` is empty (0 bytes). `DESCRIPTION` diff: the `Config/Needs/website` field only (Imports, Suggests, Version unchanged; G4 showed no missing Suggests). `tests/` and `.Rbuildignore`: no diff. |
| G2 | PASS | `Rscript dev/scripts/check_doc_deps.R` exits 0 (`g2_g3_check_doc_deps.txt`): no undeclared package; every package used by an article is in Imports, Suggests or Config/Needs/website; every package used is on CRAN, so none is non-CRAN. |
| G3 | PASS | Same output, second listing: the 7 uses of Suggests packages (other than knitr, rmarkdown, testthat) in `weightedsurv_examples.Rmd` are guarded (lines 82, 164, 286, 808, 809, 886, 887; table in section 2). `tests/` uses no Suggests package other than testthat. |
| G4 | PASS | `g4_00check.log`, `g4_testthat.Rout`: `R CMD check --as-cran`, all packages installed: `Status: 2 NOTEs`, 0 errors, 0 warnings. NOTEs are the two of Stage 3 (incoming feasibility: README URL `.../articles/articles/estimated_weights_methodology.html` status 404; HTML Tidy not recent enough). Tests `[ FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 1691 ]`. |
| G5 | PASS | `g5_00check_nosuggests.log`, `g5_testthat_nosuggests.Rout`: `Status: 1 NOTE`, 0 errors, 0 warnings. "Packages suggested but not available for checking: 'adjustedCurves', 'data.table', 'DiagrammeR', 'dplyr', 'gt', 'pammtools'". NOTE: HTML Tidy (as in G4). Vignette rebuild OK. Tests `[ FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 1691 ]`. The URL NOTE of G4 is absent: "checking CRAN incoming feasibility ... OK" (the URL check needs the `curl` package, which is not in this library). |
| G6 | PASS | `g6_compare.txt`: knit objects of `weightedsurv_examples.Rmd` at `master` against this branch, all packages installed. 87 objects in common; `identical()` with environments removed: 87 of 87 (as saved: 57 of 87, the same environment-only differences as Stage 3 G2). Objects only in the new knit: `has`, `if_has` (the guard flags and the inline helper). |
| G7 | PASS | `g7_workflows.txt`: both workflow files parse as YAML; `pkgdown.yaml` has `needs: website`; `R-CMD-check.yaml` covers macos-latest, windows-latest, ubuntu-latest (R release); `^\.github$` is in `.Rbuildignore`. |
| G8 | PASS | `g8_check_pkgdown.txt`: "No problems found". `build_site()` completed twice (`g8_build_site_live.txt`; `g8_build_site_skip_live_sim.txt` with `WS_SKIP_LIVE_SIM=1`), built from the clean clone into scratch destinations. Every article indexed in `_pkgdown.yml` has its page in both builds. With `WS_SKIP_LIVE_SIM=1` the simulations page has its 3 figures and 5 tables from the stored file and no live-run output (`g8_site_pages.txt`). `git status` shows nothing from `docs/`. |
| G9 | PASS | Section 4. Numerals added to prose: 1.96, 0.95, 95% (constants of `cumulative_rmst_bands()`), 200 and 300 (`mart_draws` settings). None of them is a result. |
| G10 | PASS | Section 5. |

## 2. Dependency inventory

From `dev/scripts/check_doc_deps.R` (R chunks and inline R of every `.Rmd` under `vignettes/`, and `tests/`;
comment lines skipped). "Guarded" counts the uses whose chunk has a non-literal `eval =`, or whose line is
a conditional, or whose inline code is wrapped in `if_has()`.

| Package | Where used | Declared | On CRAN | Guarded uses |
|---|---|---|---|---|
| adjustedCurves | examples vignette | Suggests | yes | 1 of 1 |
| data.table | examples vignette; Cox article; simulations article | Suggests | yes | 2 of 5 (vignette and Cox article guarded; simulations article not) |
| DiagrammeR | examples vignette | Suggests | yes | 1 of 1 |
| dplyr | examples vignette; simulations article | Suggests | yes | 2 of 11 (vignette guarded; simulations article not) |
| gt | examples vignette; Cox article; simulations article | Suggests | yes | 2 of 4 (vignette and Cox article guarded; simulations article not) |
| pammtools | examples vignette | Suggests | yes | 1 of 1 |
| knitr | all `.Rmd` | Suggests | yes | not guarded (vignette engine) |
| testthat | `tests/testthat.R` | Suggests | yes | not guarded (test runner) |
| doFuture, foreach, future, future.callr, parallelly, tictoc | simulations article (simulation driver, run only by the live run) | Config/Needs/website | yes | the live run chunk is guarded by `run_live`; the definitions are not |
| gsDesign2, simtrial | simulations article | Config/Needs/website | yes | `scenario-setup` chunk guarded by `live_pkgs_ok`; uses inside function definitions not |
| htmltools | simulations article (code-folding hook) | Config/Needs/website | yes | 0 of 1 |
| tibble | simulations article | Config/Needs/website | yes | 0 of 4 |
| tidyr | simulations article | Config/Needs/website | yes | 0 of 4 |
| ggplot2, survival | articles; tests | Imports | yes | not needed |
| stats, utils, grDevices | tests | base R | — | not needed |

Full per-line listing: `dev/audit/release_robust/g2_g3_check_doc_deps.txt`.

Guards added:

- `weightedsurv_examples.Rmd`: a hidden first chunk sets `has` (one `requireNamespace()` flag per Suggests
  package) and `if_has()` (returns "(needs <pkg>)" instead of evaluating its argument when the package is
  absent). Guarded: DiagrammeR flowchart; `library(dplyr)` in the setup chunk (`if`); the baseline table of
  node-positive women (two chunks); every chunk, the `fig.cap` and every inline R of the section "Compare GBSG
  experimental to Rotterdam control" (its first chunk uses `%>%`/`mutate()`, and the rest of the section
  depends on its objects); the Cox table chunk (data.table and gt); the four adjustedCurves/pammtools chunks
  and the inline R of the paragraph that reads them. Sentences added: "This table needs the dplyr package.";
  "The code of this section uses the dplyr package: its tables, figures and numbers need that package.";
  "The table needs the data.table and gt packages."; "This comparison needs the adjustedCurves and pammtools
  packages." The script detects `library()`, `require()`, `requireNamespace()` and `pkg::`, not unprefixed
  verbs such as `%>%` or `mutate()`; those uses are in chunks that are guarded.
- Simulations article: setup chunk defines `pkgs`, `live_pkgs_ok` and
  `run_live <- live_pkgs_ok && !nzchar(Sys.getenv("WS_SKIP_LIVE_SIM"))`; `run-simulations-demo` has
  `eval = run_live`; `scenario-setup` has `eval = live_pkgs_ok`. The two places that said the minimal run is
  evaluated at build time now give the condition. The results tables and figures read the stored `.rds` in
  both cases.
- `weightedcox_methodology.Rmd`: `has_tables` flag; `library(data.table)` conditional; `results-table` chunk
  `eval = has_tables`, with the sentence "The table needs the data.table and gt packages."
- Other articles use only Imports packages and knitr: no guard added.
- Tests: no change (no Suggests package used).

Renders (`render_warn2.txt`): all 7 `.Rmd` files render with `options(warn = 2)`. The simulations article
sets `options(warn = -1)` in its own setup chunk, so `warn = 2` does not apply inside that article.

## 3. Workflows

Found: `.github/workflows/R-CMD-check.yaml` (r-lib/actions v2; push to `master`, pull requests,
`workflow_dispatch`; macos-latest, windows-latest, ubuntu-latest, all R release; `needs: check`). It is not
the `check-standard` matrix (no devel, no oldrel-1). Unchanged, since a check workflow exists. No pkgdown
workflow existed.

Added: `.github/workflows/pkgdown.yaml` (r-lib/actions v2 standard): triggers push to `master`, pull
requests, published releases and `workflow_dispatch`; `setup-r-dependencies` with
`extra-packages: any::pkgdown, local::.` and `needs: website`; builds with `build_site_github_pages()`; the
deploy step (JamesIves/github-pages-deploy-action, branch `gh-pages`, folder `docs`) runs only when
`github.event_name` is `push` or `release`. `.Rbuildignore` already had `^\.github$`. A `gh-pages` branch
exists locally and on `origin`; `docs/` is not tracked by git.

On a pull request, both workflows run: R-CMD-check (three OS jobs) and pkgdown (build only, no deploy).
CC cannot see Actions results.

## 4. Documentation corrections

`R/` lines read: `R/weighted_helpers.R:83-144` (`wt.rg.S()`), `:233-270` (`get_validated_weights()`: Schemper
takes `G.pool`, XO takes `Ybar`); `R/weightedCox.R:308` (`G.pool = survG_all`, `Ybar = ybar_all`),
`:87-109` (`score_calculation()`), `:400-425` and `:520-540` (resampling, `sig_bhat_star`, `Z_bstar`,
de-biasing, `z.score`), `:355-365` (`z.score = u.zero / sqrt(sig2_score)`); `R/df_counting.R:525-530`;
`R/rmst_calculations.R:70-160` (`cumulative_rmst_bands()`); `R/kmplotting_helpers.R:451-464`
(`plotKM.band_subgroups()` arguments); `R/km_wlr_calculations_helpers.R:567` (`KM_diff()`). Line numbers
given in the task had moved; the text was located by content.

**C1 — `README.md`, weighting-scheme table.** `wt.rg.S()` uses `S_left = c(1, S[-n])` for every scheme.

| Row | Before | After |
|---|---|---|
| intro | "The package supports multiple weighting schemes for weighted log-rank tests:" | "... (S is the pooled Kaplan-Meier curve, S(t-) its left-continuous version; see `wt.rg.S()`):" |
| FH | w(t) = S(t)^ρ × (1-S(t))^γ | w(t) = S(t-)^ρ × (1-S(t-))^γ |
| MB | w(t) = 1/max(S(t), S(t*)); "Modest downweighting after cutoff time" | w(t) = 1/max(S(t-), S(t*)); "Weight rises until the cutoff time t* and is constant after it" |
| Schemper | w(t) = S(t)/G(t); "Adjust for informative censoring" | w(t) = S(t-)/G(t-), G the censoring K-M curve; "Average hazard ratio approximately free of the censoring distribution" |
| XO | w(t) = S(t)/Y(t) | w(t) = S(t-)/Y(t), Y the number at risk (use case unchanged) |
| Custom | unchanged (matches `ifelse(tpoints <= t.tau, w0.tau, w1.tau)`) | — |

**C2 — `README.md`, Dependencies.** Added after the Suggests list: "The website articles need further
packages, listed in the `Config/Needs/website` field of `DESCRIPTION`." Imports and Suggests lists unchanged
(equal to `DESCRIPTION`).

**C3 — `weightedcox_methodology.Rmd`, weighting schemes.**
- Schemper. Before: $w(t) = \hat{\bar F}(t^-)/\hat{\bar G}(t)$, "where $\hat{\bar G}$ is the pooled Kaplan-Meier
  censoring survivor function." After: $\hat{\bar F}(t^-)/\hat{\bar G}(t^-)$, "where $\hat{\bar G}(t^-)$ is the
  left-continuous pooled Kaplan-Meier censoring survivor function (the weight is 0 where
  $\hat{\bar G}(t^-) = 0$)." (code: `Scensor_left <- c(1, Scensor[-n])`; `ifelse(Scensor_left > 0, S_left / Scensor_left, 0)`.)
- XO. Before: $\hat{\bar F}(t^-)/\{\bar Y(t)/n\}$. After: $\hat{\bar F}(t^-)/\bar Y(t)$, "where $\bar Y(t)$ is the
  pooled number at risk (with subject weights, the weighted number at risk; the weight is 0 where
  $\bar Y(t) = 0$)" (code: `ifelse(Ybar > 0, S_left / Ybar, 0)`, `Ybar = ybar_all`).
- FH, MB and custom formulas agree with `wt.rg.S()`; unchanged.

**C4 — K(t).** `score_calculation()`: `K <- wt_rg * ybar1 * ybar0 / (ybar0 + ybar1)`, no
$\sqrt{(n_0+n_1)/(n_0 n_1)}$.
- `weightedcox_methodology.Rmd` (now line 443). Before: "with $K(t) = \sqrt{\frac{n_0+n_1}{n_0 n_1}}\, w(t)
  \frac{\bar Y_0(t)\bar Y_1(t)}{\bar Y(t)}$." After: "with $K(t) = w(t)\frac{\bar Y_0(t)\bar Y_1(t)}{\bar Y(t)}$,
  as `score_calculation()` computes it (without a factor $\sqrt{(n_0+n_1)/(n_0 n_1)}$). Such a constant factor
  would multiply the statistic and its standard error alike, so the z-statistic is unchanged by it."
- `martingale_resampling_applications.Rmd` (line 71). Same change, in one sentence. Lines 76-90 (martingale
  representation, resampling recipe) do not write the factor; unchanged.

**C5 — `weightedcox_methodology.Rmd`, summary.** Before: "**Simultaneous inference**: Support for
simultaneous confidence intervals needed when the weighting scheme is selected via a combination test".
After: "**Simultaneous inference**: for each scheme, `cox_rhogamma()` returns the standardized resampled values
(`fit_resamples$Z_bstar`), based on the same multipliers $G_i$ in every call; the critical value $c_{\max}$ of
simultaneous confidence intervals across schemes (needed when the weighting scheme is selected via a
combination test) is computed from them and is not returned".

**C6 — `weightedkm_methodology.Rmd`, summary.** Before: "..., with careful attention to numerical validation
against established packages." After: "The `weightedsurv` package provides an integrated implementation of these
methods. What is checked against other packages (in the `weightedsurv_examples` vignette): with weights that
are not all equal, the weighted K-M standard errors equal those of `survival::survfit()` with the same weights;
the weighted K-M curves equal those of `adjustedCurves` (`method = "iptw_km"`), whose standard errors differ
from those of `weightedsurv`."

**C7 / C8 — `weightedrmst_methodology.Rmd`, intervals and bands.** `cumulative_rmst_bands()`: pointwise
`rmst_time ± 1.96*sqrt(sig2_rmst_time)`; band critical value `quantile(sups, 0.95)`; no `alpha` argument; only
the difference `dhat` of a `KM_diff()` fit.
- Added at the start of "Pointwise Confidence Intervals": "`cumulative_rmst_bands()` computes the intervals and
  bands of this section and the next for the difference of the two arms only, $D(\tau)$ (section "Two-Sample
  Simultaneous Bands"); the one-sample formulas show the construction. It uses 95% intervals and bands (1.96
  and the 0.95 quantile) and has no `alpha` argument."
- "For any $\alpha \in (0,1)$, a two-sided $(1-\alpha)$ pointwise confidence interval ... $\pm z_{1-\alpha/2}$"
  → "A two-sided 95% pointwise confidence interval ... $\pm 1.96$".
- Simultaneous band: $c_\alpha$ → $c_{0.95}$, "$\ge 1-\alpha$" → "$\ge 0.95$", "equal-precision confidence band"
  → "equal-precision 95% confidence band".
- Two-sample band: $\tilde c_\alpha$ → $\tilde c_{0.95}$, "$\ge 1-\alpha$" → "$\ge 0.95$"; "at the $(1-\alpha)$
  level" → "at the 95% level".

**C9 — `weightedrmst_methodology.Rmd`, workflow step 1.** Before: "1. **Prepare data** via `df_counting` or
`get_dfcounting`, optionally with a `weight.name` for propensity-score-weighted analyses." After: "1. **Prepare
the data frame** with the time, event and treatment columns, and optionally a weight column (for example
stabilized propensity-score weights from `ps_weights()`). `plotKM.band_subgroups` and `cumulative_rmst_bands`
take this data frame (`df`, with `tte.name`, `event.name`, `treat.name` and `weight.name`), not a `df_counting()`
result; `cumulative_rmst_bands` also takes the `KM_diff()` fit made with `draws.band > 0` (`fit`)." Steps 2-3
and the code block agree with the argument lists; unchanged.

**C10 — `weightedcox_wald-vs-logrank_simulations.Rmd`.** Read in full at HEAD against `R/` and the stored run
(`sim_results/wald_vs_logrank_full.rds`: 5000 trials per scenario, 6 scenarios, 30,000 completed, 0 failed,
`mart_draws` 200, 12 workers, 71.8 minutes, package 0.2.0 at `ac76f82`). Statements found to contradict the code
or the stored run, and corrected:
1. Line 596 (output dictionary). Before: "z-statistics from both `simtrial` (weighted log-rank) and
   `weightedsurv` (weighted Cox / Wald)". After: "... and `weightedsurv` (the standardized score of the weighted
   Cox model at $\beta = 0$, that is the weighted log-rank z, from `cox_rhogamma()`), together with the weighted
   Cox estimates". Code: `cox_rhogamma()$z.score = u.zero / sqrt(sig2_score)`; in the stored run
   `logrankz - fh00z_mine` lies in [-2.1e-04, 8.0e-04].
2. Z-statistic table, "Notes" of `fh00z_mine`, `fh05z_mine`, `fh01z_mine`, `mb12z_mine`, `fhe1z`, `fhe2z`,
   `t601z` (7 cells). Before: "Wald z from weighted Cox". After: "Score z (weighted log-rank) from weighted Cox".
3. Z-statistic table, source note. Before: "*M* = `mart_draws` (default 300). De-biased variants use the
   martingale-residual bootstrap of Xu & O'Quigley." After: "*M* = `mart_draws` (200 in the stored run below; the
   default of `sim_fn_analysis()` is 300). De-biased variants subtract the mean of the *M* martingale-resampled
   values from the estimate (or score)." Code: `bhat_debiased <- bhat - mean(bhat_center_star)`,
   `u.zero_debiased <- u.zero - mean(score_star_null)`, multipliers `G`.
4. Companion-estimates table, `_sigbhatstar`. Before: "De-biased SE via martingale resampling". After: "SE via
   martingale resampling (SD of the resampled estimates)". Code: `sig_bhat_star <- sqrt(var(bhat_center_star))`.

Read and found consistent: scenario description and control survival (lines 471-475), sample-size design
(scenario 2, 85% power, one-sided 0.025), "more than an hour on 12 workers" (71.8 minutes), the stored-run
settings (inline), the Wald/log-rank rejection rules (line 939), the list of FH variants (line 1041), the
sign convention note, and the paired-facet statement (correlations of the pairs are 1.0000 in the stored run).
Not a prose statement and not changed: the "all-fh-data" and "simtrial-comparison-data" chunks colour
null-scenario labels with the hard-coded threshold 0.0264, while `pnull_threshold` computed in the same chunks
is 0.0293 for 5000 trials. The note "These functions will be included in a future CRAN release" (lines 74-77)
is not checkable against the code; unchanged.

`NEWS.md`, `## Documentation`, one line added: "Vignette chunks that use a suggested package run only when it
is installed; website-only packages are listed in `Config/Needs/website`."

## 5. Section 5: stored CSV, OC cell 8, RMST, month 36

Per-replicate files present: `dev/audit/oc/raw/proposed_cell08.rds` (2,000 replicates; also `cran_cell08.rds`,
`boot_cell08.rds`, 500) and `dev/audit/stage2e/raw/m7_cell31.rds` (2,000). Data checksums `$chk` of cell 8 and
cell 31 are identical; the RMST estimates at month 36 differ by at most 4.1e-14.

Cause. The two summaries measure bias and coverage against different truths:
- `dev/audit/oc/oc_sim.R:210-211` stores per replicate `R.truth.36 = R_true(at_points)[jj]`, the true RMST at
  the last grid point $t_j \le 36$ of that replicate. `dev/audit/oc/oc_summarise.R:107` passes it as
  `truth_mat`, so `pw()` (lines 37-46) computes bias = mean(estimate − R(t_j)) and coverage of R(t_j), while the
  `truth` column holds the nominal R(36) = 4.01351390580877. The bootstrap-replicate rows (line 163, estimator
  `proposed_rect`, `proposed_ew`) use the same per-replicate truth; the bootstrap rows themselves (line 151)
  use R(36).
- `dev/audit/stage2e/oc_summarise.R:100-101` (M7) uses the fixed truth R(36).
- `dev/audit/ps_by_time/ps9_summarise.R:87` and `:139` write `mean_estimate <- truth + bias`. For cell 8 this gives
  4.01351 − 0.01226 = 4.00125, which is not the mean of the estimates (3.96268). The mean of the per-replicate
  truths is 3.9749415589271, 0.0385723 below R(36): the gap between the two biases.

Whether displayed: `estimated_weights_methodology.Rmd` line 296 removes cell 8 from `main`, so Tables 1-4,
Figure 1 and the reading values built with `v()`/`vals()` do not use it. Table 5 (lines 421-433) includes the
row and shows `coverage_known_bootreps`, `coverage_estimated_bootreps`, bootstrap coverages and interval lengths;
it does not show the mean or bias. No table or inline R shows this row's mean or bias. The Table 5 reading
(line 473) is inline-computed over the Table 5 rows: "differ by at most" changes from 0.006 to 0.008; the ranges
of the length ratio (0.97 to 0.99) and of "known" minus bootstrap (0.016 to 0.036) do not change
(`sec5_table5_reading.txt`).

Correction (`sec5_cell8_R36.R` / `.txt`, `sec5_write_csv.R`): recomputed from `proposed_cell08.rds` against
`true_value` = R(36); only that CSV line changed, every other line byte-identical.

| Column | Old | New |
|---|---|---|
| mean_estimate | 4.00124978724059 | 3.96267744035892 |
| bias | -0.0122641185681805 | -0.0508364654498496 |
| bias_mcse | 0.0297311700875105 | 0.0297370073093352 |
| coverage_known | 0.98 | 0.9805 |
| coverage_known_mcse | 0.00313049516849971 | 0.00309190475273738 |
| coverage_estimated | 0.947 | 0.946 |
| coverage_estimated_mcse | 0.0050095408971282 | 0.00505390937789747 |
| coverage_known_bootreps | 0.974 | 0.976 |
| coverage_estimated_bootreps | 0.94 | 0.934 |

The new mean, bias, bias MCSE and coverages equal those of the M7 cell 31 row. `empirical_sd`, mean SEs,
interval lengths and the bootstrap columns are unchanged. The row's `note` text ("copied from
dev/audit/oc/oc_summary.csv") was left as it was.

The same per-replicate truth was found in the other OC RMST rows of the CSV (cells 5 = D-400 and 10 = F-400, all
months), which Table 2 displays (bias and coverage); for example D-400 month 36 bias −0.03777 stored against
−0.07641 with the fixed truth, F-400 month 48 −0.55203 against −0.61966, and `coverage_known` differing by up to
0.0185 (F-400 month 1: 0.6945 against 0.6760). These rows are outside section 5 and were not changed.

## 6. Check results

Tarball: built by `R CMD build` from a clean `git clone -b release-0.2.0-robust` at `8107811` (`R/`, `man/`,
`NAMESPACE`, `DESCRIPTION`, tests and vignette identical to the report commit); copied to
`~/Downloads/weightedsurv_0.2.0.tar.gz` (replacing the earlier one). 4,330,967 bytes; sha256
`a7f38c2e2fae1b156c2ece9ddddc4c95f6890c660fb48158ed5f378b91f9c0a9`.

Full check (all packages): `R CMD check --as-cran`, `Status: 2 NOTEs` (section 1, G4); tests 1691 passed, 0
failed, 0 skipped.

No-Suggests check: all packages are installed in R's system library (`.Library`), which `.libPaths()` always
includes, so a temporary library alone cannot hide them. A private copy of `R_HOME` was made in the scratch
directory (the R framework's `Resources/` without `library/`, with `bin/R` pointing to the copy) whose
`library/` holds only the base and recommended packages plus `ggplot2`, `rlang`, `survival`, `knitr`,
`rmarkdown`, `testthat` and their recursive Depends/Imports/LinkingTo (82 packages,
`g5_library_packages.txt`). `R CMD check --as-cran` with `_R_CHECK_FORCE_SUGGESTS_=false`,
`R_LIBS_USER=NULL`, `R_LIBS_SITE=NULL` on the same tarball: `Status: 1 NOTE` (HTML Tidy), 0 errors, 0 warnings;
vignette rebuilt; tests `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1691 ]`.

`cran-comments.md` unchanged: the full check result is the same as the one it records.

## 7. Branch state

- `master` HEAD at start: `672829126e59f5b608b3e4cb860f3d9af3b7f271`.
- Commits on `release-0.2.0-robust` since `master` (before this report's commit):
  - `b97d4a0` task: TASK_release_robustness_0.2.0.md
  - `94a4424` dependency guards; Config/Needs/website; `dev/scripts/check_doc_deps.R`; C3-C5, C10
  - `61a50f6` documentation corrections C1, C2, C4, C6-C9; NEWS line
  - `4c6ebab` section 5 CSV row
  - `8107811` pkgdown workflow; G2, G3, G6, G7 outputs
  - the report commit (this file and the G1, G4, G5, G8 outputs), the branch HEAD.
- `git status` before the report commit: only the untracked G1/G4/G5/G8 outputs under
  `dev/audit/release_robust/`; nothing from `docs/`. Not pushed.
