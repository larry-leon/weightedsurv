# TASK: release robustness and documentation corrections before submission — weightedsurv 0.2.0

Written 2026-11-20. Executor: Claude Code (CC). Approved by Larry: implement the practically important revisions before the CRAN submission. Version stays 0.2.0 (not yet submitted).

## 0. Rules

- `master` now contains the merged `km-fixes-0.2.0` (HEAD includes `a19e189` and the report commit after it). Record `master` HEAD. **Create branch `release-0.2.0-robust` from `master` and work there.** First action on the branch: copy this file to `dev/tasks/TASK_release_robustness_0.2.0.md` and commit it.
- **No changes under `R/`, `man/`, `NAMESPACE`.** No change in any computed result.
- Allowed paths, with limits:
  - `DESCRIPTION`: add a `Config/Needs/website:` field only. `Imports`, `Suggests`, `Version` unchanged (unless G4 shows a package used by the CRAN-built vignette or tests is missing from Suggests; then add it and record why).
  - `vignettes/` (all files): dependency guards (section 2) and the text corrections of section 4.
  - `README.md`: the corrections of section 4.
  - `tests/`: **only** `testthat::skip_if_not_installed("<pkg>")` lines where a test uses a Suggests package. No other test change.
  - `.github/workflows/`, and `.Rbuildignore` (only to add `^\.github$` if missing).
  - `NEWS.md`: one line under `## Documentation`.
  - `cran-comments.md`: update only if a check result changes.
  - `vignettes/articles/sim_results/ps_weight_results_by_time.csv`: only the correction of section 5, if needed.
  - `dev/`.
- Never push. Never merge. No tag. No submission anywhere (CRAN, win-builder, R-hub). No site deployment.
- No simulation, no bootstrap, no rerun of stored simulations. Recomputing a summary from stored per-replicate files is allowed (section 5).
- Gates (section 7) are stop-on-failure: if one fails, stop, write the report, do not ask.
- Report findings only: no recommendation, no follow-up task proposals.

## 1. Why

The simulations article uses 11 packages declared nowhere, so no build tool installs them and the site can fail to build on any machine but the one that built it. The same class of risk applies on CRAN: packages in `Suggests` must be used conditionally in vignettes and tests (Writing R Extensions; CRAN's "noSuggests" checks). If a suggested package is archived, an unconditional vignette chunk fails and CRAN asks for a fix under a deadline.

## 2. Dependency robustness

1. **Inventory.** For every `.Rmd` under `vignettes/` (including `articles/`) and every file under `tests/`, list each package used (`library()`, `require()`, `requireNamespace()`, `pkg::`). Promote `dev/audit/stage3/vignette_packages.R` to a standing script `dev/scripts/check_doc_deps.R` that prints, per package: where used, whether declared (Imports / Suggests / Config/Needs/website), whether on CRAN (`available.packages()`), and whether each use is guarded. Exit status non-zero if any package is undeclared.
2. **Declare article-only packages** in `Config/Needs/website:` (comma-separated, as pkgdown and `r-lib/actions/setup-r-dependencies` read it). Packages used only by articles under `vignettes/articles/` go here, not in Suggests.
3. **Packages not on CRAN.** Any package not on CRAN must not be needed to render: its chunks become `eval = FALSE` (code shown) or are guarded as in 4, and the article text says the package is needed to run that code.
4. **Guards.**
   - **CRAN-built vignette** (`weightedsurv_examples.Rmd`): every chunk that uses a Suggests package other than `knitr`/`rmarkdown` is evaluated only when that package is installed (`eval = requireNamespace("pkg", quietly = TRUE)`, or one setup flag per package). Where a guarded chunk's output is referred to in the text, add one sentence saying the output needs the package. Inline R that depends on a guarded object must not error when the object is absent.
   - **Simulations article:** the minimal live run is evaluated only when all its packages are installed, and can be switched off by an environment variable: `run_live <- all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)) && !nzchar(Sys.getenv("WS_SKIP_LIVE_SIM"))`. Tables and figures come from the stored file in both cases.
   - **Other articles:** the same guard for any package not in Imports.
   - **Tests:** `skip_if_not_installed()` for any Suggests package a test uses.

## 3. Continuous integration

Read `.github/workflows/`. Record what exists.

- **pkgdown workflow:** ensure `r-lib/actions/setup-r-dependencies` has `needs: website` so it installs `Config/Needs/website`. If there is no pkgdown workflow, add the standard `r-lib/actions` `pkgdown.yaml` with deployment only on push to `master` (and on release), not on pull requests.
- **R CMD check workflow:** if absent, add the standard `r-lib/actions` `check-standard.yaml` (macOS, Windows, Ubuntu release, devel, oldrel-1), triggered on push to `master` and on pull requests.
- Change nothing else in existing workflows. Validate every workflow file as YAML. `^\.github$` must be in `.Rbuildignore`.

CC cannot see Actions results. Record which workflows will run on the pull request.

## 4. Documentation corrections

Text only. Each change states what the code does, read from `R/` at HEAD. These were listed in `dev/audit/REPORT_vignette_revision.md` section 5 and left unchanged.

| # | File | Correction |
|---|---|---|
| C1 | `README.md`, weighting-scheme table | Magirr–Burman: the weight 1/max(S(t), S(t*)) rises until t* and is constant after it (not "downweighting after cutoff"). Schemper: describe it as in the Cox article (average hazard ratio approximately free of the censoring distribution). Check every other row of the table against `wt.rg.S()`. |
| C2 | `README.md` | Add one line in Dependencies or Installation noting that the website articles need further packages (listed in `DESCRIPTION` `Config/Needs/website`). The Imports and Suggests lists stay equal to `DESCRIPTION`. |
| C3 | `articles/weightedcox_methodology.Rmd` | Schemper weights: write the left-continuous censoring curve as the code uses it, F̄(t⁻)/Ḡ(t⁻). Xu–O'Quigley weights: Ȳ(t) as in the code, not Ȳ(t)/n. Check every scheme's displayed formula against `wt.rg.S()`. |
| C4 | `articles/martingale_resampling_applications.Rmd` (69, 76–85) and `articles/weightedcox_methodology.Rmd` (371) | K(t) carries a factor √{(n₀+n₁)/(n₀n₁)} that `score_calculation()` does not use. Write K(t) as coded and state that the z-statistic is unchanged by the factor. |
| C5 | `articles/weightedcox_methodology.Rmd`, summary | The bullet on simultaneous confidence intervals: say what is returned (`fit_resamples$Z_bstar` per scheme on common multipliers), as the corrected section on c_max already does. |
| C6 | `articles/weightedkm_methodology.Rmd`, summary | Replace "careful attention to numerical validation against established packages" with what is checked: weighted K-M standard errors equal `survfit()`'s; curves equal `adjustedCurves` `iptw_km`, whose standard errors differ. |
| C7 | `articles/weightedrmst_methodology.Rmd` (164–176) | Intervals and bands: `cumulative_rmst_bands()` uses 1.96 and the 0.95 quantile; no `alpha` argument. Replace the general α. |
| C8 | `articles/weightedrmst_methodology.Rmd` (162–178) | The function computes intervals and bands for the difference of two arms only; say so. |
| C9 | `articles/weightedrmst_methodology.Rmd` (277) | Workflow step 1: `plotKM.band_subgroups()` and `cumulative_rmst_bands()` take the data frame (and a `KM_diff()` fit), not a `df_counting()` result. Check against the argument lists. |
| C10 | `articles/weightedcox_wald-vs-logrank_simulations.Rmd` | Read it at HEAD for statements that contradict the code or the stored run (it was not read for accuracy before). Correct only such statements; list them. |

Do not change anything else in these files.

## 5. Stored CSV, cell 8 RMST row

In `vignettes/articles/sim_results/ps_weight_results_by_time.csv`, the row for OC cell 8, RMST, month 36 (bootstrap benchmark) has the same empirical SD, standard error and interval length as M7 cell 31, month 36, but a different mean (bias −0.0123 against −0.0508). The record says cells 8 and 31 share data and that the package values agree exactly for RMST.

1. Find the cause by reading `dev/audit/oc/oc_summarise.R`, `dev/audit/stage2e/oc_summarise.R`, and the per-replicate files `dev/audit/oc/raw/` and `dev/audit/stage2e/raw/` if present. No rerun.
2. Record whether the article displays that row's mean or bias (Table 5, reading, inline R).
3. If the stored OC summary value is wrong, write the corrected value into the CSV row (computed from per-replicate files) and record old and new values and the cause. Do not edit `dev/audit/oc/` or `dev/audit/stage2e/`. If the per-replicate files are absent, or the cause cannot be determined, change nothing and record that.

## 6. Final rebuild and checks

1. Render all vignettes and articles with `options(warn = 2)`.
2. Build the tarball from a clean `git clone` of the branch; `R CMD check --as-cran` (`DISPLAY` unset; pandoc on the path as in Stage 3).
3. **No-Suggests check:** in a fresh temporary library containing only `Imports`, `knitr`, `rmarkdown` and `testthat` (and their dependencies), run `R CMD check --as-cran` on the same tarball with `_R_CHECK_FORCE_SUGGESTS_=false`.
4. `pkgdown::check_pkgdown()`; `pkgdown::build_site()` locally, twice: once as is, once with `WS_SKIP_LIVE_SIM=1`. Do not commit `docs/`. Do not deploy.
5. Copy the tarball to `~/Downloads/weightedsurv_0.2.0.tar.gz` (replacing the earlier one).

## 7. Gates

| # | Check |
|---|---|
| G1 | `git diff --name-only <master HEAD>..HEAD` lists only allowed paths. `git diff <master HEAD>..HEAD -- R/ man/ NAMESPACE` is empty. `DESCRIPTION` diff is the `Config/Needs/website` field only (or a Suggests addition justified by G4). `tests/` diff contains only `skip_if_not_installed(` lines. |
| G2 | `Rscript dev/scripts/check_doc_deps.R` exits 0: no undeclared package; every package used by an article is in Imports, Suggests or Config/Needs/website; every non-CRAN package is off the render path. |
| G3 | Every use of a Suggests package (other than knitr, rmarkdown, testthat) in `weightedsurv_examples.Rmd` and in `tests/` is guarded (script listing, with line numbers). |
| G4 | `R CMD check --as-cran` with all packages: 0 errors, 0 warnings, NOTEs no more than the two of Stage 3; tests `FAIL 0`, `PASS 1691`. |
| G5 | No-Suggests check: 0 errors, 0 warnings. Record the NOTEs and the test counts (skips allowed). |
| G6 | Unchanged results: with all packages installed, the knitted objects of the analysis chunks of `weightedsurv_examples.Rmd` are `identical()` before and after this task (environments removed as in Stage 3 G2). |
| G7 | Workflows: every file in `.github/workflows/` parses as YAML; the pkgdown workflow has `needs: website`; an R CMD check workflow covering Windows, macOS and Linux exists; `^\.github$` in `.Rbuildignore`. |
| G8 | Site: `check_pkgdown()` no problems; `build_site()` completes both times; every indexed article has its page; with `WS_SKIP_LIVE_SIM=1` the simulations page still shows its tables and figures from the stored file. `git status` shows nothing from `docs/`. |
| G9 | Corrections C1–C10 each listed with before and after text and the `R/` lines read. No numeral in changed prose that is a result unless inline-computed. |
| G10 | Section 5 recorded (cause, whether displayed, CSV old/new values or "unchanged"). |

## 8. Report

`dev/audit/REPORT_release_robustness_0.2.0.md`: gate table; the dependency inventory (package × where used × declared field × on CRAN × guarded); workflows found and changed, and which run on a pull request; corrections C1–C10; section 5 result; check results (full and no-Suggests); tarball path, size, sha256; branch state (HEAD, commits since `master`, `git status`).

After the final commit, copy to `~/Downloads/` and print the paths: the report, the tarball, and `cran-comments.md` if changed. Do not push.
