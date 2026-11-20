# TASK: Stage 3 release prep — weightedsurv 0.2.0, branch `km-fixes-0.2.0`

Written 2026-11-20. Executor: Claude Code (CC). Approved by Larry, including D1–D4 from the vignette revision report and the compute go for the simulations article.

## 0. Rules

- Branch `km-fixes-0.2.0`; HEAD must contain `e95728a`. If not, stop.
- First action: copy this file to `dev/tasks/TASK_stage3_release_prep.md` and commit it.
- **No changes under `R/`, `man/`, `tests/`, `NAMESPACE`.** This task changes no code and no behaviour.
- Allowed paths: `vignettes/weightedsurv_examples.Rmd`, `vignettes/articles/martingale_resampling_applications.Rmd`, `README.md`, `DESCRIPTION` (`Version` and `Date` fields only), `NEWS.md`, `cran-comments.md`, `.Rbuildignore` (only if G10 requires it), `docs/`, `dev/`.
- **Never push. Never merge. No tag. No submission anywhere** (CRAN, win-builder, R-hub, GitHub release). `pkgdown::build_site()` only; no `deploy_to_branch()`.
- The forestsearch repository is read-only for this task: no edits, no commits there.
- Gates (section 6) are stop-on-failure: if one fails, stop, write the report, do not ask.
- Verify from source at the current HEAD.
- Report findings only: no recommendation, no follow-up task proposals.

**Approved compute (one run only):** the full-size knit of `articles/weightedcox_wald-vs-logrank_simulations.Rmd` inside the site build: 6 scenarios × 5,000 trials, 12 workers, about 72 minutes. Do not run it more than once. No other simulation or bootstrap.

## 1. Order of work

1. Task document committed.
2. Text rider (section 2).
3. Version, `NEWS.md`, `cran-comments.md` (section 3).
4. `R CMD check` (G3), tarball contents (G10), URLs (G9).
5. forestsearch check (section 4).
6. Site build (section 5). Last, because it is the long step.
7. Report (section 7); commit.

Commit per step; messages prefixed `release:` (rider commits `vignettes:` / `readme:`).

## 2. Text rider (approved D3, D4)

Text only. No chunk code changes except the one chunk header in R5. Locate each passage by content; line numbers are from `ac76f82`.

| # | File | Passage | Change |
|---|---|---|---|
| R1a | `weightedsurv_examples.Rmd` | Overview: "propensity-score adjustment appears to suggest a stronger impact consistent with the randomized trial" | State the three values for the restricted Rotterdam population, computed inline: weighted hazard ratio, unweighted hazard ratio, and the GBSG trial value, each with its interval. No claim that weighting supports or validates the observational comparison. No significance framing; describe where the estimates sit relative to each other. If the Overview precedes the chunks that compute these values, point to the section where they are shown instead of quoting numbers. |
| R1b | `martingale_resampling_applications.Rmd` | about 263: "broadly consistent with the randomized GBSG trial, supporting the validity of the observational comparison" | Same treatment, with the values of the analysis that passage discusses, inline. |
| R1c | `martingale_resampling_applications.Rmd` | about 408: "consistent with the GBSG trial findings" | Same treatment. |
| R2 | `martingale_resampling_applications.Rmd` | about 272: band written as the unstandardised supremum | Standardised supremum; band = estimate ± c · se, as corrected in `weightedkm_methodology.Rmd` and as coded in `KM_diff()`. |
| R3 | `README.md` | "Dependencies" | Match `DESCRIPTION` Imports and Suggests exactly; fix the unclosed parenthesis. |
| R4 | `README.md` | features / usage | Add a short mention of `ps_weights()`, `weight_variance = "estimated"` and `df_counting()$ess`, with one code snippet copied from a chunk that runs in the main vignette. Link to the estimated-weights article. |
| R5 | `weightedsurv_examples.Rmd` | "Uses resampling (bootstrap or similar)"; chunk header `messages = FALSE` | "multiplier (martingale) resampling"; `message = FALSE`. |

**Not to be touched** (recorded in the vignette report, section 5; stay in the record): the K(t) scaling factor, the Schemper and Xu–O'Quigley weight formulas, the RMST article's general-α text and single-curve description, the workflow step in the RMST article, the summary bullets of the Cox and KM articles, the README weighting-scheme table, and the simulations article (including its "NOT evaluated" comment).

## 3. Version, NEWS, cran-comments

- **`DESCRIPTION`**: `Version: 0.2.0`; update `Date` if the field exists. Nothing else. Record (do not change) whether every package used in vignette chunks is in Imports or Suggests.
- **`NEWS.md`**: a `# weightedsurv 0.2.0` section at the top covering every user-visible change since CRAN 0.1.0. Build it from the git history since the 0.1.0 release and from every `dev/audit/REPORT_km_fixes_*.md`, not from this document. One line per change, grouped:
  - New: `ps_weights()`; `weight_variance` / `weight_model`; `df_counting()$ess`; the estimated-weights article; anything else the reports show.
  - Changed results for weighted analyses: robust weighted Kaplan-Meier variance and influence-based resampling (Stage 2c); score-residual variance for the weighted log-rank and Cox (Stage 2d); weighted `logrank_results` now the `coxph` robust score test; any earlier stage on this branch.
  - Changed defaults: `plotKM.band_subgroups()` `draws.band` now 1000; any other the reports show.
  - Unchanged: state what the reports establish about unweighted results, in their words; do not claim more than they show.
  - Fixes from earlier stages on this branch, as the reports record them.
  If a `0.2.0` section already exists, complete it; do not duplicate.
- **`cran-comments.md`**: for 0.2.0: test environment(s) actually used; `R CMD check` result with every NOTE listed and explained; reverse-dependency result (section 4); a short summary of changes since 0.1.0, including that weighted standard errors and tests change and why.

## 4. forestsearch check

forestsearch uses weightedsurv. Check it against this build.

1. List CRAN reverse dependencies of weightedsurv (`tools::package_dependencies("weightedsurv", reverse = TRUE, which = "all")`); record the list.
2. Locate the forestsearch source: the local clone on this machine at its current HEAD (record path, branch, commit, and whether the working tree is clean). If there is no local clone, G5 fails.
3. Build the forestsearch tarball in a temporary directory outside both repos. Run `R CMD check` on it twice in separate temporary libraries: (a) with weightedsurv 0.1.0 from CRAN, (b) with the weightedsurv 0.2.0 tarball of G3.
4. Compare (a) and (b): errors, warnings, notes, test results (counts and names of failing tests), example and vignette failures.
5. List every weightedsurv function forestsearch calls (`grep` of its `R/`, `NAMESPACE`, vignettes and tests), and for each say whether Stage 2c–2e changed its results, defaults or arguments.

## 5. Site build

1. `pkgdown::check_pkgdown()`.
2. `pkgdown::build_site()` once, with the simulations article at full size as written (approved compute). Record wall-clock time and the article's own completion line.
3. Record for each article: rendered or not, and render time where available.
4. Commit `docs/`.

If the site build fails for a reason unrelated to the simulations article, stop (G6); do not rerun the simulations.

## 6. Gates (stop-on-failure; record PASS/FAIL with evidence)

| # | Check |
|---|---|
| G1 | `git diff --name-only e95728a..HEAD` lists nothing outside the allowed paths of section 0. `git diff e95728a..HEAD -- R/ man/ tests/ NAMESPACE` is empty. `git diff e95728a..HEAD -- DESCRIPTION` changes only `Version` and `Date`. |
| G2 | Rider is text only: in the two Rmd files, every changed line is prose, an inline-R expression, or the one chunk header of R5. The knitted objects of the main vignette's analysis chunks are `identical()` before and after the rider. |
| G3 | `R CMD check --as-cran` on the tarball built from a clean `git clone`: 0 errors, 0 warnings; every NOTE is listed in `cran-comments.md`. Tests: `FAIL 0`, `PASS 1691`. |
| G4 | `DESCRIPTION` `Version` is `0.2.0`; the first heading of `NEWS.md` is `# weightedsurv 0.2.0`; `cran-comments.md` names 0.2.0. Every exported function and argument added since 0.1.0 (`diff` of `NAMESPACE` and of the argument lists) appears in `NEWS.md`. |
| G5 | forestsearch: check (b) has no error, warning, note or failing test that check (a) does not have. |
| G6 | `pkgdown::check_pkgdown()`: no problems. `build_site()` completes. Every article in the `_pkgdown.yml` index has an HTML file under `docs/articles/`, including `estimated_weights_methodology`. The simulations article reports 30,000 completed, 0 failed. |
| G7 | Rider greps: `grep -rn "supporting the validity\|bootstrap or similar\|messages = FALSE" vignettes/` returns nothing; `README.md` contains `ps_weights(`; a script shows the README dependency lists equal `DESCRIPTION` Imports and Suggests. |
| G8 | No hard-coded result numerals in the changed prose (each numeral inline-computed, a design constant, or a literature value). |
| G9 | `urlchecker::url_check()`: no broken URL. If urlchecker is not installed, record NOT RUN (not a failure). |
| G10 | Tarball contents: no `dev/`, `docs/`, `vignettes/articles/`, `_pkgdown.yml`, `cran-comments.md`. Tarball size recorded; under 5 MB. |

## 7. Report

`dev/audit/REPORT_stage3_release_prep.md`:

1. Gate table.
2. Rider: passage, what it said, what it says now (R1a–R5).
3. `NEWS.md` 0.2.0 section in full, with the source (report or commit) of each line.
4. `R CMD check` status and NOTEs; environment.
5. forestsearch: one comparative table, check (a) vs check (b); the table of weightedsurv functions it calls and whether each changed.
6. Site build: time, per-article status, the simulations article's completion line.
7. State of the branch: HEAD, commits since `e95728a`, `git status` clean or not, commits ahead of the remote.
8. Anything found that this task does not cover: listed, not changed.

After the final commit, copy to `~/Downloads/` and print the paths: the report, `cran-comments.md`, `NEWS.md`, and the built tarball `weightedsurv_0.2.0.tar.gz`. Do not push.
