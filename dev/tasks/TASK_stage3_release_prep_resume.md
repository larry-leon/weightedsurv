# TASK: Stage 3 release prep, resume — weightedsurv 0.2.0, branch `km-fixes-0.2.0`

Written 2026-11-20. Executor: Claude Code (CC). Amends `dev/tasks/TASK_stage3_release_prep.md`, which stopped at G7 after step 2 (`dev/audit/REPORT_stage3_release_prep.md`).

## 0. Rules

- Branch `km-fixes-0.2.0`; HEAD must contain `e19df63`. If not, stop.
- First action: copy this file to `dev/tasks/TASK_stage3_release_prep_resume.md` and commit it.
- The original task document stays in force (rules, sections 3 to 7, approved compute) **except as amended below**. Where the two differ, this document wins.
- Steps 1 and 2 of the original (task document, text rider) are done and committed. **Do not redo or re-verify them.** G2 and G8 stand as reported.
- Resume at step 3 of the original order of work.
- Never push. Never merge. No tag. No submission. No site deployment (`deploy_to_branch()` pushes).
- Gates are stop-on-failure. Report findings only: no recommendation, no follow-up task proposals.

## 1. Why G7 failed

The G7 grep ran over all of `vignettes/` and matched `vignettes/.Rhistory`, a stray tracked file. No `.Rmd` matched. The gate was written too broadly; the rider is correct as committed.

## 2. Amendments

**A. `vignettes/.Rhistory`.** Remove it from the repository (`git rm vignettes/.Rhistory`) and add `.Rhistory` to `.gitignore`. Record any other tracked `.Rhistory`, `.RData` or `.DS_Store` file (`git ls-files`); remove those under `vignettes/` the same way, list the others without changing them.

**B. Site.** `docs/` is in `.gitignore` and the site is deployed from `gh-pages`. So: build the site locally with `pkgdown::build_site()` (once, with the full-size simulations article, as approved), verify it, and **do not commit `docs/`**. Do not deploy. Record that the README link to the estimated-weights article resolves only after the site is deployed, and confirm that the built site has the page at the path the README link uses (`articles/articles/estimated_weights_methodology.html`).

**C. `DESCRIPTION`.** Already `Version: 0.2.0`, no `Date` field. Make no change to `DESCRIPTION`.

**D. `NEWS.md`.** A complete-looking `# weightedsurv 0.2.0` section exists. Do not rewrite or reorder existing entries. Do two things only: (1) run the G4 comparison and add a line for any export or argument added since 0.1.0 that the section does not name; (2) add a `## Documentation` subsection with one entry: the vignettes and articles were revised for 0.2.0 (the Rotterdam example now uses postmenopausal women without chemotherapy and `ps_weights()` with `weight_variance = "estimated"`; effective sample size tables added; a new article gives the statistical details of estimated weights).

**E. `cran-comments.md`.** Does not exist: create it as the original section 3 describes. Confirm it is listed in `.Rbuildignore`; add it there if not.

**F. Allowed paths, replacing the original list:** `NEWS.md`, `cran-comments.md`, `.Rbuildignore` (only for E or G10), `.gitignore` (the `.Rhistory` line only), `vignettes/.Rhistory` (deletion only), `dev/`. Nothing else changes in this resume: no further edits to `README.md` or any `.Rmd`.

## 3. Amended gates

Gates not listed here are as in the original.

| # | Check |
|---|---|
| G1 | `git diff --name-only e95728a..HEAD` lists only: the rider files already committed (`README.md`, `vignettes/weightedsurv_examples.Rmd`, `vignettes/articles/martingale_resampling_applications.Rmd`), and the paths of amendment F. `git diff e95728a..HEAD -- R/ man/ tests/ NAMESPACE DESCRIPTION` is empty. No file under `docs/` is tracked. |
| G2, G8 | As reported; not rerun. |
| G4 | `DESCRIPTION` `Version` is `0.2.0`; the first heading of `NEWS.md` is `# weightedsurv 0.2.0`; `cran-comments.md` exists and names 0.2.0. Every exported function and argument added since 0.1.0 (`diff` of `NAMESPACE` and of the argument lists against the `v0.1.0` tag) appears in `NEWS.md`. |
| G6 | `pkgdown::check_pkgdown()`: no problems. `build_site()` completes locally. Every article in the `_pkgdown.yml` index has an HTML file in the built site, including `estimated_weights_methodology` at the path of the README link. The simulations article reports 30,000 completed, 0 failed. `git status` shows nothing from `docs/`. |
| G7 | `grep -rn --include="*.Rmd" "supporting the validity\|bootstrap or similar\|messages = FALSE" vignettes/` returns nothing. `git ls-files vignettes/ \| grep -i "rhistory\|rdata\|ds_store"` returns nothing. The README parts of G7 stand as reported. |

## 4. Report

Append to `dev/audit/REPORT_stage3_release_prep.md` a section "Resume" with: the full gate table (G1 to G10, final state), and sections 3 to 8 of the original report outline completed. Keep the first part of the report as the record of the stop.

After the final commit, copy to `~/Downloads/` and print the paths: the report, `cran-comments.md`, `NEWS.md`, and the built tarball `weightedsurv_0.2.0.tar.gz`. Do not push.
