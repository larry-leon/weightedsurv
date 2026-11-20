# REPORT: one truth convention for the RMST rows (read-off truth R(t_j)) — weightedsurv 0.2.0

Task: `dev/tasks/TASK_rmst_readoff_truth.md`. Branch `release-0.2.0-robust`, start HEAD `e48874f` (contains `8107811` and the
report commit `e48874f` after it). Executed 2026-11-20. The read-off truth R(t_j) was used for every RMST row; the fallback
(section 2.4 of the task) was not needed.

Scripts and records: `dev/audit/rmst_readoff/` (`01_recover_tj.R` / `.txt`, `02_recompute.R` / `.txt`, `03_report_table.R`,
`report_table.md`, `rmst_rows_old_new.csv`, `tj.rds`, `recompute.rds`, `render_warn2.R` / `.txt`, `check_pkgdown.txt`,
`build_article.txt`). No package function was called on any replicate; no simulation and no bootstrap were run.

## 1. Gates

| # | Check | Result |
|---|---|---|
| G1 | `git diff --name-only e48874f..HEAD`: only `dev/tasks/TASK_rmst_readoff_truth.md`, `dev/audit/REPORT_rmst_readoff_truth.md`, `dev/audit/rmst_readoff/*` (13 files), `vignettes/articles/estimated_weights_methodology.Rmd`, `vignettes/articles/sim_results/ps_weight_results_by_time.csv`. `git diff e48874f..HEAD -- . ':(exclude)vignettes/articles' ':(exclude)dev'`: empty (0 bytes). Nothing under `dev/audit/stage2e/`, `dev/audit/oc/`, `dev/audit/ps_by_time/` changed (0 files). `.Rbuildignore` holds `^dev$` and `^vignettes/articles$`. | PASS |
| G2 | Route recorded for every cell (section 2): a = OC 5, 8, 10; b = 41, 42, 31, 32, 35, 37, 38, 39. Route b: regenerated data checksum equal to the stored `$chk` in 16,000 of 16,000 replicates (8 cells × 2,000). The OC cells were also regenerated as a check: 6,000 of 6,000 equal. | PASS |
| G3 | OC cells 5 and 10, 22 rows (all 11 months): recomputed bias, `mean_estimate − mean_truth_readoff`, bias MCSE, coverage and coverage MCSE against `dev/audit/oc/oc_summary.csv` (`proposed`) and the stored CSV: 176 comparisons, largest difference 6.1e-16. | PASS |
| G4 | OC cell 8 against M7 cell 31, RMST month 36, after recomputation: `mean_estimate` 0, `bias` 1.8e-16, `bias_mcse` 2.1e-17, `coverage_known` 0, `coverage_estimated` 0. Values: 3.96267744035892, −0.0122641185681805, 0.0297311700875105, 0.98, 0.947. | PASS |
| G5 | 111 RMST rows, 1,443 values: `mean_estimate` against the mean of the stored per-replicate estimates over the evaluable replicates; `empirical_sd`, `mean_se_known`, `mean_se_estimated`, `mean_ci_length_known`, `mean_ci_length_estimated`, counts, `true_value` and the bootstrap-replicate lengths against the stored CSV: 0 differ beyond 1e-12. | PASS |
| G6 | 370 non-RMST rows (481 − 111): 370 byte-identical to before apart from the added empty field `mean_truth_readoff`, checked by script (`02_recompute.R`). In the 111 RMST rows the 26 columns not recomputed keep their text. | PASS |
| G7 | Render with `options(warn = 2)`: OK in 1.6 s (`render_warn2.txt`; weightedsurv 0.2.0 of this branch installed into a scratch library, since the default library holds 0.1.0). `pkgdown::check_pkgdown()`: "No problems found." `pkgdown::build_article("articles/estimated_weights_methodology")` into a scratch destination: completed in 4.8 s. The design paragraph and the Table 2 note state the convention (checked in both HTML outputs). `docs/` not touched. | PASS |
| G8 | Numerals in changed prose (section 5): each is inline R from the CSV, a design constant (months 36 and 48; the benchmark month 36) or formula notation. | PASS |

## 2. Route per cell, grid and read-off rule

Grid and read-off rule, as coded:

- `dev/audit/oc/oc_sim.R:195-196` and `dev/audit/stage2e/oc_sim.R:182-183` (`run_option()`, whose definitions
  `dev/audit/ps_by_time/ps9_sim.R:30-35` evaluates for cells 41, 42): `fit <- KM_diff(df, ...)` with the default
  `at_points`; `jj <- last_le(fit$at_points, tt)`.
- `oc_common.R:109`: `last_le(grid, t) = findInterval(t, grid)`, NA when 0.
- `R/km_wlr_calculations_helpers.R:567` (the same line at `d9b7558:496`, `c195953:567`, `f55001f:567`, the package
  commits of the OC, M7 and new runs): `at_points = sort(df[[tte.name]])`, every observed time of the data, events and
  censorings. The grid is not changed later (`modify_tau = FALSE`).
- `R/rmst_calculations.R`, `cumulative_rmst_bands()`: `at_points <- fit$at_points`, `rmst_time` on that grid, so
  `rm$rmst_time[jj]` (`oc_sim.R:89`, `stage2e/oc_sim.R:91`) is the RMST at t_j = `sort(df$time)[jj]`.

(The task text says the default grid is "the event times"; the code uses all observed times.)

| Cell | k | Route | Source of t_j / R(t_j) | Replicates |
|---|---|---|---|---|
| D-400 | 5 | a | `R.truth.<t>` in `dev/audit/oc/raw/proposed_cell05.rds` (`oc_sim.R:210-211`: `R_true(fit$at_points, k)[jj]`) | 2,000 |
| E-400 | 8 | a | `R.truth.<t>` in `proposed_cell08.rds` | 2,000 |
| F-400 | 10 | a | `R.truth.<t>` in `proposed_cell10.rds` | 2,000 |
| Es9-400, Es9-400-null | 41, 42 | b | data from seed 1e6·k + r with `gen9()` (verbatim copy from `ps9_sim.R`, a = −0.44, b = 0.51); checksum `sum(time) + sum(status) + sum(X)` | 2,000 each |
| Es-400, Es-400-null | 31, 32 | b | `m7_data()` of `stage2e/oc_common.R` (data of OC cells 8 and 15, seed 1e6·k + r of those cells); checksum includes `sum(w)` | 2,000 each |
| Fe-400, Fe-400-null | 35, 37 | b | `m7_data()`, `gen(400, ·, "strong")`, seed 1e6·k + r | 2,000 each |
| V-1450-null, V-1450-ph | 38, 39 | b | `m7_data()`, `gen_v(1450, ·)`, seed 1e6·k + r | 2,000 each |

Truth at t_j: M7 and new cells `R_sc()` (the closed form of `stage2e/oc_common.R`, the function of their summaries); OC
cells the stored `R.truth` (`R_true()`, the evaluation's integral). Null cells: R ≡ 0, so R(t_j) = R(t) = 0.

Checks of the rule (`01_recover_tj.txt`): in the OC cells the rule applied to the regenerated data reproduces the stored
`R.truth` exactly (22,000 of 22,000 values identical per cell, largest difference 0). Cell 31 (data of OC cell 8): t_j
identical to cell 8's in every replicate and month; `R_sc(t_j)` against cell 8's stored `R.truth`: largest difference
8.9e-15. t_j exists at every month in every replicate. No evaluable replicate lacks a finite estimate; no failed RMST
call and no non-finite interval among evaluable replicates.

## 3. Comparative table

RMST rows, months 1, 6, 12, 24, 36, 48 (all months in `rmst_rows_old_new.csv`). "Stored CSV was against": the truth of
the bias and coverage in the file before this task. Mean R(t) − R(t_j): over the replicates in the bias average.
"Bias vs R(t_j)" and "Cov … R(t_j)" are the new stored values.

| k | Cell | Month | Stored CSV was against | Mean R(t) − R(t_j) | Bias vs R(t) | Bias vs R(t_j) | Cov known R(t) | Cov known R(t_j) | Cov estimated R(t) | Cov estimated R(t_j) |
|---|---|---|---|---|---|---|---|---|---|---|
| 41 | Es9-400 | 1 | R(t) | 0.00066 | -0.00024 | 0.00042 | 0.9425 | 0.9480 | 0.9360 | 0.9415 |
| 41 | Es9-400 | 6 | R(t) | 0.00407 | -0.00080 | 0.00327 | 0.9725 | 0.9730 | 0.9540 | 0.9540 |
| 41 | Es9-400 | 12 | R(t) | 0.00903 | -0.00720 | 0.00183 | 0.9750 | 0.9750 | 0.9520 | 0.9525 |
| 41 | Es9-400 | 24 | R(t) | 0.02144 | -0.02429 | -0.00285 | 0.9750 | 0.9740 | 0.9405 | 0.9405 |
| 41 | Es9-400 | 36 | R(t) | 0.04103 | -0.03513 | 0.00590 | 0.9760 | 0.9755 | 0.9435 | 0.9435 |
| 41 | Es9-400 | 48 | R(t) | 0.06482 | -0.04558 | 0.01924 | 0.9755 | 0.9750 | 0.9445 | 0.9445 |
| 42 | Es9-400-null | 1 | R(t) | 0.00000 | -0.00050 | -0.00050 | 0.9490 | 0.9490 | 0.9425 | 0.9425 |
| 42 | Es9-400-null | 6 | R(t) | 0.00000 | -0.00277 | -0.00277 | 0.9695 | 0.9695 | 0.9475 | 0.9475 |
| 42 | Es9-400-null | 12 | R(t) | 0.00000 | -0.00266 | -0.00266 | 0.9745 | 0.9745 | 0.9500 | 0.9500 |
| 42 | Es9-400-null | 24 | R(t) | 0.00000 | 0.01076 | 0.01076 | 0.9835 | 0.9835 | 0.9515 | 0.9515 |
| 42 | Es9-400-null | 36 | R(t) | 0.00000 | 0.00253 | 0.00253 | 0.9775 | 0.9775 | 0.9445 | 0.9445 |
| 42 | Es9-400-null | 48 | R(t) | 0.00000 | -0.02164 | -0.02164 | 0.9760 | 0.9760 | 0.9409 | 0.9409 |
| 31 | Es-400 | 1 | R(t) | 0.00070 | -0.00060 | 0.00010 | 0.8825 | 0.8935 | 0.8720 | 0.8840 |
| 31 | Es-400 | 6 | R(t) | 0.00435 | -0.00154 | 0.00282 | 0.9605 | 0.9605 | 0.9380 | 0.9385 |
| 31 | Es-400 | 12 | R(t) | 0.00887 | -0.00409 | 0.00478 | 0.9785 | 0.9780 | 0.9515 | 0.9540 |
| 31 | Es-400 | 24 | R(t) | 0.02181 | -0.01670 | 0.00512 | 0.9805 | 0.9800 | 0.9520 | 0.9540 |
| 31 | Es-400 | 36 | R(t) | 0.03857 | -0.05084 | -0.01226 | 0.9805 | 0.9800 | 0.9460 | 0.9470 |
| 31 | Es-400 | 48 | R(t) | 0.06842 | -0.13222 | -0.06380 | 0.9750 | 0.9760 | 0.9415 | 0.9415 |
| 32 | Es-400-null | 1 | R(t) | 0.00000 | 0.00014 | 0.00014 | 0.9195 | 0.9195 | 0.9115 | 0.9115 |
| 32 | Es-400-null | 6 | R(t) | 0.00000 | -0.00795 | -0.00795 | 0.9585 | 0.9585 | 0.9335 | 0.9335 |
| 32 | Es-400-null | 12 | R(t) | 0.00000 | -0.02339 | -0.02339 | 0.9765 | 0.9765 | 0.9370 | 0.9370 |
| 32 | Es-400-null | 24 | R(t) | 0.00000 | -0.04409 | -0.04409 | 0.9785 | 0.9785 | 0.9475 | 0.9475 |
| 32 | Es-400-null | 36 | R(t) | 0.00000 | -0.07573 | -0.07573 | 0.9755 | 0.9755 | 0.9445 | 0.9445 |
| 32 | Es-400-null | 48 | R(t) | 0.00000 | -0.11037 | -0.11037 | 0.9686 | 0.9686 | 0.9392 | 0.9392 |
| 35 | Fe-400 | 1 | R(t) | 0.00068 | -0.00137 | -0.00069 | 0.6640 | 0.6860 | 0.6575 | 0.6780 |
| 35 | Fe-400 | 6 | R(t) | 0.00430 | -0.02546 | -0.02116 | 0.8890 | 0.8925 | 0.8540 | 0.8625 |
| 35 | Fe-400 | 12 | R(t) | 0.00914 | -0.08106 | -0.07192 | 0.9445 | 0.9450 | 0.9060 | 0.9100 |
| 35 | Fe-400 | 24 | R(t) | 0.02071 | -0.22227 | -0.20156 | 0.9595 | 0.9610 | 0.9185 | 0.9190 |
| 35 | Fe-400 | 36 | R(t) | 0.03737 | -0.40368 | -0.36631 | 0.9570 | 0.9570 | 0.9160 | 0.9155 |
| 35 | Fe-400 | 48 | R(t) | 0.06774 | -0.64226 | -0.57452 | 0.9345 | 0.9375 | 0.9010 | 0.9080 |
| 37 | Fe-400-null | 1 | R(t) | 0.00000 | -0.00182 | -0.00182 | 0.7360 | 0.7360 | 0.7160 | 0.7160 |
| 37 | Fe-400-null | 6 | R(t) | 0.00000 | -0.03043 | -0.03043 | 0.9050 | 0.9050 | 0.8635 | 0.8635 |
| 37 | Fe-400-null | 12 | R(t) | 0.00000 | -0.07929 | -0.07929 | 0.9305 | 0.9305 | 0.8920 | 0.8920 |
| 37 | Fe-400-null | 24 | R(t) | 0.00000 | -0.19521 | -0.19521 | 0.9365 | 0.9365 | 0.8895 | 0.8895 |
| 37 | Fe-400-null | 36 | R(t) | 0.00000 | -0.31998 | -0.31998 | 0.9265 | 0.9265 | 0.8834 | 0.8834 |
| 37 | Fe-400-null | 48 | R(t) | 0.00000 | -0.39370 | -0.39370 | 0.9104 | 0.9104 | 0.8729 | 0.8729 |
| 39 | V-1450-ph | 1 | R(t) | 0.00018 | -0.00028 | -0.00009 | 0.9615 | 0.9605 | 0.9405 | 0.9425 |
| 39 | V-1450-ph | 6 | R(t) | 0.00113 | -0.00136 | -0.00023 | 0.9740 | 0.9730 | 0.9440 | 0.9420 |
| 39 | V-1450-ph | 12 | R(t) | 0.00251 | -0.00874 | -0.00623 | 0.9700 | 0.9700 | 0.9410 | 0.9410 |
| 39 | V-1450-ph | 24 | R(t) | 0.00608 | -0.05292 | -0.04684 | 0.9600 | 0.9610 | 0.9310 | 0.9300 |
| 39 | V-1450-ph | 36 | R(t) | 0.01169 | -0.14057 | -0.12888 | 0.9455 | 0.9460 | 0.9260 | 0.9265 |
| 39 | V-1450-ph | 48 | R(t) | 0.02265 | -0.25586 | -0.23321 | 0.9310 | 0.9310 | 0.9120 | 0.9130 |
| 38 | V-1450-null | 1 | R(t) | 0.00000 | -0.00006 | -0.00006 | 0.9680 | 0.9680 | 0.9535 | 0.9535 |
| 38 | V-1450-null | 6 | R(t) | 0.00000 | -0.00331 | -0.00331 | 0.9745 | 0.9745 | 0.9385 | 0.9385 |
| 38 | V-1450-null | 12 | R(t) | 0.00000 | -0.01904 | -0.01904 | 0.9640 | 0.9640 | 0.9375 | 0.9375 |
| 38 | V-1450-null | 24 | R(t) | 0.00000 | -0.07267 | -0.07267 | 0.9485 | 0.9485 | 0.9225 | 0.9225 |
| 38 | V-1450-null | 36 | R(t) | 0.00000 | -0.15826 | -0.15826 | 0.9310 | 0.9310 | 0.9105 | 0.9105 |
| 38 | V-1450-null | 48 | R(t) | 0.00000 | -0.17336 | -0.17336 | 0.9143 | 0.9143 | 0.8987 | 0.8987 |
| 5 | D-400 | 1 | R(t_j) | 0.00067 | -0.00055 | 0.00012 | 0.8825 | 0.8965 | — | — |
| 5 | D-400 | 6 | R(t_j) | 0.00438 | -0.00588 | -0.00150 | 0.9370 | 0.9400 | — | — |
| 5 | D-400 | 12 | R(t_j) | 0.00890 | -0.01582 | -0.00692 | 0.9465 | 0.9465 | — | — |
| 5 | D-400 | 24 | R(t_j) | 0.02199 | -0.04162 | -0.01963 | 0.9500 | 0.9490 | — | — |
| 5 | D-400 | 36 | R(t_j) | 0.03864 | -0.07641 | -0.03777 | 0.9515 | 0.9520 | — | — |
| 5 | D-400 | 48 | R(t_j) | 0.06555 | -0.12707 | -0.06152 | 0.9475 | 0.9480 | — | — |
| 10 | F-400 | 1 | R(t_j) | 0.00069 | -0.00218 | -0.00149 | 0.6760 | 0.6945 | — | — |
| 10 | F-400 | 6 | R(t_j) | 0.00435 | -0.02593 | -0.02158 | 0.8710 | 0.8760 | — | — |
| 10 | F-400 | 12 | R(t_j) | 0.00923 | -0.07050 | -0.06127 | 0.9090 | 0.9100 | — | — |
| 10 | F-400 | 24 | R(t_j) | 0.02057 | -0.20727 | -0.18670 | 0.9190 | 0.9190 | — | — |
| 10 | F-400 | 36 | R(t_j) | 0.03692 | -0.39978 | -0.36287 | 0.9035 | 0.9050 | — | — |
| 10 | F-400 | 48 | R(t_j) | 0.06763 | -0.61966 | -0.55203 | 0.8948 | 0.8978 | — | — |
| 8 | E-400 | 36 | R(t) | 0.03857 | -0.05084 | -0.01226 | 0.9805 | 0.9800 | 0.9460 | 0.9470 |

**Reading.** In every replicate the RMST "at month t" is the RMST at a slightly earlier time t_j, on average 0.01 to 0.6
months earlier (least in the N = 1,450 cells), more at later months. Since the true RMST difference is increasing in the PH cells, R(t_j) is below R(t)
by 0.0002 to 0.068 months on average; the bias against R(t) is that against R(t_j) minus this gap (to 7e-16). In the null
cells R ≡ 0 and nothing changes. Against R(t_j):

- Near 0.9 (Es9-400) and About 0.69 (Es-400): the negative bias against R(t) that grew with t (to −0.046 and −0.132 at
  month 48) becomes small and of either sign (−0.003 to +0.019 and −0.064 to +0.009 over all months); the largest
  |bias| / MCSE over the 11 months goes from 1.4 to 1.8 (Es9-400) and from 3.2 to 1.5 (Es-400). D-400 (true weights)
  shows the same: 2.6 to 1.3.
- Strong confounding (Fe-400, F-400) and 14% treated (V-1450-ph): the bias stays negative and grows with t (Fe-400
  −0.575 at month 48, 9.1 MCSE; F-400 −0.552; V-1450-ph −0.233); the gap accounts for about a tenth of it. The null
  cells of the same regimes, where there is no gap, show bias of the same size (Fe-400-null −0.394, V-1450-null −0.173
  at month 48).
- Coverage moves up by at most 0.022 ("known") and 0.021 ("estimated"), the largest changes at month 1 (Fe-400
  0.664 → 0.686, F-400 0.676 → 0.695, Es-400 0.883 → 0.894); at months 3 and 6 up to 0.011, at months 9 to 48 up to
  0.007 (Fe-400 "estimated", month 48); in some rows coverage moves down, by at most 0.003.
- OC cells 5 and 10 were already stored against R(t_j) and do not change except `mean_estimate` (now the mean of the
  estimates instead of R(t) + bias). The rows that change are the 44 PH rows of cells 41, 31, 35, 39 and the cell 8 row.

Bootstrap benchmark columns (month 36, first 500 replicates; `02_recompute.txt`):

| Cell | known_bootreps R(t) → R(t_j) | estimated_bootreps R(t) → R(t_j) | boot normal: stored (R(36)); against R(t_j) | boot percentile: stored (R(36)); against R(t_j) |
|---|---|---|---|---|
| 8 (E-400) | 0.976 → 0.974 | 0.934 → 0.940 | 0.942; 0.938 | 0.936; 0.934 |
| 39 (V-1450-ph) | 0.936 → 0.938 | 0.920 → 0.922 | 0.920; 0.918 | 0.922; 0.920 |
| 38 (V-1450-null) | 0.926 → 0.926 | 0.900 → 0.900 | 0.898; 0.898 | 0.906; 0.906 |

`coverage_known_bootreps` and `coverage_estimated_bootreps` were recomputed against R(t_j) from the per-replicate files
(cell 8: `proposed_rect` and `proposed_ew` with all 500 replicates counted, `oc_summarise.R:148`; cells 38, 39: `known`
and `estimated` on the replicates evaluable at month 36, `stage2e/oc_summarise.R:169-172`). The per-replicate checksums
of `boot_cell*.rds` equal those of the main files (all 500, all three cells); the bootstrap-replicate lengths are
unchanged. `coverage_boot`, `coverage_boot_mcse`, `coverage_boot_pct` and `boot_se_over_sd` were left unchanged: the
bootstrap estimate `est.R36` (`dev/audit/stage2e/oc_boot.R:21-27`, `dev/audit/oc/oc_boot.R`) is
`summary(survfit(...), rmean = 36)`, the integral of each Kaplan–Meier curve over [0, 36] itself, not a value read at
t_j, so its truth is R(36); the stored values are reproduced against R(36) to 1e-12. Their R(t_j) coverages are in the
table above for the record.

## 4. Band truth

The RMST band rows (`R_band_3_42`) use the truth at the band's own grid points:

- OC: `oc_sim.R:91`, `band_stats(..., rm$at_points, rm$rmst_time_sb_lower, rm$rmst_time_sb_upper, Rgrid, ...)` with
  `Rgrid <- R_true(fit$at_points, k)` (`oc_sim.R:210`), and `rm$at_points = fit$at_points`.
- M7 and new cells: `stage2e/oc_sim.R:93` with `R_sc(fit$at_points, cl$sc)` (`stage2e/oc_sim.R:190`).
- `band_stats()` (`oc_sim.R:38-47`, `stage2e/oc_sim.R:39-48`) compares the band and the truth point by point on that grid.

The truth is R at each grid point, not the nominal-month truth applied at a different time. Band rows: unchanged.

## 5. Article changes (`vignettes/articles/estimated_weights_methodology.Rmd`, section 7.1)

Rendered text before (HEAD `e48874f`) and after; from `pandoc -t plain` of the two renders.

1. **Design paragraph**, sentences added after "A replicate counts at month $t$ … at least $t$.":
   - Before: (none)
   - After: "The RMST difference at month t is read on the grid of the KM_diff() fit (the observed times) at the last grid
     time t_j ≤ t of the replicate, so its bias and coverage are against R(t_j), the true RMST difference at that time;
     the mean of R(t) − R(t_j) over replicates is at most 0.068 months (month 48). A Kaplan–Meier curve is constant
     between grid times, so the survival difference at t_j is that at t; the bootstrap RMST of Table 5 integrates each
     curve to month 36 itself and is scored against R(36)."
   - Numerals: 0.068 and 48 inline R (`max(main$true_value - main$mean_truth_readoff)`, its month); 36 the benchmark
     month (design constant).
2. **Table 2 note**, new line after Table 2:
   - Before: (none)
   - After: "Table 2: True is R(t) at the stated month; bias and coverage are against R(t_j), t_j the last grid time ≤ t
     of each replicate (design paragraph above)."
   - Table 2 code unchanged (columns selected by name; the new CSV column is not displayed). Displayed values that
     change: Bias and both coverages of the PH rows of Es9-400, Es-400, Fe-400, V-1450-ph (D-400 and F-400 unchanged).
3. **Reading, Near 0.9** (text unchanged; inline values recomputed):
   - Before: "RMST difference: 0.936–0.954 ("estimated") and 0.943–0.984 ("known")."
   - After: "RMST difference: 0.939–0.954 ("estimated") and 0.948–0.984 ("known")."
   - The RMST band sentence ("0.967–0.971 and 0.987 for RMST") is unchanged.
4. **Reading, Strong confounding** (rewritten):
   - Before: "For the RMST difference it accumulates: -0.40 months at month 36 and -0.64 at month 48 (true values 4.01
     and 5.66)."
   - After: "For the RMST difference, against R(t_j), it accumulates: -0.37 months at month 36 and -0.57 at month 48
     (8.0 and 9.1 MCSE from zero; true values 4.01 and 5.66)."
   - "Accumulates" is kept: against R(t_j) the Fe-400 bias still grows from −0.0007 (month 1) to −0.575 (month 48) and
     the read-off gap is 0.037 and 0.068 at months 36 and 48; the existing check `v("Fe-400", "R", "bias", 48) <
     v("Fe-400", "R", "bias", 36)` holds. Numerals: inline R from the CSV (bias, bias_mcse, true_value); 36 and 48 design months.
5. **Reading, Bootstrap (Table 5)** (text unchanged; inline values over the Table 5 rows, which include the RMST rows):
   - Before: "… differ by at most 0.008, and the "estimated" interval is 0.97–0.99 times the bootstrap's length.
     "known" covers 0.016–0.036 more than the bootstrap."
   - After: "… differ by at most 0.006, and the "estimated" interval is 0.97–0.99 times the bootstrap's length. "known"
     covers 0.016–0.036 more than the bootstrap."
   - Table 5 RMST rows: E-400 0.976 / 0.934 → 0.974 / 0.940; V-1450-ph 0.936 / 0.920 → 0.938 / 0.922 (known /
     estimated); bootstrap columns unchanged.

No other RMST sentence exists in section 7.1. No statement in the article attributes the read-off gap to the estimator
after the change. Nothing else in the article changed (`git diff` of the `.Rmd`: the design paragraph, the Table 2 note
line, the Fe-400 RMST sentence).

CSV (`ps_weight_results_by_time.csv`): column `mean_truth_readoff` added after `true_value` (empty in non-RMST rows).
RMST rows, columns whose text changed (rows): `mean_estimate` 22 (OC 5, 10); `bias` 45; `bias_mcse` 45;
`coverage_known` 38; `coverage_known_mcse` 38; `coverage_estimated` 35; `coverage_estimated_mcse` 35;
`coverage_known_bootreps` 2; `coverage_estimated_bootreps` 2. A recomputed value within 1e-12 of the stored one keeps
the stored text: the stored coverage MCSEs of copied rows had been computed from coverages read back from a CSV and
differ from a full-precision recomputation in the 16th digit (for example cell 32, month 42). `true_value` (the nominal
R(t)) and the `note` texts are unchanged. `bias_over_estimated_rows_if_different` stays empty: the bias under
"estimated" is identical to that under "known" in all 89 RMST rows with an "estimated" option.

## 6. Branch state

- Branch `release-0.2.0-robust`; start HEAD `e48874f`.
- Commits since start:
  - `22df707` rmst read-off truth: task file
  - `6640fbb` rmst read-off truth: t_j recovered …, RMST rows recomputed against R(t_j), column mean_truth_readoff; G2-G6
  - `89047f9` rmst read-off truth: article section 7.1 …; render with warn = 2, check_pkgdown, build_article (G7)
  - the commit of this report (`dev/audit/REPORT_rmst_readoff_truth.md`, `dev/audit/rmst_readoff/03_report_table.R`,
    `report_table.md`)
- `git status` after the report commit: clean.
- Not pushed, not merged, no tag, no submission, no site deployment; the tarball in `~/Downloads` was not rebuilt.
