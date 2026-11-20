#!/bin/bash
# Gate G2 post-conditions (TASK_km_oc_evaluation.md, section 13). Run from the repo root.
START_SHA=d9b7558713732ef18bd2bcb6c400e057d5ff0391
MASTER_SHA=166919ec7128f165fb380847f549fa26ca1db985
O=dev/audit/oc
[ "$(git branch --show-current)" = "km-fixes-0.2.0" ] && echo "OK 1a" || echo "FAIL: wrong branch"
[ "$(git rev-parse master)" = "$MASTER_SHA" ] && echo "OK 1b" || echo "FAIL: master moved"
git diff --name-only "$START_SHA" HEAD | grep -Ev '^dev/(tasks|audit)/' && echo "FAIL: paths outside dev/tasks, dev/audit" || echo "OK 2"
for f in oc_sim.R oc_boot.R oc_summarise.R oc_summary.csv oc_profiles.pdf oc_stdout.txt .gitignore; do
  [ -f "$O/$f" ] || echo "FAIL: missing $O/$f"; done
[ -f dev/audit/REPORT_km_oc.md ] && [ -f dev/tasks/TASK_km_oc_evaluation.md ] && echo "OK 3 (if no FAIL above)" || echo "FAIL: report or task doc missing"
Rscript -e 'd <- read.csv("dev/audit/oc/oc_summary.csv")
  m <- d[d$quantity == "delta" & !is.na(d$time) & d$estimator %in% c("cran", "proposed"), ]
  stopifnot(length(unique(d$cell_id)) == 19, nrow(m) == 19 * 11 * 2, all(m$n_rep == 2000),
            all(is.finite(m$coverage[m$estimator == "proposed" & m$time <= 42])),
            all(is.finite(m$bias[m$estimator == "proposed" & m$time <= 42])))
  cat("OK 4\n")'
Rscript $O/oc_sim.R --version=proposed --pkg=. --check > /dev/null 2>&1 && echo "OK 5" || echo "FAIL: not reproducible"
git worktree list | grep -q oc_cran && echo "FAIL: 0.1.0 worktree not removed" || echo "OK 6"
[ -z "$(git status --porcelain)" ] && echo "OK 7" || echo "FAIL: uncommitted changes"
