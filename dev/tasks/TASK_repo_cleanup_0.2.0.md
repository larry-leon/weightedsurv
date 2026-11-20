# TASK: repository clean-up after 0.2.0 — merge `release-0.2.0-robust`, delete both release branches

Written 2026-11-20. Executor: Claude Code (CC). Requested by Larry ("ask CC to implement the repository clean up").

## 0. Rules

- **Pushing:** allowed in this task only if the kickoff contains the line `Push authorized for this task: merge to master and remote branch deletion.` Without that line, do sections 1–3 and stop before section 4; the report then says what is left for Larry.
- Only these remote operations are allowed: push `master`; delete the remote branches `km-fixes-0.2.0` and `release-0.2.0-robust`. **No force push. No tag. Do not touch `gh-pages` or any other branch. No submission anywhere.**
- Local branch deletion uses `git branch -d` only (never `-D`).
- Gates are stop-on-failure: if one fails, stop, write the report, do not ask. A failed gate before section 4 means nothing is pushed.
- Report findings only: no recommendation, no follow-up task proposals.

## 1. Preconditions

1. On `release-0.2.0-robust`, working tree clean (`git status --porcelain` empty). HEAD contains `89047f9` and the RMST report commit after it.
2. First action: copy this file to `dev/tasks/TASK_repo_cleanup_0.2.0.md` and commit it **on `release-0.2.0-robust`** (so it enters `master` with the merge).
3. `git fetch origin --prune`. Record `master`, `origin/master`, `km-fixes-0.2.0`, `origin/km-fixes-0.2.0`, `release-0.2.0-robust`, `origin/release-0.2.0-robust` (or "absent").

## 2. Merge

1. `git checkout master`. Local `master` must equal `origin/master`; if `origin/master` has commits that local `master` lacks, or the two have diverged, stop (G1).
2. `git merge --no-ff release-0.2.0-robust` with the message `Merge branch 'release-0.2.0-robust': release robustness and RMST read-off truth for weightedsurv 0.2.0`. If there is any conflict, `git merge --abort` and stop (G2).

## 3. Checks before pushing

- **Package content unchanged from the submitted tarball:** `git diff 8107811 master -- . ':(exclude)vignettes/articles' ':(exclude)dev'` is empty. (The tarball in `~/Downloads`, sha256 `a7f38c2e…`, was built from `8107811`; everything changed since is excluded from the package by `.Rbuildignore`.)
- **Nothing is lost by deleting the branches:**
  - `git merge-base --is-ancestor km-fixes-0.2.0 master` and, if it exists, `git merge-base --is-ancestor origin/km-fixes-0.2.0 master`: both true. (The remote branch received two pushes that no task made; this confirms it holds nothing `master` lacks.)
  - `git merge-base --is-ancestor release-0.2.0-robust master`: true; and the same for `origin/release-0.2.0-robust` if it exists.
- The merge commit has two parents.

## 4. Push and delete (only with the kickoff authorization line)

1. `git push origin master` (no `--force`). Confirm `origin/master` equals local `master`.
2. Delete the remote branches that exist: `git push origin --delete km-fixes-0.2.0` and `git push origin --delete release-0.2.0-robust`.
3. Delete the local branches: `git branch -d km-fixes-0.2.0 release-0.2.0-robust`.
4. Confirm `origin/gh-pages` is still present and unchanged by this task, and that no tag was created.

## 5. After the push (record; these are not stop gates)

The push to `master` starts the `R-CMD-check` workflow (Windows, macOS, Linux) and the `pkgdown` workflow, which deploys the site to `gh-pages`.

- If the `gh` CLI is available and authenticated, follow the two runs for this push (`gh run list --branch master`, `gh run watch <id>`; read-only), for at most 45 minutes. Record each job's result. For a failed job, save the failing step's log (`gh run view <id> --log-failed`) to `~/Downloads/actions_failed_<workflow>.log`. If `gh` is not available, record that.
- After a successful pkgdown deploy, request `https://larry-leon.github.io/weightedsurv/articles/articles/estimated_weights_methodology.html` and record the HTTP status (and where it redirects, if it does).

## 6. Gates

| # | Check |
|---|---|
| G1 | Preconditions of section 1 hold; local `master` equals `origin/master` before the merge. |
| G2 | Merge completed without conflict; merge commit has two parents. |
| G3 | `git diff 8107811 master -- . ':(exclude)vignettes/articles' ':(exclude)dev'` is empty. |
| G4 | All ancestry checks of section 3 are true. |
| G5 | (Authorized push only) `origin/master` equals local `master`; the remote branches are absent after `git fetch --prune`; the local branches are deleted with `-d`; `origin/gh-pages` present; `git tag` unchanged. |

## 7. Report

Write `~/Downloads/REPORT_repo_cleanup_0.2.0.md` (not committed: committing it would need another push). Contents: gate table; the recorded refs before and after; the merge commit id; the GitHub Actions results and the README-link status of section 5; `git branch -a` and `git status` at the end. Print its path.
