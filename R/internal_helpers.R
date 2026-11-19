# Internal helpers (not exported).

# Evaluate `expr` after set.seed(seed) and then restore the caller's RNG state:
# .Random.seed is put back exactly as it was, or removed again if it did not exist.
# `expr` is evaluated lazily in the caller's frame, after the seed is set.
#' @noRd
ws_with_seed <- function(seed, expr) {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE) else NULL
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  }, add = TRUE)
  set.seed(seed)
  expr
}
