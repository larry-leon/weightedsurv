# G2: compare the objects of two knits of the main vignette, object by object.
# usage: Rscript g2_compare.R a.rds b.rds
# identical() is applied twice: to the objects as saved, and to the objects with every environment
# removed (an environment read back from two sessions is never identical() by reference):
# closures and language objects are replaced by their deparse, environments by NULL.
a <- commandArgs(TRUE)
b <- readRDS(a[1]); f <- readRDS(a[2])
strip <- function(x) {
  if (is.environment(x)) return(NULL)
  if (is.function(x) || is.language(x)) return(deparse(x))
  at <- attributes(x)
  if (is.list(x)) x <- lapply(unclass(x), strip)
  if (!is.null(at)) attributes(x) <- lapply(at, strip)
  x
}
cat("objects in a:", length(b), " in b:", length(f), " same names:", setequal(names(b), names(f)), "\n")
res <- data.frame(object = names(b), class = vapply(b, function(x) class(x)[1], ""),
                  identical = vapply(names(b), function(n) identical(b[[n]], f[[n]]), NA),
                  identical_no_env = vapply(names(b), function(n) identical(strip(b[[n]]), strip(f[[n]])), NA),
                  row.names = NULL)
print(res, row.names = FALSE)
cat("\nidentical() as saved:", sum(res$identical), "of", nrow(res), "\n")
cat("not identical() as saved:", paste(res$object[!res$identical], collapse = ", "), "\n")
cat("identical() with environments removed:", sum(res$identical_no_env), "of", nrow(res), "\n")
cat("not identical() with environments removed:", paste(res$object[!res$identical_no_env], collapse = ", "), "\n")
