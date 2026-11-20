## G4: exports and arguments added since the v0.1.0 tag, and whether NEWS.md names each.
## Usage: Rscript dev/audit/stage3/g4_api_diff.R <dir with R/ and NAMESPACE of v0.1.0> <repo root>
args <- commandArgs(trailingOnly = TRUE)
old_dir <- args[1]; new_dir <- args[2]
load_formals <- function(dir) {
  env <- new.env()
  for (f in list.files(file.path(dir, "R"), pattern = "\\.[Rr]$", full.names = TRUE)) {
    for (e in parse(f, keep.source = FALSE)) {
      if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=") && is.name(e[[2]]) &&
          is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
        assign(as.character(e[[2]]), e[[3]][[2]], envir = env)
      }
    }
  }
  env
}
exports <- function(dir) {
  dir <- normalizePath(dir); ns <- parseNamespaceFile(basename(dir), dirname(dir), mustExist = TRUE)
  list(exports = sort(ns$exports), S3 = if (nrow(ns$S3methods)) sort(paste(ns$S3methods[, 1], ns$S3methods[, 2], sep = ".")) else character())
}
eo <- exports(old_dir); en <- exports(new_dir)
fo <- load_formals(old_dir); fn <- load_formals(new_dir)
news <- paste(readLines(file.path(new_dir, "NEWS.md")), collapse = "\n")
in_news <- function(x) grepl(x, news, fixed = TRUE)

cat("Exports at v0.1.0:", length(eo$exports), " at HEAD:", length(en$exports), "\n")
cat("Exports removed:", paste(setdiff(eo$exports, en$exports), collapse = ", "), "\n")
added <- setdiff(en$exports, eo$exports)
cat("Exports added:\n"); for (a in added) cat("  ", a, "  in NEWS.md:", in_news(a), "\n")
cat("S3 methods added:\n"); for (a in setdiff(en$S3, eo$S3)) cat("  ", a, "  in NEWS.md:", in_news(a), "\n")

cat("\nArguments of exported functions (present at both):\n")
res <- list()
for (f in intersect(eo$exports, en$exports)) {
  if (!exists(f, fo, inherits = FALSE) || !exists(f, fn, inherits = FALSE)) { cat("  ", f, ": definition not found by the parser\n"); next }
  a0 <- get(f, fo); a1 <- get(f, fn)
  n0 <- names(a0); n1 <- names(a1)
  for (a in setdiff(n1, n0)) res[[length(res) + 1]] <- data.frame(fun = f, change = "argument added", arg = a,
      detail = paste(deparse(a1[[a]]), collapse = ""), arg_in_NEWS = in_news(a))
  for (a in setdiff(n0, n1)) res[[length(res) + 1]] <- data.frame(fun = f, change = "argument removed", arg = a,
      detail = "", arg_in_NEWS = in_news(a))
  for (a in intersect(n0, n1)) {
    d0 <- paste(deparse(a0[[a]]), collapse = ""); d1 <- paste(deparse(a1[[a]]), collapse = "")
    if (!identical(d0, d1)) res[[length(res) + 1]] <- data.frame(fun = f, change = "default changed", arg = a,
        detail = paste(d0, "->", d1), arg_in_NEWS = in_news(a))
  }
}
res <- do.call(rbind, res)
print(res, row.names = FALSE, right = FALSE)
in_code <- function(x) grepl(paste0("`[^`]*\\b", x, "\\b[^`]*`"), news, perl = TRUE)
cat("\nArguments of the added exports (named inside a code span of NEWS.md):\n")
miss <- character()
for (a in added) for (g in names(get(a, fn))) {
  ok <- in_code(g); cat("  ", a, "(", g, ")  in NEWS.md:", ok, "\n"); if (!ok) miss <- c(miss, paste0(a, "(", g, ")"))
}
cat("\nAdded exports, S3 methods or arguments not named in NEWS.md:\n")
s3_added <- setdiff(en$S3, eo$S3)
miss <- c(added[!vapply(added, in_news, TRUE)], s3_added[!vapply(s3_added, in_news, logical(1))], miss)
if (!is.null(res)) { r <- res[res$change == "argument added" & !res$arg_in_NEWS, ]
  if (nrow(r)) miss <- c(miss, paste0(r$fun, "(", r$arg, ")")) }
if (length(miss)) cat(paste0("  ", miss, "\n"), sep = "") else cat("  none\n")
cat("\nDESCRIPTION Version:", read.dcf(file.path(new_dir, "DESCRIPTION"))[, "Version"], "\n")
cat("First heading of NEWS.md:", grep("^#", readLines(file.path(new_dir, "NEWS.md")), value = TRUE)[1], "\n")
cc <- file.path(new_dir, "cran-comments.md")
cat("cran-comments.md exists:", file.exists(cc), " names 0.2.0:", file.exists(cc) && any(grepl("0.2.0", readLines(cc), fixed = TRUE)), "\n")
