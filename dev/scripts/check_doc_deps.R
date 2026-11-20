## Packages used by the vignettes, the website articles and the tests, against DESCRIPTION.
## Run from the repository root:  Rscript dev/scripts/check_doc_deps.R
##
## Scans the R code of every .Rmd under vignettes/ (R chunks and inline `r ` code; comment lines
## are skipped) and every .R file under tests/ for library(), require(), requireNamespace() and pkg::.
## Prints, per package: where used (file:line), the DESCRIPTION field that declares it (Imports,
## Suggests, Config/Needs/website, base R, or the package itself), whether it is on CRAN
## (available.packages()), and whether each use is guarded. Then lists every use of a Suggests
## package (other than knitr, rmarkdown, testthat) in the CRAN-built vignette and in tests/, with its
## guard.
##
## A use is counted as guarded when
##   - in an Rmd: its chunk header has an `eval =` option that is not the literal TRUE (for example
##     `eval = has[["dplyr"]]`, `eval = run_live`, `eval = FALSE`), or the line itself is a
##     conditional (`if (...)`) or a requireNamespace() call, or the inline code is wrapped in if_has();
##   - in a test file: the test file calls skip_if_not_installed("<pkg>") before the use.
##
## Exit status: 1 if any package is undeclared, or if a package not on CRAN has an unguarded use;
## 2 if CRAN could not be reached; 0 otherwise.

d <- read.dcf("DESCRIPTION")
field <- function(f) {
  if (!f %in% colnames(d) || is.na(d[, f])) return(character(0))
  x <- trimws(strsplit(gsub("\n", " ", d[, f]), ",")[[1]])
  x <- gsub("\\s*\\(.*\\)", "", x)
  x[nzchar(x)]
}
decl <- list(Imports = field("Imports"), Suggests = field("Suggests"),
             `Config/Needs/website` = field("Config/Needs/website"))
self <- d[, "Package"]
base_pk <- c("base", "stats", "graphics", "grDevices", "utils", "methods", "datasets",
             "parallel", "tools", "grid", "splines", "stats4", "compiler", "tcltk")

pat_call <- "(?<=library\\()[A-Za-z][A-Za-z0-9.]*|(?<=require\\()[A-Za-z][A-Za-z0-9.]*|(?<=requireNamespace\\([\"'])[A-Za-z][A-Za-z0-9.]*"
pat_ns <- "(?<![A-Za-z0-9._$])[A-Za-z][A-Za-z0-9.]*(?=:::?[A-Za-z_.`])"
pkgs_in <- function(x) {
  p <- c(unlist(regmatches(x, gregexpr(pat_call, x, perl = TRUE))),
         unlist(regmatches(x, gregexpr(pat_ns, x, perl = TRUE))))
  unique(p)
}

uses <- list()
add_use <- function(file, line, pkg, guarded, how) {
  uses[[length(uses) + 1]] <<- data.frame(file = file, line = line, pkg = pkg, guarded = guarded,
                                          guard = how, stringsAsFactors = FALSE)
}

## R Markdown: R chunks and inline R
rmd <- list.files("vignettes", "\\.Rmd$", recursive = TRUE, full.names = TRUE)
for (f in rmd) {
  x <- readLines(f, warn = FALSE)
  in_chunk <- FALSE; header <- ""
  for (i in seq_along(x)) {
    l <- x[i]
    if (!in_chunk && grepl("^\\s*```\\s*\\{r", l)) { in_chunk <- TRUE; header <- l; next }
    if (in_chunk && grepl("^\\s*```\\s*$", l)) { in_chunk <- FALSE; header <- ""; next }
    if (in_chunk) {
      if (grepl("^\\s*#", l)) next
      code <- sub("\\s#[^'\"]*$", "", l)
      for (p in pkgs_in(code)) {
        ev <- regmatches(header, regexpr("eval\\s*=\\s*[^,}]+(\\([^)]*\\))?[^,}]*", header))
        hdr_guard <- length(ev) == 1 && !grepl("eval\\s*=\\s*(TRUE|T)\\s*$", ev)
        line_guard <- grepl("^\\s*if\\s*\\(", code) || grepl(sprintf("requireNamespace\\([\"']%s", p), code)
        how <- if (hdr_guard) paste0("chunk ", trimws(ev)) else if (line_guard) "conditional line" else ""
        add_use(f, i, p, hdr_guard || line_guard, how)
      }
    } else {
      inl <- regmatches(l, gregexpr("`r [^`]+`", l))[[1]]
      for (code in inl) for (p in pkgs_in(code)) {
        g <- grepl("if_has\\(", code)
        add_use(f, i, p, g, if (g) "inline if_has()" else "")
      }
    }
  }
}

## tests
tf <- list.files("tests", "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
for (f in tf) {
  x <- readLines(f, warn = FALSE)
  for (i in seq_along(x)) {
    l <- x[i]
    if (grepl("^\\s*#", l)) next
    for (p in pkgs_in(l)) {
      skips <- grep(sprintf("skip_if_not_installed\\([\"']%s[\"']", p), x[seq_len(i)])
      g <- length(skips) > 0 || grepl(sprintf("requireNamespace\\([\"']%s", p), l)
      add_use(f, i, p, g, if (length(skips)) sprintf("skip_if_not_installed() line %d", max(skips)) else "")
    }
  }
}

U <- do.call(rbind, uses)
declared_in <- function(p) {
  if (p == self) return("(this package)")
  if (p %in% base_pk) return("base R")
  hit <- names(decl)[vapply(decl, function(v) p %in% v, logical(1))]
  if (length(hit)) paste(hit, collapse = ", ") else "UNDECLARED"
}
ap <- tryCatch(rownames(utils::available.packages(repos = "https://cloud.r-project.org")),
               error = function(e) NULL)
if (is.null(ap) || length(ap) == 0) { cat("CRAN could not be reached (available.packages()).\n"); quit(status = 2) }
on_cran <- function(p) if (p %in% c(base_pk, "base")) "base R" else if (p %in% ap) "yes" else "NO"

pk <- sort(unique(U$pkg))
S <- data.frame(pkg = pk, declared = vapply(pk, declared_in, ""), on_cran = vapply(pk, on_cran, ""),
                uses = vapply(pk, function(p) sum(U$pkg == p), 0L),
                guarded = vapply(pk, function(p) sprintf("%d of %d", sum(U$guarded[U$pkg == p]), sum(U$pkg == p)), ""),
                row.names = NULL)
S$where <- vapply(pk, function(p) {
  u <- U[U$pkg == p, ]
  paste(sprintf("%s:%s", names(split(u$line, u$file)),
                vapply(split(u$line, u$file), function(l) paste(l, collapse = ","), "")), collapse = "; ")
}, "")

options(width = 250)
cat("== Packages used by vignettes/ (R chunks and inline R) and tests/\n\n")
for (i in seq_len(nrow(S))) {
  cat(sprintf("%-16s declared: %-22s on CRAN: %-6s guarded uses: %s\n    %s\n",
              S$pkg[i], S$declared[i], S$on_cran[i], S$guarded[i], S$where[i]))
}

cat("\n== Uses of Suggests packages (other than knitr, rmarkdown, testthat) in the CRAN-built vignette and in tests/\n\n")
cran_rmd <- list.files("vignettes", "\\.Rmd$", full.names = TRUE)  # top level only: articles/ is .Rbuildignore'd
sug <- setdiff(decl$Suggests, c("knitr", "rmarkdown", "testthat"))
G3 <- U[U$pkg %in% sug & (U$file %in% cran_rmd | startsWith(U$file, "tests")), ]
if (nrow(G3)) print(G3, row.names = FALSE) else cat("(none)\n")
g3_ok <- all(G3$guarded)
cat("\nG3: every such use guarded:", g3_ok, "\n")

und <- S$pkg[S$declared == "UNDECLARED"]
noncran_unguarded <- S$pkg[S$on_cran == "NO" & S$pkg != self &
                             vapply(S$pkg, function(p) !all(U$guarded[U$pkg == p]), logical(1))]
cat("\nUndeclared packages:", if (length(und)) paste(und, collapse = ", ") else "none", "\n")
cat("Packages not on CRAN:", paste(setdiff(S$pkg[S$on_cran == "NO"], self), collapse = ", "), "\n")
cat("  with an unguarded use:", if (length(noncran_unguarded)) paste(noncran_unguarded, collapse = ", ") else "none", "\n")
quit(status = if (length(und) || length(noncran_unguarded) || !g3_ok) 1 else 0)
