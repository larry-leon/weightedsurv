## Packages used in vignette and article chunks (library(), require(), requireNamespace(), pkg::)
## against DESCRIPTION Imports and Suggests. Read-only; run from the repository root.
d <- read.dcf("DESCRIPTION")
dep <- function(f) gsub("\\s*\\(.*\\)", "", trimws(strsplit(d[, f], ",")[[1]]))
declared <- c(dep("Imports"), dep("Suggests"), "weightedsurv")
base_pk <- c("base", "stats", "graphics", "grDevices", "utils", "methods", "datasets",
             "parallel", "tools", "grid", "splines", "stats4", "compiler", "tcltk")
files <- c(list.files("vignettes", "\\.Rmd$", full.names = TRUE),
           list.files("vignettes/articles", "\\.Rmd$", full.names = TRUE))
for (f in files) {
  x <- readLines(f)
  pat1 <- "(?<=library\\()[A-Za-z0-9.]+|(?<=require\\()[A-Za-z0-9.]+|(?<=requireNamespace\\([\"'])[A-Za-z0-9.]+"
  pat2 <- "[A-Za-z][A-Za-z0-9.]*(?=:::?[A-Za-z_.])"
  p <- unique(c(unlist(regmatches(x, gregexpr(pat1, x, perl = TRUE))),
                unlist(regmatches(x, gregexpr(pat2, x, perl = TRUE)))))
  p <- sort(p)
  und <- setdiff(p, declared)
  cat(sprintf("%s\n   used: %s\n   not in Imports/Suggests: %s\n", f, paste(p, collapse = ", "),
              paste(c(setdiff(und, base_pk), sprintf("%s (base R)", intersect(und, base_pk))), collapse = ", ")))
}
