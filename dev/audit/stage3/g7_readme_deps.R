# G7: the README "Dependencies" lists against DESCRIPTION Imports and Suggests. Run from the package root.
d <- read.dcf("DESCRIPTION", fields = c("Imports", "Suggests"))
dep <- function(x) trimws(strsplit(gsub("\n", " ", x), ",")[[1]])
r <- readLines("README.md")
block <- function(head) {
  i <- which(r == head) + 1; out <- character()
  while (i <= length(r) && grepl("^- ", r[i])) { out <- c(out, sub("^- ", "", r[i])); i <- i + 1 }
  out
}
for (f in c("Imports", "Suggests")) {
  a <- dep(d[, f]); b <- block(paste0("**", f, ":**"))
  cat(f, "\n  DESCRIPTION:", paste(a, collapse = ", "), "\n  README:     ", paste(b, collapse = ", "),
      "\n  identical (same entries, same order):", identical(a, b), "\n")
}
p <- paste(r, collapse = "\n")
cat("parentheses in README: open", lengths(regmatches(p, gregexpr("\\(", p))), "close", lengths(regmatches(p, gregexpr("\\)", p))), "\n")
