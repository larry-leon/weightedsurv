## G7: parse every workflow file as YAML; pkgdown needs: website; check workflow OS coverage; .Rbuildignore
for (f in list.files(".github/workflows", "\\.ya?ml$", full.names = TRUE)) {
  y <- tryCatch(yaml::read_yaml(f), error = function(e) e)
  cat(f, ": parses:", !inherits(y, "error"), "\n")
  if (inherits(y, "error")) next
  for (j in names(y$jobs)) {
    st <- y$jobs[[j]]$steps
    for (s in st) if (!is.null(s$uses) && grepl("setup-r-dependencies", s$uses)) cat("   job", j, "setup-r-dependencies needs:", s$with$needs, "\n")
    m <- y$jobs[[j]]$strategy$matrix$config
    if (!is.null(m)) cat("   job", j, "matrix:", paste(vapply(m, function(z) paste0(z$os, "/", z$r), ""), collapse = ", "), "\n")
  }
  trig <- y[["on"]]; if (is.null(trig)) trig <- y[[TRUE]]
  cat("   triggers:", paste(names(trig), collapse = ", "), "\n")
}
cat("^\\.github$ in .Rbuildignore:", "^\\.github$" %in% readLines(".Rbuildignore"), "\n")
