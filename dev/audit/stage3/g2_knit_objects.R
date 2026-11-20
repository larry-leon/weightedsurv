# G2: knit the main vignette in a fresh session and save every object left in the knit environment.
# usage: Rscript g2_knit_objects.R <lib> <dir holding weightedsurv_examples.Rmd> <out.rds>
a <- commandArgs(TRUE)
.libPaths(c(a[1], .libPaths()))
stopifnot(as.character(packageVersion("weightedsurv")) == "0.2.0")
e <- new.env(parent = globalenv())
t0 <- Sys.time()
rmarkdown::render(file.path(a[2], "weightedsurv_examples.Rmd"), output_dir = a[2], envir = e, quiet = TRUE)
cat("knit seconds:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n")
saveRDS(mget(ls(e, all.names = TRUE), envir = e), a[3])
cat("objects:", length(ls(e, all.names = TRUE)), "\n")
