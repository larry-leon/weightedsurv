a <- commandArgs(TRUE)   # input, output_dir
options(warn = 2)
t0 <- proc.time()[["elapsed"]]
rmarkdown::render(a[1], output_dir = a[2], intermediates_dir = a[2], envir = new.env(), quiet = TRUE)
cat(sprintf("render %s: OK in %.1f s (options(warn = 2); weightedsurv %s from %s)\n", basename(a[1]), proc.time()[["elapsed"]] - t0,
            packageVersion("weightedsurv"), dirname(system.file(package = "weightedsurv"))))
