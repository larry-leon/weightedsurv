## Section 5: write the recomputed fields into the one CSV row (text edit of that line only; every other
## line is left byte-identical). usage: Rscript sec5_write_csv.R <new.rds>
csv <- "vignettes/articles/sim_results/ps_weight_results_by_time.csv"
new <- readRDS(commandArgs(TRUE)[1])[c("mean_estimate", "bias", "bias_mcse", "coverage_known", "coverage_known_mcse",
                                     "coverage_estimated", "coverage_estimated_mcse",
                                     "coverage_known_bootreps", "coverage_estimated_bootreps")]
L <- readLines(csv); X <- read.csv(csv)
row <- which(X$source == "OC" & X$cell_id == 8 & X$estimand == "R" & X$time %in% 36)
li <- row + 1
f <- scan(text = L[li], what = "", sep = ",", quote = "\"", quiet = TRUE, na.strings = character(0))
stopifnot(length(f) == ncol(X))
q <- grepl('^"', strsplit(gsub('"[^"]*"', '"Q"', L[li]), ",")[[1]])
for (nm in names(new)) f[match(nm, names(X))] <- format(new[[nm]], digits = 15)
f[q] <- paste0('"', f[q], '"')
L[li] <- paste(f, collapse = ",")
writeLines(L, csv)
Y <- read.csv(csv)
stopifnot(identical(Y[-row, ], X[-row, ]), identical(Y[row, setdiff(names(X), names(new))], X[row, setdiff(names(X), names(new))]))
print(rbind(old = unlist(X[row, names(new)]), new = unlist(Y[row, names(new)])), digits = 15)
