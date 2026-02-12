library(usethis)
library(devtools)

# Step 1: Create the package structure
# This will create the template directory
# with Rproject setup
#usethis::create_package(".. /GitHub/weightedsurv")

# Copy functions into R directory
# Step 3: Add README and MIT license
usethis::use_readme_rmd(open = FALSE)
usethis::use_mit_license("Larry Leon")

# Step 4: Add dependencies to DESCRIPTION
desc::desc_set_dep("survival", file = "DESCRIPTION")

#usethis::use_package("ggplot2")

# Step 5: Generate documentation
# Also, run this if revising R files such as @importFrom
# clean up old documentation
unlink("man/*.Rd")

devtools::document()

devtools::load_all()

devtools::check()

devtools::clean_dll()

devtools::build()

# In terminal
# rm -r man/
usethis::use_git()

#usethis::use_github()

# If functions are in namespace but not directly loaded
# devtools::load_all()

# Or  access hidden files:  mypackage::my_function

# To remove file in Terminal --> git rm "file"


# Notes
# Every \item in a \describe{} block must have both a label and a description.
# Do not leave a lone \item{...} without {...} after it.
# Do not next \item inside another \item
# Use plain text, dashes, or a single paragraph for subpoints.

# Run these before CRAN submission
devtools::check(cran = TRUE)     # Full CRAN validation
#devtools::check_win_devel()       # Windows checks


# Retain dev/working locally
# Remove from Git's index only (keeps local files intact)
#git rm -r --cached dev/working/
# Then commit
#git add .gitignore

echo "dev/working/" >> .gitignore
git rm -r --cached dev/working/
git add .gitignore
git commit -m "Stop tracking dev/working/, keep locally"


# pkgdown initialization

usethis::use_pkgdown()

# update yaml via claude
#The simplest way: run one of these in R from your package root and paste the output here.
cat(readLines("NAMESPACE"), sep = "\n")

dir.create("vignettes/articles", recursive = TRUE, showWarnings = FALSE)
usethis::use_build_ignore("vignettes/articles")

pkgdown::check_pkgdown()

# Build the full site (renders all vignettes and articles)
pkgdown::build_site()

# Or build incrementally:
pkgdown::build_home()              # README → index.html
pkgdown::build_reference()         # Function reference pages
pkgdown::build_articles()          # All vignettes and articles
pkgdown::build_article("weightedsurv_examples")      # Just the main vignette
pkgdown::build_article("articles/weightedcox_methodology")     # Just the methodology article

# Preview locally
pkgdown::preview_site()



# What files are in vignettes/articles/?
list.files("vignettes/articles", recursive = TRUE)

# What does pkgdown detect as articles?
pkg <- pkgdown::as_pkgdown()
pkg$vignettes


# Keys referenced in the document
doc <- readLines("vignettes/articles/weightedcox_methodology.Rmd")
regmatches(doc, gregexpr("@[A-Za-z0-9_]+", doc)) |> unlist() |> unique()

# Keys defined in the bib file
grep("^@", readLines("vignettes/articles/weightedcox.bib"), value = TRUE)

