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


