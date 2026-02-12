# Create Baseline Characteristics Table by Treatment Arm

This function generates publication-ready baseline characteristic tables
commonly used in clinical trials and observational studies. It
calculates summary statistics, p-values, and standardized mean
differences for continuous, categorical, and binary variables.

## Usage

``` r
create_baseline_table(
  data,
  treat_var = "treat",
  vars_continuous = NULL,
  vars_categorical = NULL,
  vars_binary = NULL,
  var_labels = NULL,
  digits = 1,
  show_pvalue = TRUE,
  show_smd = TRUE,
  show_missing = TRUE
)
```

## Arguments

- data:

  Data frame containing baseline variables

- treat_var:

  Name of treatment variable (default: "treat")

- vars_continuous:

  Character vector of continuous variable names

- vars_categorical:

  Character vector of categorical variable names

- vars_binary:

  Character vector of binary variable names

- var_labels:

  Named vector for variable labels (e.g., c(age = "Age (years)"))

- digits:

  Number of decimal places for continuous variables (default: 1)

- show_pvalue:

  Logical, whether to show p-values (default: TRUE)

- show_smd:

  Logical, whether to show standardized mean differences (default: TRUE)

- show_missing:

  Logical, whether to show missing data counts (default: TRUE)

## Value

A gt table object (if gt package is available) or data frame

## Examples

``` r
# \donttest{
# Create sample data
set.seed(123)
n <- 500
sample_data <- data.frame(
  treat = rbinom(n, 1, 0.5),
  age = rnorm(n, mean = 55, sd = 10),
  stage = sample(c("I", "II", "III", "IV"), n, replace = TRUE),
  sex = rbinom(n, 1, 0.45),
  smoking = rbinom(n, 1, 0.3)
)

# Create table
table <- create_baseline_table(
  data = sample_data,
  treat_var = "treat",
  vars_continuous = "age",
  vars_categorical = "stage",
  vars_binary = c("sex", "smoking"),
  var_labels = c(
    age = "Age (years)",
    stage = "Disease Stage",
    sex = "Female",
    smoking = "Current Smoker"
  )
)

# }
```
