#' Create Baseline Characteristics Table by Treatment Arm
#'
#' This function generates publication-ready baseline characteristic tables
#' commonly used in clinical trials and observational studies. It calculates
#' summary statistics, p-values, and standardized mean differences for
#' continuous, categorical, and binary variables.
#'
#' @param data Data frame containing baseline variables
#' @param treat_var Name of treatment variable (default: "treat")
#' @param vars_continuous Character vector of continuous variable names
#' @param vars_categorical Character vector of categorical variable names
#' @param vars_binary Character vector of binary variable names
#' @param var_labels Named vector for variable labels (e.g., c(age = "Age (years)"))
#' @param digits Number of decimal places for continuous variables (default: 1)
#' @param show_pvalue Logical, whether to show p-values (default: TRUE)
#' @param show_smd Logical, whether to show standardized mean differences (default: TRUE)
#' @param show_missing Logical, whether to show missing data counts (default: TRUE)
#'
#' @return A gt table object (if gt package is available) or data frame
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' set.seed(123)
#' n <- 500
#' sample_data <- data.frame(
#'   treat = rbinom(n, 1, 0.5),
#'   age = rnorm(n, mean = 55, sd = 10),
#'   stage = sample(c("I", "II", "III", "IV"), n, replace = TRUE),
#'   sex = rbinom(n, 1, 0.45),
#'   smoking = rbinom(n, 1, 0.3)
#' )
#'
#' # Create table
#' table <- create_baseline_table(
#'   data = sample_data,
#'   treat_var = "treat",
#'   vars_continuous = "age",
#'   vars_categorical = "stage",
#'   vars_binary = c("sex", "smoking"),
#'   var_labels = c(
#'     age = "Age (years)",
#'     stage = "Disease Stage",
#'     sex = "Female",
#'     smoking = "Current Smoker"
#'   )
#' )
#'
#' }
#'
#' @export
#' @importFrom stats sd t.test chisq.test fisher.test na.omit
#'
create_baseline_table <- function(data,
                                  treat_var = "treat",
                                  vars_continuous = NULL,
                                  vars_categorical = NULL,
                                  vars_binary = NULL,
                                  var_labels = NULL,
                                  digits = 1,
                                  show_pvalue = TRUE,
                                  show_smd = TRUE,
                                  show_missing = TRUE) {

  # Check if gt package is available
  use_gt <- requireNamespace("gt", quietly = TRUE)

  # Get treatment groups
  treat_levels <- sort(unique(data[[treat_var]]))
  n_treat <- length(treat_levels)

  if(n_treat != 2) {
    stop("Currently supports only 2 treatment arms. Found ", n_treat, " levels in ", treat_var)
  }

  # Split data by treatment
  data_0 <- data[data[[treat_var]] == treat_levels[1], ]
  data_1 <- data[data[[treat_var]] == treat_levels[2], ]

  # Initialize results table
  results <- data.frame(
    Variable = character(),
    Level = character(),
    Control = character(),
    Treatment = character(),
    P_value = numeric(),
    SMD = numeric(),
    stringsAsFactors = FALSE
  )

  # =====================================
  # Process continuous variables
  # =====================================
  if(!is.null(vars_continuous)) {
    for(var in vars_continuous) {

      # Check if variable exists
      if(!var %in% names(data)) {
        warning("Variable '", var, "' not found in data")
        next
      }

      # Calculate statistics
      mean_0 <- mean(data_0[[var]], na.rm = TRUE)
      sd_0 <- stats::sd(data_0[[var]], na.rm = TRUE)
      n_0 <- sum(!is.na(data_0[[var]]))
      miss_0 <- sum(is.na(data_0[[var]]))

      mean_1 <- mean(data_1[[var]], na.rm = TRUE)
      sd_1 <- stats::sd(data_1[[var]], na.rm = TRUE)
      n_1 <- sum(!is.na(data_1[[var]]))
      miss_1 <- sum(is.na(data_1[[var]]))

      # Format as mean (SD)
      val_0 <- sprintf(paste0("%.", digits, "f (", "%.", digits, "f)"),
                       mean_0, sd_0)
      val_1 <- sprintf(paste0("%.", digits, "f (", "%.", digits, "f)"),
                       mean_1, sd_1)

      # Handle NaN values
      if(is.na(mean_0) || is.na(sd_0)) val_0 <- "N/A"
      if(is.na(mean_1) || is.na(sd_1)) val_1 <- "N/A"

      # Calculate p-value (t-test)
      if(show_pvalue && n_0 > 1 && n_1 > 1) {
        p_val <- tryCatch(
          stats::t.test(data_0[[var]], data_1[[var]])$p.value,
          error = function(e) NA
        )
      } else {
        p_val <- NA
      }

      # Calculate SMD (Cohen's d)
      if(show_smd && n_0 > 1 && n_1 > 1 && sd_0 > 0 && sd_1 > 0) {
        pooled_sd <- sqrt(((n_0 - 1) * sd_0^2 + (n_1 - 1) * sd_1^2) /
                            (n_0 + n_1 - 2))
        if(pooled_sd > 0) {
          smd <- abs(mean_1 - mean_0) / pooled_sd
        } else {
          smd <- NA
        }
      } else {
        smd <- NA
      }

      # Variable label
      var_label <- ifelse(!is.null(var_labels) && var %in% names(var_labels),
                          var_labels[var], var)

      # Add main row
      results <- rbind(results, data.frame(
        Variable = var_label,
        Level = "Mean (SD)",
        Control = val_0,
        Treatment = val_1,
        P_value = p_val,
        SMD = smd,
        stringsAsFactors = FALSE
      ))

      # Add missing data row if needed
      if(show_missing && (miss_0 > 0 || miss_1 > 0)) {
        results <- rbind(results, data.frame(
          Variable = "",
          Level = "Missing",
          Control = paste0(miss_0, " (",
                           sprintf("%.1f", 100 * miss_0/nrow(data_0)), "%)"),
          Treatment = paste0(miss_1, " (",
                             sprintf("%.1f", 100 * miss_1/nrow(data_1)), "%)"),
          P_value = NA,
          SMD = NA,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # =====================================
  # Process categorical variables
  # =====================================
  if(!is.null(vars_categorical)) {
    for(var in vars_categorical) {

      # Check if variable exists
      if(!var %in% names(data)) {
        warning("Variable '", var, "' not found in data")
        next
      }

      # Get levels
      levels_var <- sort(unique(stats::na.omit(data[[var]])))

      # Variable label
      var_label <- ifelse(!is.null(var_labels) && var %in% names(var_labels),
                          var_labels[var], var)

      # Calculate p-value (chi-square)
      if(show_pvalue && length(levels_var) > 1) {
        tab <- table(data[[var]], data[[treat_var]])
        if(all(dim(tab) > 1) && min(tab) >= 5) {
          p_val <- tryCatch(
            stats::chisq.test(tab)$p.value,
            error = function(e) NA
          )
        } else if(all(dim(tab) > 1)) {
          # Use Fisher's exact test for small cell counts
          p_val <- tryCatch(
            stats::fisher.test(tab, simulate.p.value = TRUE)$p.value,
            error = function(e) NA
          )
        } else {
          p_val <- NA
        }
      } else {
        p_val <- NA
      }

      # Calculate Cramer's V for effect size
      if(show_smd && length(levels_var) > 1) {
        tab <- table(data[[var]], data[[treat_var]])
        if(all(dim(tab) > 1) && sum(tab) > 0) {
          chi2 <- tryCatch(
            stats::chisq.test(tab)$statistic,
            error = function(e) NA
          )
          if(!is.na(chi2)) {
            n <- sum(tab)
            k <- min(nrow(tab), ncol(tab))
            cramers_v <- sqrt(chi2 / (n * (k - 1)))
          } else {
            cramers_v <- NA
          }
        } else {
          cramers_v <- NA
        }
      } else {
        cramers_v <- NA
      }

      # Add header row
      results <- rbind(results, data.frame(
        Variable = var_label,
        Level = "",
        Control = "",
        Treatment = "",
        P_value = p_val,
        SMD = cramers_v,
        stringsAsFactors = FALSE
      ))

      # Add level rows
      for(lev in levels_var) {
        n_0 <- sum(data_0[[var]] == lev, na.rm = TRUE)
        n_1 <- sum(data_1[[var]] == lev, na.rm = TRUE)

        total_0 <- sum(!is.na(data_0[[var]]))
        total_1 <- sum(!is.na(data_1[[var]]))

        pct_0 <- if(total_0 > 0) 100 * n_0 / total_0 else 0
        pct_1 <- if(total_1 > 0) 100 * n_1 / total_1 else 0

        val_0 <- sprintf("%d (%.1f%%)", n_0, pct_0)
        val_1 <- sprintf("%d (%.1f%%)", n_1, pct_1)

        results <- rbind(results, data.frame(
          Variable = "",
          Level = paste0("  ", lev),
          Control = val_0,
          Treatment = val_1,
          P_value = NA,
          SMD = NA,
          stringsAsFactors = FALSE
        ))
      }

      # Add missing row if needed
      miss_0 <- sum(is.na(data_0[[var]]))
      miss_1 <- sum(is.na(data_1[[var]]))

      if(show_missing && (miss_0 > 0 || miss_1 > 0)) {
        results <- rbind(results, data.frame(
          Variable = "",
          Level = "  Missing",
          Control = paste0(miss_0, " (",
                           sprintf("%.1f", 100 * miss_0/nrow(data_0)), "%)"),
          Treatment = paste0(miss_1, " (",
                             sprintf("%.1f", 100 * miss_1/nrow(data_1)), "%)"),
          P_value = NA,
          SMD = NA,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # =====================================
  # Process binary variables
  # =====================================
  if(!is.null(vars_binary)) {
    for(var in vars_binary) {

      # Check if variable exists
      if(!var %in% names(data)) {
        warning("Variable '", var, "' not found in data")
        next
      }

      # Calculate proportions
      n_0_yes <- sum(data_0[[var]] == 1, na.rm = TRUE)
      n_1_yes <- sum(data_1[[var]] == 1, na.rm = TRUE)

      n_0_total <- sum(!is.na(data_0[[var]]))
      n_1_total <- sum(!is.na(data_1[[var]]))

      if(n_0_total > 0 && n_1_total > 0) {
        pct_0 <- 100 * n_0_yes / n_0_total
        pct_1 <- 100 * n_1_yes / n_1_total

        val_0 <- sprintf("%d (%.1f%%)", n_0_yes, pct_0)
        val_1 <- sprintf("%d (%.1f%%)", n_1_yes, pct_1)

        # Calculate p-value (Fisher's exact or chi-square)
        if(show_pvalue) {
          tab <- table(factor(data[[var]], levels = c(0, 1)),
                       data[[treat_var]])
          if(all(dim(tab) == c(2, 2))) {
            if(min(tab) < 5) {
              p_val <- tryCatch(
                stats::fisher.test(tab)$p.value,
                error = function(e) NA
              )
            } else {
              p_val <- tryCatch(
                stats::chisq.test(tab)$p.value,
                error = function(e) NA
              )
            }
          } else {
            p_val <- NA
          }
        } else {
          p_val <- NA
        }

        # Calculate SMD for binary
        if(show_smd) {
          p0 <- n_0_yes / n_0_total
          p1 <- n_1_yes / n_1_total
          denom <- sqrt((p0 * (1 - p0) + p1 * (1 - p1)) / 2)
          if(denom > 0) {
            smd <- abs(p1 - p0) / denom
          } else {
            smd <- NA
          }
        } else {
          smd <- NA
        }
      } else {
        val_0 <- "0 (0.0%)"
        val_1 <- "0 (0.0%)"
        p_val <- NA
        smd <- NA
      }

      # Variable label
      var_label <- ifelse(!is.null(var_labels) && var %in% names(var_labels),
                          var_labels[var], var)

      results <- rbind(results, data.frame(
        Variable = var_label,
        Level = "",
        Control = val_0,
        Treatment = val_1,
        P_value = p_val,
        SMD = smd,
        stringsAsFactors = FALSE
      ))

      # Add missing row if needed
      if(show_missing) {
        miss_0 <- sum(is.na(data_0[[var]]))
        miss_1 <- sum(is.na(data_1[[var]]))

        if(miss_0 > 0 || miss_1 > 0) {
          results <- rbind(results, data.frame(
            Variable = "",
            Level = "Missing",
            Control = paste0(miss_0, " (",
                             sprintf("%.1f", 100 * miss_0/nrow(data_0)), "%)"),
            Treatment = paste0(miss_1, " (",
                               sprintf("%.1f", 100 * miss_1/nrow(data_1)), "%)"),
            P_value = NA,
            SMD = NA,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  # Format p-values
  results$P_value_formatted <- ifelse(
    is.na(results$P_value), "",
    ifelse(results$P_value < 0.001, "<0.001",
           sprintf("%.3f", results$P_value))
  )

  # Format SMD
  results$SMD_formatted <- ifelse(
    is.na(results$SMD), "",
    sprintf("%.2f", results$SMD)
  )

  # Create gt table if package is available
  if(use_gt) {
    # Select columns for display
    display_data <- results[, c("Variable", "Level", "Control", "Treatment",
                                "P_value_formatted", "SMD_formatted")]

    gt_table <- gt::gt(display_data)

    gt_table <- gt::cols_label(
      gt_table,
      Variable = "Characteristic",
      Level = "",
      Control = paste0("Control\n(n=", nrow(data_0), ")"),
      Treatment = paste0("Treatment\n(n=", nrow(data_1), ")"),
      P_value_formatted = "P-value",
      SMD_formatted = "SMD"
    )

    gt_table <- gt::tab_header(
      gt_table,
      title = "Baseline Characteristics by Treatment Arm"
    )

    # Check gt version and use appropriate function
    # sub_missing was introduced in gt v0.6.0 to replace fmt_missing
    if (utils::packageVersion("gt") >= "0.6.0") {
      gt_table <- gt::sub_missing(
        gt_table,
        columns = gt::everything(),
        missing_text = ""
      )
    } else {
      gt_table <- gt::fmt_missing(
        gt_table,
        columns = gt::everything(),
        missing_text = ""
      )
    }

    gt_table <- gt::tab_style(
      gt_table,
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(
        columns = "Variable",
        rows = display_data$Variable != ""
      )
    )

    gt_table <- gt::tab_style(
      gt_table,
      style = gt::cell_borders(
        sides = "bottom",
        color = "black",
        weight = gt::px(2)
      ),
      locations = gt::cells_column_labels()
    )

    gt_table <- gt::tab_footnote(
      gt_table,
      footnote = "SMD = Standardized mean difference",
      locations = gt::cells_column_labels(columns = "SMD_formatted")
    )

    gt_table <- gt::tab_footnote(
      gt_table,
      footnote = "P-values: t-test for continuous, chi-square/Fisher's exact for categorical/binary variables",
      locations = gt::cells_column_labels(columns = "P_value_formatted")
    )

    # Remove columns if not needed
    if(!show_pvalue) {
      gt_table <- gt::cols_hide(gt_table, columns = "P_value_formatted")
    }

    if(!show_smd) {
      gt_table <- gt::cols_hide(gt_table, columns = "SMD_formatted")
    }

    return(gt_table)

  } else {
    # Return data frame if gt not available
    message("Note: gt package not available. Returning data frame instead of formatted table.")

    # Clean up column names
    names(results)[names(results) == "Control"] <- paste0("Control (n=", nrow(data_0), ")")
    names(results)[names(results) == "Treatment"] <- paste0("Treatment (n=", nrow(data_1), ")")

    # Select columns based on options
    cols_to_keep <- c("Variable", "Level",
                      paste0("Control (n=", nrow(data_0), ")"),
                      paste0("Treatment (n=", nrow(data_1), ")"))

    if(show_pvalue) {
      cols_to_keep <- c(cols_to_keep, "P_value_formatted")
    }

    if(show_smd) {
      cols_to_keep <- c(cols_to_keep, "SMD_formatted")
    }

    return(results[, cols_to_keep])
  }
}

# =====================================
# Helper function for ForestSearch
# =====================================

#' Create ForestSearch Baseline Table
#'
#' Wrapper function specifically tailored for ForestSearch package data structure.
#' Automatically identifies variable types based on naming conventions.
#'
#' @param df_work Data frame with ForestSearch format
#' @param output_format "gt" for gt table or "data.frame" for simple data frame
#' @param save_path Optional file path to save the table (supports .html, .docx, .tex, .csv)
#'
#' @return A gt table object or data frame depending on output_format
#'
#' @export
#' @importFrom tools toTitleCase
#'
create_forestsearch_baseline_table <- function(df_work,
                                               output_format = "gt",
                                               save_path = NULL) {

  # Identify variable types based on naming convention
  vars_continuous <- c()
  vars_binary <- c()
  vars_categorical <- c()

  # Get all variable names except outcome and treatment
  all_vars <- setdiff(names(df_work), c("treat", "y", "event", "time", "status"))

  # Detect continuous variables (z-scored variables or known continuous)
  z_vars <- grep("^z_", all_vars, value = TRUE)
  continuous_keywords <- c("age", "biomarker", "bmi", "score", "value", "level")

  for(var in all_vars) {
    if(var %in% z_vars) {
      vars_continuous <- c(vars_continuous, var)
    } else if(any(sapply(continuous_keywords, function(x) grepl(x, var, ignore.case = TRUE)))) {
      # Check if actually continuous
      if(is.numeric(df_work[[var]]) && length(unique(stats::na.omit(df_work[[var]]))) > 10) {
        vars_continuous <- c(vars_continuous, var)
      }
    }
  }

  # Detect binary and categorical variables
  remaining_vars <- setdiff(all_vars, vars_continuous)

  for(var in remaining_vars) {
    unique_vals <- unique(stats::na.omit(df_work[[var]]))
    if(length(unique_vals) == 2) {
      vars_binary <- c(vars_binary, var)
    } else if(length(unique_vals) <= 10) {  # Assume categorical if ≤10 unique values
      vars_categorical <- c(vars_categorical, var)
    }
  }

  # Create nice labels
  var_labels <- c()
  for(var in c(vars_continuous, vars_categorical, vars_binary)) {
    # Remove z_ prefix for labels
    clean_name <- gsub("^z_", "", var)
    # Replace underscores with spaces and capitalize
    clean_name <- gsub("_", " ", clean_name)
    clean_name <- tools::toTitleCase(clean_name)

    # Add units for known variables
    if(grepl("age", var, ignore.case = TRUE)) clean_name <- paste0(clean_name, " (years)")
    if(grepl("bmi", var, ignore.case = TRUE)) clean_name <- paste0(clean_name, " (kg/m²)")
    if(grepl("^z_", var)) clean_name <- paste0(clean_name, " (standardized)")

    var_labels[var] <- clean_name
  }

  # Create table
  if(output_format == "gt" && !requireNamespace("gt", quietly = TRUE)) {
    message("gt package not available, switching to data.frame output")
    output_format <- "data.frame"
  }

  table <- create_baseline_table(
    data = df_work,
    treat_var = "treat",
    vars_continuous = vars_continuous,
    vars_categorical = vars_categorical,
    vars_binary = vars_binary,
    var_labels = var_labels,
    show_pvalue = TRUE,
    show_smd = TRUE
  )

  return(table)
}
