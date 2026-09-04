# 
# LE50 calculation (logistic regression + bootstrap CI)
# looped across CSV files in taxa folder
# ----------------------------------------------------------

library(dplyr)
library(boot)
library(tools)

# -------------------------------------------------------------
# Abbott's correction:
#   P_tc = 100 * [ 1 - (P_t - P_d) / (-P_d) ]
# modified, found in (Newman 1995; Williamson et al. 1999)
#
# Now also corrects each individual dark control replicate against
# the dark control mean (P_d), and tags those rows with dose = 0, so
# they can be used directly as dose=0 anchor points (e.g. for PRR-).
# -------------------------------------------------------------

abbott_correct <- function(data, response_col = "survival", group_col = "PRR",
                           dark_control_label = "dark control", dose_col = "UV.B") {
  
  control_rows <- data[data[[group_col]] == dark_control_label, ]
  control_vals <- control_rows[[response_col]]
  
  if (length(control_vals) == 0) {
    stop("No rows found matching dark_control_label = '", dark_control_label,
         "' in group_col '", group_col, "'.")
  }
  
  P_d <- mean(control_vals, na.rm = TRUE)
  
  if (P_d == 0 || is.na(P_d)) {
    stop("Dark control mean response (P_d) is 0 or NA -- cannot apply Abbott's correction.")
  }
  
  treated <- data[data[[group_col]] != dark_control_label, ]
  
  P_t <- treated[[response_col]]
  treated$response_corrected_pct <- 100 * (1 - (P_t - P_d) / (-P_d))
  treated$response_corrected_prop <- treated$response_corrected_pct / 100
  
  # ----correct each individual dark control replicate against P_d ----
  control_corrected <- control_rows
  control_corrected$response_corrected_pct <- 100 * (1 - (control_vals - P_d) / (-P_d))
  control_corrected$response_corrected_prop <- control_corrected$response_corrected_pct / 100
  control_corrected[[dose_col]] <- 0
  
  list(data_corrected = treated, dark_control_corrected = control_corrected, P_d = P_d)
}

# -------------------------------------------------------------
# Main LE50 function: fits 2-parameter logistic curve, finds LE50
# via two-point interpolation (or extrapolation), and bootstraps a CI.
# -------------------------------------------------------------

calc_LE50 <- function(dose, response,
                      n_boot = 2000,
                      conf = 0.95,
                      seed = 123) {
  
  df <- data.frame(dose = dose, response = response)
  df <- df[stats::complete.cases(df), ]
  
  if (nrow(df) < 3) stop("Fewer than 3 usable observations -- cannot fit logistic model.")
  if (length(unique(df$dose)) < 2) stop("Need at least 2 distinct dose levels to fit a curve.")
  
  # ---- compute starting point from data ---
  # (simple logit-linearization: log(p/(1-p)) ~ dose)
  p_clamped <- pmax(pmin(df$response, 0.99), 0.01)
  lin_fit <- lm(log(p_clamped / (1 - p_clamped)) ~ dose, data = df)
  start_vals <- list(a = as.numeric(coef(lin_fit)[1]), b = as.numeric(coef(lin_fit)[2]))
  
  # ---- fit logistic curve: response = 1 / (1 + exp(-(a + b*dose))) ----
  # warnOnly = FALSE: a non-converged fit throws an error instead of silently returning a bad result
  
  fit_logistic <- function(d, start) {
    nls(response ~ 1 / (1 + exp(-(a + b * dose))),
        data = d, start = start,
        control = nls.control(maxiter = 200, warnOnly = FALSE))
  }
  
  fit <- fit_logistic(df, start_vals)
  a <- coef(fit)["a"]
  b <- coef(fit)["b"]
  
  means <- aggregate(response ~ dose, data = df, FUN = mean)
  means <- means[order(means$dose), ]
  means$pred <- 1 / (1 + exp(-(a + b * means$dose)))
  
  LE50_bracket <- NA
  for (i in 1:(nrow(means) - 1)) {
    y1 <- means$pred[i]; y2 <- means$pred[i + 1]
    if ((y1 - 0.5) * (y2 - 0.5) <= 0) {
      x1 <- means$dose[i]; x2 <- means$dose[i + 1]
      LE50_bracket <- x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)
      break
    }
  }
  
  LE50_direct <- as.numeric(-a / b)
  LE50_point <- if (!is.na(LE50_bracket)) LE50_bracket else LE50_direct
  
  # ---- bootstrap CI ----
  # Reuses the starting values from the already-fitted main model
  
  boot_fn <- function(d, indices) {
    d_boot <- d[indices, ]
    if (length(unique(d_boot$dose)) < 2) return(NA)
    
    tryCatch({
      fit_b <- fit_logistic(d_boot, start = list(a = as.numeric(a), b = as.numeric(b)))
      a_b <- coef(fit_b)["a"]; b_b <- coef(fit_b)["b"]
      
      means_b <- aggregate(response ~ dose, data = d_boot, FUN = mean)
      means_b <- means_b[order(means_b$dose), ]
      means_b$pred <- 1 / (1 + exp(-(a_b + b_b * means_b$dose)))
      
      LE50_b <- NA
      for (i in 1:(nrow(means_b) - 1)) {
        y1 <- means_b$pred[i]; y2 <- means_b$pred[i + 1]
        if ((y1 - 0.5) * (y2 - 0.5) <= 0) {
          x1 <- means_b$dose[i]; x2 <- means_b$dose[i + 1]
          LE50_b <- x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)
          break
        }
      }
      if (is.na(LE50_b)) LE50_b <- as.numeric(-a_b / b_b)
      return(LE50_b)
    }, error = function(e) NA)
  }
  
  set.seed(seed)
  boot_out <- boot(df, boot_fn, R = n_boot, strata = factor(df$dose))
  boot_vals <- boot_out$t[!is.na(boot_out$t)]
  
  if (length(boot_vals) < n_boot * 0.5) {
    warning(sprintf(
      "Only %d/%d bootstrap resamples converged -- LE50 CI may be unstable for this group.",
      length(boot_vals), n_boot))
  }
  
  if (length(boot_vals) < 10) {
    stop("Bootstrap failed for most/all resamples -- model unstable for this group.")
  }
  
  alpha <- 1 - conf
  ci <- quantile(boot_vals, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  
  list(
    coefficients = c(a = as.numeric(a), b = as.numeric(b)),
    means_table = means,
    LE50_point_estimate = LE50_point,
    LE50_extrapolated = is.na(LE50_bracket),
    LE50_CI = ci,
    conf_level = conf,
    n_boot_successful = length(boot_vals),
    n_boot_requested = n_boot
  )
}

# code to loop through files in taxa folder on KS personal computer - would ned to change filepath for others
# -------------------------------------------------------------
process_taxon_file <- function(filepath,
                               response_col = "survival",
                               dose_col = "UV.B",
                               group_col = "PRR",
                               dark_control_label = "dark control",
                               n_boot = 2000,
                               conf = 0.95,
                               add_dark_control_to_PRRminus = TRUE) {
  
  notes <- character(0)
  taxon_name <- file_path_sans_ext(basename(filepath))
  
  result_row <- data.frame(
    taxon = taxon_name,
    LE50_PRRminus = NA_real_, CI_lower_PRRminus = NA_real_, CI_upper_PRRminus = NA_real_,
    LE50_PRRplus  = NA_real_, CI_lower_PRRplus  = NA_real_, CI_upper_PRRplus  = NA_real_,
    notes = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # ---- read file ----
  data <- tryCatch(read.csv(filepath, stringsAsFactors = FALSE),
                   error = function(e) {
                     notes <<- c(notes, paste("Failed to read file:", e$message))
                     NULL
                   })
  if (is.null(data)) {
    result_row$notes <- paste(notes, collapse = "; ")
    return(result_row)
  }
  
  # ---- check required columns exist ----
  required_cols <- c(dose_col, group_col, response_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    notes <- c(notes, paste0("Missing column(s): ", paste(missing_cols, collapse = ", ")))
    result_row$notes <- paste(notes, collapse = "; ")
    return(result_row)
  }
  
  # ---- Abbott correction (uses dark control group) ----
  correction <- tryCatch(
    abbott_correct(data, response_col = response_col, group_col = group_col,
                   dark_control_label = dark_control_label, dose_col = dose_col),
    error = function(e) {
      notes <<- c(notes, "Abbott correction not applied (no dark control data found)")
      NULL
    }
  )
  
  if (is.null(correction)) {
    # filter out any control data if it exists, so 0-dose
    # controls do not skew the logistic curve fit of PRR- or PRR+ experimental treatments.
    data_corrected <- data[data[[group_col]] != dark_control_label, ]
    data_corrected$response_corrected_prop <- data_corrected[[response_col]]
  } else {
    data_corrected <- correction$data_corrected
  }
  
  # ---- PRR- ----
  prr_minus <- data_corrected[data_corrected[[group_col]] == "PRR-", ]
  if (nrow(prr_minus) == 0) {
    notes <- c(notes, "No PRR- data")
  } else {
    dose_minus <- prr_minus[[dose_col]]
    resp_minus <- prr_minus$response_corrected_prop
    
    # ---- append dark control replicates as dose = 0----
    if (add_dark_control_to_PRRminus && !is.null(correction)) {
      dark_pts <- correction$dark_control_corrected
      dose_minus <- c(dose_minus, dark_pts[[dose_col]])
      resp_minus <- c(resp_minus, dark_pts$response_corrected_prop)
      notes <- c(notes, sprintf("Added %d dark control replicate(s) as dose=0 points to PRR-", nrow(dark_pts)))
    }
    
    res_minus <- tryCatch(
      calc_LE50(dose = dose_minus,
                response = resp_minus,
                n_boot = n_boot, conf = conf),
      error = function(e) {
        notes <<- c(notes, paste("PRR- model failed:", e$message))
        NULL
      }
    )
    if (!is.null(res_minus)) {
      result_row$LE50_PRRminus <- res_minus$LE50_point_estimate
      result_row$CI_lower_PRRminus <- res_minus$LE50_CI[1]
      result_row$CI_upper_PRRminus <- res_minus$LE50_CI[2]
      if (isTRUE(res_minus$LE50_extrapolated)) {
        notes <- c(notes, "PRR- LE50 extrapolated beyond tested dose range")
      }
    }
  }
  
  # ---- PRR+ ----
  # (dark control is NOT added here, per current design -- only PRR- gets
  # the dose=0 dark control anchor points)
  prr_plus <- data_corrected[data_corrected[[group_col]] == "PRR+", ]
  if (nrow(prr_plus) == 0) {
    notes <- c(notes, "No PRR+ data")
  } else {
    res_plus <- tryCatch(
      calc_LE50(dose = prr_plus[[dose_col]],
                response = prr_plus$response_corrected_prop,
                n_boot = n_boot, conf = conf),
      error = function(e) {
        notes <<- c(notes, paste("PRR+ model failed:", e$message))
        NULL
      }
    )
    if (!is.null(res_plus)) {
      result_row$LE50_PRRplus <- res_plus$LE50_point_estimate
      result_row$CI_lower_PRRplus <- res_plus$LE50_CI[1]
      result_row$CI_upper_PRRplus <- res_plus$LE50_CI[2]
      if (isTRUE(res_plus$LE50_extrapolated)) {
        notes <- c(notes, "PRR+ LE50 extrapolated beyond tested dose range")
      }
    }
  }
  
  result_row$notes <- if (length(notes) > 0) paste(notes, collapse = "; ") else ""
  result_row
}

# run it w/ data fromt taxa folder
# ------------------------------------------------------------------

data_dir <- "../data/taxa"

csv_files <- list.files(path = data_dir, pattern = "\\.csv$", full.names = TRUE)

cat(sprintf("Found %d CSV files in '%s'\n", length(csv_files), data_dir))

results_list <- lapply(csv_files, function(f) {
  cat("Processing:", basename(f), "\n")
  process_taxon_file(f)
})

results_df <- do.call(rbind, results_list)

# ---- view and save ----
print(results_df)

write.csv(results_df, "../data/LE50_summary_all_taxa2.csv", row.names = FALSE)

#  code to just rerun 1 at a time
# -------------------------------------------------------------
taxon_run <- function(filename) print(process_taxon_file(file.path(data_dir, filename)))

taxon_run("coregonus-artedi-larvae.csv")

