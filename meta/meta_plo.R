# ==============================================================================
# Integrated main-model + bootstrap meta-regression: PLO (logit-transformed prevalence)
#
# Framework:
#   1) Main model on original full data
#      - REML: coefficients / CI / I² / variance components / QM
#      - ML:   overall LRT (full vs null), per-variable LRT (full vs reduced)
#
#   2) Bootstrap sensitivity analysis (record-level bootstrap)
#      - 200 bootstrap iterations
#      - re-run VIF filtering each iteration
#      - re-fit REML/ML models
#      - summarize Prop. and median p for variable-level robustness
#
# Notes:
#   - No-intercept parameterization is retained for descriptive coefficient output
#   - Variable-level inference should rely on ML LRT, not coefficient p-values
#   - Bootstrap is used as sensitivity analysis, not the sole basis of inference
# ==============================================================================

suppressPackageStartupMessages({
  library(metafor)
  library(car)
  library(openxlsx)
  library(tidyverse)
})

# ==============================================================================
# 0) Paths and global parameters
# ==============================================================================
base_path   <- "/data/home/test01/meta"
input_file  <- file.path(base_path, "data/6_Supplementary_table_3.xlsx")
output_path <- file.path(base_path, "out/plo_integrated")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

set.seed(123)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

control <- list(
  optimizer = "optim",
  optmethod = "L-BFGS-B",
  maxiter   = 10000,
  reltol    = 1e-8,
  stepadj   = 0.5
)

BOOT_N <- 200

# ==============================================================================
# 1) Safe utility helpers
# ==============================================================================
safe_scalar <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0) return(default)
  x1 <- x[1]
  if (is.list(x1)) return(default)
  val <- suppressWarnings(as.numeric(x1))
  if (length(val) == 0 || !is.finite(val)) return(default)
  val
}

safe_char_scalar <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0) return(default)
  as.character(x[1])
}

safe_logical_scalar <- function(x, default = NA) {
  if (is.null(x) || length(x) == 0) return(default)
  as.logical(x[1])
}

safe_nrow <- function(x) {
  if (is.null(x)) return(0L)
  nrow(x)
}

safe_rbind <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (length(lst) == 0) return(NULL)
  do.call(rbind, lst)
}

`%notin%` <- function(x, y) !(x %in% y)

# ==============================================================================
# 2) Data loading and preprocessing
# ==============================================================================
df <- read.xlsx(input_file, 1)
message("Original number of rows: ", nrow(df))

dat <- df[, setdiff(names(df), c(
  "Original_habitat", "Original_rodents", "Chinese_rodent_name",
  "Specific_mixed_sample", "Testing_method", "Specific_testing_method",
  "Original_pathogens", "Standardization"
))]

data <- dat %>%
  filter(
    Health_status == "Healthy",
    !is.na(NDVI), !is.na(Time_interval),
    Habitat_type != "Unknown habitat",
    Sample_source %notin% c("Unclassified sample")
  )

message("Number of rows after filtering: ", nrow(data))

# ------------------------------------------------------------------------------
# Factor levels
# ------------------------------------------------------------------------------
define_factor_levels <- function(data, column_name, levels_order) {
  if (column_name %in% names(data)) {
    data[[column_name]] <- factor(data[[column_name]], levels = levels_order)
  }
  data
}

levels_list <- list(
  Habitat_type = c("Artificial habitat", "Natural habitat", "Mixed habitat"),
  Sample_source = c(
    "Alimentary or intestinal sample", "Blood or serum sample",
    "Body surface sample", "Brain sample", "Diaphragm sample",
    "Ear sample", "Heart sample", "Kidney sample",
    "Liver sample", "Lung or respiratory tract sample",
    "Spleen sample", "Urine sample", "Mixed sample"
  ),
  Testing_type = c(
    "Morphological examination", "Microbial culture",
    "Immunological testing", "Molecular testing"
  ),
  Rodent_family = c(
    "Caviidae", "Chinchillidae", "Cricetidae", "Dipodidae", "Hystricidae",
    "Muridae", "Myocastoridae", "Nesomyidae", "Sciuridae", "Spalacidae"
  ),
  Region = c(
    "Central China", "Eastern China", "Northern China",
    "Northeastern China", "Northwestern China",
    "Southern China", "Southwestern China"
  ),
  Time_interval = c("1950-2000", "2001-2010", "2011-2019", "2020-2023"),
  NDVI = c("Peak NDVI", "Higher NDVI", "Medium NDVI", "Low NDVI")
)

ALL_MODERATORS <- names(levels_list)

# ------------------------------------------------------------------------------
# Drop sparse levels record-wise
# ------------------------------------------------------------------------------
filter_data <- function(data, columns, min_count = 2) {
  for (column in columns) {
    if (column %in% names(data)) {
      data <- data %>%
        group_by(across(all_of(column))) %>%
        filter(n() >= min_count) %>%
        ungroup()
    }
  }
  data
}

datasets_raw <- list(
  all       = data,
  viruses   = data %>% filter(Pathogen_type == "Viruses"),
  bacteria  = data %>% filter(Pathogen_type == "Bacteria"),
  parasites = data %>% filter(Pathogen_type == "Parasites")
)

columns_to_check <- c(
  "Habitat_type", "Sample_source", "Testing_type",
  "Rodent_family", "Region", "Time_interval", "NDVI"
)

datasets_raw <- lapply(datasets_raw, function(dataset) {
  for (col_name in names(levels_list)) {
    dataset <- define_factor_levels(dataset, col_name, levels_list[[col_name]])
  }

  dataset$Study <- factor(dataset$Study)

  dataset <- filter_data(dataset, columns_to_check)
  dataset <- droplevels(dataset)

  dataset$Observation <- factor(seq_len(nrow(dataset)))
  dataset
})

# ==============================================================================
# 3) Effect size calculation (PLO)
# ==============================================================================
calc_effect_size <- function(data) {
  if (!is.data.frame(data)) data <- as.data.frame(data)

  data <- data.frame(
    data,
    escalc(
      xi = data$Events,
      ni = data$Total,
      measure = "PLO",
      add = 0.5,
      to = "only0"
    )
  )
  data$original_ratio <- plogis(data$yi)
  data
}

datasets <- lapply(datasets_raw, calc_effect_size)

# Save filtered datasets
output_excel_file <- file.path(output_path, paste0("filtered_datasets_", timestamp, ".xlsx"))
wb_data <- createWorkbook()
for (name in names(datasets)) {
  addWorksheet(wb_data, sheetName = name)
  writeData(wb_data, sheet = name, x = datasets[[name]])
}
saveWorkbook(wb_data, output_excel_file, overwrite = TRUE)
message("Filtered datasets saved to: ", output_excel_file)

# ==============================================================================
# 4) Utility functions
# ==============================================================================

# ------------------------------------------------------------------------------
# I² for multilevel model
# ------------------------------------------------------------------------------
i2 <- function(model) {
  out_na <- data.frame(I2_Overall = NA, I2_Study = NA, I2_Observation = NA)

  if (is.null(model) || !inherits(model, "rma.mv")) return(out_na)
  if (is.null(model$sigma2) || is.null(model$k) || is.null(model$p) || is.null(model$vi)) return(out_na)

  X <- tryCatch(model.matrix(model), error = function(e) NULL)
  if (is.null(X)) return(out_na)

  W <- tryCatch(diag(1 / model$vi), error = function(e) NULL)
  if (is.null(W)) return(out_na)
  if (ncol(X) == 0 || nrow(X) < ncol(X)) return(out_na)

  XWX <- tryCatch(t(X) %*% W %*% X, error = function(e) NULL)
  if (is.null(XWX)) return(out_na)

  det_val <- tryCatch(det(XWX), error = function(e) NA_real_)
  if (!is.finite(det_val) || abs(det_val) < .Machine$double.eps) return(out_na)

  P <- tryCatch(
    W - W %*% X %*% solve(XWX) %*% t(X) %*% W,
    error = function(e) NULL
  )
  if (is.null(P)) return(out_na)

  if (length(model$sigma2) == 0 || all(model$sigma2 == 0, na.rm = TRUE)) return(out_na)

  denom <- sum(model$sigma2, na.rm = TRUE) + (model$k - model$p) / sum(diag(P))
  if (!is.finite(denom) || denom <= 0) return(out_na)

  I2_overall <- 100 * sum(model$sigma2, na.rm = TRUE) / denom
  I2_levels  <- 100 * model$sigma2 / denom

  data.frame(
    I2_Overall     = round(I2_overall, 2),
    I2_Study       = round(I2_levels[1], 2),
    I2_Observation = round(if (length(I2_levels) >= 2) I2_levels[2] else NA, 2)
  )
}

# ------------------------------------------------------------------------------
# Dynamic moderator screening + VIF filtering
# ------------------------------------------------------------------------------
filter_high_vif <- function(data, threshold = 2.5) {
  candidate_vars <- c(
    "Sample_source", "Region", "Habitat_type",
    "Testing_type", "Rodent_family", "Time_interval", "NDVI"
  )

  valid_vars <- candidate_vars[sapply(candidate_vars, function(v) {
    if (!(v %in% names(data))) return(FALSE)
    x <- data[[v]]
    if (all(is.na(x))) return(FALSE)
    nlev <- dplyr::n_distinct(x, na.rm = TRUE)
    nlev >= 2
  })]

  if (length(valid_vars) == 0) {
    return(as.formula("yi ~ 1"))
  }

  full_formula <- as.formula(
    paste("yi ~", paste(valid_vars, collapse = " + "))
  )

  if (nrow(data) < (length(valid_vars) + 2)) return(NULL)

  repeat {
    lm_model <- tryCatch(
      lm(full_formula, data = data, na.action = na.omit),
      error = function(e) NULL
    )
    if (is.null(lm_model)) return(NULL)
    if (length(coef(lm_model)) < 2) break

    vif_values <- tryCatch(car::vif(lm_model), error = function(e) NULL)
    if (is.null(vif_values) || (!is.numeric(vif_values) && !is.matrix(vif_values))) break

    if (is.vector(vif_values) && !is.null(names(vif_values))) {
      gvif <- vif_values
    } else if (is.matrix(vif_values) && ncol(vif_values) >= 2) {
      gvif <- vif_values[, 1]^(1 / (2 * vif_values[, 2]))
    } else {
      break
    }

    if (all(gvif <= threshold, na.rm = TRUE)) break

    highest_vif_var <- names(gvif)[which.max(gvif)]
    if (length(highest_vif_var) == 0 || !(highest_vif_var %in% attr(terms(full_formula), "term.labels"))) break

    full_formula <- update(full_formula, paste(". ~ . -", highest_vif_var))
  }

  full_formula
}

# ------------------------------------------------------------------------------
# One-sided formula for rma.mv
# ------------------------------------------------------------------------------
make_mods_formula <- function(formula, no_intercept = TRUE) {
  rhs <- attr(terms(formula), "term.labels")
  if (length(rhs) == 0) return(~ 1)

  suffix <- if (no_intercept) " - 1" else ""
  as.formula(paste("~", paste(rhs, collapse = " + "), suffix))
}

# ------------------------------------------------------------------------------
# Drop one moderator
# ------------------------------------------------------------------------------
drop_moderator <- function(formula, moderator) {
  rhs <- setdiff(attr(terms(formula), "term.labels"), moderator)
  if (length(rhs) == 0) return(as.formula("yi ~ 1"))
  as.formula(paste("yi ~", paste(rhs, collapse = " + ")))
}

# ------------------------------------------------------------------------------
# Pseudo-R² (descriptive)
# ------------------------------------------------------------------------------
calculate_pseudoR2 <- function(full_fit, null_fit) {
  if (is.null(full_fit) || is.null(null_fit)) return(NA_real_)

  s_full <- sum(full_fit$sigma2, na.rm = TRUE)
  s_null <- sum(null_fit$sigma2, na.rm = TRUE)

  if (!is.finite(s_null) || s_null == 0) return(NA_real_)

  round(max(0, (s_null - s_full) / s_null), 4)
}

# ------------------------------------------------------------------------------
# Fit model
# ------------------------------------------------------------------------------
fit_model <- function(data, formula = NULL, method = c("REML", "ML")) {
  method <- match.arg(method)

  tryCatch({
    if (is.null(formula)) {
      rma.mv(
        yi, vi,
        random  = ~1 | Study/Observation,
        method  = method,
        data    = data,
        control = control,
        verbose = FALSE
      )
    } else {
      mods_f <- make_mods_formula(formula, no_intercept = TRUE)
      rma.mv(
        yi, vi,
        mods    = mods_f,
        random  = ~1 | Study/Observation,
        method  = method,
        data    = data,
        control = control,
        verbose = FALSE
      )
    }
  }, warning = function(w) {
    message("Model warning [", method, "]: ", conditionMessage(w))
    invokeRestart("muffleWarning")
  }, error = function(e) {
    message("Model fitting error [", method, "]: ", e$message)
    NULL
  })
}

# ------------------------------------------------------------------------------
# Reduced ML model
# ------------------------------------------------------------------------------
fit_reduced_ml <- function(data, full_formula, drop_var) {
  reduced_formula <- drop_moderator(full_formula, drop_var)
  rhs <- attr(terms(reduced_formula), "term.labels")

  if (length(rhs) == 0) {
    fit_model(data, formula = NULL, method = "ML")
  } else {
    fit_model(data, formula = reduced_formula, method = "ML")
  }
}

# ------------------------------------------------------------------------------
# Extract compare stats from metafor::anova()
# ------------------------------------------------------------------------------
extract_compare_stats <- function(cmp) {
  out <- list(LRT = NA_real_, df_diff = NA_real_, p_value = NA_real_)
  if (is.null(cmp)) return(out)

  out$LRT <- tryCatch(safe_scalar(cmp$LRT), error = function(e) NA_real_)
  out$p_value <- tryCatch(safe_scalar(cmp$pval), error = function(e) NA_real_)

  out$df_diff <- tryCatch({
    if (!is.null(cmp$parms.f) && !is.null(cmp$parms.r)) {
      pf <- safe_scalar(cmp$parms.f)
      pr <- safe_scalar(cmp$parms.r)
      if (is.finite(pf) && is.finite(pr)) pf - pr else NA_real_
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)

  out
}

# ------------------------------------------------------------------------------
# Per-variable LRT
# ------------------------------------------------------------------------------
per_variable_lrt_test <- function(data, full_formula, full_fit_ml) {
  if (is.null(full_fit_ml)) return(NULL)

  rhs <- attr(terms(full_formula), "term.labels")
  if (length(rhs) == 0) return(NULL)

  results_list <- list()

  for (mod_var in rhs) {
    reduced_ml <- tryCatch(
      fit_reduced_ml(data, full_formula, mod_var),
      error = function(e) NULL
    )

    if (is.null(reduced_ml)) {
      message("  [WARN] reduced ML fit failed for: ", mod_var)
      next
    }

    cmp <- tryCatch(
      anova(reduced_ml, full_fit_ml),
      error = function(e) NULL
    )

    if (is.null(cmp)) {
      rm(reduced_ml)
      next
    }

    stat <- extract_compare_stats(cmp)

    results_list[[mod_var]] <- data.frame(
      variable = mod_var,
      LRT      = stat$LRT,
      df_diff  = stat$df_diff,
      p_value  = stat$p_value,
      stringsAsFactors = FALSE
    )

    rm(reduced_ml)
  }

  safe_rbind(results_list)
}

# ------------------------------------------------------------------------------
# Main model runner
# ------------------------------------------------------------------------------
run_main_model <- function(data, dataset_name) {
  message("\n==============================")
  message("Running main model: ", dataset_name)
  message("==============================")

  full_formula <- filter_high_vif(data, threshold = 2.5)
  if (is.null(full_formula)) {
    warning("No valid formula after VIF filtering for: ", dataset_name)
    return(NULL)
  }

  kept_vars <- attr(terms(full_formula), "term.labels")

  full_reml <- fit_model(data, full_formula, "REML")
  null_reml <- fit_model(data, NULL, "REML")
  full_ml   <- fit_model(data, full_formula, "ML")
  null_ml   <- fit_model(data, NULL, "ML")

  if (any(sapply(list(full_reml, null_reml, full_ml, null_ml), is.null))) {
    warning("Main model fitting failed for: ", dataset_name)
    return(NULL)
  }

  overall_lrt <- tryCatch({
    cmp <- anova(null_ml, full_ml)
    stat <- extract_compare_stats(cmp)
    data.frame(
      test    = "ML_full_vs_null",
      LRT     = stat$LRT,
      df_diff = stat$df_diff,
      p_value = stat$p_value,
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)

  var_lrt <- per_variable_lrt_test(data, full_formula, full_ml)

  coef_df <- tryCatch({
    out <- as.data.frame(coef(summary(full_reml)), stringsAsFactors = FALSE)
    out$term <- rownames(out)
    rownames(out) <- NULL
    out
  }, error = function(e) {
    message("Coefficient extraction error: ", e$message)
    NULL
  })

  i2_df <- i2(full_reml)

  vc_df <- tryCatch({
    s2 <- full_reml$sigma2
    nl <- full_reml$s.nlevels
    sn <- full_reml$s.names

    if (!is.null(s2) && !is.null(nl) && !is.null(sn) &&
        length(s2) == 2 && length(nl) == 2 && length(sn) == 2) {
      data.frame(
        Factor1 = safe_char_scalar(sn[1]),
        Estim1  = safe_scalar(s2[1]),
        Sqrt1   = sqrt(safe_scalar(s2[1])),
        Nlvls1  = safe_scalar(nl[1]),
        Factor2 = safe_char_scalar(sn[2]),
        Estim2  = safe_scalar(s2[2]),
        Sqrt2   = sqrt(safe_scalar(s2[2])),
        Nlvls2  = safe_scalar(nl[2])
      )
    } else {
      NULL
    }
  }, error = function(e) NULL)

  summary_df <- data.frame(
    dataset         = dataset_name,
    n_records       = nrow(data),
    n_studies       = dplyr::n_distinct(data$Study),
    kept_variables  = if (length(kept_vars) > 0) paste(kept_vars, collapse = "; ") else NA_character_,
    QM              = safe_scalar(full_reml$QM),
    QMdf            = safe_scalar(full_reml$QMdf),
    QMp             = safe_scalar(full_reml$QMp),
    QE              = safe_scalar(full_reml$QE),
    QEdf            = if (!is.null(full_reml$k) && !is.null(full_reml$p)) safe_scalar(full_reml$k - full_reml$p) else NA_real_,
    QEp             = safe_scalar(full_reml$QEp),
    pseudo_R2       = calculate_pseudoR2(full_reml, null_reml),
    converged_REML  = safe_logical_scalar(full_reml$conv),
    converged_ML    = safe_logical_scalar(full_ml$conv),
    stringsAsFactors = FALSE
  )

  list(
    dataset_name = dataset_name,
    data         = data,
    formula      = full_formula,
    full_reml    = full_reml,
    null_reml    = null_reml,
    full_ml      = full_ml,
    null_ml      = null_ml,
    summary      = summary_df,
    coef         = coef_df,
    i2           = i2_df,
    vc           = vc_df,
    overall_lrt  = overall_lrt,
    var_lrt      = var_lrt
  )
}

# ------------------------------------------------------------------------------
# Bootstrap resampling
# ------------------------------------------------------------------------------
boot.func <- function(dataset) {
  indices <- sample(seq_len(nrow(dataset)), replace = TRUE)
  dataset[indices, , drop = FALSE]
}

process_bootstrap_data <- function(bootstrap_data, threshold = 2.5) {
  results <- lapply(bootstrap_data, function(data) {
    data <- droplevels(data)
    formula <- tryCatch(filter_high_vif(data, threshold), error = function(e) NULL)
    list(data = data, formula = formula)
  })

  filtered_results <- results[!sapply(results, function(x) is.null(x$formula))]

  list(
    data     = lapply(filtered_results, function(x) x$data),
    formulas = lapply(filtered_results, function(x) x$formula)
  )
}

progress_message <- function(i, total, label, status = "") {
  ts <- format(Sys.time(), "%H:%M:%S")
  message(sprintf("[%s] %s model %d/%d %s", ts, label, i, total, status))
}

# ------------------------------------------------------------------------------
# Extract bootstrap result from one iteration
# ------------------------------------------------------------------------------
extract_boot_results <- function(full_reml, null_reml, full_ml, null_ml, data, full_formula) {
  pseudo_r2 <- calculate_pseudoR2(full_reml, null_reml)
  i2_stats  <- i2(full_reml)

  coef_df <- tryCatch({
    out <- as.data.frame(coef(summary(full_reml)), stringsAsFactors = FALSE)
    out$term <- rownames(out)
    rownames(out) <- NULL
    out
  }, error = function(e) NULL)

  vc_df <- tryCatch({
    s2 <- full_reml$sigma2
    nl <- full_reml$s.nlevels
    sn <- full_reml$s.names

    if (!is.null(s2) && !is.null(nl) && !is.null(sn) &&
        length(s2) == 2 && length(nl) == 2 && length(sn) == 2) {
      data.frame(
        Factor1 = safe_char_scalar(sn[1]),
        Estim1  = safe_scalar(s2[1]),
        Sqrt1   = sqrt(safe_scalar(s2[1])),
        Nlvls1  = safe_scalar(nl[1]),
        Factor2 = safe_char_scalar(sn[2]),
        Estim2  = safe_scalar(s2[2]),
        Sqrt2   = sqrt(safe_scalar(s2[2])),
        Nlvls2  = safe_scalar(nl[2])
      )
    } else {
      NULL
    }
  }, error = function(e) NULL)

  overall_lrt <- tryCatch({
    cmp <- anova(null_ml, full_ml)
    stat <- extract_compare_stats(cmp)
    data.frame(
      test    = "ML_full_vs_null",
      LRT     = stat$LRT,
      df_diff = stat$df_diff,
      p_value = stat$p_value,
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)

  var_lrt <- per_variable_lrt_test(data, full_formula, full_ml)

  list(
    coef        = coef_df,
    sigma2      = full_reml$sigma2,
    vc          = vc_df,
    pseudo_r2   = pseudo_r2,
    i2          = i2_stats,
    k           = safe_scalar(full_reml$k),
    p           = safe_scalar(full_reml$p),
    QM          = safe_scalar(full_reml$QM),
    QMdf        = safe_scalar(full_reml$QMdf),
    QMp         = safe_scalar(full_reml$QMp),
    QE          = safe_scalar(full_reml$QE),
    QEdf        = if (!is.null(full_reml$k) && !is.null(full_reml$p)) safe_scalar(full_reml$k - full_reml$p) else NA_real_,
    QEp         = safe_scalar(full_reml$QEp),
    formula     = tryCatch(paste(deparse(full_formula), collapse = ""), error = function(e) NA_character_),
    converged   = safe_logical_scalar(full_reml$conv),
    overall_lrt = overall_lrt,
    var_lrt     = var_lrt
  )
}

# ------------------------------------------------------------------------------
# Bootstrap fitting loop
# ------------------------------------------------------------------------------
fit_bootstrap_list <- function(data_list, formula_list, label = "") {
  results   <- vector("list", length(data_list))
  total     <- length(data_list)
  n_success <- 0L
  n_fail    <- 0L

  for (i in seq_along(data_list)) {
    dat_i <- data_list[[i]]
    frm_i <- formula_list[[i]]

    progress_message(i, total, label, "bootstrap REML+ML fitting")

    full_reml <- fit_model(dat_i, frm_i, "REML")
    null_reml <- fit_model(dat_i, NULL,  "REML")
    full_ml   <- fit_model(dat_i, frm_i, "ML")
    null_ml   <- fit_model(dat_i, NULL,  "ML")

    all_ok <- !is.null(full_reml) && !is.null(null_reml) &&
      !is.null(full_ml) && !is.null(null_ml)

    if (all_ok) {
      results[[i]] <- extract_boot_results(
        full_reml, null_reml, full_ml, null_ml, dat_i, frm_i
      )
      n_success <- n_success + 1L
      progress_message(i, total, label, "-> OK")
    } else {
      message(sprintf(
        "  [WARN] iter %d: REML_f=%s REML_n=%s ML_f=%s ML_n=%s",
        i,
        ifelse(is.null(full_reml), "FAIL", "ok"),
        ifelse(is.null(null_reml), "FAIL", "ok"),
        ifelse(is.null(full_ml),   "FAIL", "ok"),
        ifelse(is.null(null_ml),   "FAIL", "ok")
      ))
      results[[i]] <- NULL
      n_fail <- n_fail + 1L
    }

    rm(full_reml, null_reml, full_ml, null_ml)
    if (i %% 50 == 0) gc()
  }

  message(sprintf("[%s bootstrap] Done: %d ok, %d failed / %d total", label, n_success, n_fail, total))
  results
}

# ==============================================================================
# 5) Run main models first
# ==============================================================================
main_results <- lapply(names(datasets), function(nm) {
  run_main_model(datasets[[nm]], nm)
})
names(main_results) <- names(datasets)
main_results <- main_results[!sapply(main_results, is.null)]

# ==============================================================================
# 6) Run bootstrap sensitivity analysis
# ==============================================================================
message("\n====================================")
message("Starting bootstrap sensitivity analysis")
message("====================================")

bootstrap_data_a <- lapply(seq_len(BOOT_N), function(i) boot.func(datasets[["all"]]))
bootstrap_data_v <- lapply(seq_len(BOOT_N), function(i) boot.func(datasets[["viruses"]]))
bootstrap_data_b <- lapply(seq_len(BOOT_N), function(i) boot.func(datasets[["bacteria"]]))
bootstrap_data_p <- lapply(seq_len(BOOT_N), function(i) boot.func(datasets[["parasites"]]))

message("Processing bootstrap data and VIF filtering...")
filtered_results_a <- process_bootstrap_data(bootstrap_data_a)
filtered_results_v <- process_bootstrap_data(bootstrap_data_v)
filtered_results_b <- process_bootstrap_data(bootstrap_data_b)
filtered_results_p <- process_bootstrap_data(bootstrap_data_p)
message("Bootstrap VIF filtering completed.")

fit_results_a <- fit_bootstrap_list(filtered_results_a$data, filtered_results_a$formulas, "all")
rm(filtered_results_a, bootstrap_data_a); gc()

fit_results_v <- fit_bootstrap_list(filtered_results_v$data, filtered_results_v$formulas, "viruses")
rm(filtered_results_v, bootstrap_data_v); gc()

fit_results_b <- fit_bootstrap_list(filtered_results_b$data, filtered_results_b$formulas, "bacteria")
rm(filtered_results_b, bootstrap_data_b); gc()

fit_results_p <- fit_bootstrap_list(filtered_results_p$data, filtered_results_p$formulas, "parasites")
rm(filtered_results_p, bootstrap_data_p); gc()

bootstrap_results <- list(
  all       = fit_results_a,
  viruses   = fit_results_v,
  bacteria  = fit_results_b,
  parasites = fit_results_p
)

saveRDS(
  bootstrap_results,
  file.path(output_path, paste0("bootstrap_results_", timestamp, ".RDS"))
)

saveRDS(
  main_results,
  file.path(output_path, paste0("main_model_objects_", timestamp, ".RDS"))
)

################################################################################
## 7) EXPORT MAIN MODEL RESULTS
################################################################################

# ------------------------------------------------------------------------------
# 7.1 Main summary
# ------------------------------------------------------------------------------
summary_list <- lapply(main_results, function(x) x$summary)
summary_df <- safe_rbind(summary_list)

wb_main_summary <- createWorkbook()
addWorksheet(wb_main_summary, "main_summary")
writeData(wb_main_summary, "main_summary", summary_df)
saveWorkbook(
  wb_main_summary,
  file.path(output_path, paste0("main_model_summary_", timestamp, ".xlsx")),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 7.2 Main coefficients
# ------------------------------------------------------------------------------
wb_main_coef <- createWorkbook()
for (nm in names(main_results)) {
  if (!is.null(main_results[[nm]]$coef) && nrow(main_results[[nm]]$coef) > 0) {
    addWorksheet(wb_main_coef, nm)
    writeData(wb_main_coef, nm, main_results[[nm]]$coef)
  }
}
saveWorkbook(
  wb_main_coef,
  file.path(output_path, paste0("main_model_coefficients_", timestamp, ".xlsx")),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 7.3 Main per-variable LRT
# ------------------------------------------------------------------------------
wb_main_var <- createWorkbook()
for (nm in names(main_results)) {
  if (!is.null(main_results[[nm]]$var_lrt) && nrow(main_results[[nm]]$var_lrt) > 0) {
    addWorksheet(wb_main_var, nm)
    writeData(wb_main_var, nm, main_results[[nm]]$var_lrt)
  }
}
saveWorkbook(
  wb_main_var,
  file.path(output_path, paste0("main_model_var_lrt_", timestamp, ".xlsx")),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 7.4 Main overall LRT
# ------------------------------------------------------------------------------
wb_main_overall <- createWorkbook()
for (nm in names(main_results)) {
  if (!is.null(main_results[[nm]]$overall_lrt) && nrow(main_results[[nm]]$overall_lrt) > 0) {
    addWorksheet(wb_main_overall, nm)
    writeData(wb_main_overall, nm, main_results[[nm]]$overall_lrt)
  }
}
saveWorkbook(
  wb_main_overall,
  file.path(output_path, paste0("main_model_overall_lrt_", timestamp, ".xlsx")),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 7.5 Main I²
# ------------------------------------------------------------------------------
wb_main_i2 <- createWorkbook()
for (nm in names(main_results)) {
  addWorksheet(wb_main_i2, nm)
  writeData(wb_main_i2, nm, main_results[[nm]]$i2)
}
saveWorkbook(
  wb_main_i2,
  file.path(output_path, paste0("main_model_i2_", timestamp, ".xlsx")),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 7.6 Main variance components
# ------------------------------------------------------------------------------
wb_main_vc <- createWorkbook()
for (nm in names(main_results)) {
  if (!is.null(main_results[[nm]]$vc) && nrow(main_results[[nm]]$vc) > 0) {
    addWorksheet(wb_main_vc, nm)
    writeData(wb_main_vc, nm, main_results[[nm]]$vc)
  }
}
saveWorkbook(
  wb_main_vc,
  file.path(output_path, paste0("main_model_variance_components_", timestamp, ".xlsx")),
  overwrite = TRUE
)

################################################################################
## 8) EXPORT BOOTSTRAP RESULTS
################################################################################

safe_get <- function(x, field) {
  if (!is.null(x) && !is.null(x[[field]])) x[[field]] else NULL
}

# ------------------------------------------------------------------------------
# 8.1 Bootstrap fixed effects
# ------------------------------------------------------------------------------
extract_coef_df <- function(fit_results, model_label) {
  rows <- lapply(seq_along(fit_results), function(i) {
    coef_df <- safe_get(fit_results[[i]], "coef")
    if (is.null(coef_df)) return(NULL)
    coef_df$boot_id <- i
    coef_df$model <- model_label
    coef_df
  })
  safe_rbind(rows)
}

coef_a <- extract_coef_df(fit_results_a, "all")
coef_v <- extract_coef_df(fit_results_v, "viruses")
coef_b <- extract_coef_df(fit_results_b, "bacteria")
coef_p <- extract_coef_df(fit_results_p, "parasites")

mean_coef <- function(coef_df) {
  if (is.null(coef_df) || nrow(coef_df) == 0) return(data.frame())

  coef_df %>%
    group_by(term) %>%
    summarise(
      n_boot        = n(),
      mean_estimate = mean(estimate, na.rm = TRUE),
      sd_estimate   = sd(estimate, na.rm = TRUE),
      q025_estimate = quantile(estimate, 0.025, na.rm = TRUE),
      q975_estimate = quantile(estimate, 0.975, na.rm = TRUE),
      mean_se       = mean(se, na.rm = TRUE),
      median_p      = median(pval, na.rm = TRUE),
      prop_sig      = mean(pval < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
}

mean_coef_a <- mean_coef(coef_a)
mean_coef_v <- mean_coef(coef_v)
mean_coef_b <- mean_coef(coef_b)
mean_coef_p <- mean_coef(coef_p)

boot_coef_map <- list(
  all       = list(raw = coef_a, mean = mean_coef_a),
  viruses   = list(raw = coef_v, mean = mean_coef_v),
  bacteria  = list(raw = coef_b, mean = mean_coef_b),
  parasites = list(raw = coef_p, mean = mean_coef_p)
)

wb_boot_coef <- createWorkbook()
for (nm in names(boot_coef_map)) {
  raw_df  <- boot_coef_map[[nm]]$raw
  mean_df <- boot_coef_map[[nm]]$mean

  if (!is.null(raw_df) && nrow(raw_df) > 0) {
    addWorksheet(wb_boot_coef, paste0("raw_", nm))
    writeData(wb_boot_coef, paste0("raw_", nm), raw_df)
  }
  if (!is.null(mean_df) && nrow(mean_df) > 0) {
    addWorksheet(wb_boot_coef, paste0("mean_", nm))
    writeData(wb_boot_coef, paste0("mean_", nm), mean_df)
  }
}
saveWorkbook(
  wb_boot_coef,
  file.path(output_path, paste0("bootstrap_coef_results_", timestamp, ".xlsx")),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 8.2 Bootstrap pseudo-R²
# ------------------------------------------------------------------------------
extract_pseudo_r2_df <- function(fit_results) {
  vals <- sapply(fit_results, function(x) {
    r2 <- safe_get(x, "pseudo_r2")
    if (is.null(r2)) NA else r2
  })
  data.frame(boot_id = seq_along(vals), pseudo_R2 = vals)
}

pr2_a <- extract_pseudo_r2_df(fit_results_a)
pr2_v <- extract_pseudo_r2_df(fit_results_v)
pr2_b <- extract_pseudo_r2_df(fit_results_b)
pr2_p <- extract_pseudo_r2_df(fit_results_p)

wb_boot_pr2 <- createWorkbook()
for (nm in list(
  list("all", pr2_a), list("viruses", pr2_v),
  list("bacteria", pr2_b), list("parasites", pr2_p)
)) {
  if (nrow(nm[[2]]) > 0) {
    addWorksheet(wb_boot_pr2, nm[[1]])
    writeData(wb_boot_pr2, nm[[1]], nm[[2]])
  }
}
saveWorkbook(
  wb_boot_pr2,
  file.path(output_path, paste0("bootstrap_pseudo_r2_", timestamp, ".xlsx")),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 8.3 Bootstrap I²
# ------------------------------------------------------------------------------
extract_i2_df <- function(fit_results) {
  rows <- lapply(seq_along(fit_results), function(i) {
    i2_df <- safe_get(fit_results[[i]], "i2")
    if (is.null(i2_df)) return(NULL)
    cbind(boot_id = i, i2_df)
  })
  safe_rbind(rows)
}

i2_a <- extract_i2_df(fit_results_a)
i2_v <- extract_i2_df(fit_results_v)
i2_b <- extract_i2_df(fit_results_b)
i2_p <- extract_i2_df(fit_results_p)

wb_boot_i2 <- createWorkbook()
for (nm in list(
  list("all", i2_a), list("viruses", i2_v),
  list("bacteria", i2_b), list("parasites", i2_p)
)) {
  if (!is.null(nm[[2]]) && nrow(nm[[2]]) > 0) {
    addWorksheet(wb_boot_i2, nm[[1]])
    writeData(wb_boot_i2, nm[[1]], nm[[2]])
  }
}
saveWorkbook(
  wb_boot_i2,
  file.path(output_path, paste0("bootstrap_i2_results_", timestamp, ".xlsx")),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 8.4 Bootstrap QM / QE
# ------------------------------------------------------------------------------
extract_qmqe_df <- function(fit_results) {
  rows <- lapply(seq_along(fit_results), function(i) {
    x <- fit_results[[i]]
    if (is.null(x)) return(NULL)
    data.frame(
      boot_id = i,
      QM   = safe_get(x, "QM"),
      QMdf = safe_get(x, "QMdf"),
      QMp  = safe_get(x, "QMp"),
      QE   = safe_get(x, "QE"),
      QEdf = safe_get(x, "QEdf"),
      QEp  = safe_get(x, "QEp")
    )
  })
  safe_rbind(rows)
}

qmqe_a <- extract_qmqe_df(fit_results_a)
qmqe_v <- extract_qmqe_df(fit_results_v)
qmqe_b <- extract_qmqe_df(fit_results_b)
qmqe_p <- extract_qmqe_df(fit_results_p)

wb_boot_qmqe <- createWorkbook()
for (nm in list(
  list("all", qmqe_a), list("viruses", qmqe_v),
  list("bacteria", qmqe_b), list("parasites", qmqe_p)
)) {
  if (!is.null(nm[[2]]) && nrow(nm[[2]]) > 0) {
    addWorksheet(wb_boot_qmqe, nm[[1]])
    writeData(wb_boot_qmqe, nm[[1]], nm[[2]])
  }
}
saveWorkbook(
  wb_boot_qmqe,
  file.path(output_path, paste0("bootstrap_qmqe_results_", timestamp, ".xlsx")),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 8.5 Bootstrap variance components
# ------------------------------------------------------------------------------
extract_vc_df <- function(fit_results) {
  rows <- lapply(seq_along(fit_results), function(i) {
    vc <- safe_get(fit_results[[i]], "vc")
    if (is.null(vc)) return(NULL)
    cbind(boot_id = i, vc)
  })
  safe_rbind(rows)
}

vc_a <- extract_vc_df(fit_results_a)
vc_v <- extract_vc_df(fit_results_v)
vc_b <- extract_vc_df(fit_results_b)
vc_p <- extract_vc_df(fit_results_p)

wb_boot_vc <- createWorkbook()
for (nm in list(
  list("all", vc_a), list("viruses", vc_v),
  list("bacteria", vc_b), list("parasites", vc_p)
)) {
  if (!is.null(nm[[2]]) && nrow(nm[[2]]) > 0) {
    addWorksheet(wb_boot_vc, nm[[1]])
    writeData(wb_boot_vc, nm[[1]], nm[[2]])
  }
}
saveWorkbook(
  wb_boot_vc,
  file.path(output_path, paste0("bootstrap_variance_components_", timestamp, ".xlsx")),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 8.6 Bootstrap per-variable LRT
# ------------------------------------------------------------------------------
extract_var_lrt_df <- function(fit_results, model_label) {
  rows <- lapply(seq_along(fit_results), function(i) {
    vlrt <- safe_get(fit_results[[i]], "var_lrt")
    if (is.null(vlrt)) return(NULL)
    vlrt$boot_id <- i
    vlrt$model <- model_label
    vlrt
  })
  safe_rbind(rows)
}

var_lrt_a <- extract_var_lrt_df(fit_results_a, "all")
var_lrt_v <- extract_var_lrt_df(fit_results_v, "viruses")
var_lrt_b <- extract_var_lrt_df(fit_results_b, "bacteria")
var_lrt_p <- extract_var_lrt_df(fit_results_p, "parasites")

mean_var_lrt <- function(vlrt_df) {
  if (is.null(vlrt_df) || nrow(vlrt_df) == 0) return(data.frame())

  vlrt_df %>%
    group_by(variable) %>%
    summarise(
      n_boot     = n(),
      median_LRT = median(LRT, na.rm = TRUE),
      median_p   = median(p_value, na.rm = TRUE),
      prop_sig   = mean(p_value < 0.05, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(median_p)
}

mean_var_lrt_a <- mean_var_lrt(var_lrt_a)
mean_var_lrt_v <- mean_var_lrt(var_lrt_v)
mean_var_lrt_b <- mean_var_lrt(var_lrt_b)
mean_var_lrt_p <- mean_var_lrt(var_lrt_p)

boot_varlrt_map <- list(
  all       = list(raw = var_lrt_a, mean = mean_var_lrt_a),
  viruses   = list(raw = var_lrt_v, mean = mean_var_lrt_v),
  bacteria  = list(raw = var_lrt_b, mean = mean_var_lrt_b),
  parasites = list(raw = var_lrt_p, mean = mean_var_lrt_p)
)

wb_boot_var_lrt <- createWorkbook()
for (nm in names(boot_varlrt_map)) {
  raw_df  <- boot_varlrt_map[[nm]]$raw
  mean_df <- boot_varlrt_map[[nm]]$mean

  if (!is.null(raw_df) && nrow(raw_df) > 0) {
    addWorksheet(wb_boot_var_lrt, paste0("raw_", nm))
    writeData(wb_boot_var_lrt, paste0("raw_", nm), raw_df)
  }
  if (!is.null(mean_df) && nrow(mean_df) > 0) {
    addWorksheet(wb_boot_var_lrt, paste0("mean_", nm))
    writeData(wb_boot_var_lrt, paste0("mean_", nm), mean_df)
  }
}
saveWorkbook(
  wb_boot_var_lrt,
  file.path(output_path, paste0("bootstrap_per_variable_lrt_", timestamp, ".xlsx")),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 8.7 Bootstrap overall LRT
# ------------------------------------------------------------------------------
extract_overall_lrt_df <- function(fit_results, model_label) {
  rows <- lapply(seq_along(fit_results), function(i) {
    olrt <- safe_get(fit_results[[i]], "overall_lrt")
    if (is.null(olrt)) return(NULL)
    olrt$boot_id <- i
    olrt$model <- model_label
    olrt
  })
  safe_rbind(rows)
}

overall_lrt_a <- extract_overall_lrt_df(fit_results_a, "all")
overall_lrt_v <- extract_overall_lrt_df(fit_results_v, "viruses")
overall_lrt_b <- extract_overall_lrt_df(fit_results_b, "bacteria")
overall_lrt_p <- extract_overall_lrt_df(fit_results_p, "parasites")

wb_boot_overall <- createWorkbook()
for (nm in list(
  list("all", overall_lrt_a), list("viruses", overall_lrt_v),
  list("bacteria", overall_lrt_b), list("parasites", overall_lrt_p)
)) {
  if (!is.null(nm[[2]]) && nrow(nm[[2]]) > 0) {
    addWorksheet(wb_boot_overall, nm[[1]])
    writeData(wb_boot_overall, nm[[1]], nm[[2]])
  }
}
saveWorkbook(
  wb_boot_overall,
  file.path(output_path, paste0("bootstrap_overall_lrt_", timestamp, ".xlsx")),
  overwrite = TRUE
)

# ------------------------------------------------------------------------------
# 8.8 Funnel / Effect-SE plots from bootstrap coefficients
# ------------------------------------------------------------------------------
make_plot_df <- function(coef_df, label) {
  if (is.null(coef_df) || nrow(coef_df) == 0) return(NULL)
  coef_df$Model <- label
  coef_df[, c("term", "estimate", "se", "Model")]
}

plot_list <- Filter(Negate(is.null), list(
  make_plot_df(coef_a, "All pathogens"),
  make_plot_df(coef_v, "Viruses"),
  make_plot_df(coef_b, "Bacteria"),
  make_plot_df(coef_p, "Parasites")
))

if (length(plot_list) > 0) {
  df_plot <- do.call(rbind, plot_list)
  df_plot$Model <- factor(
    df_plot$Model,
    levels = c("All pathogens", "Viruses", "Bacteria", "Parasites")
  )
  df_plot$Precision <- 1 / df_plot$se

  pal <- c(
    "All pathogens" = "#9d9fc9",
    "Viruses"       = "#dc846a",
    "Bacteria"      = "#d49e71",
    "Parasites"     = "#b395bd"
  )

  p_effect_se <- ggplot(df_plot, aes(x = estimate, y = se, color = Model)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#808080") +
    scale_color_manual(values = pal) +
    facet_wrap(~ Model) +
    theme_classic() +
    labs(x = "Effect size (yi)", y = "Standard error") +
    theme(legend.position = "bottom")

  p_funnel <- ggplot(df_plot, aes(x = estimate, y = Precision, color = Model)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#808080") +
    scale_color_manual(values = pal) +
    facet_wrap(~ Model) +
    theme_classic() +
    labs(x = "Effect size (yi)", y = "Precision (1/SE)") +
    theme(legend.position = "bottom")

  pdf(file.path(output_path, paste0("bootstrap_effect_SE_plot_", timestamp, ".pdf")), width = 10, height = 8)
  print(p_effect_se)
  dev.off()

  pdf(file.path(output_path, paste0("bootstrap_funnel_plot_", timestamp, ".pdf")), width = 10, height = 8)
  print(p_funnel)
  dev.off()
}

################################################################################
## 9) Console summary
################################################################################
message("\n==============================")
message("Main model summaries")
message("==============================")
print(summary_df)

for (nm in names(main_results)) {
  message("\n--- ", nm, ": MAIN per-variable LRT ---")
  print(main_results[[nm]]$var_lrt)
}

for (nm in c("all", "viruses", "bacteria", "parasites")) {
  mean_df <- switch(
    nm,
    all       = mean_var_lrt_a,
    viruses   = mean_var_lrt_v,
    bacteria  = mean_var_lrt_b,
    parasites = mean_var_lrt_p
  )
  message("\n--- ", nm, ": BOOTSTRAP per-variable robustness ---")
  if (!is.null(mean_df) && nrow(mean_df) > 0) print(mean_df, n = Inf)
}

message("\nIntegrated main-model + bootstrap results saved to: ", output_path)