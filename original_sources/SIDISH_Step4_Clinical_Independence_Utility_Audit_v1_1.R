#!/usr/bin/env Rscript
###############################################################################
# SIDISH Step 4: Clinical independence and clinical utility audit v1.1
#
# Purpose
# -------
# Audit whether the locked Step 2 v2.6 / Step 3 v1.4 SIDISH risk score has
# clinical independence and clinical utility.
#
# Locked-model rule:
#   - Do NOT change genes.
#   - Do NOT change coefficients.
#   - Do NOT change Step 2 development cutoff.
#   - Do NOT re-optimize cutoff in any cohort.
#
# Inputs:
#   Step 2 v2.6:
#     - 99_locked_model/SIDISH_Step2_locked_development_risk_scores.tsv
#     - 99_locked_model/SIDISH_Step2_locked_model_coefficients.tsv
#     - 99_locked_model/SIDISH_Step2_locked_model_summary.tsv
#
#   Step 3 v1.4:
#     - 99_final_external_validation/SIDISH_Step3_GSE16011_no_refit_risk_scores.tsv
#     - 99_final_external_validation/SIDISH_Step3_GSE16011_no_refit_external_validation_summary.tsv
#
#   Clinical files:
#     - TCGA_clinical_molecular.tsv
#     - CGGA325_clinical_molecular.tsv
#     - CGGA693_clinical_molecular.tsv
#     - GSE_glioma_clinical_molecular.tsv
#
# Outputs:
#   - univariable and multivariable Cox audit tables
#   - clinical covariate availability audit
#   - subgroup/stratified risk effect audit
#   - optional time-dependent binary DCA-style net benefit at 12/36/60 months
#   - optional calibration-by-risk-quantile audit
#   - publication-style figures consistent with Step 1/2/3
#
# Required packages:
#   install.packages(c("data.table", "survival", "ggplot2"))
#
# Optional packages:
#   install.packages(c("rms", "pec", "riskRegression"))
#
# Run:
#   source("G:/NEW/The2/SIDISH_Step4_Clinical_Independence_Utility_Audit_v1_1.R")
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = 3600)
set.seed(20260522)

###############################################################################
# 0. CONFIG
###############################################################################

CONFIG <- list(
  output_dir = "G:/NEW/The2/SIDISH_Step4_Clinical_Independence_Utility_Audit_v1_1",

  # Locked development model.
  step2_root = "G:/NEW/The2/SIDISH_Step2_ML175_DirectionStable_Model_Development_v2_6_strict_all3",
  step2_risk_file = "G:/NEW/The2/SIDISH_Step2_ML175_DirectionStable_Model_Development_v2_6_strict_all3/99_locked_model/SIDISH_Step2_locked_development_risk_scores.tsv",
  step2_coef_file = "G:/NEW/The2/SIDISH_Step2_ML175_DirectionStable_Model_Development_v2_6_strict_all3/99_locked_model/SIDISH_Step2_locked_model_coefficients.tsv",
  step2_summary_file = "G:/NEW/The2/SIDISH_Step2_ML175_DirectionStable_Model_Development_v2_6_strict_all3/99_locked_model/SIDISH_Step2_locked_model_summary.tsv",

  # External no-refit validation.
  step3_root = "G:/NEW/The2/SIDISH_Step3_GSE16011_no_refit_external_validation_v1_4",
  step3_risk_file = "G:/NEW/The2/SIDISH_Step3_GSE16011_no_refit_external_validation_v1_4/99_final_external_validation/SIDISH_Step3_GSE16011_no_refit_risk_scores.tsv",
  step3_summary_file = "G:/NEW/The2/SIDISH_Step3_GSE16011_no_refit_external_validation_v1_4/99_final_external_validation/SIDISH_Step3_GSE16011_no_refit_external_validation_summary.tsv",

  # Clinical files. Adjust here only if your local paths differ.
  clinical_files = list(
    TCGA = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/TCGA/TCGA_clinical_molecular.tsv",
    CGGA325 = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/CGGA_325/CGGA325_clinical_molecular.tsv",
    CGGA693 = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/CGGA_693/CGGA693_clinical_molecular.tsv",
    GSE16011 = "G:/NEW/The2/External_Cohorts/GSE16011/GSE_glioma_clinical_molecular.tsv"
  ),

  # Clinical covariate handling.
  min_complete_cases_multivariable = 50,
  min_events_multivariable = 20,
  min_subgroup_n = 30,
  min_subgroup_events = 10,
  covariate_candidates = c("Age", "WHO_grade", "IDH_status", "Sex", "1p19q_status"),

  # Utility audit time points.
  utility_time_unit = "auto", # auto/days/months
  utility_times_days = c(365, 1095, 1825),
  utility_times_months = c(12, 36, 60),
  dca_thresholds = seq(0.05, 0.80, by = 0.05),
  calibration_quantile_groups = 4,

  # Plotting.
  make_publication_figures = TRUE,
  plot_base_family = "Arial",
  plot_base_size = 8,
  plot_linewidth = 0.5,
  plot_point_size = 2.0,
  p_value_floor_for_display = 2.2e-16,

  # Safety.
  objective_reporting_note = "raw_metrics_only_no_qualitative_label"
)

###############################################################################
# 1. Utilities
###############################################################################

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

norm_path <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  gsub("\\\\", "/", x)
}

safe_dir_create <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  invisible(x)
}

OUT <- list(
  root = CONFIG$output_dir,
  audit = file.path(CONFIG$output_dir, "00_audit"),
  data = file.path(CONFIG$output_dir, "01_clinical_merged_data"),
  cox = file.path(CONFIG$output_dir, "02_clinical_independence_cox"),
  utility = file.path(CONFIG$output_dir, "03_clinical_utility"),
  fig = file.path(CONFIG$output_dir, "04_publication_figures"),
  final = file.path(CONFIG$output_dir, "99_final_clinical_audit")
)
invisible(lapply(OUT, safe_dir_create))

load_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required R package: ", pkg, "\nInstall with install.packages('", pkg, "')")
  }
  invisible(TRUE)
}
load_or_stop("data.table")
load_or_stop("survival")
load_or_stop("ggplot2")
suppressPackageStartupMessages(library(survival))
Surv <- survival::Surv

write_tsv <- function(dt, file) {
  safe_dir_create(dirname(file))
  data.table::fwrite(data.table::as.data.table(dt), file, sep = "\t", quote = FALSE, na = "NA")
}

write_text <- function(x, file) {
  safe_dir_create(dirname(file))
  writeLines(x, con = file, useBytes = TRUE)
}

clean_names <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("[.]+", "_", x)
  x
}

read_table_auto <- function(file) {
  file <- norm_path(file)
  if (!file.exists(file)) stop("File not found: ", file)
  dt <- tryCatch(
    data.table::fread(file, sep = "auto", data.table = TRUE, showProgress = FALSE, quote = "", fill = TRUE),
    error = function(e) {
      data.table::fread(file, sep = "\t", data.table = TRUE, showProgress = FALSE, quote = "", fill = TRUE)
    }
  )
  data.table::setnames(dt, clean_names(names(dt)))
  dt
}

as_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub(",", "", x)
  suppressWarnings(as.numeric(x))
}

standardize_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.", "-", x)
  x <- gsub("_", "-", x)
  x
}

first_existing_col <- function(dt, candidates, regex = FALSE, ignore.case = TRUE) {
  nms <- names(dt)
  if (!regex) {
    hit <- candidates[candidates %in% nms]
    if (length(hit) > 0) return(hit[1])
  }
  for (pat in candidates) {
    hit <- grep(pat, nms, value = TRUE, ignore.case = ignore.case)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

format_p_for_plot <- function(p, floor_value = CONFIG$p_value_floor_for_display) {
  p <- suppressWarnings(as.numeric(p))
  if (!is.finite(p)) return("p = NA")
  if (p < floor_value) return(paste0("p < ", format(floor_value, scientific = TRUE, digits = 2)))
  if (p < 0.001) return(paste0("p = ", format(p, scientific = TRUE, digits = 2)))
  paste0("p = ", sprintf("%.3f", p))
}

make_unique_factor_levels <- function(x, reverse = TRUE) {
  x <- as.character(x)
  lev <- unique(x)
  if (isTRUE(reverse)) lev <- rev(lev)
  lev
}

SIDISH_COLORS <- c(
  red = "#E64B35FF",
  blue = "#4DBBD5FF",
  green = "#00A087FF",
  navy = "#3C5488FF",
  orange = "#F39B7FFF",
  purple = "#8491B4FF",
  teal = "#91D1C2FF",
  brown = "#7E6148FF",
  grey = "#8C8C8C"
)
pal <- function(x) unname(SIDISH_COLORS[x])

pub_theme <- function(base_size = CONFIG$plot_base_size, base_family = CONFIG$plot_base_family) {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, colour = "black", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(colour = "#333333", size = base_size),
      axis.title = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      panel.grid.major = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", linewidth = CONFIG$plot_linewidth),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = CONFIG$plot_linewidth),
      legend.title = ggplot2::element_text(face = "bold", colour = "black"),
      legend.text = ggplot2::element_text(colour = "black"),
      panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = CONFIG$plot_linewidth),
      plot.margin = grid::unit(c(5.5, 12, 5.5, 5.5), "pt")
    )
}

save_pub_plot <- function(p, filename_stub, width, height) {
  safe_dir_create(OUT$fig)
  pdf_file <- file.path(OUT$fig, paste0(filename_stub, ".pdf"))
  png_file <- file.path(OUT$fig, paste0(filename_stub, ".png"))
  pdf_status <- "not_run"
  pdf_error <- NA_character_

  tryCatch({
    ggplot2::ggsave(pdf_file, plot = p, width = width, height = height, units = "in",
                    device = grDevices::cairo_pdf, family = CONFIG$plot_base_family, bg = "white")
    pdf_status <- paste0("ok_cairo_pdf_family_", CONFIG$plot_base_family)
  }, error = function(e1) {
    pdf_error <<- conditionMessage(e1)
    p2 <- p + ggplot2::theme(text = ggplot2::element_text(family = "sans"))
    ggplot2::ggsave(pdf_file, plot = p2, width = width, height = height, units = "in",
                    device = grDevices::cairo_pdf, family = "sans", bg = "white")
    pdf_status <<- "ok_cairo_pdf_family_sans_fallback"
  })

  ggplot2::ggsave(png_file, plot = p, width = width, height = height, units = "in", dpi = 600, bg = "white")

  log_file <- file.path(OUT$fig, "figure_export_status.tsv")
  row <- data.table::data.table(
    figure = filename_stub,
    pdf_file = pdf_file,
    png_file = png_file,
    requested_family = CONFIG$plot_base_family,
    pdf_status = pdf_status,
    pdf_error = pdf_error,
    width = width,
    height = height
  )
  if (file.exists(log_file)) {
    old <- data.table::fread(log_file, sep = "\t", data.table = TRUE, showProgress = FALSE)
    row <- data.table::rbindlist(list(old, row), fill = TRUE)
  }
  write_tsv(row, log_file)
}

detect_survival_time_unit <- function(time, configured = CONFIG$utility_time_unit) {
  configured <- tolower(as.character(configured))
  if (configured %in% c("days", "months")) return(configured)
  tt <- suppressWarnings(as.numeric(time))
  tt <- tt[is.finite(tt) & tt > 0]
  if (length(tt) == 0) return("unknown")
  if (stats::quantile(tt, 0.95, na.rm = TRUE) <= 240 && max(tt, na.rm = TRUE) <= 400) return("months")
  "days"
}

choose_utility_times <- function(time, unit) {
  tt <- suppressWarnings(as.numeric(time))
  tt <- tt[is.finite(tt) & tt > 0]
  if (length(tt) == 0) return(numeric(0))
  cand <- if (identical(unit, "months")) CONFIG$utility_times_months else CONFIG$utility_times_days
  cand[cand > min(tt, na.rm = TRUE) & cand < max(tt, na.rm = TRUE)]
}

###############################################################################
# 2. Read risk scores and clinical metadata
###############################################################################

msg("SIDISH Step 4 clinical independence and utility audit started.")

required_files <- c(CONFIG$step2_risk_file, CONFIG$step2_coef_file, CONFIG$step2_summary_file,
                    CONFIG$step3_risk_file, CONFIG$step3_summary_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) stop("Missing required locked-risk files:\n", paste(missing_files, collapse = "\n"))

step2_risk <- read_table_auto(CONFIG$step2_risk_file)
step3_risk <- read_table_auto(CONFIG$step3_risk_file)
step2_coef <- read_table_auto(CONFIG$step2_coef_file)
step2_summary <- read_table_auto(CONFIG$step2_summary_file)
step3_summary <- read_table_auto(CONFIG$step3_summary_file)

# Harmonize risk tables.
step2_risk[, dataset_stage := "development"]
step2_risk[, cohort_final := as.character(cohort)]
step2_risk[, sample_id := standardize_id(sample_id)]

step3_risk[, dataset_stage := "external"]
if (!"cohort" %in% names(step3_risk)) step3_risk[, cohort := if ("external_cohort" %in% names(step3_risk)) external_cohort else "GSE16011"]
step3_risk[, cohort_final := as.character(cohort)]
step3_risk[, sample_id := standardize_id(sample_id)]

risk_common_cols <- c("sample_id", "cohort_final", "dataset_stage", "OS_time", "OS_event",
                      "locked_risk_score", "locked_risk_group")
risk_all <- data.table::rbindlist(list(
  step2_risk[, ..risk_common_cols],
  step3_risk[, ..risk_common_cols]
), fill = TRUE)
risk_all[, OS_time := as_numeric_safe(OS_time)]
risk_all[, OS_event := as_numeric_safe(OS_event)]
risk_all[, locked_risk_score := as_numeric_safe(locked_risk_score)]
risk_all <- risk_all[is.finite(OS_time) & OS_time > 0 & OS_event %in% c(0, 1) & is.finite(locked_risk_score)]

write_tsv(risk_all, file.path(OUT$data, "Step4_all_locked_risk_scores_for_clinical_audit.tsv"))

###############################################################################
# 3. Clinical table standardization
###############################################################################

standardize_grade <- function(x) {
  xx <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(xx))
  out[grepl("4|IV|GBM|GLIOBLASTOMA", xx)] <- "G4"
  out[grepl("3|III", xx) & is.na(out)] <- "G3"
  out[grepl("2|II", xx) & is.na(out)] <- "G2"
  out
}

standardize_idh <- function(x) {
  xx <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(xx))
  out[grepl("MUT|MUTANT|IDH1|IDH2", xx) & !grepl("WILD|WT", xx)] <- "IDHmut"
  out[grepl("WT|WILD", xx)] <- "IDHwt"
  out[xx %in% c("1", "MUT")] <- "IDHmut"
  out[xx %in% c("0", "WT")] <- "IDHwt"
  out
}

standardize_sex <- function(x) {
  xx <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(xx))
  out[grepl("^M|MALE", xx)] <- "Male"
  out[grepl("^F|FEMALE", xx)] <- "Female"
  out
}

standardize_1p19q <- function(x) {
  xx <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(xx))
  out[grepl("CODEL|CO-DEL|CO_DEL|DELETED|LOSS", xx)] <- "codel"
  out[grepl("NON|INTACT|NO", xx)] <- "noncodel"
  out
}

standardize_clinical_table <- function(file, cohort) {
  if (!file.exists(file)) {
    return(data.table::data.table(
      cohort_final = cohort,
      sample_id = character(),
      Age = numeric(),
      WHO_grade = character(),
      IDH_status = character(),
      Sex = character(),
      `1p19q_status` = character()
    ))
  }

  dt <- read_table_auto(file)
  if (nrow(dt) == 0) stop("Clinical table empty: ", file)

  sample_col <- first_existing_col(dt, c("^sample_id$", "^sample$", "geo_accession", "GSM", "barcode", "case", "patient", "^ID$"), regex = TRUE)
  if (is.na(sample_col)) sample_col <- names(dt)[1]

  age_col <- first_existing_col(dt, c("^Age$", "^age$"), regex = TRUE)
  grade_col <- first_existing_col(dt, c("WHO_grade", "grade", "histology"), regex = TRUE)
  idh_col <- first_existing_col(dt, c("IDH_status", "IDH", "idh"), regex = TRUE)
  sex_col <- first_existing_col(dt, c("^Sex$", "^sex$", "gender"), regex = TRUE)
  p19q_col <- first_existing_col(dt, c("1p19q", "1p_19q", "codeletion", "co_deletion", "codel"), regex = TRUE)

  out <- data.table::data.table(
    cohort_final = cohort,
    sample_id_raw = as.character(dt[[sample_col]]),
    sample_id = standardize_id(dt[[sample_col]]),
    Age = if (!is.na(age_col)) as_numeric_safe(dt[[age_col]]) else NA_real_,
    WHO_grade_raw = if (!is.na(grade_col)) as.character(dt[[grade_col]]) else NA_character_,
    IDH_status_raw = if (!is.na(idh_col)) as.character(dt[[idh_col]]) else NA_character_,
    Sex_raw = if (!is.na(sex_col)) as.character(dt[[sex_col]]) else NA_character_,
    `1p19q_status_raw` = if (!is.na(p19q_col)) as.character(dt[[p19q_col]]) else NA_character_
  )
  out[, WHO_grade := standardize_grade(WHO_grade_raw)]
  out[, IDH_status := standardize_idh(IDH_status_raw)]
  out[, Sex := standardize_sex(Sex_raw)]
  out[, `1p19q_status` := standardize_1p19q(`1p19q_status_raw`)]
  unique(out, by = "sample_id")
}

clinical_list <- list()
clinical_status <- list()
i <- 1L
for (coh in names(CONFIG$clinical_files)) {
  file <- CONFIG$clinical_files[[coh]]
  exists <- file.exists(file)
  clinical_status[[i]] <- data.table::data.table(cohort_final = coh, clinical_file = file, exists = exists)
  i <- i + 1L
  if (exists) clinical_list[[coh]] <- standardize_clinical_table(file, coh)
}
clinical_all <- data.table::rbindlist(clinical_list, fill = TRUE)
write_tsv(data.table::rbindlist(clinical_status, fill = TRUE), file.path(OUT$audit, "Step4_clinical_file_status.tsv"))
write_tsv(clinical_all, file.path(OUT$data, "Step4_standardized_clinical_metadata.tsv"))

merged <- merge(
  risk_all,
  clinical_all[, .(cohort_final, sample_id, Age, WHO_grade, IDH_status, Sex, `1p19q_status`)],
  by = c("cohort_final", "sample_id"),
  all.x = TRUE
)

merged[, risk_z := as.numeric(scale(locked_risk_score))]
merged[, risk_group01 := ifelse(locked_risk_group == "High", 1, ifelse(locked_risk_group == "Low", 0, NA_real_))]
merged[, grade_binary_G4 := ifelse(WHO_grade == "G4", "G4", ifelse(!is.na(WHO_grade), "G2_G3", NA_character_))]
merged[, IDH_binary := IDH_status]
merged[, Age_z := as.numeric(scale(Age))]

write_tsv(merged, file.path(OUT$data, "Step4_clinical_merged_locked_risk_dataset.tsv"))

cov_avail <- merged[, .(
  n = .N,
  events = sum(OS_event == 1, na.rm = TRUE),
  Age_nonmissing = sum(!is.na(Age)),
  WHO_grade_nonmissing = sum(!is.na(WHO_grade)),
  IDH_status_nonmissing = sum(!is.na(IDH_status)),
  Sex_nonmissing = sum(!is.na(Sex)),
  `1p19q_status_nonmissing` = sum(!is.na(`1p19q_status`)),
  risk_nonmissing = sum(is.finite(locked_risk_score))
), by = .(dataset_stage, cohort_final)]
write_tsv(cov_avail, file.path(OUT$audit, "Step4_clinical_covariate_availability_by_cohort.tsv"))

###############################################################################
# 4. Cox model helpers
###############################################################################

cox_term_table <- function(fit, model_label, cohort_label, dataset_stage, n, events) {
  if (is.null(fit)) return(data.table::data.table())
  sm <- summary(fit)
  cf <- as.data.frame(sm$coefficients)
  ci <- as.data.frame(sm$conf.int)
  terms <- rownames(cf)
  data.table::data.table(
    dataset_stage = dataset_stage,
    cohort = cohort_label,
    model = model_label,
    term = terms,
    n = n,
    events = events,
    beta = cf[, "coef"],
    HR = ci[, "exp(coef)"],
    lower95 = ci[, "lower .95"],
    upper95 = ci[, "upper .95"],
    p = cf[, "Pr(>|z|)"],
    concordance = as.numeric(sm$concordance[1])
  )
}

fit_cox_safe <- function(dt, formula, model_label, cohort_label, dataset_stage) {
  vars <- all.vars(formula)
  dd <- as.data.frame(dt)
  dd <- dd[stats::complete.cases(dd[, vars, drop = FALSE]), , drop = FALSE]
  n <- nrow(dd)
  events <- sum(dd$OS_event == 1, na.rm = TRUE)
  if (n < CONFIG$min_complete_cases_multivariable || events < CONFIG$min_events_multivariable) {
    return(list(
      fit = NULL,
      table = data.table::data.table(
        dataset_stage = dataset_stage,
        cohort = cohort_label,
        model = model_label,
        term = NA_character_,
        n = n,
        events = events,
        beta = NA_real_,
        HR = NA_real_,
        lower95 = NA_real_,
        upper95 = NA_real_,
        p = NA_real_,
        concordance = NA_real_,
        status = "skipped_insufficient_complete_cases_or_events"
      )
    ))
  }

  fit <- tryCatch(survival::coxph(formula, data = dd), error = function(e) {
    attr(e, "sidish_error") <- conditionMessage(e)
    e
  })
  if (inherits(fit, "error")) {
    return(list(
      fit = NULL,
      table = data.table::data.table(
        dataset_stage = dataset_stage,
        cohort = cohort_label,
        model = model_label,
        term = NA_character_,
        n = n,
        events = events,
        beta = NA_real_,
        HR = NA_real_,
        lower95 = NA_real_,
        upper95 = NA_real_,
        p = NA_real_,
        concordance = NA_real_,
        status = paste0("cox_error: ", attr(fit, "sidish_error"))
      )
    ))
  }
  tab <- cox_term_table(fit, model_label, cohort_label, dataset_stage, n, events)
  tab[, status := "ok"]
  list(fit = fit, table = tab)
}

make_formula_set <- function(dt, pooled = FALSE) {
  forms <- list()
  forms[["risk_only_continuous"]] <- survival::Surv(OS_time, OS_event) ~ risk_z
  forms[["risk_only_high_low"]] <- survival::Surv(OS_time, OS_event) ~ risk_group01

  if (sum(!is.na(dt$Age_z)) >= CONFIG$min_complete_cases_multivariable) {
    forms[["risk_plus_age"]] <- survival::Surv(OS_time, OS_event) ~ risk_z + Age_z
  }
  if (sum(!is.na(dt$grade_binary_G4)) >= CONFIG$min_complete_cases_multivariable &&
      data.table::uniqueN(dt$grade_binary_G4[!is.na(dt$grade_binary_G4)]) >= 2) {
    forms[["risk_plus_grade"]] <- survival::Surv(OS_time, OS_event) ~ risk_z + grade_binary_G4
  }
  if (sum(!is.na(dt$IDH_binary)) >= CONFIG$min_complete_cases_multivariable &&
      data.table::uniqueN(dt$IDH_binary[!is.na(dt$IDH_binary)]) >= 2) {
    forms[["risk_plus_IDH"]] <- survival::Surv(OS_time, OS_event) ~ risk_z + IDH_binary
  }

  covs <- c()
  if (sum(!is.na(dt$Age_z)) >= CONFIG$min_complete_cases_multivariable) covs <- c(covs, "Age_z")
  if (sum(!is.na(dt$grade_binary_G4)) >= CONFIG$min_complete_cases_multivariable &&
      data.table::uniqueN(dt$grade_binary_G4[!is.na(dt$grade_binary_G4)]) >= 2) covs <- c(covs, "grade_binary_G4")
  if (sum(!is.na(dt$IDH_binary)) >= CONFIG$min_complete_cases_multivariable &&
      data.table::uniqueN(dt$IDH_binary[!is.na(dt$IDH_binary)]) >= 2) covs <- c(covs, "IDH_binary")
  if (length(covs) > 0) {
    rhs <- paste(c("risk_z", covs), collapse = " + ")
    forms[["risk_plus_available_clinical_covariates"]] <- stats::as.formula(paste0("survival::Surv(OS_time, OS_event) ~ ", rhs))
  }

  if (isTRUE(pooled)) {
    forms2 <- forms
    for (nm in names(forms)) {
      forms2[[paste0(nm, "_strata_cohort")]] <- stats::as.formula(paste0(deparse(forms[[nm]]), " + strata(cohort_final)"))
    }
    forms <- forms2
  }

  forms
}

###############################################################################
# 5. Clinical independence Cox audit
###############################################################################

cox_results <- list()
idx <- 1L

# Per cohort.
for (coh in unique(merged$cohort_final)) {
  dt <- merged[cohort_final == coh]
  fset <- make_formula_set(dt, pooled = FALSE)
  for (nm in names(fset)) {
    res <- fit_cox_safe(dt, fset[[nm]], nm, coh, unique(dt$dataset_stage)[1])
    cox_results[[idx]] <- res$table
    idx <- idx + 1L
  }
}

# Pooled development with cohort strata.
dt_dev <- merged[dataset_stage == "development"]
fset_dev <- make_formula_set(dt_dev, pooled = TRUE)
for (nm in names(fset_dev)) {
  res <- fit_cox_safe(dt_dev, fset_dev[[nm]], nm, "DEVELOPMENT_POOLED", "development")
  cox_results[[idx]] <- res$table
  idx <- idx + 1L
}

# All cohorts pooled with cohort strata, for audit only.
dt_all <- merged
fset_all <- make_formula_set(dt_all, pooled = TRUE)
for (nm in names(fset_all)) {
  res <- fit_cox_safe(dt_all, fset_all[[nm]], nm, "ALL_POOLED_WITH_EXTERNAL", "development_plus_external")
  cox_results[[idx]] <- res$table
  idx <- idx + 1L
}

cox_table <- data.table::rbindlist(cox_results, fill = TRUE)
write_tsv(cox_table, file.path(OUT$cox, "Step4_clinical_independence_Cox_all_models.tsv"))

risk_term_table <- cox_table[grepl("^risk_z$|^risk_group01$", term)]
risk_term_table[, reporting_note := CONFIG$objective_reporting_note]
write_tsv(risk_term_table, file.path(OUT$cox, "Step4_clinical_independence_risk_term_summary.tsv"))

###############################################################################
# 6. Subgroup audit
###############################################################################

subgroup_results <- list()
si <- 1L
subgroup_vars <- c("WHO_grade", "grade_binary_G4", "IDH_binary", "Sex", "1p19q_status")
for (sg in subgroup_vars) {
  if (!sg %in% names(merged)) next
  vals <- sort(unique(na.omit(as.character(merged[[sg]]))))
  for (vv in vals) {
    dd <- merged[as.character(get(sg)) == vv]
    if (nrow(dd) < CONFIG$min_subgroup_n || sum(dd$OS_event == 1, na.rm = TRUE) < CONFIG$min_subgroup_events) next
    res <- fit_cox_safe(dd, survival::Surv(OS_time, OS_event) ~ risk_z, paste0("risk_within_", sg), paste0(sg, "=", vv), "subgroup")
    subgroup_results[[si]] <- res$table
    subgroup_results[[si]][, subgroup_variable := sg]
    subgroup_results[[si]][, subgroup_level := vv]
    si <- si + 1L
  }
}
subgroup_table <- data.table::rbindlist(subgroup_results, fill = TRUE)
write_tsv(subgroup_table, file.path(OUT$cox, "Step4_subgroup_risk_effect_Cox.tsv"))

###############################################################################
# 7. Calibration-by-risk-quantile and DCA-style utility audit
###############################################################################

predict_survival_from_cox <- function(fit, newdata, times) {
  sf <- tryCatch(survival::survfit(fit, newdata = newdata), error = function(e) NULL)
  if (is.null(sf)) return(matrix(NA_real_, nrow = nrow(newdata), ncol = length(times)))
  ss <- summary(sf, times = times, extend = TRUE)
  # summary.survfit with newdata returns matrix-like survival; handle both forms.
  surv <- ss$surv
  if (is.null(dim(surv))) {
    surv <- matrix(surv, ncol = length(times), byrow = FALSE)
  }
  if (nrow(surv) != nrow(newdata) && ncol(surv) == nrow(newdata)) surv <- t(surv)
  surv
}

binary_outcome_at_time <- function(time, event, t0) {
  out <- rep(NA_integer_, length(time))
  out[event == 1 & time <= t0] <- 1L
  out[time > t0] <- 0L
  out
}

empirical_net_benefit <- function(prob, y, thresholds) {
  rows <- list()
  for (i in seq_along(thresholds)) {
    pt <- thresholds[i]
    pred <- prob >= pt
    keep <- y %in% c(0L, 1L) & is.finite(prob)
    yy <- y[keep]
    pp <- pred[keep]
    n <- length(yy)
    if (n == 0) {
      rows[[i]] <- data.table::data.table(threshold = pt, net_benefit = NA_real_, n = 0, TP = NA_integer_, FP = NA_integer_)
      next
    }
    TP <- sum(pp & yy == 1L)
    FP <- sum(pp & yy == 0L)
    nb <- TP / n - FP / n * (pt / (1 - pt))
    rows[[i]] <- data.table::data.table(threshold = pt, net_benefit = nb, n = n, TP = TP, FP = FP)
  }
  data.table::rbindlist(rows)
}

utility_rows <- list()
calibration_rows <- list()
dca_rows <- list()
ui <- 1L
ci <- 1L
di <- 1L

for (coh in unique(merged$cohort_final)) {
  dd <- merged[cohort_final == coh]
  unit <- detect_survival_time_unit(dd$OS_time)
  times <- choose_utility_times(dd$OS_time, unit)
  if (length(times) == 0) next

  # Fit risk-only Cox inside each cohort for calibration/utility of the fixed score.
  # This does not change the gene model; it calibrates survival probability from fixed risk score.
  fit <- tryCatch(survival::coxph(survival::Surv(OS_time, OS_event) ~ risk_z, data = as.data.frame(dd)), error = function(e) NULL)
  if (is.null(fit)) next

  for (t0 in times) {
    y_t <- binary_outcome_at_time(dd$OS_time, dd$OS_event, t0)
    pred_surv <- predict_survival_from_cox(fit, as.data.frame(dd), t0)
    pred_risk <- 1 - as.numeric(pred_surv[, 1])
    valid <- y_t %in% c(0L, 1L) & is.finite(pred_risk)

    utility_rows[[ui]] <- data.table::data.table(
      cohort = coh,
      time_unit = unit,
      time = t0,
      n_evaluable = sum(valid),
      n_events_by_time = sum(y_t[valid] == 1L),
      predicted_risk_median = median(pred_risk[valid], na.rm = TRUE),
      predicted_risk_IQR = stats::IQR(pred_risk[valid], na.rm = TRUE)
    )
    ui <- ui + 1L

    # Calibration by predicted-risk quantile.
    dd_cal <- data.table::data.table(pred_risk = pred_risk, y = y_t)
    dd_cal <- dd_cal[y %in% c(0L, 1L) & is.finite(pred_risk)]
    if (nrow(dd_cal) >= 40 && data.table::uniqueN(dd_cal$y) >= 2) {
      dd_cal[, risk_bin := data.table::frank(pred_risk, ties.method = "average")]
      dd_cal[, risk_bin := cut(risk_bin,
                               breaks = unique(stats::quantile(risk_bin, probs = seq(0, 1, length.out = CONFIG$calibration_quantile_groups + 1), na.rm = TRUE)),
                               include.lowest = TRUE)]
      cal <- dd_cal[, .(
        n = .N,
        observed_event_rate = mean(y == 1L),
        mean_predicted_risk = mean(pred_risk, na.rm = TRUE)
      ), by = risk_bin]
      cal[, cohort := coh]
      cal[, time_unit := unit]
      cal[, time := t0]
      calibration_rows[[ci]] <- cal
      ci <- ci + 1L
    }

    # DCA-style net benefit; censored before/at t are excluded.
    nb_model <- empirical_net_benefit(pred_risk, y_t, CONFIG$dca_thresholds)
    nb_model[, strategy := "SIDISH_locked_score"]
    prevalence <- mean(y_t[valid] == 1L)
    nb_all <- data.table::data.table(
      threshold = CONFIG$dca_thresholds,
      net_benefit = prevalence - (1 - prevalence) * (CONFIG$dca_thresholds / (1 - CONFIG$dca_thresholds)),
      n = sum(valid),
      TP = sum(y_t[valid] == 1L),
      FP = sum(y_t[valid] == 0L),
      strategy = "treat_all"
    )
    nb_none <- data.table::data.table(
      threshold = CONFIG$dca_thresholds,
      net_benefit = 0,
      n = sum(valid),
      TP = 0L,
      FP = 0L,
      strategy = "treat_none"
    )
    nb <- data.table::rbindlist(list(nb_model, nb_all, nb_none), fill = TRUE)
    nb[, cohort := coh]
    nb[, time_unit := unit]
    nb[, time := t0]
    nb[, method_note := "binary_timepoint_DCA_style_net_benefit_excluding_censored_before_or_at_time"]
    dca_rows[[di]] <- nb
    di <- di + 1L
  }
}

utility_table <- data.table::rbindlist(utility_rows, fill = TRUE)
calibration_table <- data.table::rbindlist(calibration_rows, fill = TRUE)
dca_table <- data.table::rbindlist(dca_rows, fill = TRUE)

write_tsv(utility_table, file.path(OUT$utility, "Step4_timepoint_utility_audit_summary.tsv"))
write_tsv(calibration_table, file.path(OUT$utility, "Step4_calibration_by_predicted_risk_quantile.tsv"))
write_tsv(dca_table, file.path(OUT$utility, "Step4_DCA_style_net_benefit.tsv"))

###############################################################################
# 8. Figures
###############################################################################

if (isTRUE(CONFIG$make_publication_figures)) {
  msg("Generating Step 4 figures...")

  # Figure 1: Clinical covariate availability.
  fig1_dt <- data.table::melt(
    cov_avail,
    id.vars = c("dataset_stage", "cohort_final", "n", "events"),
    measure.vars = c("Age_nonmissing", "WHO_grade_nonmissing", "IDH_status_nonmissing", "Sex_nonmissing", "1p19q_status_nonmissing"),
    variable.name = "covariate",
    value.name = "nonmissing"
  )
  fig1_dt[, fraction := nonmissing / n]
  fig1_dt[, covariate := gsub("_nonmissing", "", covariate)]
  p1 <- ggplot2::ggplot(fig1_dt, ggplot2::aes(x = cohort_final, y = fraction, fill = covariate)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.66, color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::scale_fill_manual(values = c(
      Age = pal("navy"), WHO_grade = pal("orange"), IDH_status = pal("green"), Sex = pal("blue"), `1p19q_status` = pal("purple")
    )) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = "Step 4 | Clinical covariate availability",
                  subtitle = "Availability determines which adjusted Cox models are interpretable",
                  x = NULL, y = "Non-missing fraction", fill = "Covariate") +
    pub_theme()
  save_pub_plot(p1, "Figure_1_clinical_covariate_availability", width = 5.4, height = 3.4)

  # Figure 2: Risk HR in multivariable Cox models.
  fig2_dt <- risk_term_table[term == "risk_z" & status == "ok" & is.finite(HR)]
  fig2_dt <- fig2_dt[grepl("risk_only_continuous|available_clinical_covariates|risk_plus_age|risk_plus_grade|risk_plus_IDH", model)]
  fig2_dt[, label := paste(cohort, model, sep = " | ")]
  fig2_dt[, label := make.unique(as.character(label), sep = " #")]
  fig2_dt[, label := factor(label, levels = make_unique_factor_levels(label, reverse = TRUE))]
  write_tsv(fig2_dt, file.path(OUT$fig, "Figure_2_clinical_independence_risk_HR_forest_plot_data.tsv"))
  p2 <- ggplot2::ggplot(fig2_dt, ggplot2::aes(x = HR, y = label, color = dataset_stage)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower95, xmax = upper95), height = 0.18, linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_point(size = CONFIG$plot_point_size) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c(development = pal("navy"), external = pal("red"), development_plus_external = pal("green"))) +
    ggplot2::labs(title = "Step 4 | Risk-score Cox effect before and after clinical adjustment",
                  subtitle = "HR per 1 SD locked SIDISH risk score; raw estimates only",
                  x = "Hazard ratio, log scale", y = NULL, color = "Dataset") +
    pub_theme()
  save_pub_plot(p2, "Figure_2_clinical_independence_risk_HR_forest", width = 7.0, height = max(4.0, 0.18 * nrow(fig2_dt) + 1.8))

  # Figure 3: C-index by model.
  fig3_dt <- cox_table[status == "ok" & is.finite(concordance)]
  fig3_dt <- fig3_dt[grepl("risk_only_continuous|available_clinical_covariates|strata_cohort", model)]
  fig3_dt[, label := paste(cohort, model, sep = " | ")]
  fig3_dt[, label := make.unique(as.character(label), sep = " #")]
  fig3_dt[, label := factor(label, levels = make_unique_factor_levels(label, reverse = TRUE))]
  write_tsv(fig3_dt, file.path(OUT$fig, "Figure_3_objective_Cindex_audit_plot_data.tsv"))
  p3 <- ggplot2::ggplot(fig3_dt, ggplot2::aes(x = concordance, y = label, fill = dataset_stage)) +
    ggplot2::geom_col(width = 0.68, color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::scale_fill_manual(values = c(development = pal("navy"), external = pal("red"), development_plus_external = pal("green"))) +
    ggplot2::coord_cartesian(xlim = c(0.45, max(0.85, max(fig3_dt$concordance, na.rm = TRUE) + 0.03))) +
    ggplot2::labs(title = "Step 4 | Objective C-index audit",
                  subtitle = "Raw concordance values; no qualitative labels",
                  x = "C-index", y = NULL, fill = "Dataset") +
    pub_theme()
  save_pub_plot(p3, "Figure_3_objective_Cindex_audit", width = 7.0, height = max(4.0, 0.18 * nrow(fig3_dt) + 1.8))

  # Figure 4: subgroup HR.
  fig4_dt <- subgroup_table[term == "risk_z" & status == "ok" & is.finite(HR)]
  if (nrow(fig4_dt) > 0) {
    fig4_dt[, label := paste(cohort, subgroup_variable, subgroup_level, sep = " | ")]
    fig4_dt[, label := make.unique(as.character(label), sep = " #")]
    fig4_dt[, label := factor(label, levels = make_unique_factor_levels(label, reverse = TRUE))]
    write_tsv(fig4_dt, file.path(OUT$fig, "Figure_4_subgroup_risk_effect_forest_plot_data.tsv"))
    p4 <- ggplot2::ggplot(fig4_dt, ggplot2::aes(x = HR, y = label, color = subgroup_variable)) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = CONFIG$plot_linewidth) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower95, xmax = upper95), height = 0.18, linewidth = CONFIG$plot_linewidth) +
      ggplot2::geom_point(size = CONFIG$plot_point_size) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "Step 4 | Subgroup risk-effect audit",
                    subtitle = "HR per 1 SD locked SIDISH risk score within clinical strata",
                    x = "Hazard ratio, log scale", y = NULL, color = "Subgroup") +
      pub_theme()
    save_pub_plot(p4, "Figure_4_subgroup_risk_effect_forest", width = 5.8, height = max(3.8, 0.22 * nrow(fig4_dt) + 1.4))
  }

  # Figure 5: calibration by risk quantile.
  if (nrow(calibration_table) > 0) {
    p5 <- ggplot2::ggplot(calibration_table, ggplot2::aes(x = mean_predicted_risk, y = observed_event_rate, color = cohort)) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = CONFIG$plot_linewidth, color = "black") +
      ggplot2::geom_point(size = CONFIG$plot_point_size) +
      ggplot2::geom_line(linewidth = CONFIG$plot_linewidth) +
      ggplot2::facet_wrap(~ paste0(time, " ", time_unit), scales = "free") +
      ggplot2::labs(title = "Step 4 | Calibration by predicted-risk quantile",
                    subtitle = "Binary timepoint calibration; censored before/at time excluded",
                    x = "Mean predicted risk", y = "Observed event rate", color = "Cohort") +
      pub_theme()
    save_pub_plot(p5, "Figure_5_calibration_by_risk_quantile", width = 6.2, height = 4.2)
  }

  # Figure 6: DCA-style net benefit.
  if (nrow(dca_table) > 0) {
    fig6_dt <- dca_table[cohort %in% c("GSE16011", "TCGA", "CGGA325", "CGGA693")]
    write_tsv(fig6_dt, file.path(OUT$fig, "Figure_6_DCA_style_net_benefit_plot_data.tsv"))
    if (nrow(fig6_dt) > 0) {
    p6 <- ggplot2::ggplot(fig6_dt, ggplot2::aes(x = threshold, y = net_benefit, color = strategy, linetype = strategy)) +
      ggplot2::geom_line(linewidth = CONFIG$plot_linewidth) +
      ggplot2::facet_grid(cohort ~ paste0(time, " ", time_unit), scales = "free_y") +
      ggplot2::scale_color_manual(values = c(SIDISH_locked_score = pal("red"), treat_all = pal("grey"), treat_none = "black")) +
      ggplot2::labs(title = "Step 4 | DCA-style net-benefit audit",
                    subtitle = "Binary timepoint approximation; not a replacement for full survival DCA",
                    x = "Threshold probability", y = "Net benefit", color = "Strategy", linetype = "Strategy") +
      pub_theme()
    save_pub_plot(p6, "Figure_6_DCA_style_net_benefit", width = 7.2, height = 6.0)
    }
  }
}

###############################################################################
# 9. Final summary and README
###############################################################################

main_risk_independence <- risk_term_table[
  term == "risk_z" &
    grepl("available_clinical_covariates|risk_plus_age|risk_plus_grade|risk_plus_IDH|risk_only_continuous", model) &
    status == "ok",
  .(
    dataset_stage, cohort, model, n, events, HR, lower95, upper95, p, concordance
  )
]
write_tsv(main_risk_independence, file.path(OUT$final, "SIDISH_Step4_main_clinical_independence_summary.tsv"))

final_summary <- data.table::data.table(
  metric = c(
    "n_total_samples_in_clinical_audit",
    "n_development_samples",
    "n_external_samples",
    "n_locked_model_genes",
    "n_cox_rows_total",
    "n_cox_rows_ok",
    "n_subgroup_rows_ok",
    "n_calibration_rows",
    "n_dca_rows",
    "reporting_policy"
  ),
  value = c(
    nrow(merged),
    nrow(merged[dataset_stage == "development"]),
    nrow(merged[dataset_stage == "external"]),
    nrow(step2_coef),
    nrow(cox_table),
    nrow(cox_table[status == "ok"]),
    nrow(subgroup_table[status == "ok"]),
    nrow(calibration_table),
    nrow(dca_table),
    CONFIG$objective_reporting_note
  )
)
write_tsv(final_summary, file.path(OUT$final, "SIDISH_Step4_clinical_audit_final_summary.tsv"))

readme <- c(
  "SIDISH Step 4 Clinical independence and utility audit v1.1",
  "",
  "This step does not change the Step 2 v2.6 locked gene model.",
  "It audits the fixed SIDISH risk score against available clinical variables.
Version v1.1 fixes duplicated factor-level labels in forest/C-index/subgroup figures and exports figure-level plot-data TSV files for easier debugging.",
  "",
  "Key outputs:",
  file.path(OUT$cox, "Step4_clinical_independence_Cox_all_models.tsv"),
  file.path(OUT$cox, "Step4_clinical_independence_risk_term_summary.tsv"),
  file.path(OUT$cox, "Step4_subgroup_risk_effect_Cox.tsv"),
  file.path(OUT$utility, "Step4_calibration_by_predicted_risk_quantile.tsv"),
  file.path(OUT$utility, "Step4_DCA_style_net_benefit.tsv"),
  file.path(OUT$final, "SIDISH_Step4_main_clinical_independence_summary.tsv"),
  "",
  "Reporting policy:",
  "Report raw HR, C-index, calibration and DCA-style values objectively.",
  "Do not overstate C-index. Use DCA-style net benefit as an audit, not as a definitive clinical decision analysis unless a full survival-DCA method is later added."
)
write_text(readme, file.path(OUT$root, "README_SIDISH_Step4_Clinical_Independence_Utility_Audit_v1_1.txt"))
write_text(capture.output(sessionInfo()), file.path(OUT$audit, "sessionInfo.txt"))

msg("SIDISH Step 4 clinical independence and utility audit complete.")
msg("Output directory: ", OUT$root)
