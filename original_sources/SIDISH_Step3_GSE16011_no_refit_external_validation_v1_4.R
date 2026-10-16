#!/usr/bin/env Rscript
###############################################################################
# SIDISH Step 3: GSE16011 no-refit external validation v1.4
#
# Purpose
# -------
# Externally validate the Step 2 v2.6 direction-stable locked model in GSE16011.
#
# Strict no-refit rule:
#   - Do NOT select new genes.
#   - Do NOT refit coefficients.
#   - Do NOT re-optimize cutoff.
#   - Do NOT use external outcome to tune anything.
#
# Only locked Step 2 v2.6 artifacts are allowed:
#   1) SIDISH_Step2_locked_model_coefficients.tsv
#   2) SIDISH_Step2_locked_model_summary.tsv
#
# External validation workflow:
#   locked 9 genes + locked coefficients + locked development cutoff
#   → map GSE16011 expression and clinical OS
#   → cohort-wise z-score external expression for locked genes
#   → compute risk score using locked coefficients
#   → classify High/Low using locked development cutoff
#   → report C-index, Cox HR, KM/log-rank, optional time-dependent AUC
#
# Required packages:
#   install.packages(c("data.table", "survival", "ggplot2"))
#
# Optional:
#   install.packages("timeROC")
#
# Run:
#   source("G:/NEW/The2/SIDISH_Step3_GSE16011_no_refit_external_validation_v1_4.R")
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
  step2_root = "G:/NEW/The2/SIDISH_Step2_ML175_DirectionStable_Model_Development_v2_6_strict_all3",
  locked_coef_file = "G:/NEW/The2/SIDISH_Step2_ML175_DirectionStable_Model_Development_v2_6_strict_all3/99_locked_model/SIDISH_Step2_locked_model_coefficients.tsv",
  locked_summary_file = "G:/NEW/The2/SIDISH_Step2_ML175_DirectionStable_Model_Development_v2_6_strict_all3/99_locked_model/SIDISH_Step2_locked_model_summary.tsv",

  external_cohort = "GSE16011",
  output_dir = "G:/NEW/The2/SIDISH_Step3_GSE16011_no_refit_external_validation_v1_4",

  # Preferred paths. If missing, script recursively searches external_search_roots.
  external_expr_file = "G:/NEW/The2/External_Cohorts/GSE16011/GSE16011_expression_matrix.tsv",
  external_clin_file = "G:/NEW/The2/External_Cohorts/GSE16011/GSE_glioma_clinical_molecular.tsv",

  external_search_roots = c(
    "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/external",
    "G:/NEW/The2/External_Cohorts",
    "G:/NEW/转化/External_Cohorts",
    "G:/NEW/转化/GSE16011",
    "G:/NEW/The2/GSE16011",
    "G:/NEW/The2"
  ),

  # No-refit validation.
  require_all_locked_genes = TRUE,
  external_expression_scale = "within_external_cohort_zscore",
  primary_cutoff_source = "locked_development_cutoff",

  # v1.1 external file validation.
  # A file name is not enough. Expression/clinical files must pass minimum
  # dimensionality checks before they can be selected.
  min_external_expr_rows = 50,
  min_external_expr_cols = 20,
  min_external_clin_rows = 20,
  min_external_clin_cols = 2,

  # Optional time-dependent AUC.
  # v1.3: never write an empty timeROC file. Time points are chosen according
  # to the detected OS_time unit. GSE16011 commonly uses months, so using
  # 365/1095/1825 on a month-scale OS_time would be invalid.
  compute_timeROC_if_available = TRUE,
  timeROC_time_unit = "auto",  # "auto", "days", or "months"
  timeROC_times_days = c(365, 1095, 1825),
  timeROC_times_months = c(12, 36, 60),

  # Plotting, consistent with Step 1/Step 2.
  make_publication_figures = TRUE,
  plot_base_family = "Arial",
  plot_base_size = 8,
  plot_linewidth = 0.5,
  plot_point_size = 2.0,
  km_time_unit_label = "Time",
  km_palette_high = "#E64B35FF",
  km_palette_low = "#4DBBD5FF",
  p_value_floor_for_display = 2.2e-16,

  # Safety.
  forbidden_development_patterns = c("TCGA", "CGGA325", "CGGA_325", "CGGA693", "CGGA_693", "Bulk_Cohorts"),
  required_external_pattern = "GSE16011"
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
  data = file.path(CONFIG$output_dir, "01_external_validation_data"),
  perf = file.path(CONFIG$output_dir, "02_external_validation_performance"),
  fig = file.path(CONFIG$output_dir, "03_publication_figures"),
  final = file.path(CONFIG$output_dir, "99_final_external_validation")
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

# v1.4: timeROC may call Surv() from the search path in some Windows/R setups.
# Explicitly attach survival and expose Surv in the global environment.
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

detect_delim <- function(file) {
  con <- file(file, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  ln <- readLines(con, n = 1, warn = FALSE)
  if (length(ln) == 0) return("\t")
  n_tab <- lengths(regmatches(ln, gregexpr("\t", ln, fixed = TRUE)))
  n_comma <- lengths(regmatches(ln, gregexpr(",", ln, fixed = TRUE)))
  ifelse(n_tab > n_comma, "\t", ",")
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
    data.table::fread(file, sep = "auto", data.table = TRUE, showProgress = FALSE, quote = ""),
    error = function(e1) {
      ext <- tolower(tools::file_ext(file))
      sep <- if (ext %in% c("tsv", "txt")) "\t" else detect_delim(file)
      data.table::fread(file, sep = sep, data.table = TRUE, showProgress = FALSE, quote = "")
    }
  )
  data.table::setnames(dt, clean_names(names(dt)))
  dt
}

read_table_robust_dimensions <- function(file, nrows = Inf) {
  file <- norm_path(file)
  if (!file.exists(file)) stop("File not found: ", file)

  seps <- c("auto", "\t", ",", ";", " ")
  tries <- list()
  k <- 1L
  for (sp in seps) {
    dt <- tryCatch(
      data.table::fread(file, sep = sp, nrows = nrows, data.table = TRUE,
                        showProgress = FALSE, quote = "", fill = TRUE),
      error = function(e) NULL
    )
    if (!is.null(dt)) {
      data.table::setnames(dt, clean_names(names(dt)))
      tries[[k]] <- list(sep = sp, dt = dt, nrow = nrow(dt), ncol = ncol(dt))
      k <- k + 1L
    }
  }
  if (length(tries) == 0) stop("Could not read table with common delimiters: ", file)

  score <- vapply(tries, function(z) {
    # Prefer real matrices: many columns and rows. Penalize one-column parses.
    z$nrow * max(0, z$ncol - 1)
  }, numeric(1))
  best <- tries[[which.max(score)]]
  best$dt
}

normalize_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("^['\"]|['\"]$", "", x)
  x <- toupper(x)
  x[x == "" | x == "NA" | is.na(x)] <- NA_character_
  x
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

zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 1, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv <= 0] <- 1
  out <- sweep(mat, 1, mu, "-")
  out <- sweep(out, 1, sdv, "/")
  out[!is.finite(out)] <- 0
  out
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

km_summary_table <- function(time, event, group) {
  dd <- data.frame(time = as.numeric(time), event = as.numeric(event), group = as.character(group))
  dd <- dd[is.finite(dd$time) & dd$time > 0 & dd$event %in% c(0,1) & !is.na(dd$group), , drop = FALSE]
  if (nrow(dd) == 0 || length(unique(dd$group)) < 2) return(NULL)
  sf <- tryCatch(survival::survfit(survival::Surv(time, event) ~ group, data = dd), error = function(e) NULL)
  if (is.null(sf)) return(NULL)
  ss <- summary(sf)
  if (length(ss$time) == 0) return(NULL)
  data.table::data.table(
    time = ss$time,
    surv = ss$surv,
    strata = gsub("^group=", "", ss$strata),
    n_risk = ss$n.risk,
    n_event = ss$n.event
  )
}

km_logrank_p <- function(time, event, group) {
  dd <- data.frame(time = as.numeric(time), event = as.numeric(event), group = as.character(group))
  dd <- dd[is.finite(dd$time) & dd$time > 0 & dd$event %in% c(0,1) & !is.na(dd$group), , drop = FALSE]
  if (nrow(dd) == 0 || length(unique(dd$group)) < 2) return(NA_real_)
  sdif <- tryCatch(survival::survdiff(survival::Surv(time, event) ~ group, data = dd), error = function(e) NULL)
  if (is.null(sdif)) return(NA_real_)
  1 - stats::pchisq(sdif$chisq, df = max(1, length(sdif$n) - 1))
}

format_p_for_plot <- function(p, floor_value = CONFIG$p_value_floor_for_display) {
  p <- suppressWarnings(as.numeric(p))
  if (!is.finite(p)) return("p = NA")
  if (p < floor_value) return(paste0("p < ", format(floor_value, scientific = TRUE, digits = 2)))
  if (p < 0.001) return(paste0("p = ", format(p, scientific = TRUE, digits = 2)))
  paste0("p = ", sprintf("%.3f", p))
}

eval_risk <- function(time, event, risk) {
  dat <- data.frame(time = as.numeric(time), event = as.numeric(event), risk = as.numeric(risk))
  dat <- dat[is.finite(dat$time) & dat$time > 0 & dat$event %in% c(0, 1) & is.finite(dat$risk), , drop = FALSE]
  out <- c(n = nrow(dat), events = sum(dat$event == 1), cindex = NA_real_, beta = NA_real_, HR = NA_real_, p = NA_real_)
  if (nrow(dat) < 20 || sum(dat$event == 1) < 5 || stats::sd(dat$risk, na.rm = TRUE) <= 0) return(out)
  dat$risk_z <- as.numeric(scale(dat$risk))
  fit <- tryCatch(survival::coxph(survival::Surv(time, event) ~ risk_z, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(out)
  sm <- summary(fit)
  out["beta"] <- as.numeric(stats::coef(fit)["risk_z"])
  out["HR"] <- exp(out["beta"])
  out["p"] <- sm$coefficients["risk_z", "Pr(>|z|)"]
  out["cindex"] <- as.numeric(sm$concordance[1])
  out
}

detect_survival_time_unit <- function(time, configured = CONFIG$timeROC_time_unit) {
  configured <- tolower(as.character(configured))
  if (configured %in% c("days", "months")) return(configured)

  tt <- suppressWarnings(as.numeric(time))
  tt <- tt[is.finite(tt) & tt > 0]
  if (length(tt) == 0) return("unknown")

  # Heuristic:
  # - If values look like clinical months, max is usually < 300 and median often < 100.
  # - If values are days, max often exceeds 1000.
  if (stats::quantile(tt, 0.95, na.rm = TRUE) <= 240 && max(tt, na.rm = TRUE) <= 400) {
    return("months")
  }
  return("days")
}

choose_timeROC_times <- function(time, unit) {
  tt <- suppressWarnings(as.numeric(time))
  tt <- tt[is.finite(tt) & tt > 0]
  if (length(tt) == 0) return(numeric(0))

  candidate_times <- if (identical(unit, "months")) CONFIG$timeROC_times_months else CONFIG$timeROC_times_days

  # timeROC requires meaningful time points within the observed follow-up range.
  max_t <- max(tt, na.rm = TRUE)
  min_t <- min(tt, na.rm = TRUE)
  valid_times <- candidate_times[candidate_times > min_t & candidate_times < max_t]

  # If predefined times are outside range, use robust quantile-based audit times.
  if (length(valid_times) == 0) {
    q_times <- as.numeric(stats::quantile(tt, probs = c(0.25, 0.50, 0.75), na.rm = TRUE, names = FALSE))
    q_times <- unique(round(q_times, 6))
    valid_times <- q_times[q_times > min_t & q_times < max_t]
  }

  unique(valid_times)
}

make_timeROC_status_table <- function(times, unit, status, reason, auc = NULL) {
  if (length(times) == 0) {
    return(data.table::data.table(
      external_cohort = CONFIG$external_cohort,
      time_unit = unit,
      time = NA_real_,
      AUC = NA_real_,
      status = status,
      reason = reason
    ))
  }

  if (is.null(auc)) auc <- rep(NA_real_, length(times))
  data.table::data.table(
    external_cohort = CONFIG$external_cohort,
    time_unit = unit,
    time = times,
    AUC = as.numeric(auc),
    status = status,
    reason = reason
  )
}

empirical_auc_cases_controls <- function(marker, case) {
  marker <- as.numeric(marker)
  case <- as.integer(case)
  keep <- is.finite(marker) & case %in% c(0L, 1L)
  marker <- marker[keep]
  case <- case[keep]
  n_case <- sum(case == 1L)
  n_ctrl <- sum(case == 0L)
  if (n_case < 1 || n_ctrl < 1) return(NA_real_)

  # AUC = P(marker_case > marker_control) + 0.5 * P(tie)
  ranks <- rank(marker, ties.method = "average")
  auc <- (sum(ranks[case == 1L]) - n_case * (n_case + 1) / 2) / (n_case * n_ctrl)
  as.numeric(auc)
}

fallback_dynamic_auc_table <- function(time, event, marker, times, unit, fail_reason) {
  time <- as.numeric(time)
  event <- as.integer(event)
  marker <- as.numeric(marker)

  if (length(times) == 0) {
    return(data.table::data.table(
      external_cohort = CONFIG$external_cohort,
      time_unit = unit,
      time = NA_real_,
      AUC = NA_real_,
      status = "failed_no_fallback",
      reason = paste0(fail_reason, "; no_valid_time_points_for_fallback"),
      method = "none",
      n_case = NA_integer_,
      n_control = NA_integer_,
      n_excluded_censored_before_or_at_time = NA_integer_
    ))
  }

  rows <- vector("list", length(times))
  for (i in seq_along(times)) {
    tt <- times[i]
    case <- rep(NA_integer_, length(time))
    case[event == 1L & time <= tt] <- 1L
    case[time > tt] <- 0L
    # Censored before/at t are unknown for status at t and are excluded.
    excluded <- sum(event == 0L & time <= tt, na.rm = TRUE)
    auc <- empirical_auc_cases_controls(marker, case)

    rows[[i]] <- data.table::data.table(
      external_cohort = CONFIG$external_cohort,
      time_unit = unit,
      time = tt,
      AUC = auc,
      status = ifelse(is.finite(auc), "fallback_ok", "fallback_failed"),
      reason = ifelse(
        is.finite(auc),
        paste0(fail_reason, "; fallback_empirical_dynamic_AUC_excluding_censored_before_or_at_time"),
        paste0(fail_reason, "; fallback_failed_insufficient_cases_or_controls")
      ),
      method = "empirical_dynamic_AUC_excluding_censored_before_or_at_time_not_IPCW",
      n_case = sum(case == 1L, na.rm = TRUE),
      n_control = sum(case == 0L, na.rm = TRUE),
      n_excluded_censored_before_or_at_time = excluded
    )
  }
  data.table::rbindlist(rows, fill = TRUE)
}

###############################################################################
# 2. Safety and locked model input
###############################################################################

msg("SIDISH Step 3 GSE16011 no-refit external validation started.")

if (!file.exists(CONFIG$locked_coef_file)) stop("Locked coefficient file not found: ", CONFIG$locked_coef_file)
if (!file.exists(CONFIG$locked_summary_file)) stop("Locked summary file not found: ", CONFIG$locked_summary_file)

locked_coef <- read_table_auto(CONFIG$locked_coef_file)
locked_summary <- read_table_auto(CONFIG$locked_summary_file)

if (!all(c("gene", "locked_coefficient") %in% names(locked_coef))) {
  stop("Locked coefficient file must contain gene and locked_coefficient columns.")
}
locked_coef[, gene := normalize_gene(gene)]
locked_coef[, locked_coefficient := as_numeric_safe(locked_coefficient)]
locked_coef <- locked_coef[!is.na(gene) & gene != "" & is.finite(locked_coefficient)]
locked_coef <- locked_coef[order(model_rank)]

if (nrow(locked_coef) < 3) stop("Too few locked genes found: ", nrow(locked_coef))

cutoff_row <- locked_summary[metric == "locked_development_cutoff"]
if (nrow(cutoff_row) == 0) stop("locked_development_cutoff not found in locked summary.")
locked_cutoff <- as_numeric_safe(cutoff_row$value[1])
if (!is.finite(locked_cutoff)) stop("Invalid locked development cutoff: ", cutoff_row$value[1])

write_tsv(locked_coef, file.path(OUT$audit, "Step3_locked_coefficients_used.tsv"))
write_tsv(data.table::data.table(
  metric = c("external_cohort", "n_locked_genes", "locked_development_cutoff", "no_refit_rule"),
  value = c(CONFIG$external_cohort, nrow(locked_coef), locked_cutoff,
            "genes_coefficients_and_cutoff_imported_from_Step2_v2_6_without_modification")
), file.path(OUT$audit, "Step3_no_refit_locking_audit.tsv"))

###############################################################################
# 3. External file discovery
###############################################################################

is_safe_external_file <- function(file) {
  file <- norm_path(file)
  if (file == "" || !file.exists(file)) return(FALSE)
  has_external <- grepl(CONFIG$required_external_pattern, file, ignore.case = TRUE)
  has_dev <- grepl(paste(CONFIG$forbidden_development_patterns, collapse = "|"), file, ignore.case = TRUE)
  has_external && !has_dev
}

list_external_files <- function(roots) {
  roots <- norm_path(roots)
  roots <- roots[roots != "" & dir.exists(roots)]
  if (length(roots) == 0) return(data.table::data.table())
  out <- data.table::rbindlist(lapply(roots, function(r) {
    files <- list.files(r, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
    files <- norm_path(files[file.exists(files)])
    if (length(files) == 0) return(NULL)
    info <- file.info(files)
    data.table::data.table(
      file = files,
      basename = basename(files),
      ext = tolower(tools::file_ext(files)),
      size_bytes = as.numeric(info$size),
      modified_time = as.character(info$mtime)
    )
  }), fill = TRUE)
  if (is.null(out)) data.table::data.table() else out
}

quick_table_probe <- function(file, nrows = 100) {
  file <- norm_path(file)
  if (!file.exists(file)) {
    return(list(ok = FALSE, nrow = 0L, ncol = 0L, reason = "file_not_found"))
  }
  dt <- tryCatch(read_table_robust_dimensions(file, nrows = nrows), error = function(e) NULL)
  if (is.null(dt)) {
    return(list(ok = FALSE, nrow = 0L, ncol = 0L, reason = "fread_failed_all_common_delimiters"))
  }
  list(ok = TRUE, nrow = nrow(dt), ncol = ncol(dt), reason = "readable")
}

validate_candidate_file <- function(file, type = c("expr", "clin")) {
  type <- match.arg(type)
  file <- norm_path(file)
  if (!is_safe_external_file(file)) {
    return(data.table::data.table(
      file = file, type = type, valid = FALSE, probe_nrow = 0L, probe_ncol = 0L,
      reason = "not_safe_external_or_missing"
    ))
  }

  pr <- quick_table_probe(file, nrows = 200)
  if (!isTRUE(pr$ok)) {
    return(data.table::data.table(
      file = file, type = type, valid = FALSE, probe_nrow = pr$nrow, probe_ncol = pr$ncol,
      reason = pr$reason
    ))
  }

  if (type == "expr") {
    valid <- pr$nrow >= CONFIG$min_external_expr_rows && pr$ncol >= CONFIG$min_external_expr_cols
    reason <- if (valid) "valid_expression_candidate" else paste0(
      "expression_table_too_small_probe_nrow_", pr$nrow, "_ncol_", pr$ncol,
      "_required_nrow_", CONFIG$min_external_expr_rows, "_ncol_", CONFIG$min_external_expr_cols
    )
  } else {
    valid <- pr$nrow >= CONFIG$min_external_clin_rows && pr$ncol >= CONFIG$min_external_clin_cols
    reason <- if (valid) "valid_clinical_candidate" else paste0(
      "clinical_table_too_small_probe_nrow_", pr$nrow, "_ncol_", pr$ncol,
      "_required_nrow_", CONFIG$min_external_clin_rows, "_ncol_", CONFIG$min_external_clin_cols
    )
  }

  data.table::data.table(
    file = file, type = type, valid = valid, probe_nrow = pr$nrow, probe_ncol = pr$ncol,
    reason = reason
  )
}

external_inventory <- list_external_files(CONFIG$external_search_roots)
if (nrow(external_inventory) > 0) {
  external_inventory[, external_pattern := grepl(CONFIG$required_external_pattern, file, ignore.case = TRUE)]
  external_inventory[, development_pattern := grepl(paste(CONFIG$forbidden_development_patterns, collapse = "|"), file, ignore.case = TRUE)]
  external_inventory <- external_inventory[external_pattern == TRUE & development_pattern == FALSE]
}
write_tsv(external_inventory, file.path(OUT$audit, "Step3_external_candidate_file_inventory.tsv"))

score_file_candidate <- function(files_dt, include_patterns, exclude_patterns = character(0),
                                 allowed_ext = c("tsv", "txt", "csv", "gz")) {
  if (nrow(files_dt) == 0) return(data.table::data.table())
  dt <- data.table::copy(files_dt[ext %in% allowed_ext])
  if (nrow(dt) == 0) return(data.table::data.table())
  dt[, score := 0]
  for (pat in include_patterns) {
    dt[grepl(pat, basename, ignore.case = TRUE), score := score + 3]
    dt[grepl(pat, file, ignore.case = TRUE), score := score + 1]
  }
  for (pat in exclude_patterns) {
    dt[grepl(pat, basename, ignore.case = TRUE) | grepl(pat, file, ignore.case = TRUE), score := score - 6]
  }
  dt <- dt[score > 0]
  if (nrow(dt) == 0) return(data.table::data.table())
  dt[order(-score, -size_bytes)]
}

discover_external_file <- function(type = c("expr", "clin")) {
  type <- match.arg(type)
  exact <- if (type == "expr") CONFIG$external_expr_file else CONFIG$external_clin_file
  exact <- norm_path(exact)

  validation_rows <- list()
  rr <- 1L

  # Try exact configured file first, but only accept if it passes matrix checks.
  if (exact != "") {
    v <- validate_candidate_file(exact, type)
    v[, source := "configured_exact_path"]
    validation_rows[[rr]] <- v
    rr <- rr + 1L
    if (isTRUE(v$valid[1])) {
      write_tsv(data.table::rbindlist(validation_rows, fill = TRUE),
                file.path(OUT$audit, paste0("Step3_", type, "_candidate_validation_report.tsv")))
      return(exact)
    }
  }

  if (type == "expr") {
    include <- c("GSE16011", "expr|expression|matrix|mRNA|RNAseq|microarray|normalized|TPM|FPKM|RSEM|count")
    exclude <- c("clinical|clin|survival|phenotype|annotation|meta|platform")
  } else {
    include <- c("GSE16011", "GSE_glioma", "clinical|clin|survival|phenotype|molecular|meta|annotation")
    exclude <- c("expr|expression|matrix|mRNA|RNAseq|TPM|FPKM|RSEM|count")
  }

  cand_dt <- score_file_candidate(external_inventory, include, exclude)
  if (nrow(cand_dt) > 0) {
    for (ii in seq_len(nrow(cand_dt))) {
      f <- cand_dt$file[ii]
      v <- validate_candidate_file(f, type)
      v[, source := "inventory_search"]
      v[, search_score := cand_dt$score[ii]]
      validation_rows[[rr]] <- v
      rr <- rr + 1L
      if (isTRUE(v$valid[1])) {
        write_tsv(data.table::rbindlist(validation_rows, fill = TRUE),
                  file.path(OUT$audit, paste0("Step3_", type, "_candidate_validation_report.tsv")))
        return(f)
      }
    }
  }

  if (length(validation_rows) > 0) {
    write_tsv(data.table::rbindlist(validation_rows, fill = TRUE),
              file.path(OUT$audit, paste0("Step3_", type, "_candidate_validation_report.tsv")))
  }
  ""
}

expr_file <- discover_external_file("expr")
clin_file <- discover_external_file("clin")

selected_external_files <- data.table::data.table(
  external_cohort = CONFIG$external_cohort,
  expr_file = expr_file,
  clin_file = clin_file,
  expr_exists = expr_file != "" & file.exists(expr_file),
  clin_exists = clin_file != "" & file.exists(clin_file)
)
write_tsv(selected_external_files, file.path(OUT$audit, "Step3_selected_external_input_files.tsv"))

if (expr_file == "" || clin_file == "" || !file.exists(expr_file) || !file.exists(clin_file)) {
  stop(
    "Valid external expression/clinical files not found. ",
    "Inspect Step3_external_candidate_file_inventory.tsv, ",
    "Step3_expr_candidate_validation_report.tsv and Step3_clin_candidate_validation_report.tsv. ",
    "The previous file may exist but be too small or not a matrix. ",
    "Set CONFIG$external_expr_file / CONFIG$external_clin_file manually to the real GSE16011 expression matrix and clinical table."
  )
}

###############################################################################
# 4. Read and standardize GSE16011 expression/clinical
###############################################################################

standardize_external_expr <- function(file) {
  dt <- read_table_robust_dimensions(file)
  if (nrow(dt) == 0 || ncol(dt) < 2) stop("Expression table too small after robust delimiter reading: ", file)

  # Case A: conventional gene x sample matrix with a gene column.
  gene_col <- first_existing_col(dt, c("^gene$", "gene_symbol", "^symbol$", "hgnc", "feature", "probe", "ID_REF", "ensembl"), regex = TRUE)

  # Case B: transposed sample x gene matrix. If locked genes appear in column names,
  # transpose the locked-gene columns into gene x sample form.
  norm_col_names <- normalize_gene(names(dt))
  locked_hits_in_columns <- intersect(locked_coef$gene, norm_col_names)
  if (is.na(gene_col) && length(locked_hits_in_columns) >= min(3, nrow(locked_coef))) {
    sample_col <- names(dt)[1]
    gene_col_idx <- which(norm_col_names %in% locked_coef$gene)
    sample_ids <- standardize_id(dt[[sample_col]])
    mat <- as.matrix(dt[, ..gene_col_idx])
    storage.mode(mat) <- "numeric"
    mat <- t(mat)
    rownames(mat) <- norm_col_names[gene_col_idx]
    colnames(mat) <- sample_ids
    mat[!is.finite(mat)] <- NA_real_

    sample_vals <- as.numeric(mat[seq_len(min(nrow(mat), 100)), seq_len(min(ncol(mat), 100)), drop = FALSE])
    if (length(sample_vals) > 0 && max(sample_vals, na.rm = TRUE) > 50) {
      mat <- log2(mat + 1)
    }
    return(mat)
  }

  if (is.na(gene_col)) {
    first_col <- names(dt)[1]
    first_num_frac <- mean(!is.na(as_numeric_safe(dt[[first_col]])))
    if (first_num_frac < 0.5) gene_col <- first_col
  }
  if (is.na(gene_col)) stop("Cannot identify gene/probe column in external expression file: ", file)

  data.table::setnames(dt, gene_col, "gene")
  dt[, gene := normalize_gene(gene)]
  dt <- dt[!is.na(gene) & gene != ""]
  sample_cols <- setdiff(names(dt), "gene")
  for (cc in sample_cols) dt[[cc]] <- as_numeric_safe(dt[[cc]])

  if (anyDuplicated(dt$gene)) {
    dt <- dt[, lapply(.SD, function(z) mean(as.numeric(z), na.rm = TRUE)), by = gene, .SDcols = sample_cols]
  }

  mat <- as.matrix(dt[, -"gene"])
  rownames(mat) <- dt$gene
  colnames(mat) <- standardize_id(colnames(mat))
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- NA_real_

  sample_vals <- as.numeric(mat[seq_len(min(nrow(mat), 100)), seq_len(min(ncol(mat), 100)), drop = FALSE])
  if (length(sample_vals) > 0 && max(sample_vals, na.rm = TRUE) > 50) {
    mat <- log2(mat + 1)
  }
  mat
}

standardize_external_clinical <- function(file) {
  dt <- read_table_robust_dimensions(file)
  if (nrow(dt) == 0) stop("External clinical table empty: ", file)

  sample_col <- first_existing_col(dt, c("^sample_id$", "^sample$", "geo_accession", "GSM", "barcode", "case", "patient", "^ID$"), regex = TRUE)
  if (is.na(sample_col)) sample_col <- names(dt)[1]

  time_col <- first_existing_col(dt, c("^OS_time$", "^os_time$", "survival.*time", "follow.*time", "days", "months", "time"), regex = TRUE)
  event_col <- first_existing_col(dt, c("^OS_event$", "^os_event$", "survival.*event", "status", "death", "event", "vital"), regex = TRUE)

  age_col <- first_existing_col(dt, c("^Age$", "^age$"), regex = TRUE)
  grade_col <- first_existing_col(dt, c("WHO_grade", "grade", "histology"), regex = TRUE)
  idh_col <- first_existing_col(dt, c("IDH_status", "IDH", "idh"), regex = TRUE)

  if (is.na(time_col) || is.na(event_col)) stop("External clinical survival columns not found: ", file)

  out <- data.table::data.table(
    sample_id_raw = as.character(dt[[sample_col]]),
    sample_id = standardize_id(dt[[sample_col]]),
    OS_time = as_numeric_safe(dt[[time_col]]),
    OS_event_raw = dt[[event_col]],
    Age = if (!is.na(age_col)) as_numeric_safe(dt[[age_col]]) else NA_real_,
    WHO_grade = if (!is.na(grade_col)) as.character(dt[[grade_col]]) else NA_character_,
    IDH_status = if (!is.na(idh_col)) as.character(dt[[idh_col]]) else NA_character_
  )

  ev <- toupper(as.character(out$OS_event_raw))
  out[, OS_event := data.table::fcase(
    ev %in% c("1", "DEAD", "DECEASED", "EVENT", "TRUE", "YES", "1:DECEASED", "DIED"), 1,
    ev %in% c("0", "ALIVE", "CENSORED", "FALSE", "NO", "LIVING", "0:LIVING"), 0,
    default = suppressWarnings(as.numeric(ev))
  )]
  out <- out[!is.na(sample_id) & sample_id != "" & is.finite(OS_time) & OS_time > 0 & OS_event %in% c(0, 1)]
  unique(out, by = "sample_id")
}

match_external_samples <- function(expr_samples, clin_dt) {
  expr_dt <- data.table::data.table(expr_sample = expr_samples)
  expr_dt[, exact := expr_sample]
  expr_dt[, clean := standardize_id(expr_sample)]

  clin_map <- data.table::copy(clin_dt)
  clin_map[, exact := sample_id]
  clin_map[, clean := standardize_id(sample_id)]

  modes <- c("exact", "clean")
  reports <- list()
  best <- NULL
  best_n <- -1
  best_mode <- NA_character_

  for (mode in modes) {
    e <- expr_dt[, .(expr_sample, key = get(mode))]
    c <- clin_map[, .(sample_id, key = get(mode))]
    e <- e[!is.na(key) & key != ""]
    c <- c[!is.na(key) & key != ""]
    c <- unique(c, by = "key")
    m <- merge(e, c, by = "key", all = FALSE)
    m <- unique(m, by = "expr_sample")
    reports[[mode]] <- data.table::data.table(match_mode = mode, n_matched = nrow(m))
    if (nrow(m) > best_n) {
      best <- m
      best_n <- nrow(m)
      best_mode <- mode
    }
  }

  list(matches = best, mode = best_mode, report = data.table::rbindlist(reports))
}

msg("Reading external expression and clinical data...")
expr <- standardize_external_expr(expr_file)
clin <- standardize_external_clinical(clin_file)
sample_match <- match_external_samples(colnames(expr), clin)
write_tsv(sample_match$report, file.path(OUT$audit, "Step3_external_sample_match_mode_report.tsv"))

if (nrow(sample_match$matches) < 20) {
  stop("Too few matched external samples: ", nrow(sample_match$matches))
}

expr <- expr[, sample_match$matches$expr_sample, drop = FALSE]
clin <- clin[match(sample_match$matches$sample_id, sample_id)]

###############################################################################
# 5. Locked gene mapping and no-refit score calculation
###############################################################################

locked_genes <- locked_coef$gene
present_locked <- intersect(locked_genes, rownames(expr))
missing_locked <- setdiff(locked_genes, rownames(expr))

gene_mapping_audit <- data.table::data.table(
  gene = locked_genes,
  present_in_external_expr = locked_genes %in% rownames(expr),
  locked_coefficient = locked_coef$locked_coefficient[match(locked_genes, locked_coef$gene)]
)
write_tsv(gene_mapping_audit, file.path(OUT$data, "Step3_locked_gene_external_mapping_audit.tsv"))

if (length(missing_locked) > 0 && isTRUE(CONFIG$require_all_locked_genes)) {
  stop("Missing locked genes in external expression: ", paste(missing_locked, collapse = ", "),
       ". No-refit validation requires all locked genes.")
}

expr_locked <- expr[locked_genes, , drop = FALSE]
expr_locked_z <- zscore_rows(expr_locked)

external_matrix <- data.table::as.data.table(t(expr_locked_z), keep.rownames = "expr_sample")
data.table::setnames(external_matrix, "expr_sample", "sample_id_expr")
external_matrix[, sample_id := clin$sample_id]
external_matrix[, OS_time := clin$OS_time]
external_matrix[, OS_event := clin$OS_event]
write_tsv(external_matrix, file.path(OUT$data, "Step3_external_locked_gene_matrix_zscored.tsv"))

coef_vec <- locked_coef$locked_coefficient
names(coef_vec) <- locked_coef$gene
score_mat <- as.matrix(external_matrix[, locked_genes, with = FALSE])
risk_score <- as.numeric(score_mat[, names(coef_vec), drop = FALSE] %*% coef_vec)

risk_dt <- data.table::data.table(
  sample_id = external_matrix$sample_id,
  external_cohort = CONFIG$external_cohort,
  OS_time = external_matrix$OS_time,
  OS_event = external_matrix$OS_event,
  locked_risk_score = risk_score,
  locked_development_cutoff = locked_cutoff,
  locked_risk_group = ifelse(risk_score >= locked_cutoff, "High", "Low")
)

write_tsv(risk_dt, file.path(OUT$final, "SIDISH_Step3_GSE16011_no_refit_risk_scores.tsv"))

risk_group_count <- risk_dt[, .(
  n = .N,
  events = sum(OS_event == 1, na.rm = TRUE),
  median_OS_time = median(OS_time, na.rm = TRUE),
  median_risk_score = median(locked_risk_score, na.rm = TRUE)
), by = locked_risk_group]
write_tsv(risk_group_count, file.path(OUT$perf, "SIDISH_Step3_GSE16011_locked_risk_group_counts.tsv"))

if (data.table::uniqueN(risk_dt$locked_risk_group) < 2) {
  warning("Locked development cutoff produced only one external risk group. No cutoff optimization will be performed.")
}

###############################################################################
# 6. External performance
###############################################################################

cont_ev <- eval_risk(risk_dt$OS_time, risk_dt$OS_event, risk_dt$locked_risk_score)

high_low_perf <- data.table::data.table()
if (data.table::uniqueN(risk_dt$locked_risk_group) >= 2) {
  dd <- data.table::copy(risk_dt)
  dd[, group01 := ifelse(locked_risk_group == "High", 1, 0)]
  fit_hl <- tryCatch(survival::coxph(survival::Surv(OS_time, OS_event) ~ group01, data = dd), error = function(e) NULL)
  if (!is.null(fit_hl)) {
    sm <- summary(fit_hl)
    high_low_perf <- data.table::data.table(
      external_cohort = CONFIG$external_cohort,
      n = nrow(dd),
      events = sum(dd$OS_event == 1),
      cutoff = locked_cutoff,
      high_n = sum(dd$locked_risk_group == "High"),
      low_n = sum(dd$locked_risk_group == "Low"),
      HR_high_vs_low = exp(as.numeric(stats::coef(fit_hl)[1])),
      p = sm$coefficients[1, "Pr(>|z|)"]
    )
  }
}

continuous_perf <- data.table::data.table(
  external_cohort = CONFIG$external_cohort,
  n = cont_ev["n"],
  events = cont_ev["events"],
  cindex = cont_ev["cindex"],
  beta_per_1SD = cont_ev["beta"],
  HR_per_1SD = cont_ev["HR"],
  p = cont_ev["p"],
  reporting_note = "raw_external_C_index_values_only_no_qualitative_label"
)
write_tsv(continuous_perf, file.path(OUT$perf, "SIDISH_Step3_GSE16011_continuous_score_Cox_Cindex.tsv"))
write_tsv(high_low_perf, file.path(OUT$perf, "SIDISH_Step3_GSE16011_locked_cutoff_high_low_Cox.tsv"))

logrank_p <- NA_real_
if (data.table::uniqueN(risk_dt$locked_risk_group) >= 2) {
  logrank_p <- km_logrank_p(risk_dt$OS_time, risk_dt$OS_event, risk_dt$locked_risk_group)
}
km_stats <- data.table::data.table(
  external_cohort = CONFIG$external_cohort,
  n = nrow(risk_dt),
  events = sum(risk_dt$OS_event == 1),
  high_n = sum(risk_dt$locked_risk_group == "High"),
  low_n = sum(risk_dt$locked_risk_group == "Low"),
  logrank_p = logrank_p,
  HR_high_vs_low = if (nrow(high_low_perf) > 0) high_low_perf$HR_high_vs_low[1] else NA_real_
)
write_tsv(km_stats, file.path(OUT$perf, "SIDISH_Step3_GSE16011_KM_statistics.tsv"))

# Optional time-dependent AUC.
# v1.4: this table is never empty and preferentially reports IPCW timeROC AUC.
# If timeROC fails, it writes the failure reason and also reports an empirical
# fallback AUC excluding subjects censored before/at the time point. The fallback
# is explicitly labelled as non-IPCW and should be treated as sensitivity only.
time_unit_detected <- detect_survival_time_unit(risk_dt$OS_time)
timeROC_times_used <- choose_timeROC_times(risk_dt$OS_time, time_unit_detected)

timeROC_result <- make_timeROC_status_table(
  times = timeROC_times_used,
  unit = time_unit_detected,
  status = "not_run",
  reason = "initial_state"
)
timeROC_method_used <- "not_run"

if (!isTRUE(CONFIG$compute_timeROC_if_available)) {
  timeROC_result <- make_timeROC_status_table(
    times = timeROC_times_used,
    unit = time_unit_detected,
    status = "not_run",
    reason = "CONFIG_compute_timeROC_if_available_FALSE"
  )
  timeROC_result[, method := "none"]
  timeROC_method_used <- "none"
} else if (!requireNamespace("timeROC", quietly = TRUE)) {
  timeROC_result <- fallback_dynamic_auc_table(
    time = risk_dt$OS_time,
    event = risk_dt$OS_event,
    marker = risk_dt$locked_risk_score,
    times = timeROC_times_used,
    unit = time_unit_detected,
    fail_reason = "timeROC_package_not_installed"
  )
  timeROC_method_used <- "empirical_fallback_no_timeROC_package"
} else if (length(timeROC_times_used) == 0) {
  timeROC_result <- make_timeROC_status_table(
    times = numeric(0),
    unit = time_unit_detected,
    status = "not_run",
    reason = "no_valid_time_points_within_observed_followup_range"
  )
  timeROC_result[, method := "none"]
  timeROC_method_used <- "none"
} else if (data.table::uniqueN(risk_dt$OS_event) < 2 || sum(risk_dt$OS_event == 1, na.rm = TRUE) < 5) {
  timeROC_result <- make_timeROC_status_table(
    times = timeROC_times_used,
    unit = time_unit_detected,
    status = "not_run",
    reason = "insufficient_events_for_timeROC"
  )
  timeROC_result[, method := "none"]
  timeROC_method_used <- "none"
} else {
  # Explicitly attach survival/timeROC before calling timeROC.
  suppressPackageStartupMessages(library(survival))
  Surv <- survival::Surv
  suppressPackageStartupMessages(library(timeROC))

  troc <- tryCatch(timeROC::timeROC(
    T = risk_dt$OS_time,
    delta = risk_dt$OS_event,
    marker = risk_dt$locked_risk_score,
    cause = 1,
    times = timeROC_times_used,
    iid = TRUE
  ), error = function(e) {
    attr(e, "sidish_timeROC_error") <- conditionMessage(e)
    e
  })

  if (inherits(troc, "error")) {
    fail_reason <- paste0("timeROC_error: ", attr(troc, "sidish_timeROC_error"))
    timeROC_result <- fallback_dynamic_auc_table(
      time = risk_dt$OS_time,
      event = risk_dt$OS_event,
      marker = risk_dt$locked_risk_score,
      times = timeROC_times_used,
      unit = time_unit_detected,
      fail_reason = fail_reason
    )
    timeROC_method_used <- "empirical_fallback_after_timeROC_error"
    msg("timeROC failed; fallback empirical AUC written: ", attr(troc, "sidish_timeROC_error"))
  } else {
    auc_vec <- suppressWarnings(as.numeric(troc$AUC))
    if (length(auc_vec) == 0 || all(is.na(auc_vec))) {
      fail_reason <- "timeROC_returned_no_finite_AUC"
      timeROC_result <- fallback_dynamic_auc_table(
        time = risk_dt$OS_time,
        event = risk_dt$OS_event,
        marker = risk_dt$locked_risk_score,
        times = timeROC_times_used,
        unit = time_unit_detected,
        fail_reason = fail_reason
      )
      timeROC_method_used <- "empirical_fallback_after_empty_timeROC_AUC"
    } else {
      timeROC_result <- make_timeROC_status_table(
        times = troc$times,
        unit = time_unit_detected,
        status = "ok",
        reason = "timeROC_completed",
        auc = auc_vec
      )
      timeROC_result[, method := "IPCW_timeROC"]
      timeROC_result[, n_case := NA_integer_]
      timeROC_result[, n_control := NA_integer_]
      timeROC_result[, n_excluded_censored_before_or_at_time := NA_integer_]
      timeROC_method_used <- "IPCW_timeROC"
    }
  }
}

write_tsv(timeROC_result, file.path(OUT$perf, "SIDISH_Step3_GSE16011_time_dependent_AUC_optional.tsv"))
write_tsv(data.table::data.table(
  metric = c(
    "detected_OS_time_unit",
    "timeROC_times_used",
    "timeROC_status",
    "timeROC_reason",
    "timeROC_method_used",
    "timeROC_method_used"
  ),
  value = c(
    time_unit_detected,
    paste(timeROC_times_used, collapse = ";"),
    paste(unique(timeROC_result$status), collapse = ";"),
    paste(unique(timeROC_result$reason), collapse = ";"),
    timeROC_method_used,
    timeROC_method_used
  )
), file.path(OUT$perf, "SIDISH_Step3_timeROC_audit.tsv"))

summary_dt <- data.table::data.table(
  metric = c(
    "external_cohort",
    "n_external_samples_matched",
    "n_external_events",
    "n_locked_genes_required",
    "n_locked_genes_present",
    "external_expression_scale",
    "locked_development_cutoff",
    "high_group_n_by_locked_cutoff",
    "low_group_n_by_locked_cutoff",
    "external_Cindex_continuous_score",
    "external_HR_per_1SD_continuous_score",
    "external_p_continuous_score",
    "external_HR_high_vs_low_locked_cutoff",
    "external_p_high_vs_low_locked_cutoff",
    "external_logrank_p_locked_cutoff",
    "timeROC_detected_OS_time_unit",
    "timeROC_status",
    "timeROC_reason",
    "no_refit_validation_status"
  ),
  value = c(
    CONFIG$external_cohort,
    nrow(risk_dt),
    sum(risk_dt$OS_event == 1),
    length(locked_genes),
    length(present_locked),
    CONFIG$external_expression_scale,
    locked_cutoff,
    sum(risk_dt$locked_risk_group == "High"),
    sum(risk_dt$locked_risk_group == "Low"),
    continuous_perf$cindex,
    continuous_perf$HR_per_1SD,
    continuous_perf$p,
    if (nrow(high_low_perf) > 0) high_low_perf$HR_high_vs_low[1] else NA_real_,
    if (nrow(high_low_perf) > 0) high_low_perf$p[1] else NA_real_,
    logrank_p,
    time_unit_detected,
    paste(unique(timeROC_result$status), collapse = ";"),
    paste(unique(timeROC_result$reason), collapse = ";"),
    "completed_without_gene_selection_refit_or_cutoff_optimization"
  )
)
write_tsv(summary_dt, file.path(OUT$final, "SIDISH_Step3_GSE16011_no_refit_external_validation_summary.tsv"))

###############################################################################
# 7. Publication-style figures
###############################################################################

if (isTRUE(CONFIG$make_publication_figures)) {
  msg("Generating Step 3 external validation figures...")

  # Figure 1: locked coefficients used.
  fig1_dt <- data.table::copy(locked_coef)
  fig1_dt[, direction := ifelse(locked_coefficient >= 0, "Risk-positive", "Protective/negative")]
  fig1_dt <- fig1_dt[order(locked_coefficient)]
  fig1_dt[, gene := factor(gene, levels = gene)]
  p1 <- ggplot2::ggplot(fig1_dt, ggplot2::aes(x = locked_coefficient, y = gene, color = direction)) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = locked_coefficient, yend = gene), linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_point(size = CONFIG$plot_point_size + 0.25) +
    ggplot2::scale_color_manual(values = c("Risk-positive" = pal("red"), "Protective/negative" = pal("blue"))) +
    ggplot2::labs(title = "Step 3 | Locked model coefficients applied to GSE16011",
                  subtitle = "Coefficients imported from Step 2 v2.6 without refit",
                  x = "Locked coefficient", y = NULL, color = "Direction") +
    pub_theme()
  save_pub_plot(p1, "Figure_1_locked_coefficients_applied_to_GSE16011", width = 4.8, height = 3.6)

  # Figure 2: external risk distribution.
  p2 <- ggplot2::ggplot(risk_dt, ggplot2::aes(x = locked_risk_group, y = locked_risk_score, fill = locked_risk_group)) +
    ggplot2::geom_violin(width = 0.78, alpha = 0.75, color = "black", linewidth = CONFIG$plot_linewidth, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.16, outlier.size = 0.4, linewidth = CONFIG$plot_linewidth, color = "black", fill = "white") +
    ggplot2::geom_hline(yintercept = locked_cutoff, linetype = "dashed", color = pal("red"), linewidth = CONFIG$plot_linewidth) +
    ggplot2::scale_fill_manual(values = c(Low = CONFIG$km_palette_low, High = CONFIG$km_palette_high), guide = "none") +
    ggplot2::labs(title = "Step 3 | GSE16011 locked risk-score distribution",
                  subtitle = "Dashed line: Step 2 locked development cutoff",
                  x = "Locked risk group", y = "Locked risk score") +
    pub_theme()
  save_pub_plot(p2, "Figure_2_GSE16011_locked_risk_score_distribution", width = 4.2, height = 3.4)

  # Figure 3: KM curve.
  if (data.table::uniqueN(risk_dt$locked_risk_group) >= 2) {
    km_dt <- km_summary_table(risk_dt$OS_time, risk_dt$OS_event, risk_dt$locked_risk_group)
    if (!is.null(km_dt) && nrow(km_dt) > 0) {
      km_dt[, strata := factor(strata, levels = c("Low", "High"))]
      hr_text <- if (nrow(high_low_perf) > 0 && is.finite(high_low_perf$HR_high_vs_low[1])) {
        paste0("HR = ", sprintf("%.2f", high_low_perf$HR_high_vs_low[1]))
      } else {
        "HR = NA"
      }
      p_text <- paste0("log-rank ", format_p_for_plot(logrank_p))

      p3 <- ggplot2::ggplot(km_dt, ggplot2::aes(x = time, y = surv, color = strata)) +
        ggplot2::geom_step(linewidth = 0.7) +
        ggplot2::scale_color_manual(values = c(Low = CONFIG$km_palette_low, High = CONFIG$km_palette_high)) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(title = "Step 3 | GSE16011 no-refit Kaplan–Meier validation",
                      subtitle = paste0("Locked cutoff | ", hr_text, " | ", p_text),
                      x = CONFIG$km_time_unit_label,
                      y = "Survival probability",
                      color = "Risk group") +
        pub_theme()
      save_pub_plot(p3, "Figure_3_GSE16011_no_refit_KM_locked_cutoff", width = 4.6, height = 3.6)
    }
  }

  # Figure 4: objective external metrics.
  fig4_dt <- data.table::data.table(
    metric = c("C-index", "HR per 1SD", "HR high vs low"),
    value = c(
      continuous_perf$cindex,
      continuous_perf$HR_per_1SD,
      if (nrow(high_low_perf) > 0) high_low_perf$HR_high_vs_low[1] else NA_real_
    )
  )
  fig4_dt <- fig4_dt[is.finite(value)]
  p4 <- ggplot2::ggplot(fig4_dt, ggplot2::aes(x = metric, y = value, fill = metric)) +
    ggplot2::geom_col(width = 0.68, color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3g", value)), vjust = -0.25,
                       family = CONFIG$plot_base_family, size = 2.7) +
    ggplot2::scale_fill_manual(values = c("C-index" = pal("navy"), "HR per 1SD" = pal("orange"), "HR high vs low" = pal("red")), guide = "none") +
    ggplot2::expand_limits(y = max(fig4_dt$value, na.rm = TRUE) * 1.18) +
    ggplot2::labs(title = "Step 3 | Objective external validation metrics",
                  subtitle = "Raw values only; no model refit or cutoff optimization",
                  x = NULL, y = "Value") +
    pub_theme()
  save_pub_plot(p4, "Figure_4_GSE16011_objective_external_metrics", width = 4.4, height = 3.2)

  # Figure 5: optional time-dependent AUC.
  if (nrow(timeROC_result) > 0 && any(is.finite(timeROC_result$AUC))) {
    p5 <- ggplot2::ggplot(timeROC_result[is.finite(AUC)], ggplot2::aes(x = factor(time), y = AUC, group = 1)) +
      ggplot2::geom_line(color = pal("navy"), linewidth = CONFIG$plot_linewidth) +
      ggplot2::geom_point(color = pal("red"), size = CONFIG$plot_point_size) +
      ggplot2::coord_cartesian(ylim = c(0.45, max(0.80, max(timeROC_result$AUC, na.rm = TRUE) + 0.03))) +
      ggplot2::labs(title = "Step 3 | Optional time-dependent AUC",
                    subtitle = paste0("AUC method: ", timeROC_method_used, " | OS time unit: ", time_unit_detected),
                    x = paste0("Time (", time_unit_detected, ")"), y = "AUC") +
      pub_theme()
    save_pub_plot(p5, "Figure_5_GSE16011_time_dependent_AUC_optional", width = 4.4, height = 3.2)
  }
}

###############################################################################
# 8. README and session info
###############################################################################

readme <- c(
  "SIDISH Step 3 GSE16011 no-refit external validation v1.4",
  "",
  "Input locked model: Step 2 v2.6 direction-stable locked model.",
  "External cohort: GSE16011.
Version v1.4 fixes time-dependent AUC calculation more robustly: survival is explicitly attached for timeROC, the optional AUC TSV is never empty, OS_time units are auto-detected, and an explicitly labelled empirical non-IPCW fallback AUC is written if timeROC fails.",
  "",
  "Strict validation rule:",
  "- no new gene selection",
  "- no coefficient refit",
  "- no cutoff optimization",
  "- locked Step 2 development cutoff used as primary grouping cutoff",
  "",
  "Key final files:",
  file.path(OUT$final, "SIDISH_Step3_GSE16011_no_refit_risk_scores.tsv"),
  file.path(OUT$final, "SIDISH_Step3_GSE16011_no_refit_external_validation_summary.tsv"),
  file.path(OUT$perf, "SIDISH_Step3_GSE16011_continuous_score_Cox_Cindex.tsv"),
  file.path(OUT$perf, "SIDISH_Step3_GSE16011_locked_cutoff_high_low_Cox.tsv"),
  "",
  "Interpretation policy:",
  "C-index and AUC values are reported objectively. Do not label them as strong/weak without context."
)
write_text(readme, file.path(OUT$root, "README_SIDISH_Step3_GSE16011_no_refit_external_validation_v1_4.txt"))
write_text(capture.output(sessionInfo()), file.path(OUT$audit, "sessionInfo.txt"))

msg("SIDISH Step 3 GSE16011 no-refit external validation complete.")
msg("Summary: ", file.path(OUT$final, "SIDISH_Step3_GSE16011_no_refit_external_validation_summary.tsv"))
