#!/usr/bin/env Rscript
###############################################################################
# SIDISH Step 2: ML175-adapted direction-stable model development v2.6
#
# Purpose
# -------
# Use the locked Step 1 v3.5 Class A signed core-marker pool as the ONLY feature
# candidate set for model development.
#
# Input:
#   Step 1 v3.5 Class A pool:
#   G:/NEW/The2/SIDISH_Step1_Core_marker_refinement_v3_5_flat_inputs/
#     99_final_outputs/SIDISH_Step1_ClassA_signed_core_marker_pool_for_Step2_ML.tsv
#
# Development cohorts:
#   TCGA + CGGA325 + CGGA693 only.
#
# Strict held-out rule:
#   GSE16011 / GSE6011 must NOT be used in feature selection, model selection,
#   alpha/lambda tuning, bootstrap stability, coefficient locking, or cutoff
#   optimization. They are reserved for Step 3 external no-refit validation.
#
# Main outputs:
#   - LOCO alpha/lambda tuning table
#   - bootstrap gene stability
#   - final full elastic-net coefficients
#   - compact locked model coefficients
#   - locked development risk scores
#   - cohort-wise model performance
#   - publication-style audit figures
#
# Required packages:
#   install.packages(c("data.table", "survival", "glmnet", "ggplot2"))
#
# Run:
#   source("G:/NEW/The2/SIDISH_Step2_ML175_DirectionStable_Model_Development_v2_6_strict_all3.R")
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
  step1_classA_file = "G:/NEW/The2/SIDISH_Step1_Core_marker_refinement_v3_5_flat_inputs/99_final_outputs/SIDISH_Step1_ClassA_signed_core_marker_pool_for_Step2_ML.tsv",

  input_root = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs",
  output_dir = "G:/NEW/The2/SIDISH_Step2_ML175_DirectionStable_Model_Development_v2_6_strict_all3",

  # Exact preferred cohort files. If missing, the script searches bulk_search_roots.
  bulk_cohorts = list(
    TCGA = list(
      expr = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/TCGA/TCGA_expression_matrix.tsv",
      clin = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/TCGA/TCGA_clinical_molecular.tsv"
    ),
    CGGA325 = list(
      expr = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/CGGA_325/CGGA325_expression_matrix.tsv",
      clin = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/CGGA_325/CGGA325_clinical_molecular.tsv"
    ),
    CGGA693 = list(
      expr = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/CGGA_693/CGGA693_expression_matrix.tsv",
      clin = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/CGGA_693/CGGA693_clinical_molecular.tsv"
    )
  ),

  bulk_search_roots = c(
    "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk",
    "G:/NEW/The2/Bulk_Cohorts",
    "G:/NEW/转化/Bulk_Cohorts",
    "G:/NEW/转化",
    "G:/NEW/The2"
  ),

  # Modeling.
  expected_development_cohorts = c("TCGA", "CGGA325", "CGGA693"),

  # Strict modelability rule:
  # A gene can enter Step 2 only if it is measurable in ALL THREE development cohorts.
  # No cross-cohort missing-gene imputation is allowed.
  min_gene_available_cohorts = 3,
  require_all_development_cohorts_for_modeling = TRUE,
  allow_missing_gene_imputation = FALSE,
  min_development_genes = 5,

  # v2.1: ML175-inspired two-stage combination framework.
  # The uploaded reference script is a binary-classification AUC workflow.
  # Here it is survival-adapted: selector + learner combinations are evaluated
  # by LOCO C-index / Cox HR under strict all-three-cohort development.
  # Important: LOCO prediction learners must not include strata(cohort), because
  # the held-out cohort is absent from the training strata levels.
  feature_selector_grid = c(
    "simple_all",
    "univ_p",
    "univ_abs_beta",
    "class_priority_abs_beta",
    "cohort_consensus",
    "hybrid_rank",
    "positive_beta"
  ),
  feature_size_grid = c(5, 8, 10, 15, 33),
  learner_config_grid = c(
    "CoxRidge[alpha=0;lambda=lambda.min]",
    "CoxEnet[alpha=0.05;lambda=lambda.min]",
    "CoxEnet[alpha=0.25;lambda=lambda.min]",
    "CoxEnet[alpha=0.50;lambda=lambda.min]",
    "CoxLasso[alpha=1;lambda=lambda.min]",
    "CoxPH[nostrata]",
    "StepCoxPH[direction=both;nostrata]"
  ),
  nfolds_cv = 5,

  # Model selection rule:
  # first find the maximum LOCO mean C-index; among models within this tolerance,
  # choose the most compact model with acceptable worst-cohort performance.
  loco_mean_cindex_tolerance = 0.01,
  loco_min_cindex_floor = 0.60,
  min_locked_genes = 3,

  # v2.4 parsimony lock:
  # the ML screen is still broad, but final model locking prefers <=12 genes and
  # does not allow >15 genes unless the user explicitly changes this setting.
  parsimony_prefer_max_locked_genes = 12,
  parsimony_hard_max_locked_genes = 15,
  bootstrap_lock_min_frequency = 0.70,

  # v2.6 direction-stability gate.
  # Genes with mixed univariable Cox beta signs across TCGA/CGGA325/CGGA693
  # are excluded from the modelable feature space before ML screening.
  direction_min_valid_cohorts = 3,
  require_no_mixed_univariable_cox_sign = TRUE,
  require_step1_sign_compatibility = TRUE,
  bootstrap_sign_consistency_min = 0.80,
  final_refit_require_direction_consistency = TRUE,
  bootstrap_n = 100,

  # Reporting.
  p_value_floor_for_display = 2.2e-16,

  # Risk score/cutoff locking.
  development_cutoff_method = "median",

  # Plotting.
  make_publication_figures = TRUE,
  plot_base_family = "Arial",
  plot_base_size = 8,
  plot_linewidth = 0.5,
  plot_point_size = 2.0,
  km_time_unit_label = "Time",
  km_palette_high = "#E64B35FF",
  km_palette_low = "#4DBBD5FF",

  # Safety.
  forbidden_external_patterns = c("GSE16011", "GSE6011"),
  forbidden_old_model_patterns = c(
    "Step3",
    "external_validation",
    "GSE16011",
    "GSE6011"
  )
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
  data = file.path(CONFIG$output_dir, "01_development_data"),
  model = file.path(CONFIG$output_dir, "02_model_development"),
  perf = file.path(CONFIG$output_dir, "03_model_performance"),
  fig = file.path(CONFIG$output_dir, "04_publication_figures"),
  final = file.path(CONFIG$output_dir, "99_locked_model")
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
load_or_stop("glmnet")
load_or_stop("ggplot2")

write_tsv <- function(dt, file) {
  safe_dir_create(dirname(file))
  data.table::fwrite(data.table::as.data.table(dt), file, sep = "\t", quote = FALSE, na = "NA")
}

write_text <- function(x, file) {
  safe_dir_create(dirname(file))
  writeLines(x, con = file, useBytes = TRUE)
}

assert_safe_input <- function(file, label = "input") {
  file <- norm_path(file)
  if (file == "") return(invisible(TRUE))
  for (pat in CONFIG$forbidden_external_patterns) {
    if (grepl(pat, file, ignore.case = TRUE)) {
      stop("Forbidden external validation file detected in ", label, ": ", file)
    }
  }
  invisible(TRUE)
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
  assert_safe_input(file, "read_table_auto")
  if (!file.exists(file)) stop("File not found: ", file)
  ext <- tolower(tools::file_ext(file))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else detect_delim(file)
  dt <- data.table::fread(file, sep = sep, data.table = TRUE, showProgress = FALSE)
  data.table::setnames(dt, clean_names(names(dt)))
  dt
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

tcga_patient_id <- function(x) {
  x <- standardize_id(x)
  out <- ifelse(grepl("^TCGA-", x), sub("^(([^-]+-){2}[^-]+).*$", "\\1", x), x)
  out
}

tcga_sample_id <- function(x) {
  x <- standardize_id(x)
  out <- ifelse(grepl("^TCGA-", x), sub("^(([^-]+-){3}[^-]+).*$", "\\1", x), x)
  out
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

###############################################################################
# 2. Input and bulk discovery
###############################################################################

msg("SIDISH Step 2 v2.6 ML175-adapted direction-stable strict-all3 model development started.")

if (!file.exists(CONFIG$step1_classA_file)) stop("Step 1 Class A file not found: ", CONFIG$step1_classA_file)
assert_safe_input(CONFIG$step1_classA_file, "step1_classA_file")

classA <- read_table_auto(CONFIG$step1_classA_file)
gene_col <- first_existing_col(classA, c("^gene$|gene_symbol|symbol"), regex = TRUE)
if (is.na(gene_col)) stop("Cannot identify gene column in Step 1 Class A file.")

classA[, gene := normalize_gene(get(gene_col))]
classA <- classA[!is.na(gene) & gene != ""]
classA <- unique(classA, by = "gene")
candidate_genes_all <- classA$gene

if (length(candidate_genes_all) < CONFIG$min_development_genes) {
  stop("Too few Step 1 Class A genes: ", length(candidate_genes_all))
}

write_tsv(classA, file.path(OUT$audit, "Step2_input_ClassA_gene_pool.tsv"))

list_bulk_files <- function(roots) {
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

bulk_inventory <- list_bulk_files(CONFIG$bulk_search_roots)
if (nrow(bulk_inventory) > 0) {
  bulk_inventory[, forbidden_external := grepl(paste(CONFIG$forbidden_external_patterns, collapse = "|"), file, ignore.case = TRUE)]
  bulk_inventory <- bulk_inventory[forbidden_external == FALSE]
}
write_tsv(bulk_inventory, file.path(OUT$audit, "bulk_candidate_file_inventory.tsv"))

score_file_candidate <- function(files_dt, include_patterns, exclude_patterns = character(0),
                                 allowed_ext = c("tsv", "txt", "csv")) {
  if (nrow(files_dt) == 0) return("")
  dt <- data.table::copy(files_dt[ext %in% allowed_ext])
  if (nrow(dt) == 0) return("")
  dt[, score := 0]
  for (pat in include_patterns) {
    dt[grepl(pat, basename, ignore.case = TRUE), score := score + 3]
    dt[grepl(pat, file, ignore.case = TRUE), score := score + 1]
  }
  for (pat in exclude_patterns) {
    dt[grepl(pat, basename, ignore.case = TRUE) | grepl(pat, file, ignore.case = TRUE), score := score - 6]
  }
  dt <- dt[score > 0]
  if (nrow(dt) == 0) return("")
  dt <- dt[order(-score, size_bytes)]
  dt$file[1]
}

is_valid_cohort_file <- function(file, cohort) {
  file <- norm_path(file)
  if (file == "" || !file.exists(file)) return(FALSE)
  b <- basename(file)
  d <- basename(dirname(file))

  if (cohort == "TCGA") {
    return(grepl("TCGA", b, ignore.case = TRUE) || grepl("TCGA", d, ignore.case = TRUE))
  }
  if (cohort == "CGGA325") {
    return(grepl("CGGA.*325|CGGA325|CGGA_325", b, ignore.case = TRUE) ||
             grepl("CGGA.*325|CGGA325|CGGA_325", d, ignore.case = TRUE))
  }
  if (cohort == "CGGA693") {
    return(grepl("CGGA.*693|CGGA693|CGGA_693", b, ignore.case = TRUE) ||
             grepl("CGGA.*693|CGGA693|CGGA_693", d, ignore.case = TRUE))
  }
  FALSE
}

discover_bulk_file <- function(cohort, type = c("expr", "clin")) {
  type <- match.arg(type)
  exact <- CONFIG$bulk_cohorts[[cohort]][[type]]
  exact <- norm_path(exact)
  if (exact != "" && file.exists(exact)) {
    assert_safe_input(exact, paste0(cohort, "_", type))
    if (!is_valid_cohort_file(exact, cohort)) {
      stop("Configured ", cohort, " ", type, " file does not look cohort-specific: ", exact)
    }
    return(exact)
  }

  # v1.3: search within the cohort-specific folder first, then global inventory.
  cohort_dir_patterns <- switch(cohort,
    TCGA = "TCGA",
    CGGA325 = "CGGA.*325|CGGA325|CGGA_325",
    CGGA693 = "CGGA.*693|CGGA693|CGGA_693"
  )

  inv <- bulk_inventory
  if (nrow(inv) > 0) {
    inv[, cohort_path_flag := grepl(cohort_dir_patterns, file, ignore.case = TRUE)]
    inv <- inv[cohort_path_flag == TRUE]
  }

  if (type == "expr") {
    include <- switch(cohort,
      TCGA = c("TCGA.*expr|TCGA.*expression|TCGA.*matrix", "expr|expression|matrix|mRNA|RNAseq|TPM|FPKM|RSEM|count"),
      CGGA325 = c("CGGA.*325.*expr|CGGA325.*expr|CGGA_325.*expr|CGGA.*325.*matrix", "expr|expression|matrix|mRNA|RNAseq|TPM|FPKM|RSEM|count"),
      CGGA693 = c("CGGA.*693.*expr|CGGA693.*expr|CGGA_693.*expr|CGGA.*693.*matrix", "expr|expression|matrix|mRNA|RNAseq|TPM|FPKM|RSEM|count")
    )
    exclude <- c("clinical|clin|survival|phenotype|annotation|meta")
  } else {
    include <- switch(cohort,
      TCGA = c("TCGA.*clinical|TCGA.*clin|TCGA.*molecular", "clinical|clin|survival|phenotype|molecular|meta"),
      CGGA325 = c("CGGA.*325.*clinical|CGGA325.*clinical|CGGA_325.*clinical|CGGA.*325.*molecular", "clinical|clin|survival|phenotype|molecular|meta"),
      CGGA693 = c("CGGA.*693.*clinical|CGGA693.*clinical|CGGA_693.*clinical|CGGA.*693.*molecular", "clinical|clin|survival|phenotype|molecular|meta")
    )
    exclude <- c("expr|expression|matrix|mRNA|RNAseq|TPM|FPKM|RSEM|count")
  }

  cand <- score_file_candidate(inv, include, exclude)
  if (cand != "" && is_valid_cohort_file(cand, cohort)) return(cand)

  ""
}

selected_bulk <- data.table::rbindlist(lapply(names(CONFIG$bulk_cohorts), function(coh) {
  data.table::data.table(
    cohort = coh,
    expr_file = discover_bulk_file(coh, "expr"),
    clin_file = discover_bulk_file(coh, "clin")
  )
}), fill = TRUE)
selected_bulk[, expr_exists := expr_file != "" & file.exists(expr_file)]
selected_bulk[, clin_exists := clin_file != "" & file.exists(clin_file)]
write_tsv(selected_bulk, file.path(OUT$audit, "selected_bulk_input_files.tsv"))

###############################################################################
# 3. Bulk read, sample matching, and development matrix
###############################################################################

standardize_bulk_expr <- function(file) {
  dt <- read_table_auto(file)
  if (nrow(dt) == 0 || ncol(dt) < 2) stop("Expression table too small: ", file)

  gene_col <- first_existing_col(dt, c("^gene$", "gene_symbol", "^symbol$", "hgnc", "feature", "ensembl"), regex = TRUE)
  if (is.na(gene_col)) {
    first_col <- names(dt)[1]
    first_num_frac <- mean(!is.na(as_numeric_safe(dt[[first_col]])))
    if (first_num_frac < 0.5) gene_col <- first_col
  }
  if (is.na(gene_col)) stop("Cannot identify gene column in expression file: ", file)

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
  mat
}

standardize_clinical <- function(file) {
  dt <- read_table_auto(file)
  if (nrow(dt) == 0) stop("Clinical table empty: ", file)

  sample_col <- first_existing_col(dt, c("^sample_id$", "^sample$", "barcode", "case", "patient", "^ID$"), regex = TRUE)
  if (is.na(sample_col)) sample_col <- names(dt)[1]

  time_col <- first_existing_col(dt, c("^OS_time$", "^os_time$", "survival.*time", "follow.*time", "days", "months"), regex = TRUE)
  event_col <- first_existing_col(dt, c("^OS_event$", "^os_event$", "survival.*event", "status", "death", "event"), regex = TRUE)
  age_col <- first_existing_col(dt, c("^Age$", "^age$"), regex = TRUE)
  grade_col <- first_existing_col(dt, c("WHO_grade", "grade", "histology"), regex = TRUE)
  idh_col <- first_existing_col(dt, c("IDH_status", "IDH", "idh"), regex = TRUE)

  if (is.na(time_col) || is.na(event_col)) stop("Clinical survival columns not found: ", file)

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
    ev %in% c("1", "DEAD", "DECEASED", "EVENT", "TRUE", "YES", "1:DECEASED"), 1,
    ev %in% c("0", "ALIVE", "CENSORED", "FALSE", "NO", "LIVING", "0:LIVING"), 0,
    default = suppressWarnings(as.numeric(ev))
  )]
  out <- out[!is.na(sample_id) & sample_id != "" & is.finite(OS_time) & OS_time > 0 & OS_event %in% c(0, 1)]
  unique(out, by = "sample_id")
}

match_samples <- function(expr_samples, clin_dt, cohort) {
  expr_dt <- data.table::data.table(expr_sample = expr_samples)
  expr_dt[, exact := expr_sample]
  expr_dt[, clean := standardize_id(expr_sample)]
  expr_dt[, tcga_sample := tcga_sample_id(expr_sample)]
  expr_dt[, tcga_patient := tcga_patient_id(expr_sample)]

  clin_map <- data.table::copy(clin_dt)
  clin_map[, exact := sample_id]
  clin_map[, clean := standardize_id(sample_id)]
  clin_map[, tcga_sample := tcga_sample_id(sample_id)]
  clin_map[, tcga_patient := tcga_patient_id(sample_id)]

  modes <- c("exact", "clean", "tcga_sample", "tcga_patient")
  reports <- list()
  best <- NULL
  best_n <- -1

  for (mode in modes) {
    e <- expr_dt[, .(expr_sample, key = get(mode))]
    c <- clin_map[, .(sample_id, key = get(mode))]
    e <- e[!is.na(key) & key != ""]
    c <- c[!is.na(key) & key != ""]
    c <- unique(c, by = "key")
    m <- merge(e, c, by = "key", all = FALSE)
    m <- unique(m, by = "expr_sample")
    reports[[mode]] <- data.table::data.table(cohort = cohort, match_mode = mode, n_matched = nrow(m))
    if (nrow(m) > best_n) {
      best <- m
      best_n <- nrow(m)
      best_mode <- mode
    }
  }

  list(matches = best, mode = best_mode, report = data.table::rbindlist(reports))
}

cohort_objects <- list()
cohort_status <- data.table::data.table()
gene_avail <- data.table::data.table()

for (i in seq_len(nrow(selected_bulk))) {
  coh <- selected_bulk$cohort[i]
  expr_file <- selected_bulk$expr_file[i]
  clin_file <- selected_bulk$clin_file[i]

  if (expr_file == "" || clin_file == "" || !file.exists(expr_file) || !file.exists(clin_file)) {
    cohort_status <- rbind(cohort_status, data.table::data.table(
      cohort = coh, status = "missing_expr_or_clin", n_samples = 0, n_genes_present = 0, match_mode = NA_character_
    ))
    next
  }

  msg("Reading cohort: ", coh)
  expr <- tryCatch(standardize_bulk_expr(expr_file), error = function(e) {
    cohort_status <<- rbind(cohort_status, data.table::data.table(cohort = coh, status = paste0("expr_error:", conditionMessage(e)), n_samples = 0, n_genes_present = 0, match_mode = NA_character_))
    NULL
  })
  clin <- tryCatch(standardize_clinical(clin_file), error = function(e) {
    cohort_status <<- rbind(cohort_status, data.table::data.table(cohort = coh, status = paste0("clin_error:", conditionMessage(e)), n_samples = 0, n_genes_present = 0, match_mode = NA_character_))
    NULL
  })
  if (is.null(expr) || is.null(clin)) next

  matches <- match_samples(colnames(expr), clin, coh)
  write_tsv(matches$report, file.path(OUT$audit, paste0("Step2_", coh, "_sample_match_mode_report.tsv")))

  if (nrow(matches$matches) < 20) {
    cohort_status <- rbind(cohort_status, data.table::data.table(
      cohort = coh, status = "too_few_matched_samples", n_samples = nrow(matches$matches), n_genes_present = 0, match_mode = matches$mode
    ))
    next
  }

  common_expr_samples <- matches$matches$expr_sample
  common_clin_samples <- matches$matches$sample_id
  expr <- expr[, common_expr_samples, drop = FALSE]
  clin <- clin[match(common_clin_samples, sample_id)]

  # Optional log transform if values are count-like or large scale.
  sample_vals <- as.numeric(expr[seq_len(min(nrow(expr), 100)), seq_len(min(ncol(expr), 100)), drop = FALSE])
  if (length(sample_vals) > 0 && max(sample_vals, na.rm = TRUE) > 50) {
    expr <- log2(expr + 1)
  }

  genes_present <- intersect(candidate_genes_all, rownames(expr))
  gene_avail <- rbind(gene_avail, data.table::data.table(cohort = coh, gene = candidate_genes_all, present = candidate_genes_all %in% rownames(expr)))

  if (length(genes_present) < CONFIG$min_development_genes) {
    cohort_status <- rbind(cohort_status, data.table::data.table(
      cohort = coh, status = "too_few_candidate_genes_present", n_samples = ncol(expr), n_genes_present = length(genes_present), match_mode = matches$mode
    ))
    next
  }

  cohort_objects[[coh]] <- list(expr = expr, clin = clin, samples = common_expr_samples, genes_present = genes_present)
  cohort_status <- rbind(cohort_status, data.table::data.table(
    cohort = coh, status = "ok", n_samples = ncol(expr), n_genes_present = length(genes_present), match_mode = matches$mode
  ))
}

write_tsv(cohort_status, file.path(OUT$audit, "Step2_bulk_cohort_read_status.tsv"))
write_tsv(gene_avail, file.path(OUT$audit, "Step2_ClassA_gene_availability_long.tsv"))

# v1.3 hard gate: strict-all3 means all expected development cohorts must be present
# and successfully matched. A two-cohort run must not continue.
missing_expected <- setdiff(CONFIG$expected_development_cohorts, names(cohort_objects))
not_ok_expected <- cohort_status[cohort %in% CONFIG$expected_development_cohorts & status != "ok", cohort]
if (length(missing_expected) > 0 || length(not_ok_expected) > 0) {
  stop(
    "Strict all-three-cohort model development cannot proceed. Missing/not-ok cohorts: ",
    paste(unique(c(missing_expected, not_ok_expected)), collapse = ", "),
    ". Inspect 00_audit/selected_bulk_input_files.tsv and Step2_*_sample_match_mode_report.tsv. ",
    "TCGA must use TCGA_clinical_molecular.tsv, not a CGGA clinical file."
  )
}
if (length(cohort_objects) != length(CONFIG$expected_development_cohorts)) {
  stop("Internal strict-all3 failure: expected ", length(CONFIG$expected_development_cohorts),
       " cohort objects, got ", length(cohort_objects))
}

gene_avail_summary <- gene_avail[, .(
  n_cohorts_present = sum(present, na.rm = TRUE),
  cohorts_present = paste(cohort[present == TRUE], collapse = ";"),
  cohorts_absent = paste(cohort[present == FALSE], collapse = ";")
), by = gene]
gene_avail_summary[, present_in_all_development_cohorts := n_cohorts_present == length(CONFIG$expected_development_cohorts)]
gene_avail_summary[, included_in_modeling := present_in_all_development_cohorts == TRUE]
gene_avail_summary <- gene_avail_summary[order(-included_in_modeling, -n_cohorts_present, gene)]
write_tsv(gene_avail_summary, file.path(OUT$audit, "Step2_ClassA_gene_availability_summary.tsv"))

gene_presence_wide <- data.table::dcast(
  gene_avail,
  gene ~ cohort,
  value.var = "present",
  fill = FALSE
)
gene_presence_wide <- merge(
  gene_avail_summary[, .(gene, n_cohorts_present, cohorts_present, cohorts_absent, included_in_modeling)],
  gene_presence_wide,
  by = "gene",
  all.x = TRUE
)
write_tsv(gene_presence_wide, file.path(OUT$audit, "Step2_ClassA_gene_presence_matrix_strict_all3.tsv"))

# Export expression-gene examples to diagnose gene-symbol mismatch when all-three overlap is unexpectedly small.
expr_symbol_examples <- data.table::rbindlist(lapply(names(cohort_objects), function(coh) {
  rr <- rownames(cohort_objects[[coh]]$expr)
  data.table::data.table(
    cohort = coh,
    example_rank = seq_len(min(100, length(rr))),
    example_expression_gene_symbol = rr[seq_len(min(100, length(rr)))]
  )
}), fill = TRUE)
write_tsv(expr_symbol_examples, file.path(OUT$audit, "Step2_expression_gene_symbol_examples_first100.tsv"))

model_genes <- gene_avail_summary[included_in_modeling == TRUE, gene]

if (length(model_genes) < CONFIG$min_development_genes) {
  stop(
    "Too few Class A genes are measurable in all development cohorts: ", length(model_genes),
    ". Strict all-three-cohort modeling is not valid with the current gene pool. ",
    "Inspect 00_audit/Step2_ClassA_gene_availability_summary.tsv, ",
    "00_audit/Step2_ClassA_gene_presence_matrix_strict_all3.tsv, and ",
    "00_audit/Step2_expression_gene_symbol_examples_first100.tsv. ",
    "Do NOT use missing-gene imputation for the final model."
  )
}

# Construct combined development dataset.
# Strict rule: no missing-gene imputation. All model_genes must exist in every
# development cohort expression matrix.
dev_list <- list()
strict_matrix_audit <- list()
audit_i <- 1L

for (coh in names(cohort_objects)) {
  obj <- cohort_objects[[coh]]
  missing_here <- setdiff(model_genes, rownames(obj$expr))
  if (length(missing_here) > 0) {
    stop("Internal strict-rule failure: cohort ", coh, " is missing model genes: ",
         paste(missing_here, collapse = ", "))
  }

  expr <- obj$expr[model_genes, , drop = FALSE]
  expr_z <- zscore_rows(expr)

  strict_matrix_audit[[audit_i]] <- data.table::data.table(
    cohort = coh,
    gene = model_genes,
    present_in_expression_matrix = TRUE,
    imputed = FALSE
  )
  audit_i <- audit_i + 1L

  x_dt <- data.table::as.data.table(t(expr_z), keep.rownames = "expr_sample")
  data.table::setnames(x_dt, "expr_sample", "sample_id_expr")
  clin <- data.table::copy(obj$clin)
  x_dt[, sample_id := clin$sample_id]
  x_dt[, cohort := coh]
  x_dt[, OS_time := clin$OS_time]
  x_dt[, OS_event := clin$OS_event]
  dev_list[[coh]] <- x_dt

  msg("Strict development matrix cohort ", coh, ": genes=", length(model_genes),
      "; samples=", ncol(expr), "; no imputation")
}

strict_matrix_audit <- data.table::rbindlist(strict_matrix_audit, fill = TRUE)
write_tsv(strict_matrix_audit, file.path(OUT$data, "Step2_strict_all3_model_gene_matrix_audit.tsv"))
dev_dt <- data.table::rbindlist(dev_list, fill = TRUE)
dev_dt <- dev_dt[is.finite(OS_time) & OS_time > 0 & OS_event %in% c(0, 1)]
write_tsv(dev_dt, file.path(OUT$data, "Step2_development_matrix_zscored.tsv"))

x <- as.matrix(dev_dt[, model_genes, with = FALSE])
storage.mode(x) <- "numeric"
x[!is.finite(x)] <- 0
y <- survival::Surv(dev_dt$OS_time, dev_dt$OS_event)

write_tsv(data.table::data.table(
  metric = c(
    "n_samples",
    "n_events",
    "n_cohorts",
    "n_input_ClassA_genes",
    "n_modeling_genes_present_in_all_development_cohorts_before_direction_gate",
    "required_gene_available_cohorts",
    "missing_gene_imputation_allowed"
  ),
  value = c(
    nrow(dev_dt),
    sum(dev_dt$OS_event == 1),
    data.table::uniqueN(dev_dt$cohort),
    length(candidate_genes_all),
    length(model_genes),
    CONFIG$min_gene_available_cohorts,
    CONFIG$allow_missing_gene_imputation
  )
), file.path(OUT$data, "Step2_development_matrix_summary.tsv"))


###############################################################################
# 3.5 Direction-stability audit before ML screening
###############################################################################

msg("Running v2.6 direction-stability audit before ML screening...")

model_genes_all3 <- as.character(model_genes)

per_cohort_gene_cox <- function(gene_vec) {
  out <- list()
  kk <- 1L
  for (g in gene_vec) {
    for (coh in CONFIG$expected_development_cohorts) {
      idx <- which(dev_dt$cohort == coh)
      dat <- data.frame(
        OS_time = dev_dt$OS_time[idx],
        OS_event = dev_dt$OS_event[idx],
        x = as.numeric(x[idx, g])
      )
      dat <- dat[is.finite(dat$OS_time) & dat$OS_time > 0 &
                   dat$OS_event %in% c(0, 1) & is.finite(dat$x), , drop = FALSE]
      beta <- HR <- p <- cindex <- NA_real_
      n <- nrow(dat)
      events <- sum(dat$OS_event == 1, na.rm = TRUE)
      if (n >= 20 && events >= 5 && stats::sd(dat$x, na.rm = TRUE) > 0) {
        dat$x <- as.numeric(scale(dat$x))
        fit <- tryCatch(survival::coxph(survival::Surv(OS_time, OS_event) ~ x, data = dat), error = function(e) NULL)
        if (!is.null(fit)) {
          sm <- summary(fit)
          beta <- as.numeric(stats::coef(fit)["x"])
          HR <- exp(beta)
          p <- sm$coefficients["x", "Pr(>|z|)"]
          cindex <- as.numeric(sm$concordance[1])
        }
      }
      out[[kk]] <- data.table::data.table(
        gene = g,
        cohort = coh,
        n = n,
        events = events,
        beta = beta,
        HR = HR,
        p = p,
        cindex = cindex,
        beta_sign = data.table::fcase(
          is.na(beta), NA_integer_,
          beta > 0, 1L,
          beta < 0, -1L,
          beta == 0, 0L
        )
      )
      kk <- kk + 1L
    }
  }
  data.table::rbindlist(out, fill = TRUE)
}

direction_long <- per_cohort_gene_cox(model_genes_all3)
write_tsv(direction_long, file.path(OUT$model, "Step2_v2_6_per_cohort_univariable_Cox_direction_long.tsv"))

direction_summary <- direction_long[, .(
  n_valid_cohorts = sum(!is.na(beta_sign) & beta_sign != 0),
  n_beta_positive = sum(beta_sign == 1L, na.rm = TRUE),
  n_beta_negative = sum(beta_sign == -1L, na.rm = TRUE),
  n_beta_zero = sum(beta_sign == 0L, na.rm = TRUE),
  beta_sign_pattern = paste(ifelse(is.na(beta_sign), "NA", beta_sign), collapse = ";"),
  mean_beta = mean(beta, na.rm = TRUE),
  median_beta = median(beta, na.rm = TRUE),
  min_beta = suppressWarnings(min(beta, na.rm = TRUE)),
  max_beta = suppressWarnings(max(beta, na.rm = TRUE)),
  min_p = suppressWarnings(min(p, na.rm = TRUE)),
  mean_cindex_univariable = mean(cindex, na.rm = TRUE)
), by = gene]
direction_summary[is.infinite(min_beta), min_beta := NA_real_]
direction_summary[is.infinite(max_beta), max_beta := NA_real_]
direction_summary[is.infinite(min_p), min_p := NA_real_]

direction_summary[, univariable_consensus_sign := data.table::fcase(
  n_valid_cohorts >= CONFIG$direction_min_valid_cohorts & n_beta_positive == n_valid_cohorts, 1L,
  n_valid_cohorts >= CONFIG$direction_min_valid_cohorts & n_beta_negative == n_valid_cohorts, -1L,
  default = NA_integer_
)]
direction_summary[, no_mixed_univariable_cox_sign := !is.na(univariable_consensus_sign)]

direction_summary <- merge(
  direction_summary,
  classA[, .(gene, clinical_model_class, Step2_feature_sign, bulk_direction, mean_cox_beta, median_cox_beta)],
  by = "gene",
  all.x = TRUE
)

direction_summary[, step1_expected_sign := data.table::fcase(
  Step2_feature_sign == "risk_positive", 1L,
  Step2_feature_sign == "protective_or_subtype_negative", -1L,
  default = NA_integer_
)]
direction_summary[, step1_sign_compatible := data.table::fcase(
  is.na(step1_expected_sign), TRUE,
  !is.na(univariable_consensus_sign) & univariable_consensus_sign == step1_expected_sign, TRUE,
  default = FALSE
)]

direction_summary[, direction_stable_for_modeling := (
  n_valid_cohorts >= CONFIG$direction_min_valid_cohorts &
    (!isTRUE(CONFIG$require_no_mixed_univariable_cox_sign) | no_mixed_univariable_cox_sign == TRUE) &
    (!isTRUE(CONFIG$require_step1_sign_compatibility) | step1_sign_compatible == TRUE)
)]

direction_summary[, direction_exclusion_reason := data.table::fcase(
  direction_stable_for_modeling == TRUE, "kept_direction_stable",
  n_valid_cohorts < CONFIG$direction_min_valid_cohorts, "excluded_insufficient_valid_cohort_cox",
  no_mixed_univariable_cox_sign == FALSE, "excluded_mixed_per_cohort_univariable_cox_sign",
  step1_sign_compatible == FALSE, "excluded_step1_sign_incompatible",
  default = "excluded_other_direction_instability"
)]

direction_summary <- direction_summary[order(-direction_stable_for_modeling, direction_exclusion_reason, gene)]
write_tsv(direction_summary, file.path(OUT$model, "Step2_v2_6_gene_direction_stability_summary.tsv"))
write_tsv(direction_summary[direction_stable_for_modeling == FALSE],
          file.path(OUT$model, "Step2_v2_6_direction_unstable_excluded_genes.tsv"))

model_genes_before_direction_gate <- model_genes_all3
model_genes <- direction_summary[direction_stable_for_modeling == TRUE, gene]

if (length(model_genes) < CONFIG$min_development_genes) {
  stop(
    "Too few direction-stable genes remain after v2.6 Cox direction gate: ",
    length(model_genes), ". Inspect Step2_v2_6_gene_direction_stability_summary.tsv. ",
    "Do not force unstable genes into the final model."
  )
}

# Restrict the ML matrix to direction-stable genes only.
x <- x[, model_genes, drop = FALSE]

write_tsv(data.table::data.table(
  metric = c(
    "n_modeling_genes_before_direction_gate",
    "n_direction_stable_modeling_genes",
    "n_direction_unstable_excluded_genes",
    "direction_min_valid_cohorts",
    "require_no_mixed_univariable_cox_sign",
    "require_step1_sign_compatibility"
  ),
  value = c(
    length(model_genes_before_direction_gate),
    length(model_genes),
    length(setdiff(model_genes_before_direction_gate, model_genes)),
    CONFIG$direction_min_valid_cohorts,
    CONFIG$require_no_mixed_univariable_cox_sign,
    CONFIG$require_step1_sign_compatibility
  )
), file.path(OUT$model, "Step2_v2_6_direction_gate_summary.tsv"))


###############################################################################
# 4. Large multi-configuration LOCO model screening
###############################################################################

msg("Building large candidate model library and running LOCO screening on direction-stable genes only...")

model_genes <- as.character(model_genes)
cohorts <- unique(dev_dt$cohort)
expected_n_cohorts <- length(CONFIG$expected_development_cohorts)

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

safe_univ_gene_stats <- function(train_idx, genes) {
  rows <- vector("list", length(genes))
  for (ii in seq_along(genes)) {
    g <- genes[ii]
    dat <- data.frame(
      OS_time = dev_dt$OS_time[train_idx],
      OS_event = dev_dt$OS_event[train_idx],
      cohort = dev_dt$cohort[train_idx],
      x = as.numeric(x[train_idx, g])
    )
    dat <- dat[is.finite(dat$OS_time) & dat$OS_time > 0 & dat$OS_event %in% c(0, 1) & is.finite(dat$x), , drop = FALSE]
    beta <- HR <- p <- cindex <- NA_real_
    if (nrow(dat) >= 20 && sum(dat$OS_event == 1) >= 5 && stats::sd(dat$x, na.rm = TRUE) > 0) {
      dat$x <- as.numeric(scale(dat$x))
      fit <- tryCatch(survival::coxph(survival::Surv(OS_time, OS_event) ~ x + strata(cohort), data = dat), error = function(e) NULL)
      if (!is.null(fit)) {
        sm <- summary(fit)
        beta <- as.numeric(stats::coef(fit)["x"])
        HR <- exp(beta)
        p <- sm$coefficients["x", "Pr(>|z|)"]
        cindex <- as.numeric(sm$concordance[1])
      }
    }
    rows[[ii]] <- data.table::data.table(gene = g, beta = beta, HR = HR, p = p, cindex = cindex)
  }
  data.table::rbindlist(rows, fill = TRUE)
}

cohort_consensus_stats <- function(train_idx, genes) {
  train_cohorts <- unique(dev_dt$cohort[train_idx])
  rows <- list()
  rr <- 1L
  for (g in genes) {
    betas <- numeric(0)
    for (coh in train_cohorts) {
      idx <- train_idx[dev_dt$cohort[train_idx] == coh]
      dat <- data.frame(
        OS_time = dev_dt$OS_time[idx],
        OS_event = dev_dt$OS_event[idx],
        x = as.numeric(x[idx, g])
      )
      dat <- dat[is.finite(dat$OS_time) & dat$OS_time > 0 & dat$OS_event %in% c(0, 1) & is.finite(dat$x), , drop = FALSE]
      if (nrow(dat) < 20 || sum(dat$OS_event == 1) < 5 || stats::sd(dat$x, na.rm = TRUE) <= 0) next
      dat$x <- as.numeric(scale(dat$x))
      fit <- tryCatch(survival::coxph(survival::Surv(OS_time, OS_event) ~ x, data = dat), error = function(e) NULL)
      if (!is.null(fit)) betas <- c(betas, as.numeric(stats::coef(fit)["x"]))
    }
    mean_beta <- if (length(betas) > 0) mean(betas, na.rm = TRUE) else NA_real_
    frac_same <- if (length(betas) > 0 && is.finite(mean_beta) && mean_beta != 0) {
      mean(sign(betas) == sign(mean_beta), na.rm = TRUE)
    } else 0
    rows[[rr]] <- data.table::data.table(
      gene = g,
      cohort_mean_beta = mean_beta,
      cohort_abs_mean_beta = abs(mean_beta),
      cohort_direction_fraction = frac_same,
      cohort_valid_n = length(betas),
      cohort_consensus_score = abs(mean_beta) * frac_same * length(betas)
    )
    rr <- rr + 1L
  }
  data.table::rbindlist(rows, fill = TRUE)
}

rank_genes_for_training <- function(train_idx, selector, genes, feature_size) {
  feature_size <- min(feature_size, length(genes))

  u <- safe_univ_gene_stats(train_idx, genes)
  cc <- cohort_consensus_stats(train_idx, genes)

  rr <- merge(u, cc, by = "gene", all = TRUE)
  rr <- merge(
    rr,
    classA[, .(gene, clinical_model_class, Step2_feature_sign)],
    by = "gene",
    all.x = TRUE
  )
  rr[, class_priority := data.table::fcase(
    grepl("Class_A1", clinical_model_class), 1,
    clinical_model_class == "Class_A2_borderline_context_dependent_candidate", 2,
    grepl("Class_A2", clinical_model_class), 3,
    default = 4
  )]
  rr[, rank_abs_beta := data.table::frank(-abs(beta), ties.method = "average", na.last = "keep")]
  rr[, rank_p := data.table::frank(p, ties.method = "average", na.last = "keep")]
  rr[, rank_class := data.table::frank(class_priority, ties.method = "average", na.last = "keep")]
  rr[, rank_consensus := data.table::frank(-cohort_consensus_score, ties.method = "average", na.last = "keep")]
  rr[, hybrid_score := rank_abs_beta + rank_p + rank_class + rank_consensus]

  if (selector == "simple_all") {
    rr <- rr[order(gene)]
  } else if (selector == "univ_abs_beta") {
    rr <- rr[order(-abs(beta), p, gene)]
  } else if (selector == "univ_p") {
    rr <- rr[order(p, -abs(beta), gene)]
  } else if (selector == "positive_beta") {
    rr <- rr[order(-beta, p, gene)]
  } else if (selector == "class_priority_abs_beta") {
    rr <- rr[order(class_priority, -abs(beta), p, gene)]
  } else if (selector == "cohort_consensus") {
    rr <- rr[order(-cohort_consensus_score, p, gene)]
  } else if (selector == "hybrid_rank") {
    rr <- rr[order(hybrid_score, p, gene)]
  } else {
    stop("Unknown selector: ", selector)
  }

  selected <- rr[!is.na(gene), gene][seq_len(min(feature_size, nrow(rr)))]
  list(genes = selected, ranking_table = rr)
}

parse_learner_config <- function(x) {
  learner <- sub("\\[.*$", "", x)
  inside <- ifelse(grepl("\\[", x), sub("^.*\\[(.*)\\]$", "\\1", x), "")
  out <- list(
    learner_config = x,
    learner = learner,
    alpha = NA_real_,
    lambda_policy = NA_character_,
    step_direction = NA_character_,
    use_strata = grepl("strata=cohort", inside)
  )
  if (grepl("alpha=", inside)) {
    out$alpha <- as.numeric(sub("^.*alpha=([^;]+).*$", "\\1", inside))
  }
  if (grepl("lambda=", inside)) {
    out$lambda_policy <- sub("^.*lambda=([^;]+).*$", "\\1", inside)
  }
  if (grepl("direction=", inside)) {
    out$step_direction <- sub("^.*direction=([^;]+).*$", "\\1", inside)
  }
  out
}

safe_gene_terms <- function(genes) {
  paste0("`", gsub("`", "", genes), "`")
}


fit_predict_candidate <- function(train_idx, test_idx, selected_genes, learner_config) {
  selected_genes <- intersect(selected_genes, colnames(x))
  if (length(selected_genes) < CONFIG$min_locked_genes) return(NULL)

  lc <- parse_learner_config(learner_config)

  if (isTRUE(lc$use_strata)) {
    train_cohorts <- unique(as.character(dev_dt$cohort[train_idx]))
    test_cohorts <- unique(as.character(dev_dt$cohort[test_idx]))
    if (length(setdiff(test_cohorts, train_cohorts)) > 0) {
      # LOCO holdout cohort is by definition absent from training; strata(cohort)
      # cannot extrapolate to an unseen stratum. Disable strata for this fit.
      lc$use_strata <- FALSE
    }
  }

  if (lc$learner %in% c("CoxRidge", "CoxEnet", "CoxLasso")) {
    cvfit <- tryCatch(glmnet::cv.glmnet(
      x = x[train_idx, selected_genes, drop = FALSE],
      y = y[train_idx],
      family = "cox",
      alpha = lc$alpha,
      nfolds = CONFIG$nfolds_cv,
      standardize = FALSE,
      cox.ties = "breslow"
    ), error = function(e) NULL)
    if (is.null(cvfit)) return(NULL)

    risk <- as.numeric(stats::predict(cvfit, newx = x[test_idx, selected_genes, drop = FALSE], s = lc$lambda_policy, type = "link"))
    cc <- as.matrix(stats::coef(cvfit, s = lc$lambda_policy))
    coef_dt <- data.table::data.table(gene = rownames(cc), coefficient = as.numeric(cc[, 1]))
    coef_dt <- coef_dt[coefficient != 0 & !is.na(coefficient)]
    return(list(risk = risk, coef = coef_dt, model_object = cvfit, n_nonzero = nrow(coef_dt)))
  }

  if (lc$learner == "CoxPH") {
    train_dt <- as.data.frame(dev_dt[train_idx, c("OS_time", "OS_event", "cohort", selected_genes), with = FALSE])
    rhs <- paste(safe_gene_terms(selected_genes), collapse = " + ")
    if (isTRUE(lc$use_strata)) rhs <- paste0(rhs, " + strata(cohort)")
    form <- stats::as.formula(paste0("survival::Surv(OS_time, OS_event) ~ ", rhs))
    fit <- tryCatch(survival::coxph(form, data = train_dt), error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    test_dt <- as.data.frame(dev_dt[test_idx, c("OS_time", "OS_event", "cohort", selected_genes), with = FALSE])
    risk <- as.numeric(stats::predict(fit, newdata = test_dt, type = "lp"))
    coef_dt <- data.table::data.table(gene = gsub("`", "", names(stats::coef(fit))), coefficient = as.numeric(stats::coef(fit)))
    coef_dt <- coef_dt[coefficient != 0 & !is.na(coefficient)]
    return(list(risk = risk, coef = coef_dt, model_object = fit, n_nonzero = nrow(coef_dt)))
  }

  if (lc$learner == "StepCoxPH") {
    train_dt <- as.data.frame(dev_dt[train_idx, c("OS_time", "OS_event", "cohort", selected_genes), with = FALSE])
    rhs <- paste(safe_gene_terms(selected_genes), collapse = " + ")
    if (isTRUE(lc$use_strata)) rhs <- paste0(rhs, " + strata(cohort)")
    full_form <- stats::as.formula(paste0("survival::Surv(OS_time, OS_event) ~ ", rhs))
    full_fit <- tryCatch(survival::coxph(full_form, data = train_dt), error = function(e) NULL)
    if (is.null(full_fit)) return(NULL)

    fit <- tryCatch(stats::step(full_fit, direction = ifelse(is.na(lc$step_direction), "both", lc$step_direction), trace = 0), error = function(e) full_fit)
    if (is.null(fit)) return(NULL)

    test_dt <- as.data.frame(dev_dt[test_idx, c("OS_time", "OS_event", "cohort", selected_genes), with = FALSE])
    risk <- as.numeric(stats::predict(fit, newdata = test_dt, type = "lp"))
    coef_dt <- data.table::data.table(gene = gsub("`", "", names(stats::coef(fit))), coefficient = as.numeric(stats::coef(fit)))
    coef_dt <- coef_dt[coefficient != 0 & !is.na(coefficient)]
    return(list(risk = risk, coef = coef_dt, model_object = fit, n_nonzero = nrow(coef_dt)))
  }

  stop("Unknown learner_config: ", learner_config)
}

# Create ML175-inspired candidate library:
# 7 selectors × 5 feature sizes × 7 LOCO-compatible survival learners = 245 configurations.
# This preserves the uploaded script's two-stage combination idea while using
# survival-safe learners and LOCO-C-index evaluation.
candidate_grid <- data.table::CJ(
  selector = CONFIG$feature_selector_grid,
  feature_size = CONFIG$feature_size_grid,
  learner_config = CONFIG$learner_config_grid,
  unique = TRUE
)
candidate_grid[, learner := vapply(learner_config, function(z) parse_learner_config(z)$learner, character(1))]
candidate_grid[, alpha := vapply(learner_config, function(z) parse_learner_config(z)$alpha, numeric(1))]
candidate_grid[, lambda_policy := vapply(learner_config, function(z) parse_learner_config(z)$lambda_policy, character(1))]
candidate_grid[, candidate_id := sprintf("ML%04d", seq_len(.N))]
write_tsv(candidate_grid, file.path(OUT$model, "Step2_ML175_adapted_candidate_library.tsv"))
write_tsv(candidate_grid, file.path(OUT$model, "Step2_model_candidate_library.tsv"))

loco_rows <- list()
loco_coef_rows <- list()
rank_rows <- list()
row_i <- 1L
coef_i <- 1L
rank_i <- 1L

for (ci in seq_len(nrow(candidate_grid))) {
  cand <- candidate_grid[ci]
  if (ci %% 25 == 0 || ci == 1) msg("LOCO candidate ", ci, "/", nrow(candidate_grid), ": ", cand$candidate_id)

  for (holdout in CONFIG$expected_development_cohorts) {
    train_idx <- which(dev_dt$cohort != holdout)
    test_idx <- which(dev_dt$cohort == holdout)
    if (length(test_idx) < 20 || length(unique(dev_dt$cohort[train_idx])) < 2) next

    ranked <- rank_genes_for_training(
      train_idx = train_idx,
      selector = cand$selector,
      genes = model_genes,
      feature_size = cand$feature_size
    )
    selected_genes <- ranked$genes

    # Store only selected-rank rows to avoid massive files.
    rank_tmp <- ranked$ranking_table[gene %in% selected_genes]
    rank_tmp[, candidate_id := cand$candidate_id]
    rank_tmp[, holdout_cohort := holdout]
    rank_tmp[, selector := cand$selector]
    rank_tmp[, feature_size := cand$feature_size]
    rank_rows[[rank_i]] <- rank_tmp
    rank_i <- rank_i + 1L

    fitpred <- fit_predict_candidate(
      train_idx = train_idx,
      test_idx = test_idx,
      selected_genes = selected_genes,
      learner_config = cand$learner_config
    )
    if (is.null(fitpred)) next

    ev <- eval_risk(dev_dt$OS_time[test_idx], dev_dt$OS_event[test_idx], fitpred$risk)
    loco_rows[[row_i]] <- data.table::data.table(
      candidate_id = cand$candidate_id,
      selector = cand$selector,
      feature_size = cand$feature_size,
      learner = cand$learner,
      learner_config = cand$learner_config,
      alpha = cand$alpha,
      lambda_policy = cand$lambda_policy,
      holdout_cohort = holdout,
      selected_gene_set = paste(selected_genes, collapse = ";"),
      n_selected_genes_before_model = length(selected_genes),
      n_nonzero_genes = fitpred$n_nonzero,
      n_test = ev["n"],
      events_test = ev["events"],
      cindex = ev["cindex"],
      beta_per_1SD = ev["beta"],
      HR_per_1SD = ev["HR"],
      p = ev["p"]
    )
    row_i <- row_i + 1L

    if (nrow(fitpred$coef) > 0) {
      coef_dt <- fitpred$coef
      coef_dt[, candidate_id := cand$candidate_id]
      coef_dt[, holdout_cohort := holdout]
      coef_dt[, selector := cand$selector]
      coef_dt[, feature_size := cand$feature_size]
      coef_dt[, learner := cand$learner]
      coef_dt[, learner_config := cand$learner_config]
      coef_dt[, alpha := cand$alpha]
      coef_dt[, lambda_policy := cand$lambda_policy]
      loco_coef_rows[[coef_i]] <- coef_dt
      coef_i <- coef_i + 1L
    }
  }
}

loco_results <- data.table::rbindlist(loco_rows, fill = TRUE)
loco_coef <- data.table::rbindlist(loco_coef_rows, fill = TRUE)
loco_ranked_features <- data.table::rbindlist(rank_rows, fill = TRUE)

write_tsv(loco_results, file.path(OUT$model, "Step2_multiconfig_LOCO_results.tsv"))
write_tsv(loco_coef, file.path(OUT$model, "Step2_multiconfig_LOCO_coefficients_long.tsv"))
write_tsv(loco_ranked_features, file.path(OUT$model, "Step2_multiconfig_LOCO_selected_feature_rankings.tsv"))

if (nrow(loco_results) == 0) stop("Multi-configuration LOCO screening produced no results.")

# v2.3 integrity check: summary grouping fields must exist even for non-glmnet learners.
required_loco_cols <- c(
  "candidate_id", "selector", "feature_size", "learner", "learner_config",
  "alpha", "lambda_policy", "holdout_cohort", "cindex", "HR_per_1SD",
  "n_nonzero_genes", "selected_gene_set"
)
missing_loco_cols <- setdiff(required_loco_cols, names(loco_results))
if (length(missing_loco_cols) > 0) {
  stop("Internal LOCO table schema error. Missing columns: ", paste(missing_loco_cols, collapse = ", "))
}

model_screen_summary <- loco_results[, .(
  n_loco = .N,
  mean_cindex = mean(cindex, na.rm = TRUE),
  sd_cindex = stats::sd(cindex, na.rm = TRUE),
  min_cindex = min(cindex, na.rm = TRUE),
  mean_HR_per_1SD = mean(HR_per_1SD, na.rm = TRUE),
  median_nonzero_genes = median(n_nonzero_genes, na.rm = TRUE),
  mean_nonzero_genes = mean(n_nonzero_genes, na.rm = TRUE),
  selected_gene_sets = paste(unique(selected_gene_set), collapse = " || ")
), by = .(candidate_id, selector, feature_size, learner, learner_config, alpha, lambda_policy)]

model_screen_summary <- model_screen_summary[n_loco == expected_n_cohorts]
best_mean <- max(model_screen_summary$mean_cindex, na.rm = TRUE)
model_screen_summary[, within_best_band := mean_cindex >= best_mean - CONFIG$loco_mean_cindex_tolerance]
model_screen_summary[, min_cindex_pass := min_cindex >= CONFIG$loco_min_cindex_floor]
model_screen_summary[, parsimony_prefer_pass := median_nonzero_genes <= CONFIG$parsimony_prefer_max_locked_genes]
model_screen_summary[, parsimony_hard_pass := median_nonzero_genes <= CONFIG$parsimony_hard_max_locked_genes]
model_screen_summary[, selection_eligible := within_best_band == TRUE &
                       min_cindex_pass == TRUE &
                       median_nonzero_genes >= CONFIG$min_locked_genes &
                       parsimony_hard_pass == TRUE]

# v2.4 selection order:
# 1) within-best-band and worst-cohort-pass models
# 2) hard parsimony pass (<=15)
# 3) preferred parsimony pass (<=12)
# 4) smaller median nonzero genes
# 5) stronger worst-cohort and mean LOCO C-index
model_screen_summary <- model_screen_summary[order(
  -selection_eligible,
  -within_best_band,
  -min_cindex_pass,
  -parsimony_hard_pass,
  -parsimony_prefer_pass,
  median_nonzero_genes,
  -min_cindex,
  -mean_cindex,
  candidate_id
)]

# v2.5: robust selected-config fallback.
# data.table's DT[condition][1] can return a one-row NA table when the filtered
# table is empty. This caused selected_config$selector to become NA in v2.4.
model_screen_summary[, finite_screen := is.finite(mean_cindex) & is.finite(min_cindex) & !is.na(selector) & !is.na(candidate_id)]

selection_tiers <- list(
  tier1_eligible_parsimony_locked = model_screen_summary[
    finite_screen == TRUE & selection_eligible == TRUE
  ],
  tier2_bestband_worstcohort_pass_any_size = model_screen_summary[
    finite_screen == TRUE & within_best_band == TRUE & min_cindex_pass == TRUE
  ][order(parsimony_hard_pass == FALSE, parsimony_prefer_pass == FALSE, median_nonzero_genes, -min_cindex, -mean_cindex, candidate_id)],
  tier3_bestband_any_worstcohort = model_screen_summary[
    finite_screen == TRUE & within_best_band == TRUE
  ][order(parsimony_hard_pass == FALSE, parsimony_prefer_pass == FALSE, median_nonzero_genes, -min_cindex, -mean_cindex, candidate_id)],
  tier4_best_finite_model = model_screen_summary[
    finite_screen == TRUE
  ][order(parsimony_hard_pass == FALSE, parsimony_prefer_pass == FALSE, median_nonzero_genes, -min_cindex, -mean_cindex, candidate_id)]
)

selected_config <- data.table::data.table()
selection_tier_used <- NA_character_
for (nm in names(selection_tiers)) {
  cand_tier <- selection_tiers[[nm]]
  if (nrow(cand_tier) > 0) {
    selected_config <- cand_tier[1]
    selection_tier_used <- nm
    break
  }
}

if (nrow(selected_config) == 0 || is.na(selected_config$candidate_id[1]) || is.na(selected_config$selector[1])) {
  write_tsv(model_screen_summary, file.path(OUT$model, "Step2_multiconfig_model_screen_summary.tsv"))
  stop(
    "No valid selected model configuration could be identified. ",
    "Inspect Step2_multiconfig_model_screen_summary.tsv and Step2_multiconfig_LOCO_results.tsv."
  )
}

selected_config[, selection_tier_used := selection_tier_used]

# Reorder final summary after adding diagnostics.
model_screen_summary <- model_screen_summary[order(
  -selection_eligible,
  -within_best_band,
  -min_cindex_pass,
  -parsimony_hard_pass,
  -parsimony_prefer_pass,
  median_nonzero_genes,
  -min_cindex,
  -mean_cindex,
  candidate_id
)]
write_tsv(model_screen_summary, file.path(OUT$model, "Step2_multiconfig_model_screen_summary.tsv"))
write_tsv(selected_config, file.path(OUT$model, "Step2_selected_model_configuration.tsv"))

selection_diagnostic <- data.table::data.table(
  tier = names(selection_tiers),
  n_candidates = vapply(selection_tiers, nrow, integer(1))
)
selection_diagnostic[, selected := tier == selection_tier_used]
write_tsv(selection_diagnostic, file.path(OUT$model, "Step2_selected_model_configuration_tier_diagnostic.tsv"))

msg("Selected model configuration: ", selected_config$candidate_id,
    "; selection_tier=", selection_tier_used,
    "; selector=", selected_config$selector,
    "; feature_size=", selected_config$feature_size,
    "; learner=", selected_config$learner,
    "; learner_config=", selected_config$learner_config,
    "; alpha=", selected_config$alpha,
    "; lambda=", selected_config$lambda_policy,
    "; median_nonzero=", selected_config$median_nonzero_genes)

###############################################################################
# 5. Refit selected configuration on all development cohorts and bootstrap
###############################################################################

msg("Refitting selected configuration on all development data...")

all_idx <- seq_len(nrow(dev_dt))
final_ranked <- rank_genes_for_training(
  train_idx = all_idx,
  selector = selected_config$selector,
  genes = model_genes,
  feature_size = selected_config$feature_size
)
final_pre_model_genes <- final_ranked$genes
write_tsv(final_ranked$ranking_table, file.path(OUT$model, "Step2_final_selected_configuration_gene_ranking.tsv"))

final_fit <- fit_predict_candidate(
  train_idx = all_idx,
  test_idx = all_idx,
  selected_genes = final_pre_model_genes,
  learner_config = selected_config$learner_config
)
if (is.null(final_fit) || nrow(final_fit$coef) < CONFIG$min_locked_genes) {
  stop("Final selected model refit failed or selected too few genes.")
}

final_coef <- final_fit$coef
data.table::setnames(final_coef, "coefficient", "locked_coefficient")
final_coef <- final_coef[order(-abs(locked_coefficient), gene)]
final_coef[, locked_rank := seq_len(.N)]
write_tsv(final_coef, file.path(OUT$model, "Step2_final_selected_model_coefficients.tsv"))

msg("Running bootstrap stability for the selected configuration: n=", CONFIG$bootstrap_n)

bootstrap_rows <- list()
boot_i <- 1L

for (b in seq_len(CONFIG$bootstrap_n)) {
  boot_idx <- integer(0)
  for (coh in CONFIG$expected_development_cohorts) {
    idx <- which(dev_dt$cohort == coh)
    boot_idx <- c(boot_idx, sample(idx, size = length(idx), replace = TRUE))
  }

  ranked_b <- rank_genes_for_training(
    train_idx = boot_idx,
    selector = selected_config$selector,
    genes = model_genes,
    feature_size = selected_config$feature_size
  )
  fit_b <- fit_predict_candidate(
    train_idx = boot_idx,
    test_idx = boot_idx,
    selected_genes = ranked_b$genes,
    learner_config = selected_config$learner_config
  )
  if (!is.null(fit_b) && nrow(fit_b$coef) > 0) {
    coef_b <- fit_b$coef
    coef_b[, bootstrap_id := b]
    coef_b[, pre_model_gene_set := paste(ranked_b$genes, collapse = ";")]
    bootstrap_rows[[boot_i]] <- coef_b
    boot_i <- boot_i + 1L
  }

  if (b %% 20 == 0) msg("Bootstrap completed: ", b, "/", CONFIG$bootstrap_n)
}

bootstrap_long <- data.table::rbindlist(bootstrap_rows, fill = TRUE)
write_tsv(bootstrap_long, file.path(OUT$model, "Step2_selected_config_bootstrap_coefficients_long.tsv"))

if (nrow(bootstrap_long) > 0) {
  bootstrap_stability <- bootstrap_long[, .(
    bootstrap_selected_n = .N,
    selection_frequency = .N / CONFIG$bootstrap_n,
    n_bootstrap_coef_positive = sum(coefficient > 0, na.rm = TRUE),
    n_bootstrap_coef_negative = sum(coefficient < 0, na.rm = TRUE),
    n_bootstrap_coef_zero = sum(coefficient == 0, na.rm = TRUE),
    bootstrap_positive_fraction = mean(coefficient > 0, na.rm = TRUE),
    bootstrap_negative_fraction = mean(coefficient < 0, na.rm = TRUE),
    bootstrap_sign_consistency = max(
      mean(coefficient > 0, na.rm = TRUE),
      mean(coefficient < 0, na.rm = TRUE)
    ),
    bootstrap_consensus_sign = data.table::fcase(
      mean(coefficient > 0, na.rm = TRUE) >= mean(coefficient < 0, na.rm = TRUE), 1L,
      mean(coefficient < 0, na.rm = TRUE) > mean(coefficient > 0, na.rm = TRUE), -1L,
      default = NA_integer_
    ),
    median_coefficient = median(coefficient, na.rm = TRUE),
    mean_coefficient = mean(coefficient, na.rm = TRUE),
    median_abs_coefficient = median(abs(coefficient), na.rm = TRUE)
  ), by = gene]
} else {
  bootstrap_stability <- data.table::data.table()
}

bootstrap_stability <- merge(data.table::data.table(gene = model_genes), bootstrap_stability, by = "gene", all.x = TRUE)
bootstrap_stability[is.na(bootstrap_selected_n), bootstrap_selected_n := 0]
bootstrap_stability[is.na(selection_frequency), selection_frequency := 0]
bootstrap_stability <- merge(
  bootstrap_stability,
  final_coef[, .(gene, locked_coefficient, locked_rank)],
  by = "gene",
  all.x = TRUE
)
bootstrap_stability[, final_locked_gene := !is.na(locked_coefficient)]
bootstrap_stability[is.na(locked_coefficient), locked_coefficient := 0]
bootstrap_stability <- bootstrap_stability[order(-final_locked_gene, locked_rank, -selection_frequency, gene)]
write_tsv(bootstrap_stability, file.path(OUT$model, "Step2_selected_config_bootstrap_gene_stability_summary.tsv"))

###############################################################################
# 6. Direction-stable parsimony-locked model, risk score and performance
###############################################################################

msg("Locking final direction-stable parsimony-constrained model...")

# Provisional final model from the selected ML configuration.
provisional_coef <- data.table::copy(final_coef)
data.table::setnames(provisional_coef, "locked_rank", "provisional_rank")
data.table::setnames(provisional_coef, "locked_coefficient", "provisional_coefficient")

parsimony_audit <- merge(
  provisional_coef[, .(gene, provisional_rank, provisional_coefficient)],
  bootstrap_stability[, .(
    gene, selection_frequency, bootstrap_selected_n, median_coefficient, mean_coefficient,
    bootstrap_sign_consistency, bootstrap_consensus_sign,
    bootstrap_positive_fraction, bootstrap_negative_fraction
  )],
  by = "gene",
  all.x = TRUE
)
parsimony_audit[is.na(selection_frequency), selection_frequency := 0]
parsimony_audit[is.na(bootstrap_sign_consistency), bootstrap_sign_consistency := 0]
parsimony_audit[, abs_provisional_coefficient := abs(provisional_coefficient)]
parsimony_audit[, bootstrap_pass := selection_frequency >= CONFIG$bootstrap_lock_min_frequency]
parsimony_audit[, bootstrap_sign_pass := bootstrap_sign_consistency >= CONFIG$bootstrap_sign_consistency_min]

parsimony_audit <- merge(
  parsimony_audit,
  direction_summary[, .(
    gene, n_valid_cohorts, n_beta_positive, n_beta_negative,
    univariable_consensus_sign, no_mixed_univariable_cox_sign,
    step1_sign_compatible, direction_stable_for_modeling,
    direction_exclusion_reason, beta_sign_pattern
  )],
  by = "gene",
  all.x = TRUE
)
parsimony_audit <- merge(
  parsimony_audit,
  classA[, .(gene, clinical_model_class, Step2_feature_sign, bulk_direction, mean_cox_beta, median_cox_beta)],
  by = "gene",
  all.x = TRUE
)
parsimony_audit[, provisional_sign := data.table::fcase(
  provisional_coefficient > 0, 1L,
  provisional_coefficient < 0, -1L,
  default = NA_integer_
)]
parsimony_audit[, provisional_direction_consistent := (
  !is.na(provisional_sign) &
    !is.na(univariable_consensus_sign) &
    provisional_sign == univariable_consensus_sign
)]
parsimony_audit[, bootstrap_direction_consistent := (
  !is.na(bootstrap_consensus_sign) &
    !is.na(univariable_consensus_sign) &
    bootstrap_consensus_sign == univariable_consensus_sign
)]
parsimony_audit[, direction_lock_pass := (
  direction_stable_for_modeling == TRUE &
    bootstrap_pass == TRUE &
    bootstrap_sign_pass == TRUE &
    provisional_direction_consistent == TRUE &
    bootstrap_direction_consistent == TRUE
)]

parsimony_audit <- parsimony_audit[order(
  -direction_lock_pass,
  -bootstrap_pass,
  -bootstrap_sign_pass,
  -selection_frequency,
  provisional_rank,
  -abs_provisional_coefficient,
  gene
)]

selected_lock_genes <- parsimony_audit[direction_lock_pass == TRUE, gene]
parsimony_lock_reason <- paste0(
  "direction_stable_univariable_cox_plus_bootstrap_frequency_ge_",
  CONFIG$bootstrap_lock_min_frequency,
  "_plus_bootstrap_sign_consistency_ge_",
  CONFIG$bootstrap_sign_consistency_min
)

if (length(selected_lock_genes) > CONFIG$parsimony_prefer_max_locked_genes) {
  selected_lock_genes <- selected_lock_genes[seq_len(CONFIG$parsimony_prefer_max_locked_genes)]
  parsimony_lock_reason <- paste0(parsimony_lock_reason, "_then_truncated_to_prefer_max_", CONFIG$parsimony_prefer_max_locked_genes)
}
if (length(selected_lock_genes) > CONFIG$parsimony_hard_max_locked_genes) {
  selected_lock_genes <- selected_lock_genes[seq_len(CONFIG$parsimony_hard_max_locked_genes)]
  parsimony_lock_reason <- paste0(parsimony_lock_reason, "_then_truncated_to_hard_max_", CONFIG$parsimony_hard_max_locked_genes)
}
if (length(selected_lock_genes) < CONFIG$min_locked_genes) {
  write_tsv(parsimony_audit, file.path(OUT$model, "Step2_direction_stable_parsimony_lock_gene_audit.tsv"))
  stop(
    "Too few genes passed the final v2.6 direction-stable lock: ",
    length(selected_lock_genes),
    ". Inspect Step2_direction_stable_parsimony_lock_gene_audit.tsv. ",
    "No direction-unstable fallback is allowed in v2.6."
  )
}

refit_direction_stable_cox <- function(lock_genes) {
  current <- unique(lock_genes)
  removal_log <- list()
  rr <- 1L

  repeat {
    if (length(current) < CONFIG$min_locked_genes) {
      stop("Too few genes remain during direction-stable refit: ", length(current))
    }

    fit_dt <- as.data.frame(dev_dt[, c("OS_time", "OS_event", current), with = FALSE])
    form <- stats::as.formula(
      paste0("survival::Surv(OS_time, OS_event) ~ ", paste(safe_gene_terms(current), collapse = " + "))
    )
    fit <- tryCatch(survival::coxph(form, data = fit_dt), error = function(e) {
      msg("Direction-stable Cox refit failed: ", conditionMessage(e))
      NULL
    })
    if (is.null(fit)) stop("Direction-stable Cox refit failed.")

    coef_dt <- data.table::data.table(
      gene = gsub("`", "", names(stats::coef(fit))),
      locked_coefficient = as.numeric(stats::coef(fit))
    )
    coef_dt <- merge(
      coef_dt,
      direction_summary[, .(gene, univariable_consensus_sign)],
      by = "gene",
      all.x = TRUE
    )
    coef_dt[, locked_coefficient_sign := data.table::fcase(
      locked_coefficient > 0, 1L,
      locked_coefficient < 0, -1L,
      default = NA_integer_
    )]
    coef_dt[, refit_direction_consistent := (
      !is.na(locked_coefficient_sign) &
        !is.na(univariable_consensus_sign) &
        locked_coefficient_sign == univariable_consensus_sign
    )]

    bad <- coef_dt[refit_direction_consistent != TRUE | is.na(refit_direction_consistent), gene]
    if (length(bad) == 0) {
      return(list(fit = fit, coef = coef_dt, removal_log = data.table::rbindlist(removal_log, fill = TRUE)))
    }

    removal_log[[rr]] <- data.table::data.table(
      refit_round = rr,
      removed_gene = bad,
      removal_reason = "multivariable_refit_coefficient_sign_inconsistent_with_per_cohort_univariable_consensus"
    )
    rr <- rr + 1L
    current <- setdiff(current, bad)

    if (!isTRUE(CONFIG$final_refit_require_direction_consistency)) {
      return(list(fit = fit, coef = coef_dt, removal_log = data.table::rbindlist(removal_log, fill = TRUE)))
    }
  }
}

refit_obj <- refit_direction_stable_cox(selected_lock_genes)
compact_coef <- refit_obj$coef[, .(gene, locked_coefficient, locked_coefficient_sign, refit_direction_consistent)]
if (!is.null(refit_obj$removal_log) && nrow(refit_obj$removal_log) > 0) {
  write_tsv(refit_obj$removal_log, file.path(OUT$model, "Step2_v2_6_direction_refit_removed_genes.tsv"))
} else {
  write_tsv(data.table::data.table(), file.path(OUT$model, "Step2_v2_6_direction_refit_removed_genes.tsv"))
}

locked_model_type <- "direction_stable_parsimony_refit_multivariable_cox_no_strata"
selected_lock_genes <- compact_coef$gene

parsimony_audit[, parsimony_locked_pre_refit := gene %in% selected_lock_genes]
parsimony_audit <- merge(
  parsimony_audit,
  compact_coef[, .(gene, locked_coefficient, locked_coefficient_sign, refit_direction_consistent)],
  by = "gene",
  all.x = TRUE
)
parsimony_audit[, parsimony_locked := !is.na(locked_coefficient)]
parsimony_audit[, parsimony_rank := ifelse(parsimony_locked, match(gene, selected_lock_genes), NA_integer_)]
write_tsv(parsimony_audit, file.path(OUT$model, "Step2_direction_stable_parsimony_lock_gene_audit.tsv"))
# Back-compatible file name.
write_tsv(parsimony_audit, file.path(OUT$model, "Step2_parsimony_lock_gene_audit.tsv"))

compact_coef <- merge(
  data.table::data.table(gene = selected_lock_genes, model_rank = seq_along(selected_lock_genes)),
  compact_coef,
  by = "gene",
  all.x = TRUE
)
compact_coef <- merge(
  compact_coef,
  parsimony_audit[, .(
    gene, provisional_rank, provisional_coefficient, selection_frequency,
    bootstrap_selected_n, bootstrap_pass, bootstrap_sign_consistency, bootstrap_sign_pass,
    bootstrap_consensus_sign, n_valid_cohorts, n_beta_positive, n_beta_negative,
    univariable_consensus_sign, beta_sign_pattern, clinical_model_class, Step2_feature_sign,
    bulk_direction, mean_cox_beta, median_cox_beta
  )],
  by = "gene",
  all.x = TRUE
)
compact_coef <- compact_coef[order(model_rank)]
write_tsv(compact_coef, file.path(OUT$final, "SIDISH_Step2_locked_model_coefficients.tsv"))

score_mat <- as.matrix(dev_dt[, compact_coef$gene, with = FALSE])
score_mat[!is.finite(score_mat)] <- 0
coef_vec <- compact_coef$locked_coefficient
names(coef_vec) <- compact_coef$gene
risk_score <- as.numeric(score_mat[, names(coef_vec), drop = FALSE] %*% coef_vec)

risk_dt <- data.table::data.table(
  sample_id = dev_dt$sample_id,
  cohort = dev_dt$cohort,
  OS_time = dev_dt$OS_time,
  OS_event = dev_dt$OS_event,
  locked_risk_score = risk_score
)

if (CONFIG$development_cutoff_method == "median") {
  locked_cutoff <- median(risk_dt$locked_risk_score, na.rm = TRUE)
} else {
  stop("Unsupported cutoff method: ", CONFIG$development_cutoff_method)
}
risk_dt[, locked_risk_group := ifelse(locked_risk_score >= locked_cutoff, "High", "Low")]
write_tsv(risk_dt, file.path(OUT$final, "SIDISH_Step2_locked_development_risk_scores.tsv"))

perf_rows <- list()
i <- 1L
for (coh in c("ALL", CONFIG$expected_development_cohorts)) {
  dd <- if (coh == "ALL") risk_dt else risk_dt[cohort == coh]
  ev <- eval_risk(dd$OS_time, dd$OS_event, dd$locked_risk_score)
  perf_rows[[i]] <- data.table::data.table(
    cohort = coh,
    n = ev["n"],
    events = ev["events"],
    cindex = ev["cindex"],
    beta_per_1SD = ev["beta"],
    HR_per_1SD = ev["HR"],
    p = ev["p"]
  )
  i <- i + 1L
}
perf <- data.table::rbindlist(perf_rows)
write_tsv(perf, file.path(OUT$perf, "SIDISH_Step2_locked_model_performance_by_cohort.tsv"))

# Objective C-index table: raw metrics only, no qualitative labels.
objective_cindex_report <- perf[, .(
  cohort,
  n,
  events,
  cindex,
  HR_per_1SD,
  p
)]
objective_cindex_report[, reporting_note := "raw_C_index_values_only_no_qualitative_label"]
write_tsv(objective_cindex_report, file.path(OUT$perf, "SIDISH_Step2_objective_Cindex_report.tsv"))

hl_rows <- list()
i <- 1L
for (coh in c("ALL", CONFIG$expected_development_cohorts)) {
  dd <- if (coh == "ALL") data.table::copy(risk_dt) else data.table::copy(risk_dt[cohort == coh])
  dd[, group01 := ifelse(locked_risk_group == "High", 1, 0)]
  fit <- tryCatch(survival::coxph(survival::Surv(OS_time, OS_event) ~ group01, data = dd), error = function(e) NULL)
  if (!is.null(fit)) {
    sm <- summary(fit)
    hl_rows[[i]] <- data.table::data.table(
      cohort = coh,
      n = nrow(dd),
      events = sum(dd$OS_event == 1),
      cutoff = locked_cutoff,
      high_n = sum(dd$locked_risk_group == "High"),
      low_n = sum(dd$locked_risk_group == "Low"),
      HR_high_vs_low = exp(as.numeric(stats::coef(fit)[1])),
      p = sm$coefficients[1, "Pr(>|z|)"]
    )
    i <- i + 1L
  }
}
hl_perf <- data.table::rbindlist(hl_rows, fill = TRUE)
write_tsv(hl_perf, file.path(OUT$perf, "SIDISH_Step2_locked_model_high_low_performance.tsv"))

locked_summary <- data.table::data.table(
  metric = c(
    "n_input_ClassA_genes",
    "n_modeling_genes_present_in_all_development_cohorts_before_direction_gate",
    "n_direction_stable_modeling_genes",
    "n_direction_unstable_excluded_genes",
    "n_candidate_model_configurations",
    "selected_candidate_id",
    "selected_config_selection_tier",
    "selected_selector",
    "selected_feature_size_before_model",
    "selected_learner",
    "selected_learner_config",
    "selected_alpha",
    "selected_lambda_policy",
    "selected_LOCO_mean_cindex",
    "selected_LOCO_min_cindex",
    "n_provisional_final_model_genes",
    "bootstrap_lock_min_frequency",
    "bootstrap_sign_consistency_min",
    "parsimony_prefer_max_locked_genes",
    "parsimony_hard_max_locked_genes",
    "parsimony_lock_reason",
    "locked_model_type",
    "n_locked_genes",
    "locked_cutoff_method",
    "locked_development_cutoff",
    "development_n",
    "development_events",
    "development_cindex_ALL",
    "development_HR_per_1SD_ALL",
    "development_p_ALL",
    "cindex_reporting_policy"
  ),
  value = c(
    length(candidate_genes_all),
    length(model_genes_before_direction_gate),
    length(model_genes),
    length(setdiff(model_genes_before_direction_gate, model_genes)),
    nrow(candidate_grid),
    selected_config$candidate_id,
    selection_tier_used,
    selected_config$selector,
    selected_config$feature_size,
    selected_config$learner,
    selected_config$learner_config,
    selected_config$alpha,
    selected_config$lambda_policy,
    selected_config$mean_cindex,
    selected_config$min_cindex,
    nrow(provisional_coef),
    CONFIG$bootstrap_lock_min_frequency,
    CONFIG$bootstrap_sign_consistency_min,
    CONFIG$parsimony_prefer_max_locked_genes,
    CONFIG$parsimony_hard_max_locked_genes,
    parsimony_lock_reason,
    locked_model_type,
    nrow(compact_coef),
    CONFIG$development_cutoff_method,
    locked_cutoff,
    nrow(risk_dt),
    sum(risk_dt$OS_event == 1),
    perf[cohort == "ALL", cindex],
    perf[cohort == "ALL", HR_per_1SD],
    perf[cohort == "ALL", p],
    "raw_C_index_values_only_no_qualitative_label"
  )
)
write_tsv(locked_summary, file.path(OUT$final, "SIDISH_Step2_locked_model_summary.tsv"))

###############################################################################
# 7. Publication-style figures matching Step 1
###############################################################################

if (isTRUE(CONFIG$make_publication_figures)) {
  msg("Generating Step 2 v2.6 publication-style figures...")

  # Figure 0: Cox direction-stability gate.
  fig0_dt <- direction_summary[, .N, by = direction_exclusion_reason][order(-N)]
  fig0_dt[, reason_label := gsub("_", " ", direction_exclusion_reason)]
  fig0_cols <- c(
    kept_direction_stable = pal("green"),
    excluded_mixed_per_cohort_univariable_cox_sign = pal("red"),
    excluded_step1_sign_incompatible = pal("orange"),
    excluded_insufficient_valid_cohort_cox = pal("grey"),
    excluded_other_direction_instability = pal("purple")
  )
  p0 <- ggplot2::ggplot(fig0_dt, ggplot2::aes(x = stats::reorder(reason_label, N), y = N, fill = direction_exclusion_reason)) +
    ggplot2::geom_col(width = 0.70, color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_text(ggplot2::aes(label = N), hjust = -0.12, family = CONFIG$plot_base_family, size = 2.6) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::expand_limits(y = max(fig0_dt$N) * 1.15) +
    ggplot2::scale_fill_manual(values = fig0_cols, guide = "none") +
    ggplot2::labs(title = "Step 2 | Cox direction-stability gate",
                  subtitle = "Mixed per-cohort Cox directions are excluded before ML screening",
                  x = NULL, y = "Number of Class A genes") +
    pub_theme()
  save_pub_plot(p0, "Figure_0_Cox_direction_stability_gate", width = 5.0, height = 3.2)

  # Figure 1: model library overview.
  fig1_dt <- data.table::copy(model_screen_summary)
  fig1_dt[, learner_label := ifelse(learner == "glmnet", paste0("EN ", alpha, " / ", lambda_policy), "Cox-strata")]
  p1 <- ggplot2::ggplot(fig1_dt, ggplot2::aes(x = median_nonzero_genes, y = mean_cindex, color = learner, shape = selection_eligible)) +
    ggplot2::geom_hline(yintercept = max(fig1_dt$mean_cindex, na.rm = TRUE) - CONFIG$loco_mean_cindex_tolerance,
                        linetype = "dashed", linewidth = CONFIG$plot_linewidth, color = pal("grey")) +
    ggplot2::geom_point(size = 2.0, alpha = 0.75) +
    ggplot2::geom_point(data = fig1_dt[candidate_id == selected_config$candidate_id],
                        size = 3.3, shape = 21, fill = pal("red"), color = "black", stroke = 0.5) +
    ggplot2::scale_color_manual(values = c(
      CoxRidge = pal("navy"),
      CoxEnet = pal("blue"),
      CoxLasso = pal("red"),
      CoxPH = pal("green"),
      StepCoxPH = pal("orange")
    )) +
    ggplot2::labs(title = "Step 2 | Direction-stable model screening",
                  subtitle = paste0(nrow(candidate_grid), " ML175-adapted survival configurations; red point = locked configuration"),
                  x = "Median non-zero genes across LOCO folds", y = "Mean LOCO C-index",
                  color = "Learner", shape = "Eligible") +
    pub_theme()
  save_pub_plot(p1, "Figure_1_direction_stable_multiconfig_model_screening", width = 5.2, height = 3.8)

  # Figure 2: selected configuration LOCO performance.
  fig2_dt <- loco_results[candidate_id == selected_config$candidate_id]
  p2 <- ggplot2::ggplot(fig2_dt, ggplot2::aes(x = holdout_cohort, y = cindex, fill = holdout_cohort)) +
    ggplot2::geom_col(width = 0.70, color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = CONFIG$plot_linewidth, color = "black") +
    ggplot2::scale_fill_manual(values = c(TCGA = pal("red"), CGGA325 = pal("blue"), CGGA693 = pal("green")), guide = "none") +
    ggplot2::coord_cartesian(ylim = c(0.45, max(0.78, max(fig2_dt$cindex, na.rm = TRUE) + 0.03))) +
    ggplot2::labs(title = "Step 2 | Selected model LOCO performance",
                  subtitle = paste0(selected_config$selector, " | ", selected_config$learner, " | final size selected by data"),
                  x = NULL, y = "Held-out C-index") +
    pub_theme()
  save_pub_plot(p2, "Figure_2_selected_model_LOCO_performance", width = 4.4, height = 3.2)

  # Figure 3: final locked coefficients.
  fig3_dt <- compact_coef[order(locked_coefficient)]
  fig3_dt[, gene := factor(gene, levels = gene)]
  fig3_dt[, direction := ifelse(locked_coefficient >= 0, "Risk-positive", "Protective/negative")]
  p3 <- ggplot2::ggplot(fig3_dt, ggplot2::aes(x = locked_coefficient, y = gene, color = direction)) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = locked_coefficient, yend = gene), linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_point(size = CONFIG$plot_point_size + 0.25) +
    ggplot2::scale_color_manual(values = c("Risk-positive" = pal("red"), "Protective/negative" = pal("blue"))) +
    ggplot2::labs(title = "Step 2 | Direction-stable locked model coefficients",
                  subtitle = paste0(nrow(compact_coef), "-gene model after direction-stability and bootstrap-parsimony locking"),
                  x = "Locked coefficient", y = NULL, color = "Direction") +
    pub_theme()
  save_pub_plot(p3, "Figure_3_locked_model_coefficients", width = 4.8, height = max(3.4, 0.18 * nrow(fig3_dt) + 1.3))

  # Figure 4: bootstrap stability.
  fig4_dt <- bootstrap_stability[order(-final_locked_gene, -selection_frequency)][1:min(25, .N)]
  fig4_dt[, gene := factor(gene, levels = rev(gene))]
  p4 <- ggplot2::ggplot(fig4_dt, ggplot2::aes(x = selection_frequency, y = gene, fill = final_locked_gene)) +
    ggplot2::geom_col(width = 0.70, color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = pal("red"), `FALSE` = pal("grey"))) +
    ggplot2::labs(title = "Step 2 | Bootstrap stability of the selected configuration",
                  subtitle = paste0(CONFIG$bootstrap_n, " cohort-stratified bootstraps"),
                  x = "Selection frequency", y = NULL, fill = "Locked") +
    pub_theme()
  save_pub_plot(p4, "Figure_4_selected_config_bootstrap_stability", width = 4.8, height = 4.4)

  # Figure 5: Class A to locked model mapping.
  fig5_dt <- merge(
    classA[, .(gene, clinical_model_class, Step2_feature_sign)],
    compact_coef[, .(gene, locked_coefficient)],
    by = "gene",
    all.x = TRUE
  )
  fig5_dt[, locked_model_status := ifelse(!is.na(locked_coefficient), "Selected", "Not selected")]
  fig5_dt[, gene := factor(gene, levels = rev(fig5_dt$gene))]
  p5 <- ggplot2::ggplot(fig5_dt, ggplot2::aes(x = locked_model_status, y = gene, fill = Step2_feature_sign)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.20, width = 0.88, height = 0.88) +
    ggplot2::scale_fill_manual(values = c(
      risk_positive = pal("red"),
      protective_or_subtype_negative = pal("blue"),
      signed_context_dependent = pal("orange"),
      signed_context_dependent_borderline = pal("purple")
    )) +
    ggplot2::labs(title = "Step 2 | Class A pool to direction-stable locked model",
                  subtitle = "Final gene count is constrained by Cox direction stability, bootstrap stability and parsimony",
                  x = NULL, y = NULL, fill = "Step 1 sign") +
    pub_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  save_pub_plot(p5, "Figure_5_ClassA_to_locked_model_mapping", width = 5.2, height = max(4.4, 0.12 * nrow(fig5_dt) + 1.7))

  # Figure 6: risk-score distribution.
  p6 <- ggplot2::ggplot(risk_dt, ggplot2::aes(x = cohort, y = locked_risk_score, fill = cohort)) +
    ggplot2::geom_violin(width = 0.78, alpha = 0.75, color = "black", linewidth = CONFIG$plot_linewidth, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.16, outlier.size = 0.4, linewidth = CONFIG$plot_linewidth, color = "black", fill = "white") +
    ggplot2::geom_hline(yintercept = locked_cutoff, linetype = "dashed", color = pal("red"), linewidth = CONFIG$plot_linewidth) +
    ggplot2::scale_fill_manual(values = c(TCGA = pal("red"), CGGA325 = pal("blue"), CGGA693 = pal("green")), guide = "none") +
    ggplot2::labs(title = "Step 2 | Locked risk-score distribution",
                  subtitle = "Dashed line: development median cutoff",
                  x = NULL, y = "Locked risk score") +
    pub_theme()
  save_pub_plot(p6, "Figure_6_locked_risk_score_distribution_by_cohort", width = 4.7, height = 3.4)

  # Figure 7: cohort-wise performance.
  fig7_dt <- perf[cohort != "ALL"]
  p7 <- ggplot2::ggplot(fig7_dt, ggplot2::aes(x = cohort, y = cindex, fill = cohort)) +
    ggplot2::geom_col(width = 0.70, color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::scale_fill_manual(values = c(TCGA = pal("red"), CGGA325 = pal("blue"), CGGA693 = pal("green")), guide = "none") +
    ggplot2::coord_cartesian(ylim = c(0.45, max(0.78, max(fig7_dt$cindex, na.rm = TRUE) + 0.03))) +
    ggplot2::labs(title = "Step 2 | Locked model performance",
                  subtitle = "C-index within development cohorts",
                  x = NULL, y = "C-index") +
    pub_theme()
  save_pub_plot(p7, "Figure_7_locked_model_cindex_by_cohort", width = 4.4, height = 3.2)

  # KM curves: ALL and each cohort.
  km_sets <- c("ALL", CONFIG$expected_development_cohorts)
  km_stats <- list()
  ks_i <- 1L
  for (coh in km_sets) {
    dd <- if (coh == "ALL") data.table::copy(risk_dt) else data.table::copy(risk_dt[cohort == coh])
    dd <- dd[!is.na(locked_risk_group)]
    if (nrow(dd) < 20 || data.table::uniqueN(dd$locked_risk_group) < 2) next

    km_dt <- km_summary_table(dd$OS_time, dd$OS_event, dd$locked_risk_group)
    if (is.null(km_dt) || nrow(km_dt) == 0) next
    p_lr <- km_logrank_p(dd$OS_time, dd$OS_event, dd$locked_risk_group)
    hr_row <- hl_perf[cohort == coh]
    hr_text <- if (nrow(hr_row) > 0 && is.finite(hr_row$HR_high_vs_low[1])) paste0("HR = ", sprintf("%.2f", hr_row$HR_high_vs_low[1])) else "HR = NA"
    p_text <- paste0("log-rank ", format_p_for_plot(p_lr))

    km_dt[, strata := factor(strata, levels = c("Low", "High"))]
    pkm <- ggplot2::ggplot(km_dt, ggplot2::aes(x = time, y = surv, color = strata)) +
      ggplot2::geom_step(linewidth = 0.7) +
      ggplot2::scale_color_manual(values = c(Low = CONFIG$km_palette_low, High = CONFIG$km_palette_high)) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = paste0("Step 2 | Kaplan–Meier survival: ", coh),
        subtitle = paste0("Locked cutoff | ", hr_text, " | ", p_text),
        x = CONFIG$km_time_unit_label,
        y = "Survival probability",
        color = "Risk group"
      ) +
      pub_theme()
    save_pub_plot(pkm, paste0("Figure_KM_", coh), width = 4.6, height = 3.6)

    km_stats[[ks_i]] <- data.table::data.table(
      cohort = coh,
      n = nrow(dd),
      events = sum(dd$OS_event == 1),
      high_n = sum(dd$locked_risk_group == "High"),
      low_n = sum(dd$locked_risk_group == "Low"),
      logrank_p = p_lr,
      HR_high_vs_low = if (nrow(hr_row) > 0) hr_row$HR_high_vs_low[1] else NA_real_
    )
    ks_i <- ks_i + 1L
  }
  if (length(km_stats) > 0) {
    write_tsv(data.table::rbindlist(km_stats, fill = TRUE), file.path(OUT$perf, "SIDISH_Step2_KM_curve_statistics.tsv"))
  }
}

###############################################################################
# 9. README and session info
###############################################################################

readme <- c(
  "SIDISH Step 2 ML175-adapted direction-stable model development v2.6 strict-all3",
  "",
  "Input: Step 1 v3.5 Class A signed core-marker pool. The modelable subset is restricted to genes present in all three development cohorts.",
  "Development cohorts: TCGA + CGGA325 + CGGA693.
Strict v1.3 rule: TCGA, CGGA325 and CGGA693 must all be status==ok, and a gene can enter model development only if measurable in all three cohorts. No missing-gene imputation is allowed. TCGA clinical files are restricted to TCGA-specific clinical tables.",
  "Held-out rule: GSE16011/GSE6011 are not used in Step 2.
This v2.0 script fully re-runs Step 2 after fixing the cohort-file matching root cause. It screens ML175-adapted survival configurations across feature selectors, feature sizes, learners and elastic-net settings. The locked model size is selected by LOCO performance and parsimony rather than being fixed at 10 genes. It exports a Step 1-consistent publication-style figure set.",
  "",
  "Key final files:",
  file.path(OUT$final, "SIDISH_Step2_locked_model_coefficients.tsv"),
  file.path(OUT$final, "SIDISH_Step2_locked_development_risk_scores.tsv"),
  file.path(OUT$final, "SIDISH_Step2_locked_model_summary.tsv"),
  "",
  "Next step:",
  "Use locked coefficients and locked development cutoff in Step 3 external no-refit validation."
)
write_text(readme, file.path(OUT$root, "README_SIDISH_Step2_ML175_DirectionStable_Model_Development_v2_6_strict_all3.txt"))
write_text(capture.output(sessionInfo()), file.path(OUT$audit, "sessionInfo.txt"))

msg("SIDISH Step 2 v2.6 ML175-adapted direction-stable strict-all3 model development complete.")
msg("Locked model coefficients: ", file.path(OUT$final, "SIDISH_Step2_locked_compact_model_coefficients.tsv"))
