#!/usr/bin/env Rscript
###############################################################################
# SIDISH Step 1 Core Marker Refinement v3.5
# Fresh reanalysis WITHOUT h5ad parsing
#
# Why this rewrite
# ----------------
# Previous versions repeatedly failed because R tried to parse .h5ad directly.
# The uploaded Python workflow shows that the stable authorstyle metadata workflow
# is built from flat files:
#
#   inputs/scRNA_raw.mtx
#   inputs/scRNA_genes.tsv
#   inputs/scRNA_cells.tsv
#   inputs/scRNA_metadata.tsv
#   runs/os/exports/cell_scores.tsv
#   runs/os/marker/hr_markers.tsv
#   runs/os/marker/celltype_hr_markers_all_threshold.tsv
#   prepared/bulk_expression_dedup.csv + prepared/clinical_aligned.tsv
#
# Therefore this R script reads these explicit flat files and recomputes:
#   Step 1.1 Marker source audit
#   Step 1.2 Patient-level HR/BG pseudobulk robustness
#   Step 1.3 De novo NMF metaprogram construction
#   Step 1.4 TCGA/CGGA bulk transferability and clinical-model eligibility
#
# It does NOT read:
#   - h5ad
#   - previous Step1 outputs
#   - precomputed Gate2 robustness
#   - precomputed NMF/metaprogram results
#   - ClassA old results
#   - GSE16011 / GSE6011 external validation
#
# Required R packages:
#   install.packages(c("data.table", "Matrix", "survival"))
#
# Run:
#   source("G:/NEW/The2/SIDISH_Step1_Core_marker_refinement_reanalysis_v3_5_flat_inputs.R")
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
  # Put all files in this folder structure. See final assistant message.
  input_root = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs",
  output_dir = "G:/NEW/The2/SIDISH_Step1_Core_marker_refinement_v3_5_flat_inputs",

  # Single-cell raw files exported/used by Python 02_train_sidish_authorstyle_metadata.py.
  sc_mtx_file      = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/singlecell/scRNA_raw.mtx",
  sc_genes_file    = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/singlecell/scRNA_genes.tsv",
  sc_cells_file    = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/singlecell/scRNA_cells.tsv",
  sc_metadata_file = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/singlecell/scRNA_metadata.tsv",

  # Cell-level SIDISH labels exported by Python 02 export_cell_scores().
  cell_scores_file = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/sidish_exports/cell_scores.tsv",

  # Marker files from Python marker step / your SIDISH result extraction.
  hr_marker_file = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/markers/hr_markers.tsv",
  celltype_marker_file = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/markers/celltype_hr_markers_all_threshold.tsv",

  # Bulk cohorts. GSE16011/GSE6011 must not be placed here.
  bulk_cohorts = list(
    TCGA = list(
      expr = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/TCGA/bulk_expression_dedup.csv",
      clin = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/TCGA/clinical_aligned.tsv"
    ),
    CGGA325 = list(
      expr = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/CGGA325/bulk_expression_dedup.csv",
      clin = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/CGGA325/clinical_aligned.tsv"
    ),
    CGGA693 = list(
      expr = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/CGGA693/bulk_expression_dedup.csv",
      clin = "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk/CGGA693/clinical_aligned.tsv"
    )
  ),

  # v3.1: if the exact files above do not exist, recursively search these roots.
  # This fixes the common case where files are named TCGA_clinical_TCGA_format.tsv,
  # CGGA.mRNAseq_325_clinical..., or cohort-specific expression names rather than
  # exactly bulk_expression_dedup.csv / clinical_aligned.tsv.
  auto_discover_bulk_files = TRUE,
  bulk_search_roots = c(
    "G:/NEW/The2/SIDISH_core_marker_reanalysis_inputs/bulk",
    "G:/NEW/The2/Bulk_Cohorts",
    "G:/NEW/转化/Bulk_Cohorts",
    "G:/NEW/转化",
    "G:/NEW/The2"
  ),

  # Column preferences.
  cell_id_cols = c("cell_id", "cell_barcode", "cells", "original_cell_id"),
  patient_cols = c("patient_id_sc", "patient_id", "sample_id_sc", "sample_id"),
  celltype_cols = c("celltype_major", "celltype_level3", "celltype_level2", "celltype_level1", "Cell_Type", "cell_type"),
  sidish_value_cols = c("SIDISH_value", "sidish_value"),
  sidish_label_cols = c("sidish_label", "sidish_label_plot", "SIDISH", "SIDISH_"),

  # scRNA normalization.
  normalize_scRNA = TRUE,
  normalize_target_sum = 1e4,

  # Minimum cells per patient per group for paired pseudobulk.
  min_cells_per_patient_group = 10,
  min_paired_patients = 5,

  # Patient-level robustness thresholds.
  strict_median_delta = 0.03,
  strict_support_fraction = 0.65,
  supportive_median_delta = 0.01,
  supportive_support_fraction = 0.55,
  zero_fraction_cutoff = 0.60,

  # NMF.
  k_range = 2:6,
  selected_k_override = NA_integer_,
  nmf_modest_gain_threshold = 0.12,
  nmf_prefer_k_cap = 4,
  nmf_nrun = 50,
  nmf_max_iter = 1000,
  nmf_tol = 1e-5,
  top_n_genes_per_mp = 30,
  min_genes_per_mp = 5,

  # Publication-style plotting (Nature-level clean style)
  make_publication_figures = TRUE,
  plot_base_family = "Arial",
  plot_base_size = 8,
  plot_linewidth = 0.5,
  plot_point_size = 2.0,
  figure_width_single = 105/25.4,
  figure_width_double = 183/25.4,
  figure_height_unit = 2.9,

  # Bulk.
  endpoint = "os",
  time_col = "os_time",
  event_col = "os_event",
  sample_id_col = "sample_id",
  bulk_min_detected_fraction = 0.20,
  bulk_min_variance = 1e-8,
  bulk_min_available_cohorts = 2,

  # v3.5 eligibility hardening.
  # A1 requires both majority-positive cohort direction and positive mean/median beta.
  # This avoids borderline cases such as SEMA5A being called risk-positive when mean beta is ~0 or negative.
  A1_require_mean_beta_positive = TRUE,
  A1_require_median_beta_positive = TRUE,
  A1_min_mean_cox_beta = 0,
  A1_min_median_cox_beta = 0,
  A1_min_positive_cohort_fraction = 2/3,
  borderline_A2_label = "Class_A2_borderline_context_dependent_candidate",

  expected_marker_n = 245,

  # Safety.
  external_validation_forbidden_patterns = c("GSE16011", "GSE6011"),
  forbidden_old_output_patterns = c(
    "SIDISH_Step1_Core_marker_refinement_v1",
    "SIDISH_Step1_Core_marker_refinement_v2",
    "99_final_outputs",
    "ClassA_signed_core_marker_pool",
    "metaprogram_gene_loadings",
    "Gate2",
    "Gate3"
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
  step11 = file.path(CONFIG$output_dir, "01_1_marker_source_audit"),
  step12 = file.path(CONFIG$output_dir, "01_2_patient_level_pseudobulk_robustness"),
  step13 = file.path(CONFIG$output_dir, "01_3_de_novo_metaprogram_construction"),
  step14 = file.path(CONFIG$output_dir, "01_4_bulk_transferability_eligibility"),
  fig = file.path(CONFIG$output_dir, "02_publication_figures"),
  final = file.path(CONFIG$output_dir, "99_final_outputs")
)
invisible(lapply(OUT, safe_dir_create))

load_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required R package: ", pkg, "\nInstall with install.packages('", pkg, "')")
  }
  invisible(TRUE)
}
load_or_stop("data.table")
load_or_stop("Matrix")
load_or_stop("survival")
load_or_stop("ggplot2")

# Unified Nature-like palette (NPG-inspired), consistent across all figures.
SIDISH_COLORS <- c(
  red = "#E64B35FF",
  blue = "#4DBBD5FF",
  green = "#00A087FF",
  navy = "#3C5488FF",
  orange = "#F39B7FFF",
  purple = "#8491B4FF",
  teal = "#91D1C2FF",
  brown = "#7E6148FF",
  beige = "#B09C85FF",
  yellow = "#F7DC6FFF",
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
      axis.ticks.length = grid::unit(1.5, "mm"),
      legend.title = ggplot2::element_text(face = "bold", colour = "black"),
      legend.text = ggplot2::element_text(colour = "black"),
      legend.key.size = grid::unit(3.5, "mm"),
      panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = CONFIG$plot_linewidth),
      strip.background = ggplot2::element_rect(fill = "#F7F7F7", colour = "black", linewidth = CONFIG$plot_linewidth),
      strip.text = ggplot2::element_text(face = "bold", colour = "black"),
      plot.margin = grid::unit(c(5.5, 12, 5.5, 5.5), "pt")
    )
}

wrap_axis_labels <- function(x, width = 24) {
  vapply(as.character(x), function(xx) paste(strwrap(xx, width = width), collapse = "\n"), character(1))
}

save_pub_plot <- function(p, filename_stub, width, height) {
  safe_dir_create(OUT$fig)
  pdf_file <- file.path(OUT$fig, paste0(filename_stub, ".pdf"))
  png_file <- file.path(OUT$fig, paste0(filename_stub, ".png"))

  pdf_status <- "not_run"
  pdf_error <- NA_character_
  png_status <- "not_run"
  png_error <- NA_character_

  # Primary PDF export: cairo_pdf can usually use system fonts such as Arial on Windows.
  # Fallback: if Arial is not registered/available, export with base sans instead of
  # stopping the whole analysis. The requested plot theme remains Arial where supported.
  tryCatch({
    ggplot2::ggsave(
      filename = pdf_file,
      plot = p,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf,
      family = CONFIG$plot_base_family,
      bg = "white"
    )
    pdf_status <- paste0("ok_cairo_pdf_family_", CONFIG$plot_base_family)
  }, error = function(e1) {
    pdf_error <<- conditionMessage(e1)
    p_fallback <- p + ggplot2::theme(text = ggplot2::element_text(family = "sans"))
    tryCatch({
      ggplot2::ggsave(
        filename = pdf_file,
        plot = p_fallback,
        width = width,
        height = height,
        units = "in",
        device = grDevices::cairo_pdf,
        family = "sans",
        bg = "white"
      )
      pdf_status <<- "ok_cairo_pdf_family_sans_fallback"
    }, error = function(e2) {
      pdf_error <<- paste(pdf_error, "|| fallback:", conditionMessage(e2))
      tryCatch({
        grDevices::pdf(pdf_file, width = width, height = height, family = "sans", useDingbats = FALSE)
        print(p_fallback)
        grDevices::dev.off()
        pdf_status <<- "ok_pdf_family_sans_fallback"
      }, error = function(e3) {
        pdf_error <<- paste(pdf_error, "|| pdf fallback:", conditionMessage(e3))
        pdf_status <<- "failed"
        if (grDevices::dev.cur() > 1) grDevices::dev.off()
      })
    })
  })

  tryCatch({
    ggplot2::ggsave(
      filename = png_file,
      plot = p,
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      bg = "white"
    )
    png_status <- "ok_png_600dpi"
  }, error = function(e1) {
    png_error <<- conditionMessage(e1)
    p_fallback <- p + ggplot2::theme(text = ggplot2::element_text(family = "sans"))
    ggplot2::ggsave(
      filename = png_file,
      plot = p_fallback,
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      bg = "white"
    )
    png_status <<- "ok_png_600dpi_sans_fallback"
  })

  log_file <- file.path(OUT$fig, "figure_export_status.tsv")
  log_row <- data.table::data.table(
    figure = filename_stub,
    pdf_file = pdf_file,
    png_file = png_file,
    requested_family = CONFIG$plot_base_family,
    pdf_status = pdf_status,
    pdf_error = pdf_error,
    png_status = png_status,
    png_error = png_error,
    width = width,
    height = height
  )
  if (file.exists(log_file)) {
    old_log <- data.table::fread(log_file, sep = "\t", data.table = TRUE, showProgress = FALSE)
    log_row <- data.table::rbindlist(list(old_log, log_row), fill = TRUE)
  }
  write_tsv(log_row, log_file)

  if (identical(pdf_status, "failed")) {
    warning("PDF export failed for ", filename_stub, ". PNG may still be available. See figure_export_status.tsv")
  }
}

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
  if (stringi_available()) {
    # no-op; kept only for future extension
  }
  n_tab <- lengths(regmatches(ln, gregexpr("\t", ln, fixed = TRUE)))
  n_comma <- lengths(regmatches(ln, gregexpr(",", ln, fixed = TRUE)))
  ifelse(n_tab > n_comma, "\t", ",")
}

stringi_available <- function() FALSE

read_table_auto <- function(file) {
  file <- norm_path(file)
  assert_safe_input(file)
  if (!file.exists(file)) stop("File not found: ", file)
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("tsv", "txt")) {
    sep <- "\t"
  } else if (ext %in% c("csv")) {
    sep <- detect_delim(file)
  } else {
    sep <- detect_delim(file)
  }
  dt <- data.table::fread(file, sep = sep, data.table = TRUE, showProgress = FALSE)
  data.table::setnames(dt, clean_names(names(dt)))
  dt
}

clean_names <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("[.]+", "_", x)
  x
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

collapse_unique <- function(x, sep = ";") {
  x <- unique(as.character(x[!is.na(x) & x != ""]))
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = sep)
}

assert_safe_input <- function(file) {
  file <- norm_path(file)
  if (file == "") return(invisible(TRUE))
  for (pat in CONFIG$external_validation_forbidden_patterns) {
    if (grepl(pat, file, ignore.case = TRUE)) {
      stop("Forbidden external validation input path detected: ", file)
    }
  }
  for (pat in CONFIG$forbidden_old_output_patterns) {
    if (grepl(pat, file, ignore.case = TRUE)) {
      stop("Forbidden old/intermediate Step1 output input path detected: ", file)
    }
  }
  invisible(TRUE)
}

require_input_file <- function(file, label) {
  file <- norm_path(file)
  assert_safe_input(file)
  if (!file.exists(file)) stop("Missing ", label, ": ", file)
  invisible(file)
}

read_one_col <- function(file, label) {
  dt <- data.table::fread(file, sep = "\t", header = FALSE, data.table = FALSE, showProgress = FALSE)
  if (ncol(dt) < 1 || nrow(dt) < 1) stop(label, " is empty: ", file)
  as.character(dt[[1]])
}

standardize_sidish_label <- function(scores_dt) {
  value_col <- first_existing_col(scores_dt, CONFIG$sidish_value_cols, regex = FALSE)
  if (!is.na(value_col)) {
    v <- as_numeric_safe(scores_dt[[value_col]])
    uniq <- sort(unique(v[is.finite(v)]))
    if (length(uniq) > 0 && all(uniq %in% c(0, 1))) {
      return(ifelse(v == 1, "HR", ifelse(v == 0, "BG", NA_character_)))
    }
  }

  label_col <- first_existing_col(scores_dt, CONFIG$sidish_label_cols, regex = FALSE)
  if (!is.na(label_col)) {
    x <- toupper(trimws(as.character(scores_dt[[label_col]])))
    out <- rep(NA_character_, length(x))
    out[x %in% c("HR", "H", "HIGH", "HIGH-RISK", "HIGH_RISK", "1", "1.0", "TRUE")] <- "HR"
    out[x %in% c("BG", "B", "BACKGROUND", "LOW", "LOW-RISK", "LOW_RISK", "0", "0.0", "FALSE")] <- "BG"
    out[is.na(out) & grepl("HIGH|RISK_HIGH|SIDISH_HIGH", x)] <- "HR"
    out[is.na(out) & grepl("BACKGROUND|LOW|CONTROL|REF", x)] <- "BG"
    return(out)
  }

  stop("Cannot infer SIDISH HR/BG labels from cell_scores. Expected SIDISH_value or sidish_label/SIDISH.")
}


###############################################################################
# 1b. Bulk file discovery and sample matching helpers
###############################################################################

list_candidate_bulk_files <- function(roots) {
  roots <- norm_path(roots)
  roots <- roots[roots != "" & dir.exists(roots)]
  if (length(roots) == 0) {
    return(data.table::data.table(
      file = character(), basename = character(), ext = character(),
      size_bytes = numeric(), modified_time = character()
    ))
  }
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
  if (is.null(out)) {
    out <- data.table::data.table(
      file = character(), basename = character(), ext = character(),
      size_bytes = numeric(), modified_time = character()
    )
  }
  out[, forbidden_external := grepl(paste(CONFIG$external_validation_forbidden_patterns, collapse = "|"), file, ignore.case = TRUE)]
  out[, forbidden_old_output := grepl(paste(CONFIG$forbidden_old_output_patterns, collapse = "|"), file, ignore.case = TRUE)]
  out <- out[forbidden_external == FALSE & forbidden_old_output == FALSE]
  out[ext %in% c("csv", "tsv", "txt")]
}

score_bulk_candidate <- function(inv, cohort, type = c("expr", "clin")) {
  type <- match.arg(type)
  if (nrow(inv) == 0) return(data.table::data.table())

  dt <- data.table::copy(inv)
  dt[, score := 0]

  cohort_patterns <- switch(
    cohort,
    TCGA = c("TCGA"),
    CGGA325 = c("CGGA.*325", "325.*CGGA", "CGGA_325", "mRNAseq_325"),
    CGGA693 = c("CGGA.*693", "693.*CGGA", "CGGA_693", "mRNAseq_693"),
    c(cohort)
  )

  for (pat in cohort_patterns) {
    dt[grepl(pat, file, ignore.case = TRUE), score := score + 20]
    dt[grepl(pat, basename, ignore.case = TRUE), score := score + 10]
  }

  # Prefer cohort subfolders.
  dt[grepl(paste0("/|\\\\", cohort, "(/|\\\\)"), file, ignore.case = TRUE), score := score + 20]

  if (type == "expr") {
    include <- c(
      "bulk_expression_dedup", "expression", "expr", "mRNAseq", "RNAseq",
      "RSEM", "FPKM", "TPM", "count", "counts", "matrix"
    )
    exclude <- c(
      "clinical", "clin", "survival", "phenotype", "pheno", "meta", "metadata",
      "sample_mapping", "report", "summary", "marker", "cell_scores"
    )
  } else {
    include <- c(
      "clinical_aligned", "clinical", "clin", "survival", "phenotype",
      "pheno", "TCGA_format", "sample"
    )
    exclude <- c(
      "bulk_expression", "expression", "expr", "mRNAseq", "RNAseq", "RSEM",
      "FPKM", "TPM", "count", "counts", "matrix", "marker", "cell_scores"
    )
  }

  for (pat in include) {
    dt[grepl(pat, basename, ignore.case = TRUE), score := score + 8]
    dt[grepl(pat, file, ignore.case = TRUE), score := score + 2]
  }
  for (pat in exclude) {
    dt[grepl(pat, basename, ignore.case = TRUE), score := score - 25]
    dt[grepl(pat, file, ignore.case = TRUE), score := score - 5]
  }

  # Keep only plausible candidates.
  dt <- dt[score > 0]
  if (nrow(dt) == 0) return(dt)
  dt[order(-score, -size_bytes)]
}

discover_one_bulk_file <- function(cohort, type = c("expr", "clin"), inv) {
  type <- match.arg(type)
  exact <- CONFIG$bulk_cohorts[[cohort]][[type]]
  exact <- norm_path(exact)
  if (file.exists(exact)) {
    return(list(file = exact, mode = "exact_CONFIG_path", score = NA_real_))
  }
  if (!isTRUE(CONFIG$auto_discover_bulk_files)) {
    return(list(file = exact, mode = "missing_exact_CONFIG_path", score = NA_real_))
  }

  cand <- score_bulk_candidate(inv, cohort, type)
  if (nrow(cand) == 0) {
    return(list(file = "", mode = "not_found_by_recursive_discovery", score = NA_real_))
  }
  list(file = cand$file[1], mode = "auto_discovered_recursive", score = cand$score[1])
}

tcga_sample_barcode <- function(x) {
  x <- standardize_id(x)
  parts <- strsplit(x, "-", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) >= 4 && p[1] == "TCGA") paste(p[1:4], collapse = "-") else paste(p, collapse = "-")
  }, character(1))
}

tcga_patient_barcode <- function(x) {
  x <- standardize_id(x)
  parts <- strsplit(x, "-", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) >= 3 && p[1] == "TCGA") paste(p[1:3], collapse = "-") else paste(p, collapse = "-")
  }, character(1))
}

sample_match_key <- function(x, mode) {
  if (mode == "exact") return(as.character(x))
  if (mode == "clean") return(standardize_id(x))
  if (mode == "tcga_sample") return(tcga_sample_barcode(x))
  if (mode == "tcga_patient") return(tcga_patient_barcode(x))
  stop("Unknown sample match mode: ", mode)
}

evaluate_sample_match_modes <- function(expr_samples, clin_samples) {
  modes <- c("exact", "clean", "tcga_sample", "tcga_patient")
  data.table::rbindlist(lapply(modes, function(md) {
    ek <- sample_match_key(expr_samples, md)
    ck <- sample_match_key(clin_samples, md)
    data.table::data.table(
      mode = md,
      n_expr_keys = data.table::uniqueN(ek),
      n_clin_keys = data.table::uniqueN(ck),
      n_overlap = length(intersect(unique(ek), unique(ck)))
    )
  }))
}

align_bulk_expr_clin <- function(expr, clin, cohort, out_dir) {
  expr_samples <- colnames(expr)
  clin_samples <- clin$sample_id
  mode_report <- evaluate_sample_match_modes(expr_samples, clin_samples)
  mode_report[, cohort := cohort]
  write_tsv(mode_report, file.path(out_dir, paste0("SIDISH_Step1_4_", cohort, "_sample_match_mode_report.tsv")))

  safe_modes <- c("exact", "clean", "tcga_sample")
  safe_report <- mode_report[mode %in% safe_modes][order(-n_overlap)]
  if (nrow(safe_report) > 0 && safe_report$n_overlap[1] > 0) {
    selected_mode <- safe_report$mode[1]
  } else {
    patient_n <- mode_report[mode == "tcga_patient", n_overlap]
    if (length(patient_n) > 0 && patient_n > 0) {
      selected_mode <- "tcga_patient"
    } else {
      return(list(expr = NULL, clin = NULL, selected_mode = "none", n_common = 0L, mapping = data.table::data.table()))
    }
  }

  expr_dt <- data.table::data.table(
    expr_col = expr_samples,
    expr_key = sample_match_key(expr_samples, selected_mode)
  )
  clin_dt <- data.table::data.table(
    sample_id = clin_samples,
    clin_key = sample_match_key(clin_samples, selected_mode),
    clin_row = seq_len(nrow(clin))
  )

  # Drop duplicate clinical keys conservatively: keep first, report all duplicates.
  clin_dup <- clin_dt[duplicated(clin_key) | duplicated(clin_key, fromLast = TRUE)]
  if (nrow(clin_dup) > 0) {
    write_tsv(clin_dup, file.path(out_dir, paste0("SIDISH_Step1_4_", cohort, "_duplicate_clinical_match_keys.tsv")))
    clin_dt <- clin_dt[!duplicated(clin_key)]
  }

  common_keys <- intersect(unique(expr_dt$expr_key), unique(clin_dt$clin_key))
  if (length(common_keys) == 0) {
    return(list(expr = NULL, clin = NULL, selected_mode = selected_mode, n_common = 0L, mapping = data.table::data.table()))
  }

  # Aggregate duplicate expression columns that map to the same clinical key.
  expr_aligned <- matrix(NA_real_, nrow = nrow(expr), ncol = length(common_keys))
  rownames(expr_aligned) <- rownames(expr)
  colnames(expr_aligned) <- common_keys

  mapping_rows <- list()
  for (i in seq_along(common_keys)) {
    key <- common_keys[i]
    cols <- expr_dt[expr_key == key, expr_col]
    mat <- expr[, cols, drop = FALSE]
    if (length(cols) == 1) {
      expr_aligned[, i] <- as.numeric(mat[, 1])
    } else {
      expr_aligned[, i] <- rowMeans(mat, na.rm = TRUE)
    }
    mapping_rows[[i]] <- data.table::data.table(
      cohort = cohort,
      selected_match_mode = selected_mode,
      match_key = key,
      clinical_sample_id = clin_dt[clin_key == key, sample_id][1],
      n_expr_columns_mapped = length(cols),
      expr_columns_mapped = paste(cols, collapse = ";")
    )
  }

  mapping <- data.table::rbindlist(mapping_rows, fill = TRUE)
  write_tsv(mapping, file.path(out_dir, paste0("SIDISH_Step1_4_", cohort, "_sample_mapping.tsv")))

  clin_aligned <- clin[match(mapping$clinical_sample_id, clin$sample_id)]
  colnames(expr_aligned) <- mapping$clinical_sample_id
  list(expr = expr_aligned, clin = clin_aligned, selected_mode = selected_mode, n_common = ncol(expr_aligned), mapping = mapping)
}


###############################################################################
# 2. Input validation
###############################################################################

msg("SIDISH Step 1 v3.5 flat-input reanalysis started.")
msg("No h5ad will be read.")

mandatory <- c(
  CONFIG$sc_mtx_file,
  CONFIG$sc_genes_file,
  CONFIG$sc_cells_file,
  CONFIG$sc_metadata_file,
  CONFIG$cell_scores_file,
  CONFIG$hr_marker_file,
  CONFIG$celltype_marker_file
)
invisible(mapply(require_input_file, mandatory, c(
  "scRNA_raw.mtx",
  "scRNA_genes.tsv",
  "scRNA_cells.tsv",
  "scRNA_metadata.tsv",
  "cell_scores.tsv",
  "hr_markers.tsv",
  "celltype_hr_markers_all_threshold.tsv"
)))

input_manifest <- data.table::data.table(
  role = c("sc_mtx", "sc_genes", "sc_cells", "sc_metadata", "cell_scores", "hr_marker", "celltype_marker"),
  file = mandatory,
  exists = file.exists(mandatory),
  size_bytes = as.numeric(file.info(mandatory)$size)
)
write_tsv(input_manifest, file.path(OUT$audit, "selected_flat_input_files.tsv"))

bulk_inventory <- list_candidate_bulk_files(CONFIG$bulk_search_roots)
write_tsv(bulk_inventory, file.path(OUT$audit, "bulk_candidate_file_inventory.tsv"))

bulk_manifest <- data.table::rbindlist(lapply(names(CONFIG$bulk_cohorts), function(coh) {
  expr_hit <- discover_one_bulk_file(coh, "expr", bulk_inventory)
  clin_hit <- discover_one_bulk_file(coh, "clin", bulk_inventory)
  data.table::data.table(
    cohort = coh,
    expr_file = expr_hit$file,
    clin_file = clin_hit$file,
    expr_discovery_mode = expr_hit$mode,
    clin_discovery_mode = clin_hit$mode,
    expr_score = expr_hit$score,
    clin_score = clin_hit$score,
    expr_exists = expr_hit$file != "" & file.exists(expr_hit$file),
    clin_exists = clin_hit$file != "" & file.exists(clin_hit$file)
  )
}), fill = TRUE)

write_tsv(bulk_manifest, file.path(OUT$audit, "selected_bulk_input_files.tsv"))

if (sum(bulk_manifest$expr_exists & bulk_manifest$clin_exists, na.rm = TRUE) == 0) {
  stop("No usable bulk cohort expression+clinical pair was found. Inspect: ",
       file.path(OUT$audit, "bulk_candidate_file_inventory.tsv"), " and ",
       file.path(OUT$audit, "selected_bulk_input_files.tsv"))
}

###############################################################################
# 3. Step 1.1 Marker source audit
###############################################################################

msg("Step 1.1 Marker source audit...")

standardize_marker_table <- function(dt) {
  gene_col <- first_existing_col(dt, c("^gene$", "gene_symbol", "symbol", "names", "feature", "marker"), regex = TRUE)
  if (is.na(gene_col)) gene_col <- names(dt)[1]
  out <- data.table::data.table(gene = normalize_gene(dt[[gene_col]]))

  for (nm in c("logfoldchanges", "log2FC", "logFC", "scores", "pvals", "pvals_adj", "padj", "p_adj")) {
    col <- first_existing_col(dt, c(nm), regex = FALSE)
    if (!is.na(col)) out[[nm]] <- dt[[col]]
  }

  out <- out[!is.na(gene) & gene != ""]
  out[, duplicated_in_marker_file := duplicated(gene) | duplicated(gene, fromLast = TRUE)]
  unique(out, by = "gene")
}

standardize_celltype_marker_table <- function(dt) {
  gene_col <- first_existing_col(dt, c("^gene$", "gene_symbol", "symbol", "names", "feature", "marker"), regex = TRUE)
  if (is.na(gene_col)) stop("Cannot identify gene column in celltype marker file.")
  cell_col <- first_existing_col(dt, c("celltype_major", "celltype", "cell_type", "cell.type", "annotation", "cluster", "level"), regex = TRUE)
  if (is.na(cell_col)) stop("Cannot identify celltype column in celltype marker file.")

  out <- data.table::data.table(
    gene = normalize_gene(dt[[gene_col]]),
    celltype = as.character(dt[[cell_col]])
  )
  for (nm in c("logfoldchanges", "log2FC", "logFC", "avg_log2FC", "avg_logFC", "scores", "pvals", "pvals_adj", "padj", "p_adj", "n_hr", "n_bg")) {
    col <- first_existing_col(dt, c(nm), regex = FALSE)
    if (!is.na(col)) out[[nm]] <- dt[[col]]
  }
  out <- out[!is.na(gene) & gene != ""]
  out[, celltype := ifelse(is.na(celltype) | celltype == "", "celltype_unknown", celltype)]
  unique(out)
}

hr_raw <- read_table_auto(CONFIG$hr_marker_file)
global_markers <- standardize_marker_table(hr_raw)

if (nrow(global_markers) != CONFIG$expected_marker_n) {
  stop("hr_markers.tsv marker count mismatch. Expected ", CONFIG$expected_marker_n, ", got ", nrow(global_markers))
}
global_markers[, in_245_marker_universe := TRUE]

ct_raw <- read_table_auto(CONFIG$celltype_marker_file)
celltype_markers <- standardize_celltype_marker_table(ct_raw)

celltype_attr <- celltype_markers[, .(
  n_celltypes_with_HR_support = data.table::uniqueN(celltype),
  celltype_HR_support = collapse_unique(celltype),
  n_gene_celltype_evidence = .N
), by = gene]

marker_source_master <- merge(global_markers, celltype_attr, by = "gene", all.x = TRUE)
marker_source_master[, has_celltype_HR_support := !is.na(n_celltypes_with_HR_support) & n_celltypes_with_HR_support > 0]
marker_source_master[, attribution_category := data.table::fcase(
  is.na(n_celltypes_with_HR_support) | n_celltypes_with_HR_support == 0, "unassigned_no_celltype_HR_support",
  n_celltypes_with_HR_support == 1, "single_celltype_associated",
  n_celltypes_with_HR_support > 1, "shared_multi_celltype_associated",
  default = "unclassified"
)]

write_tsv(global_markers, file.path(OUT$step11, "SIDISH_Step1_1_global_marker_universe.tsv"))
write_tsv(celltype_markers, file.path(OUT$step11, "SIDISH_Step1_1_celltype_marker_evidence_long.tsv"))
write_tsv(marker_source_master, file.path(OUT$step11, "SIDISH_Step1_1_marker_source_audit_master.tsv"))

step11_summary <- data.table::data.table(
  metric = c(
    "n_global_markers",
    "n_markers_with_celltype_HR_support",
    "n_single_celltype_associated",
    "n_shared_multi_celltype_associated",
    "n_unassigned_no_celltype_HR_support",
    "n_gene_celltype_HR_evidence_pairs"
  ),
  value = c(
    nrow(marker_source_master),
    sum(marker_source_master$has_celltype_HR_support, na.rm = TRUE),
    sum(marker_source_master$attribution_category == "single_celltype_associated", na.rm = TRUE),
    sum(marker_source_master$attribution_category == "shared_multi_celltype_associated", na.rm = TRUE),
    sum(marker_source_master$attribution_category == "unassigned_no_celltype_HR_support", na.rm = TRUE),
    nrow(celltype_markers)
  )
)
write_tsv(step11_summary, file.path(OUT$step11, "SIDISH_Step1_1_marker_source_summary.tsv"))

###############################################################################
# 4. Step 1.2 Patient-level pseudobulk robustness from MTX + cell_scores
###############################################################################

msg("Step 1.2 Patient-level pseudobulk robustness from MTX + cell_scores...")

genes <- read_one_col(CONFIG$sc_genes_file, "scRNA_genes.tsv")
cells <- read_one_col(CONFIG$sc_cells_file, "scRNA_cells.tsv")

genes_norm <- normalize_gene(genes)
cells_std <- as.character(cells)

msg("Reading scRNA_raw.mtx. This may take time for large data.")
x <- Matrix::readMM(CONFIG$sc_mtx_file)
x <- methods::as(x, "dgCMatrix")

if (nrow(x) != length(genes) || ncol(x) != length(cells)) {
  # Some exported MTX may be cells x genes; fix if exactly transposed.
  if (nrow(x) == length(cells) && ncol(x) == length(genes)) {
    msg("MTX appears to be cells x genes; transposing to genes x cells.")
    x <- Matrix::t(x)
  } else {
    stop("MTX dimension mismatch. MTX=", paste(dim(x), collapse = "x"),
         "; genes=", length(genes), "; cells=", length(cells))
  }
}
rownames(x) <- genes_norm
colnames(x) <- cells_std

if (anyDuplicated(rownames(x))) {
  msg("Duplicated gene symbols in scRNA matrix; collapsing duplicated gene rows by sum.")
  grp <- split(seq_len(nrow(x)), rownames(x))
  collapsed <- lapply(names(grp), function(g) Matrix::colSums(x[grp[[g]], , drop = FALSE]))
  x <- do.call(rbind, collapsed)
  x <- Matrix::Matrix(x, sparse = TRUE)
  rownames(x) <- names(grp)
  colnames(x) <- cells_std
}

metadata <- read_table_auto(CONFIG$sc_metadata_file)
cell_scores <- read_table_auto(CONFIG$cell_scores_file)

# Identify cell IDs.
meta_cell_col <- first_existing_col(metadata, CONFIG$cell_id_cols, regex = FALSE)
if (is.na(meta_cell_col)) meta_cell_col <- first_existing_col(metadata, c("cell", "barcode"), regex = TRUE)
score_cell_col <- first_existing_col(cell_scores, CONFIG$cell_id_cols, regex = FALSE)
if (is.na(score_cell_col)) score_cell_col <- first_existing_col(cell_scores, c("cell", "barcode"), regex = TRUE)

if (is.na(meta_cell_col)) stop("Cannot identify cell ID column in scRNA_metadata.tsv")
if (is.na(score_cell_col)) stop("Cannot identify cell ID column in cell_scores.tsv")

metadata[, cell_key := as.character(get(meta_cell_col))]
cell_scores[, cell_key := as.character(get(score_cell_col))]
metadata[, cell_key_std := standardize_id(cell_key)]
cell_scores[, cell_key_std := standardize_id(cell_key)]
cell_order_dt <- data.table::data.table(cell_key = cells_std, mtx_col_index = seq_along(cells_std))
cell_order_dt[, cell_key_std := standardize_id(cell_key)]

# Merge mtx cell order with metadata and scores.
cell_annot <- merge(cell_order_dt, metadata, by = "cell_key_std", all.x = TRUE, suffixes = c("", "_meta"))
cell_annot <- merge(cell_annot, cell_scores, by = "cell_key_std", all.x = TRUE, suffixes = c("", "_score"))

# Detect patient column after merge.
patient_col <- first_existing_col(cell_annot, CONFIG$patient_cols, regex = FALSE)
if (is.na(patient_col)) patient_col <- first_existing_col(cell_annot, c("patient", "sample", "donor", "case"), regex = TRUE)
if (is.na(patient_col)) stop("Cannot identify patient/sample column after merging metadata and cell_scores.")

celltype_col <- first_existing_col(cell_annot, CONFIG$celltype_cols, regex = FALSE)
sidish_label <- standardize_sidish_label(cell_annot)

cell_annot[, patient_id := standardize_id(get(patient_col))]
cell_annot[, sidish_label := sidish_label]
cell_annot[, celltype_for_audit := if (!is.na(celltype_col)) as.character(get(celltype_col)) else NA_character_]

merge_qc <- data.table::data.table(
  metric = c(
    "n_mtx_cells",
    "n_cells_with_metadata",
    "n_cells_with_scores",
    "patient_col_used",
    "celltype_col_used",
    "n_HR_cells",
    "n_BG_cells",
    "n_UNK_or_missing_sidish_cells"
  ),
  value = c(
    ncol(x),
    sum(!is.na(cell_annot[[meta_cell_col]])),
    sum(!is.na(cell_annot[[score_cell_col]])),
    patient_col,
    ifelse(is.na(celltype_col), "NA", celltype_col),
    sum(cell_annot$sidish_label == "HR", na.rm = TRUE),
    sum(cell_annot$sidish_label == "BG", na.rm = TRUE),
    sum(!cell_annot$sidish_label %in% c("HR", "BG") | is.na(cell_annot$sidish_label))
  )
)
write_tsv(merge_qc, file.path(OUT$step12, "SIDISH_Step1_2_cell_metadata_score_merge_QC.tsv"))
write_tsv(cell_annot[, .(cell_key, mtx_col_index, patient_id, sidish_label, celltype_for_audit)],
          file.path(OUT$step12, "SIDISH_Step1_2_cell_annotation_used.tsv"))

# Select marker genes.
marker_genes <- marker_source_master$gene
present_marker_genes <- intersect(marker_genes, rownames(x))
missing_marker_genes <- setdiff(marker_genes, rownames(x))
write_tsv(data.table::data.table(gene = missing_marker_genes), file.path(OUT$step12, "SIDISH_Step1_2_missing_marker_genes_in_scRNA.tsv"))
if (length(present_marker_genes) < 20) {
  stop("Too few marker genes found in scRNA matrix: ", length(present_marker_genes), "/", length(marker_genes))
}

# Valid cells: have patient_id and HR/BG.
valid_cell_dt <- cell_annot[
  !is.na(patient_id) & patient_id != "" & sidish_label %in% c("HR", "BG"),
  .(cell_key, mtx_col_index, patient_id, sidish_label)
]

# Patients with enough HR and BG cells.
cell_count_wide <- data.table::dcast(
  valid_cell_dt[, .N, by = .(patient_id, sidish_label)],
  patient_id ~ sidish_label,
  value.var = "N",
  fill = 0
)
if (!all(c("HR", "BG") %in% names(cell_count_wide))) {
  stop("Both HR and BG cells are required. Detected groups: ", paste(names(cell_count_wide), collapse = ", "))
}
valid_patients <- cell_count_wide[
  HR >= CONFIG$min_cells_per_patient_group & BG >= CONFIG$min_cells_per_patient_group,
  patient_id
]
write_tsv(cell_count_wide, file.path(OUT$step12, "SIDISH_Step1_2_cell_counts_by_patient_HR_BG.tsv"))

if (length(valid_patients) < CONFIG$min_paired_patients) {
  stop("Too few patients with both HR and BG cells. Required ", CONFIG$min_paired_patients, ", got ", length(valid_patients))
}

valid_cell_dt <- valid_cell_dt[patient_id %in% valid_patients]
valid_col_idx <- valid_cell_dt$mtx_col_index

# Normalize scRNA counts for selected marker genes.
marker_counts <- x[present_marker_genes, valid_col_idx, drop = FALSE]
colnames(marker_counts) <- valid_cell_dt$cell_key
rownames(marker_counts) <- present_marker_genes

lib_size <- Matrix::colSums(x[, valid_col_idx, drop = FALSE])
lib_size[!is.finite(lib_size) | lib_size <= 0] <- NA_real_

if (CONFIG$normalize_scRNA) {
  scale_factor <- CONFIG$normalize_target_sum / lib_size
  scale_factor[!is.finite(scale_factor)] <- 0
  expr_sub <- marker_counts %*% Matrix::Diagonal(x = scale_factor)
  expr_sub <- log1p(expr_sub)
} else {
  expr_sub <- marker_counts
}

expr_qc <- data.table::data.table(
  metric = c("n_marker_genes_present", "n_valid_cells", "n_valid_patients", "nonzero_fraction", "min_lib_size", "median_lib_size", "max_lib_size"),
  value = c(
    nrow(expr_sub),
    ncol(expr_sub),
    length(valid_patients),
    Matrix::nnzero(expr_sub) / (nrow(expr_sub) * ncol(expr_sub)),
    min(lib_size, na.rm = TRUE),
    median(lib_size, na.rm = TRUE),
    max(lib_size, na.rm = TRUE)
  )
)
write_tsv(expr_qc, file.path(OUT$step12, "SIDISH_Step1_2_marker_expression_matrix_QC.tsv"))

# Pseudobulk means.
pb_list <- vector("list", length(valid_patients) * 2L)
idx <- 1L
for (pid in valid_patients) {
  for (grp in c("HR", "BG")) {
    cols <- which(valid_cell_dt$patient_id == pid & valid_cell_dt$sidish_label == grp)
    if (length(cols) < CONFIG$min_cells_per_patient_group) next
    mu <- Matrix::rowMeans(expr_sub[, cols, drop = FALSE])
    pb_list[[idx]] <- data.table::data.table(
      gene = rownames(expr_sub),
      patient_id = pid,
      HR_BG_group = grp,
      pseudobulk_mean = as.numeric(mu),
      n_cells = length(cols)
    )
    idx <- idx + 1L
  }
}
pseudobulk_long <- data.table::rbindlist(pb_list, fill = TRUE)
write_tsv(pseudobulk_long, file.path(OUT$step12, "SIDISH_Step1_2_de_novo_patient_HR_BG_pseudobulk_long.tsv"))

pb_wide <- data.table::dcast(pseudobulk_long, gene + patient_id ~ HR_BG_group, value.var = "pseudobulk_mean")
n_wide <- data.table::dcast(pseudobulk_long, gene + patient_id ~ HR_BG_group, value.var = "n_cells")
data.table::setnames(n_wide, intersect(c("HR", "BG"), names(n_wide)), paste0("n_cells_", intersect(c("HR", "BG"), names(n_wide))))

paired <- merge(pb_wide, n_wide, by = c("gene", "patient_id"), all = FALSE)
paired <- paired[!is.na(HR) & !is.na(BG)]
paired[, patient_delta_HR_minus_BG := HR - BG]
paired[, support_HR_up := patient_delta_HR_minus_BG > 0]
paired[, zero_like := HR <= 0 | BG <= 0]
write_tsv(paired, file.path(OUT$step12, "SIDISH_Step1_2_de_novo_patient_paired_HR_BG_delta.tsv"))

robustness <- paired[, .(
  n_paired_patients = .N,
  patient_median_delta = median(patient_delta_HR_minus_BG, na.rm = TRUE),
  patient_mean_delta = mean(patient_delta_HR_minus_BG, na.rm = TRUE),
  patient_support_fraction = mean(support_HR_up, na.rm = TRUE),
  zero_fraction = mean(zero_like, na.rm = TRUE),
  median_HR_pseudobulk = median(HR, na.rm = TRUE),
  median_BG_pseudobulk = median(BG, na.rm = TRUE)
), by = gene]

# Direction diagnostic: 245 marker universe should be HR-up globally.
direction_audit <- data.table::data.table(
  metric = c(
    "global_median_patient_delta",
    "fraction_genes_median_delta_positive",
    "fraction_genes_median_delta_negative",
    "n_genes_positive_median_delta",
    "n_genes_negative_median_delta"
  ),
  value = c(
    median(robustness$patient_median_delta, na.rm = TRUE),
    mean(robustness$patient_median_delta > 0, na.rm = TRUE),
    mean(robustness$patient_median_delta < 0, na.rm = TRUE),
    sum(robustness$patient_median_delta > 0, na.rm = TRUE),
    sum(robustness$patient_median_delta < 0, na.rm = TRUE)
  )
)
write_tsv(direction_audit, file.path(OUT$step12, "SIDISH_Step1_2_HR_BG_direction_audit.tsv"))

robustness[, patient_direction := data.table::fcase(
  patient_median_delta > 0, "HR_up",
  patient_median_delta < 0, "HR_down",
  patient_median_delta == 0, "neutral",
  default = "unknown"
)]

robustness[, robustness_tier := data.table::fcase(
  n_paired_patients < CONFIG$min_paired_patients,
  "Step1_2D_hold_exclude_insufficient_paired_patients",

  patient_median_delta >= CONFIG$strict_median_delta &
    patient_support_fraction >= CONFIG$strict_support_fraction &
    (is.na(zero_fraction) | zero_fraction < CONFIG$zero_fraction_cutoff),
  "Step1_2A_strict_patient_robust",

  patient_median_delta > CONFIG$supportive_median_delta &
    patient_support_fraction >= CONFIG$supportive_support_fraction &
    (is.na(zero_fraction) | zero_fraction < CONFIG$zero_fraction_cutoff),
  "Step1_2B_direction_consistent_supportive",

  patient_median_delta > CONFIG$supportive_median_delta &
    patient_support_fraction >= CONFIG$supportive_support_fraction &
    !is.na(zero_fraction) & zero_fraction >= CONFIG$zero_fraction_cutoff,
  "Step1_2C_zero_inflated_supportive",

  default = "Step1_2D_hold_exclude_not_patient_robust"
)]

robustness[, primary_MP_input_flag := robustness_tier %in% c(
  "Step1_2A_strict_patient_robust",
  "Step1_2B_direction_consistent_supportive"
)]
robustness[, sensitivity_MP_input_flag := robustness_tier %in% c(
  "Step1_2A_strict_patient_robust",
  "Step1_2B_direction_consistent_supportive",
  "Step1_2C_zero_inflated_supportive"
)]

robustness_master <- merge(
  marker_source_master[, .(gene, attribution_category, has_celltype_HR_support, celltype_HR_support)],
  robustness,
  by = "gene",
  all.x = TRUE
)
robustness_master[is.na(robustness_tier), robustness_tier := "Step1_2D_hold_exclude_not_matched_or_not_paired"]
robustness_master[is.na(primary_MP_input_flag), primary_MP_input_flag := FALSE]
robustness_master[is.na(sensitivity_MP_input_flag), sensitivity_MP_input_flag := FALSE]

write_tsv(robustness_master, file.path(OUT$step12, "SIDISH_Step1_2_patient_level_robustness_master.tsv"))

tier_counts <- robustness_master[, .(
  n_genes = .N,
  n_primary_MP_input = sum(primary_MP_input_flag, na.rm = TRUE),
  n_sensitivity_MP_input = sum(sensitivity_MP_input_flag, na.rm = TRUE)
), by = robustness_tier][order(robustness_tier)]
write_tsv(tier_counts, file.path(OUT$step12, "SIDISH_Step1_2_tier_counts.tsv"))

if (sum(robustness_master$primary_MP_input_flag, na.rm = TRUE) < CONFIG$min_genes_per_mp) {
  stop(
    "Too few primary MP input genes after flat-input patient-level robustness: ",
    sum(robustness_master$primary_MP_input_flag, na.rm = TRUE),
    ". Inspect Step1_2_HR_BG_direction_audit.tsv, marker_expression_matrix_QC.tsv, and robustness_master.tsv."
  )
}

###############################################################################
# 5. Step 1.3 De novo NMF metaprogram construction
###############################################################################

msg("Step 1.3 De novo NMF metaprogram construction...")

primary_genes <- robustness_master[primary_MP_input_flag == TRUE, gene]
hr_pb <- pseudobulk_long[HR_BG_group == "HR" & gene %in% primary_genes]
hr_mat_dt <- hr_pb[, .(value = mean(pseudobulk_mean, na.rm = TRUE)), by = .(gene, patient_id)]
hr_wide <- data.table::dcast(hr_mat_dt, gene ~ patient_id, value.var = "value")
nmf_genes <- hr_wide$gene
nmf_mat <- as.matrix(hr_wide[, -"gene"])
rownames(nmf_mat) <- nmf_genes
storage.mode(nmf_mat) <- "numeric"

keep_g <- rowSums(is.finite(nmf_mat)) >= max(3, floor(0.2 * ncol(nmf_mat)))
keep_s <- colSums(is.finite(nmf_mat)) >= max(3, floor(0.2 * nrow(nmf_mat)))
nmf_mat <- nmf_mat[keep_g, keep_s, drop = FALSE]

for (i in seq_len(nrow(nmf_mat))) {
  med <- median(nmf_mat[i, ], na.rm = TRUE)
  nmf_mat[i, !is.finite(nmf_mat[i, ])] <- med
}
nmf_scaled <- t(scale(t(nmf_mat)))
nmf_scaled[!is.finite(nmf_scaled)] <- 0
nmf_input <- nmf_scaled - min(nmf_scaled, na.rm = TRUE) + 1e-6

write_tsv(data.table::as.data.table(nmf_input, keep.rownames = "gene"),
          file.path(OUT$step13, "SIDISH_Step1_3_NMF_input_HR_pseudobulk_scaled_nonnegative.tsv"))

nmf_mu_single <- function(V, k, max_iter = 1000, tol = 1e-5) {
  n <- nrow(V); m <- ncol(V)
  W <- matrix(runif(n * k, 0.01, 1), nrow = n, ncol = k)
  H <- matrix(runif(k * m, 0.01, 1), nrow = k, ncol = m)
  eps <- 1e-9
  prev <- Inf
  for (iter in seq_len(max_iter)) {
    H <- H * ((t(W) %*% V) / (t(W) %*% W %*% H + eps))
    W <- W * ((V %*% t(H)) / (W %*% H %*% t(H) + eps))
    sc <- colSums(W) + eps
    W <- sweep(W, 2, sc, "/")
    H <- sweep(H, 1, sc, "*")
    if (iter %% 25 == 0 || iter == max_iter) {
      rss <- sum((V - W %*% H)^2)
      if (is.finite(prev) && abs(prev - rss) / (prev + eps) < tol) break
      prev <- rss
    }
  }
  list(W = W, H = H, rss = sum((V - W %*% H)^2), iter = iter)
}

nmf_best <- function(V, k, nrun) {
  best <- NULL
  for (r in seq_len(nrun)) {
    fit <- nmf_mu_single(V, k, CONFIG$nmf_max_iter, CONFIG$nmf_tol)
    if (is.null(best) || fit$rss < best$rss) best <- fit
  }
  best
}

rank_metrics <- data.table::data.table()
fits <- list()

for (k in CONFIG$k_range) {
  if (k >= min(nrow(nmf_input), ncol(nmf_input))) next
  msg("Running NMF k=", k)
  fit <- nmf_best(nmf_input, k, CONFIG$nmf_nrun)
  fits[[as.character(k)]] <- fit
  W <- fit$W
  rownames(W) <- rownames(nmf_input)
  colnames(W) <- paste0("MP", seq_len(k))
  assn <- apply(W, 1, which.max)
  comp_counts <- tabulate(assn, nbins = k)
  rank_metrics <- rbind(rank_metrics, data.table::data.table(
    k = k,
    rss = fit$rss,
    n_genes = nrow(nmf_input),
    n_patients = ncol(nmf_input),
    min_assigned_genes_per_MP = min(comp_counts),
    all_MP_pass_min_genes = all(comp_counts >= CONFIG$min_genes_per_mp)
  ))
}

if (nrow(rank_metrics) == 0) stop("NMF rank scan failed.")
rank_metrics <- rank_metrics[order(k)]
rank_metrics[, rss_relative := rss / max(rss, na.rm = TRUE)]
rank_metrics[, rss_improvement_from_previous := c(NA_real_, head(rss_relative, -1) - tail(rss_relative, -1))]
write_tsv(rank_metrics, file.path(OUT$step13, "SIDISH_Step1_3_NMF_rank_metrics.tsv"))

if (!is.na(CONFIG$selected_k_override)) {
  selected_k <- CONFIG$selected_k_override
  selection_reason <- "manual_selected_k_override"
} else {
  candidates <- rank_metrics[all_MP_pass_min_genes == TRUE]
  if (nrow(candidates) == 0) candidates <- rank_metrics
  candidates[, modest_gain := !is.na(rss_improvement_from_previous) &
               rss_improvement_from_previous < CONFIG$nmf_modest_gain_threshold]
  small_k_candidates <- candidates[k <= CONFIG$nmf_prefer_k_cap]
  if (nrow(small_k_candidates) > 0 && any(small_k_candidates$modest_gain)) {
    selected_k <- small_k_candidates[modest_gain == TRUE, k][1]
    selection_reason <- "automatic_balanced_first_modest_gain_with_small_k_preference"
  } else if (any(candidates$modest_gain)) {
    selected_k <- candidates[modest_gain == TRUE, k][1]
    selection_reason <- "automatic_first_modest_rss_gain"
  } else if (nrow(small_k_candidates) > 0) {
    selected_k <- small_k_candidates[which.min(rss), k]
    selection_reason <- "automatic_small_k_preference_lowest_rss"
  } else {
    selected_k <- candidates[which.min(rss), k]
    selection_reason <- "automatic_lowest_rss"
  }
}

fit <- fits[[as.character(selected_k)]]
W <- fit$W
H <- fit$H
rownames(W) <- rownames(nmf_input)
colnames(W) <- paste0("MP", seq_len(selected_k))
rownames(H) <- paste0("MP", seq_len(selected_k))
colnames(H) <- colnames(nmf_input)

loading_wide <- data.table::as.data.table(W, keep.rownames = "gene")
metaprogram_loadings <- data.table::melt(loading_wide, id.vars = "gene", variable.name = "metaprogram", value.name = "loading")
metaprogram_loadings[, loading_rank_within_MP := data.table::frank(-loading, ties.method = "first"), by = metaprogram]
metaprogram_loadings[, top_gene_for_MP := loading_rank_within_MP <= CONFIG$top_n_genes_per_mp]

assign_dt <- data.table::data.table(
  gene = rownames(W),
  assigned_MP = colnames(W)[apply(W, 1, which.max)],
  assigned_MP_loading = apply(W, 1, max)
)
metaprogram_loadings <- merge(metaprogram_loadings, assign_dt, by = "gene", all.x = TRUE)

write_tsv(metaprogram_loadings, file.path(OUT$step13, "SIDISH_Step1_3_metaprogram_gene_loadings.tsv"))
write_tsv(metaprogram_loadings[top_gene_for_MP == TRUE][order(metaprogram, loading_rank_within_MP)],
          file.path(OUT$step13, "SIDISH_Step1_3_top_genes_per_metaprogram.tsv"))
write_tsv(data.table::as.data.table(t(H), keep.rownames = "patient_id"),
          file.path(OUT$step13, "SIDISH_Step1_3_patient_MP_scores_from_NMF.tsv"))

selected_rank_summary <- data.table::data.table(
  selected_k = selected_k,
  selection_reason = selection_reason,
  n_primary_input_genes = length(primary_genes),
  n_NMF_genes = nrow(nmf_input),
  n_NMF_patients = ncol(nmf_input),
  k_range = paste(CONFIG$k_range, collapse = ","),
  nmf_nrun = CONFIG$nmf_nrun
)
write_tsv(selected_rank_summary, file.path(OUT$step13, "SIDISH_Step1_3_selected_rank_summary.tsv"))

###############################################################################
# 6. Step 1.4 Bulk transferability and eligibility
###############################################################################

msg("Step 1.4 Bulk transferability and clinical-model eligibility...")

standardize_bulk_expr <- function(file) {
  dt <- read_table_auto(file)
  gene_col <- names(dt)[1]
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

standardize_clin <- function(file) {
  dt <- read_table_auto(file)
  if (!CONFIG$sample_id_col %in% names(dt)) {
    sample_col <- first_existing_col(dt, c("sample_id", "sample", "case", "patient", "ID"), regex = TRUE)
    if (is.na(sample_col)) sample_col <- names(dt)[1]
  } else {
    sample_col <- CONFIG$sample_id_col
  }
  time_col <- if (CONFIG$time_col %in% names(dt)) CONFIG$time_col else first_existing_col(dt, c("os_time", "OS_time", "survival.*time", "follow.*time", "days", "months"), regex = TRUE)
  event_col <- if (CONFIG$event_col %in% names(dt)) CONFIG$event_col else first_existing_col(dt, c("os_event", "OS_event", "survival.*event", "status", "death", "event"), regex = TRUE)
  if (is.na(time_col) || is.na(event_col)) stop("Clinical missing survival columns: ", file)

  out <- data.table::data.table(
    sample_id = standardize_id(dt[[sample_col]]),
    OS_time = as_numeric_safe(dt[[time_col]]),
    OS_event_raw = dt[[event_col]]
  )
  ev <- toupper(as.character(out$OS_event_raw))
  out[, OS_event := data.table::fcase(
    ev %in% c("1", "DEAD", "DECEASED", "EVENT", "TRUE", "YES"), 1,
    ev %in% c("0", "ALIVE", "CENSORED", "FALSE", "NO", "LIVING"), 0,
    default = suppressWarnings(as.numeric(ev))
  )]
  out <- out[!is.na(sample_id) & sample_id != "" & is.finite(OS_time) & OS_time > 0 & OS_event %in% c(0, 1)]
  unique(out, by = "sample_id")
}

safe_cox <- function(x, clin) {
  out <- c(beta = NA_real_, HR = NA_real_, p = NA_real_, n = NA_real_, events = NA_real_)
  dat <- data.frame(x = as.numeric(x), OS_time = clin$OS_time, OS_event = clin$OS_event)
  dat <- dat[is.finite(dat$x) & is.finite(dat$OS_time) & dat$OS_time > 0 & dat$OS_event %in% c(0,1), ]
  out["n"] <- nrow(dat)
  out["events"] <- sum(dat$OS_event == 1)
  if (nrow(dat) < 20 || sum(dat$OS_event == 1) < 5 || sd(dat$x, na.rm = TRUE) <= 0) return(out)
  dat$x <- as.numeric(scale(dat$x))
  fit <- tryCatch(survival::coxph(survival::Surv(OS_time, OS_event) ~ x, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(out)
  sm <- summary(fit)
  out["beta"] <- as.numeric(stats::coef(fit)[1])
  out["HR"] <- exp(out["beta"])
  out["p"] <- sm$coefficients[1, "Pr(>|z|)"]
  out
}

bulk_gene_long <- data.table::data.table()
bulk_mp_scores <- data.table::data.table()
bulk_read_status <- data.table::data.table()

top_weights <- metaprogram_loadings[top_gene_for_MP == TRUE & !is.na(loading)]

for (ii in seq_len(nrow(bulk_manifest))) {
  coh <- bulk_manifest$cohort[ii]
  expr_file <- bulk_manifest$expr_file[ii]
  clin_file <- bulk_manifest$clin_file[ii]

  status_row <- data.table::data.table(
    cohort = coh,
    expr_file = expr_file,
    clin_file = clin_file,
    expr_exists = bulk_manifest$expr_exists[ii],
    clin_exists = bulk_manifest$clin_exists[ii],
    read_status = NA_character_,
    selected_sample_match_mode = NA_character_,
    n_matched_samples = NA_integer_,
    n_marker_genes_in_bulk = NA_integer_,
    error_message = NA_character_
  )

  if (!isTRUE(bulk_manifest$expr_exists[ii]) || !isTRUE(bulk_manifest$clin_exists[ii])) {
    msg("Skipping ", coh, ": missing expr/clin after discovery.")
    status_row$read_status <- "missing_expr_or_clin_after_discovery"
    bulk_read_status <- rbind(bulk_read_status, status_row, fill = TRUE)
    next
  }

  require_input_file(expr_file, paste0(coh, " bulk expr"))
  require_input_file(clin_file, paste0(coh, " clinical"))

  msg("Reading bulk cohort: ", coh)
  expr <- tryCatch(standardize_bulk_expr(expr_file), error = function(e) {
    status_row$read_status <<- "failed_to_read_expression"
    status_row$error_message <<- conditionMessage(e)
    NULL
  })
  if (is.null(expr)) {
    bulk_read_status <- rbind(bulk_read_status, status_row, fill = TRUE)
    next
  }

  clin <- tryCatch(standardize_clin(clin_file), error = function(e) {
    status_row$read_status <<- "failed_to_read_clinical"
    status_row$error_message <<- conditionMessage(e)
    NULL
  })
  if (is.null(clin)) {
    bulk_read_status <- rbind(bulk_read_status, status_row, fill = TRUE)
    next
  }

  aligned <- align_bulk_expr_clin(expr, clin, cohort = coh, out_dir = OUT$step14)
  status_row$selected_sample_match_mode <- aligned$selected_mode
  status_row$n_matched_samples <- aligned$n_common

  if (is.null(aligned$expr) || is.null(aligned$clin) || aligned$n_common < 20) {
    msg("Skipping ", coh, ": too few matched samples after auto sample matching: ", aligned$n_common,
        " (mode=", aligned$selected_mode, ")")
    status_row$read_status <- "too_few_matched_samples_after_auto_matching"
    bulk_read_status <- rbind(bulk_read_status, status_row, fill = TRUE)
    next
  }

  expr <- aligned$expr
  clin <- aligned$clin

  common_genes <- intersect(marker_genes, rownames(expr))
  status_row$n_marker_genes_in_bulk <- length(common_genes)

  if (length(common_genes) < 5) {
    msg("Skipping ", coh, ": too few marker genes in bulk: ", length(common_genes))
    status_row$read_status <- "too_few_marker_genes_in_bulk"
    bulk_read_status <- rbind(bulk_read_status, status_row, fill = TRUE)
    next
  }

  for (g in common_genes) {
    v <- as.numeric(expr[g, ])
    cox <- safe_cox(v, clin)
    bulk_gene_long <- rbind(bulk_gene_long, data.table::data.table(
      cohort = coh,
      gene = g,
      detected_fraction = mean(is.finite(v) & v > 0, na.rm = TRUE),
      expression_variance = var(v, na.rm = TRUE),
      bulk_detectable = mean(is.finite(v) & v > 0, na.rm = TRUE) >= CONFIG$bulk_min_detected_fraction & var(v, na.rm = TRUE) > CONFIG$bulk_min_variance,
      cox_beta = cox["beta"],
      cox_HR = cox["HR"],
      cox_p = cox["p"],
      cox_n = cox["n"],
      cox_events = cox["events"]
    ), fill = TRUE)
  }

  for (mp in unique(top_weights$metaprogram)) {
    ww <- top_weights[metaprogram == mp]
    gg <- intersect(ww$gene, rownames(expr))
    if (length(gg) < CONFIG$min_genes_per_mp) next
    w <- ww[match(gg, gene), loading]
    mat <- expr[gg, , drop = FALSE]
    mat <- t(scale(t(mat)))
    mat[!is.finite(mat)] <- 0
    score <- as.numeric(crossprod(w / sum(abs(w)), mat))
    bulk_mp_scores <- rbind(bulk_mp_scores, data.table::data.table(
      cohort = coh,
      sample_id = colnames(expr),
      metaprogram = mp,
      MP_score = score,
      n_genes_used = length(gg)
    ), fill = TRUE)
  }

  status_row$read_status <- "success"
  bulk_read_status <- rbind(bulk_read_status, status_row, fill = TRUE)
}

write_tsv(bulk_read_status, file.path(OUT$step14, "SIDISH_Step1_4_bulk_cohort_read_status.tsv"))

write_tsv(bulk_gene_long, file.path(OUT$step14, "SIDISH_Step1_4_bulk_gene_transferability_long.tsv"))
write_tsv(bulk_mp_scores, file.path(OUT$step14, "SIDISH_Step1_4_bulk_MP_scores_by_sample.tsv"))

if (nrow(bulk_gene_long) == 0) {
  stop("No bulk transferability results were generated. Inspect these audit files:\n",
       "1) ", file.path(OUT$audit, "selected_bulk_input_files.tsv"), "\n",
       "2) ", file.path(OUT$audit, "bulk_candidate_file_inventory.tsv"), "\n",
       "3) ", file.path(OUT$step14, "SIDISH_Step1_4_bulk_cohort_read_status.tsv"), "\n",
       "4) cohort-specific sample_match_mode_report.tsv files under ", OUT$step14, "\n",
       "Most likely causes: expression/clinical files were not placed under the expected bulk folders, ",
       "or sample IDs require TCGA barcode trimming / CGGA ID normalization.")
}

bulk_gene_summary <- bulk_gene_long[, .(
  n_bulk_cohorts_tested = data.table::uniqueN(cohort),
  n_bulk_cohorts_detectable = sum(bulk_detectable, na.rm = TRUE),
  mean_detected_fraction = mean(detected_fraction, na.rm = TRUE),
  min_detected_fraction = min(detected_fraction, na.rm = TRUE),
  mean_expression_variance = mean(expression_variance, na.rm = TRUE),
  n_cohorts_with_valid_cox = sum(!is.na(cox_beta)),
  n_cohorts_beta_positive = sum(!is.na(cox_beta) & cox_beta > 0),
  n_cohorts_beta_negative = sum(!is.na(cox_beta) & cox_beta < 0),
  mean_cox_beta = mean(cox_beta, na.rm = TRUE),
  median_cox_beta = median(cox_beta, na.rm = TRUE),
  min_cox_beta = suppressWarnings(min(cox_beta, na.rm = TRUE)),
  max_cox_beta = suppressWarnings(max(cox_beta, na.rm = TRUE)),
  min_cox_p = suppressWarnings(min(cox_p, na.rm = TRUE))
), by = gene]
bulk_gene_summary[is.infinite(min_cox_p), min_cox_p := NA_real_]
bulk_gene_summary[is.infinite(min_cox_beta), min_cox_beta := NA_real_]
bulk_gene_summary[is.infinite(max_cox_beta), max_cox_beta := NA_real_]
bulk_gene_summary[, positive_cohort_fraction := data.table::fifelse(
  n_cohorts_with_valid_cox > 0,
  n_cohorts_beta_positive / n_cohorts_with_valid_cox,
  NA_real_
)]
bulk_gene_summary[, bulk_transferable := n_bulk_cohorts_detectable >= CONFIG$bulk_min_available_cohorts]
bulk_gene_summary[, A1_bulk_direction_pass := (
  n_cohorts_with_valid_cox >= 2 &
    n_cohorts_beta_positive >= 2 &
    positive_cohort_fraction >= CONFIG$A1_min_positive_cohort_fraction &
    (!isTRUE(CONFIG$A1_require_mean_beta_positive) | (!is.na(mean_cox_beta) & mean_cox_beta > CONFIG$A1_min_mean_cox_beta)) &
    (!isTRUE(CONFIG$A1_require_median_beta_positive) | (!is.na(median_cox_beta) & median_cox_beta > CONFIG$A1_min_median_cox_beta))
)]
bulk_gene_summary[, A1_borderline_direction_flag := (
  n_cohorts_with_valid_cox >= 2 &
    n_cohorts_beta_positive >= 2 &
    positive_cohort_fraction >= CONFIG$A1_min_positive_cohort_fraction &
    !A1_bulk_direction_pass
)]
bulk_gene_summary[, bulk_direction := data.table::fcase(
  A1_bulk_direction_pass == TRUE, "bulk_risk_positive_consistent_strict",
  A1_borderline_direction_flag == TRUE, "bulk_borderline_majority_positive_but_mean_or_median_not_positive",
  n_cohorts_with_valid_cox >= 2 & n_cohorts_beta_negative >= 2, "bulk_negative_consistent",
  n_cohorts_with_valid_cox >= 1, "bulk_mixed_or_single_cohort",
  default = "bulk_survival_not_evaluable"
)]
write_tsv(bulk_gene_summary, file.path(OUT$step14, "SIDISH_Step1_4_bulk_gene_transferability_summary.tsv"))

# Eligibility.
mp_gene_best <- metaprogram_loadings[, .SD[which.max(abs(loading))], by = gene]
mp_gene_best <- mp_gene_best[, .(gene, primary_metaprogram = metaprogram, primary_MP_loading = loading, assigned_MP, assigned_MP_loading)]

eligibility <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), list(
  marker_source_master,
  robustness_master[, .(gene, robustness_tier, patient_median_delta, patient_support_fraction, n_paired_patients, zero_fraction, patient_direction, primary_MP_input_flag, sensitivity_MP_input_flag)],
  mp_gene_best,
  bulk_gene_summary
))
eligibility[is.na(primary_MP_input_flag), primary_MP_input_flag := FALSE]
eligibility[is.na(sensitivity_MP_input_flag), sensitivity_MP_input_flag := FALSE]
eligibility[is.na(bulk_transferable), bulk_transferable := FALSE]
eligibility[is.na(n_bulk_cohorts_detectable), n_bulk_cohorts_detectable := 0L]
eligibility[is.na(n_cohorts_with_valid_cox), n_cohorts_with_valid_cox := 0L]
eligibility[is.na(bulk_direction), bulk_direction := "bulk_not_evaluable"]

eligibility[, clinical_model_class := data.table::fcase(
  primary_MP_input_flag == TRUE & bulk_transferable == TRUE & A1_bulk_direction_pass == TRUE,
  "Class_A1_risk_directional_translational_candidate",

  primary_MP_input_flag == TRUE & bulk_transferable == TRUE & A1_borderline_direction_flag == TRUE,
  CONFIG$borderline_A2_label,

  primary_MP_input_flag == TRUE & bulk_transferable == TRUE & bulk_direction %in% c(
    "bulk_negative_consistent",
    "bulk_mixed_or_single_cohort",
    "bulk_survival_not_evaluable",
    "bulk_borderline_majority_positive_but_mean_or_median_not_positive"
  ),
  "Class_A2_signed_protective_or_subtype_associated_candidate",

  sensitivity_MP_input_flag == TRUE & !is.na(primary_metaprogram) &
    (bulk_transferable == FALSE | bulk_direction %in% c("bulk_mixed_or_single_cohort", "bulk_borderline_majority_positive_but_mean_or_median_not_positive")),
  "Class_B_biological_HR_program_mixed_bulk_direction",

  bulk_transferable == TRUE & primary_MP_input_flag == FALSE,
  "Class_C_bulk_supportive_not_primary_SIDISH_core",

  default = "Class_D_excluded_or_annotation_only"
)]

eligibility[, enter_Step2_ML := clinical_model_class %in% c(
  "Class_A1_risk_directional_translational_candidate",
  "Class_A2_signed_protective_or_subtype_associated_candidate",
  CONFIG$borderline_A2_label
)]
eligibility[, Step2_feature_sign := data.table::fcase(
  clinical_model_class == "Class_A1_risk_directional_translational_candidate", "risk_positive",
  clinical_model_class == CONFIG$borderline_A2_label, "signed_context_dependent_borderline",
  clinical_model_class == "Class_A2_signed_protective_or_subtype_associated_candidate" & bulk_direction == "bulk_negative_consistent", "protective_or_subtype_negative",
  clinical_model_class == "Class_A2_signed_protective_or_subtype_associated_candidate", "signed_context_dependent",
  default = "not_for_primary_ML"
)]
eligibility[, eligibility_rationale := paste(
  "source=", attribution_category,
  "; robustness=", robustness_tier,
  "; MP=", ifelse(is.na(primary_metaprogram), "NA", primary_metaprogram),
  "; bulk_transferable=", bulk_transferable,
  "; bulk_direction=", bulk_direction,
  "; n_pos=", ifelse(is.na(n_cohorts_beta_positive), "NA", n_cohorts_beta_positive),
  "/", ifelse(is.na(n_cohorts_with_valid_cox), "NA", n_cohorts_with_valid_cox),
  "; mean_beta=", ifelse(is.na(mean_cox_beta), "NA", sprintf("%.4g", mean_cox_beta)),
  "; median_beta=", ifelse(is.na(median_cox_beta), "NA", sprintf("%.4g", median_cox_beta)),
  sep = ""
)]

write_tsv(eligibility, file.path(OUT$step14, "SIDISH_Step1_4_clinical_model_eligibility_master.tsv"))

class_counts <- eligibility[, .(
  n_genes = .N,
  n_enter_Step2_ML = sum(enter_Step2_ML, na.rm = TRUE)
), by = clinical_model_class][order(clinical_model_class)]
write_tsv(class_counts, file.path(OUT$step14, "SIDISH_Step1_4_class_counts.tsv"))

borderline_audit <- eligibility[
  clinical_model_class == CONFIG$borderline_A2_label |
    A1_borderline_direction_flag == TRUE,
  .(
    gene,
    clinical_model_class,
    Step2_feature_sign,
    n_cohorts_with_valid_cox,
    n_cohorts_beta_positive,
    n_cohorts_beta_negative,
    positive_cohort_fraction,
    mean_cox_beta,
    median_cox_beta,
    min_cox_beta,
    max_cox_beta,
    bulk_direction,
    eligibility_rationale
  )
][order(gene)]
write_tsv(borderline_audit, file.path(OUT$step14, "SIDISH_Step1_4_A1_borderline_reclassification_audit.tsv"))

classA <- eligibility[enter_Step2_ML == TRUE][order(clinical_model_class, gene)]
write_tsv(classA, file.path(OUT$final, "SIDISH_Step1_ClassA_signed_core_marker_pool_for_Step2_ML.tsv"))

master <- eligibility[order(clinical_model_class, gene)]
write_tsv(master, file.path(OUT$final, "SIDISH_Step1_Core_marker_refinement_master_table.tsv"))

summary <- data.table::data.table(
  section = c("Input", "Step1.1", "Step1.2", "Step1.2", "Step1.3", "Step1.4", "Step1.4", "Final"),
  metric = c(
    "input_mode",
    "n_marker_universe",
    "n_valid_patients_for_pseudobulk",
    "n_primary_MP_input_genes",
    "selected_k",
    "n_bulk_cohorts_with_results",
    "n_A1_borderline_reclassified_to_A2",
    "n_ClassA_signed_core_marker_pool"
  ),
  value = c(
    "flat_files_no_h5ad",
    nrow(marker_source_master),
    length(valid_patients),
    sum(robustness_master$primary_MP_input_flag, na.rm = TRUE),
    selected_k,
    data.table::uniqueN(bulk_gene_long$cohort),
    nrow(borderline_audit),
    nrow(classA)
  )
)
write_tsv(summary, file.path(OUT$final, "SIDISH_Step1_Core_marker_refinement_final_summary.tsv"))

if (nrow(classA) == 0) {
  stop("Class A pool is empty. This is not accepted as a completed result. Inspect Step1.2 and Step1.4 audit files.")
}


###############################################################################
# 7. Publication-quality figures (Nature-style consistent aesthetics, v3.4)
###############################################################################

if (isTRUE(CONFIG$make_publication_figures)) {
  msg("Generating publication-quality figures...")

  # Compact label dictionaries for publication figures.
  attribution_label_map <- c(
    single_celltype_associated = "Single-celltype\nassociated",
    shared_multi_celltype_associated = "Shared / multi-celltype\nassociated",
    unassigned_no_celltype_HR_support = "No celltype\nsupport"
  )
  tier_label_map <- c(
    Step1_2A_strict_patient_robust = "Strict patient-\nrobust",
    Step1_2B_direction_consistent_supportive = "Direction-consistent\nsupportive",
    Step1_2C_zero_inflated_supportive = "Zero-inflated\nsupportive",
    Step1_2D_hold_exclude_not_patient_robust = "Hold / exclude:\nnot robust",
    Step1_2D_hold_exclude_not_matched_or_not_paired = "Hold / exclude:\nnot paired",
    Step1_2D_hold_exclude_insufficient_paired_patients = "Hold / exclude:\ninsufficient patients"
  )
  class_label_map <- c(
    Class_A1_risk_directional_translational_candidate = "A1 risk-positive\ntranslational",
    Class_A2_signed_protective_or_subtype_associated_candidate = "A2 signed / protective\nor subtype-associated",
    Class_A2_borderline_context_dependent_candidate = "A2 borderline\ncontext-dependent",
    Class_B_biological_HR_program_mixed_bulk_direction = "B biological HR\nmixed bulk",
    Class_C_bulk_supportive_not_primary_SIDISH_core = "C bulk-supportive\nnot primary",
    Class_D_excluded_or_annotation_only = "D excluded /\nannotation-only"
  )

  # Fig 1A: marker attribution categories.
  fig1a_dt <- marker_source_master[, .N, by = attribution_category][order(-N)]
  fig1a_dt[, label := attribution_label_map[attribution_category]]
  fig1a_dt[is.na(label), label := attribution_category]
  fig1a_cols <- c(
    "Single-celltype\nassociated" = pal("blue"),
    "Shared / multi-celltype\nassociated" = pal("green"),
    "No celltype\nsupport" = pal("grey")
  )
  p_fig1a <- ggplot2::ggplot(fig1a_dt, ggplot2::aes(x = stats::reorder(label, N), y = N, fill = label)) +
    ggplot2::geom_col(width = 0.70, color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_text(ggplot2::aes(label = N), hjust = -0.15, family = CONFIG$plot_base_family, size = 2.7) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::expand_limits(y = max(fig1a_dt$N) * 1.14) +
    ggplot2::scale_fill_manual(values = fig1a_cols, guide = "none") +
    ggplot2::labs(title = "Marker provenance", subtitle = "245 SIDISH-HR marker universe",
                  x = NULL, y = "Number of markers") +
    pub_theme()
  save_pub_plot(p_fig1a, "Figure_1A_marker_attribution_categories", width = 4.5, height = 3.0)

  # Fig 1B: top cell types supporting markers.
  fig1b_dt <- celltype_markers[, .(n_markers = data.table::uniqueN(gene)), by = celltype][order(-n_markers)]
  fig1b_dt <- fig1b_dt[1:min(15, .N)]
  fig1b_dt[, celltype_label := wrap_axis_labels(celltype, 18)]
  p_fig1b <- ggplot2::ggplot(fig1b_dt, ggplot2::aes(x = stats::reorder(celltype_label, n_markers), y = n_markers)) +
    ggplot2::geom_col(width = 0.70, fill = pal("navy"), color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_text(ggplot2::aes(label = n_markers), hjust = -0.12, family = CONFIG$plot_base_family, size = 2.6) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::expand_limits(y = max(fig1b_dt$n_markers) * 1.12) +
    ggplot2::labs(title = "Cell-type support", subtitle = "Top contributors to HR-marker evidence",
                  x = NULL, y = "Supported markers") +
    pub_theme()
  save_pub_plot(p_fig1b, "Figure_1B_top_celltypes_supporting_markers", width = 4.8, height = 4.2)

  # Fig 2A: robustness tier counts.
  fig2a_dt <- data.table::copy(tier_counts)
  fig2a_dt[, label := tier_label_map[robustness_tier]]
  fig2a_dt[is.na(label), label := robustness_tier]
  fig2a_cols <- c(
    "Strict patient-\nrobust" = pal("red"),
    "Direction-consistent\nsupportive" = pal("orange"),
    "Zero-inflated\nsupportive" = pal("yellow"),
    "Hold / exclude:\nnot robust" = pal("grey"),
    "Hold / exclude:\nnot paired" = pal("beige"),
    "Hold / exclude:\ninsufficient patients" = pal("brown")
  )
  p_fig2a <- ggplot2::ggplot(fig2a_dt, ggplot2::aes(x = stats::reorder(label, n_genes), y = n_genes, fill = label)) +
    ggplot2::geom_col(width = 0.70, color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_text(ggplot2::aes(label = n_genes), hjust = -0.12, family = CONFIG$plot_base_family, size = 2.6) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::expand_limits(y = max(fig2a_dt$n_genes) * 1.12) +
    ggplot2::scale_fill_manual(values = fig2a_cols, guide = "none") +
    ggplot2::labs(title = "Patient-level robustness", subtitle = "De novo HR/BG pseudobulk audit",
                  x = NULL, y = "Number of markers") +
    pub_theme()
  save_pub_plot(p_fig2a, "Figure_2A_patient_level_robustness_tiers", width = 4.6, height = 3.4)

  # Fig 2B: top primary genes by patient median delta.
  fig2b_dt <- robustness_master[primary_MP_input_flag == TRUE & !is.na(patient_median_delta),
                                .(gene, patient_median_delta, patient_support_fraction, robustness_tier)]
  if (nrow(fig2b_dt) > 0) {
    fig2b_dt <- fig2b_dt[order(-patient_median_delta)][1:min(20, .N)]
    fig2b_dt[, tier_short := ifelse(robustness_tier == "Step1_2A_strict_patient_robust", "Strict", "Supportive")]
    p_fig2b <- ggplot2::ggplot(fig2b_dt, ggplot2::aes(x = stats::reorder(gene, patient_median_delta), y = patient_median_delta,
                                                       fill = patient_support_fraction)) +
      ggplot2::geom_col(width = 0.70, color = "black", linewidth = CONFIG$plot_linewidth) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient(low = pal("teal"), high = pal("red")) +
      ggplot2::labs(title = "Primary pseudobulk-supported markers", subtitle = "Top 20 by median HR-BG delta",
                    x = NULL, y = "Median HR-BG delta", fill = "Support\nfraction") +
      pub_theme()
    save_pub_plot(p_fig2b, "Figure_2B_top_primary_markers_by_patient_delta", width = 4.7, height = 4.3)
  }

  # Fig 3A: NMF rank scan.
  p_fig3a <- ggplot2::ggplot(rank_metrics, ggplot2::aes(x = k, y = rss_relative)) +
    ggplot2::geom_line(color = pal("navy"), linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_point(size = CONFIG$plot_point_size + 0.2, shape = 21, stroke = CONFIG$plot_linewidth,
                        fill = ifelse(rank_metrics$k == selected_k, pal("red"), "white"),
                        color = pal("navy")) +
    ggplot2::geom_vline(xintercept = selected_k, linetype = "dashed", linewidth = CONFIG$plot_linewidth, color = pal("red")) +
    ggplot2::scale_x_continuous(breaks = rank_metrics$k) +
    ggplot2::labs(title = "NMF rank selection", subtitle = paste0("Selected k = ", selected_k, " | ", selection_reason),
                  x = "NMF rank (k)", y = "Relative RSS") +
    pub_theme()
  save_pub_plot(p_fig3a, "Figure_3A_NMF_rank_scan", width = 4.4, height = 3.2)

  # Fig 3B: top-gene loadings across metaprograms.
  fig3b_dt <- metaprogram_loadings[top_gene_for_MP == TRUE][order(metaprogram, loading_rank_within_MP)]
  if (nrow(fig3b_dt) > 0) {
    fig3b_dt <- fig3b_dt[loading_rank_within_MP <= min(8, CONFIG$top_n_genes_per_mp)]
    fig3b_dt[, gene_mp := paste0(gene, "  ")]
    fig3b_dt[, gene_mp := factor(gene_mp, levels = rev(unique(gene_mp)))]
    p_fig3b <- ggplot2::ggplot(fig3b_dt, ggplot2::aes(x = metaprogram, y = gene_mp, fill = loading)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.18) +
      ggplot2::scale_fill_gradientn(colors = c("#F7F7F7", pal("orange"), pal("red"))) +
      ggplot2::labs(title = "Metaprogram loading structure", subtitle = "Top loading genes per selected MP",
                    x = NULL, y = NULL, fill = "Loading") +
      pub_theme() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    save_pub_plot(p_fig3b, "Figure_3B_metaprogram_top_gene_heatmap", width = 5.4, height = max(4.6, 0.16 * nrow(fig3b_dt) + 1.5))
  }

  # Fig 4A: clinical-model eligibility classes.
  fig4a_dt <- data.table::copy(class_counts)
  fig4a_dt[, label := class_label_map[clinical_model_class]]
  fig4a_dt[is.na(label), label := clinical_model_class]
  fig4a_cols <- c(
    "A1 risk-positive\ntranslational" = pal("red"),
    "A2 signed / protective\nor subtype-associated" = pal("blue"),
    "A2 borderline\ncontext-dependent" = pal("purple"),
    "B biological HR\nmixed bulk" = pal("orange"),
    "C bulk-supportive\nnot primary" = pal("green"),
    "D excluded /\nannotation-only" = pal("grey")
  )
  p_fig4a <- ggplot2::ggplot(fig4a_dt, ggplot2::aes(x = stats::reorder(label, n_genes), y = n_genes, fill = label)) +
    ggplot2::geom_col(width = 0.70, color = "black", linewidth = CONFIG$plot_linewidth) +
    ggplot2::geom_text(ggplot2::aes(label = n_genes), hjust = -0.12, family = CONFIG$plot_base_family, size = 2.6) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::expand_limits(y = max(fig4a_dt$n_genes) * 1.14) +
    ggplot2::scale_fill_manual(values = fig4a_cols, guide = "none") +
    ggplot2::labs(title = "Clinical-model eligibility", subtitle = "Only signed Class A genes enter Step 2 ML",
                  x = NULL, y = "Number of markers") +
    pub_theme()
  save_pub_plot(p_fig4a, "Figure_4A_clinical_model_eligibility_classes", width = 4.8, height = 4.0)

  # Fig 4B: Class A bulk Cox direction as a lollipop plot.
  fig4b_dt <- eligibility[enter_Step2_ML == TRUE & !is.na(mean_cox_beta),
                          .(gene, clinical_model_class, mean_detected_fraction, mean_cox_beta)]
  if (nrow(fig4b_dt) > 0) {
    fig4b_dt[, class_short := data.table::fcase(grepl("Class_A1", clinical_model_class), "A1 risk-positive",
                              clinical_model_class == CONFIG$borderline_A2_label, "A2 borderline",
                              default = "A2 signed-negative")]
    fig4b_dt[, gene := factor(gene, levels = gene[order(mean_cox_beta)])]
    p_fig4b <- ggplot2::ggplot(fig4b_dt, ggplot2::aes(y = gene, x = mean_cox_beta, color = class_short)) +
      ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = CONFIG$plot_linewidth) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = mean_cox_beta, yend = gene), linewidth = CONFIG$plot_linewidth) +
      ggplot2::geom_point(size = CONFIG$plot_point_size + 0.2) +
      ggplot2::scale_color_manual(values = c("A1 risk-positive" = pal("red"), "A2 signed-negative" = pal("blue"), "A2 borderline" = pal("purple"))) +
      ggplot2::labs(title = "Bulk survival direction of the Class A pool",
                    subtitle = "Mean standardized Cox beta across available bulk cohorts",
                    x = "Mean Cox beta", y = NULL, color = "Class") +
      pub_theme()
    save_pub_plot(p_fig4b, "Figure_4B_ClassA_bulk_cox_direction_lollipop", width = 5.2, height = max(4.8, 0.13 * nrow(fig4b_dt) + 1.6))
  }

  # Fig 4C: compact Class A marker-pool map.
  fig4c_dt <- classA[, .(gene, clinical_model_class, Step2_feature_sign, robustness_tier)]
  if (nrow(fig4c_dt) > 0) {
    fig4c_dt[, class_short := data.table::fcase(grepl("Class_A1", clinical_model_class), "A1 risk-positive",
                              clinical_model_class == CONFIG$borderline_A2_label, "A2 borderline",
                              default = "A2 signed-negative")]
    fig4c_dt[, tier_short := ifelse(robustness_tier == "Step1_2A_strict_patient_robust", "Strict", "Supportive")]
    fig4c_dt[, gene := factor(gene, levels = rev(gene))]
    p_fig4c <- ggplot2::ggplot(fig4c_dt, ggplot2::aes(x = class_short, y = gene, fill = Step2_feature_sign)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.20, width = 0.86, height = 0.86) +
      ggplot2::scale_fill_manual(values = c(
        risk_positive = pal("red"),
        protective_or_subtype_negative = pal("blue"),
        signed_context_dependent = pal("orange"),
        signed_context_dependent_borderline = pal("purple"),
        not_for_primary_ML = pal("grey")
      )) +
      ggplot2::labs(title = "Signed Class A core-marker pool", subtitle = "Input feature pool for Step 2 cohort-aware ML",
                    x = NULL, y = NULL, fill = "Feature sign") +
      pub_theme() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    save_pub_plot(p_fig4c, "Figure_4C_signed_ClassA_core_marker_pool", width = 4.9, height = max(5.0, 0.13 * nrow(fig4c_dt) + 1.8))
  }
}

readme <- c(
  "SIDISH Step 1 Core Marker Refinement v3.5",
  "",
  "This version does not read h5ad.
It also generates publication-style figures using an Arial-first theme, linewidth 0.5, and a unified Nature-like palette. PDF export uses cairo_pdf first and falls back to sans if Arial is unavailable in the local R graphics device.",
  "It reads explicit flat files from Python SIDISH workflow: MTX/genes/cells/metadata/cell_scores/marker/bulk clinical files.",
  "",
  "v3.5 hardens A1 assignment: majority-positive Cox direction must also have positive mean and median beta. Borderline majority-positive genes are retained in Class A but reclassified as A2-borderline/context-dependent.

Key outputs:",
  file.path(OUT$final, "SIDISH_Step1_Core_marker_refinement_master_table.tsv"),
  file.path(OUT$final, "SIDISH_Step1_ClassA_signed_core_marker_pool_for_Step2_ML.tsv")
)
write_text(readme, file.path(OUT$root, "README_SIDISH_Step1_Core_marker_refinement_v3_5_flat_inputs.txt"))
write_text(capture.output(sessionInfo()), file.path(OUT$audit, "sessionInfo.txt"))

msg("SIDISH Step 1 v3.5 flat-input reanalysis completed.")
msg("Final output: ", file.path(OUT$final, "SIDISH_Step1_ClassA_signed_core_marker_pool_for_Step2_ML.tsv"))
