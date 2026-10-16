#!/usr/bin/env Rscript
###############################################################################
# SIDISH Step 4.5: Clinical incremental value and translational utility audit v1.4
#
# Goal: quantify whether the locked SIDISH risk score adds prognostic value beyond
# standard clinical-molecular variables. This script never changes the locked
# 9-gene model, coefficients, or cutoff.
#
# Core comparison:
#   Clinical-only      = Age + WHO grade + IDH + 1p19q + Sex, as available
#   SIDISH-only        = locked SIDISH risk score
#   Clinical + SIDISH  = Clinical-only + locked SIDISH risk score
#
# Main outputs:
#   ΔC-index, LRT, AIC/BIC, time-dependent AUC, Brier/IBS-like error,
#   NRI/IDI, clinical intermediate-risk rescue, DCA-style net benefit,
#   and clinical impact curve.
###############################################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = 3600)
set.seed(20260522)

CONFIG <- list(
  output_dir = "G:/NEW/The2/SIDISH_Step4_5_Clinical_Incremental_Value_Audit_v1_4",
  step4_merged_file = "G:/NEW/The2/SIDISH_Step4_Clinical_Independence_Utility_Audit_v1_1/01_clinical_merged_data/Step4_clinical_merged_locked_risk_dataset.tsv",
  step2_coef_file = "G:/NEW/The2/SIDISH_Step2_ML175_DirectionStable_Model_Development_v2_6_strict_all3/99_locked_model/SIDISH_Step2_locked_model_coefficients.tsv",
  clinical_covariate_priority = c("Age_z", "grade_binary_G4", "IDH_binary", "1p19q_status", "Sex"),
  cohort_order = c("TCGA", "CGGA325", "CGGA693", "GSE16011"),
  include_development_pooled = TRUE,
  include_all_pooled = TRUE,
  run_development_trained_external_application = TRUE,
  min_complete_cases = 50,
  min_events = 20,
  min_factor_level_n = 10,
  time_unit = "auto",
  eval_times_days = c(365, 1095, 1825),
  eval_times_months = c(12, 36, 60),
  run_bootstrap_delta_cindex = TRUE,
  bootstrap_n = 200,
  intermediate_lower_quantile = 1/3,
  intermediate_upper_quantile = 2/3,
  dca_thresholds = seq(0.05, 0.80, by = 0.05),
  make_publication_figures = TRUE,
  plot_base_family = "Arial",
  plot_base_size = 8,
  plot_linewidth = 0.5,
  plot_point_size = 2.0,
  p_value_floor_for_display = 2.2e-16,
  objective_reporting_note = "raw_incremental_value_metrics_no_qualitative_inflation"
)

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
norm_path <- function(x) gsub("\\\\", "/", as.character(x))
safe_dir_create <- function(x) { if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE); invisible(x) }

OUT <- list(
  root = CONFIG$output_dir,
  audit = file.path(CONFIG$output_dir, "00_audit"),
  data = file.path(CONFIG$output_dir, "01_model_comparison_data"),
  discrimination = file.path(CONFIG$output_dir, "02_incremental_discrimination"),
  prediction_error = file.path(CONFIG$output_dir, "03_prediction_error_calibration"),
  reclassification = file.path(CONFIG$output_dir, "04_reclassification_intermediate_rescue"),
  utility = file.path(CONFIG$output_dir, "05_clinical_utility"),
  fig = file.path(CONFIG$output_dir, "06_publication_figures"),
  final = file.path(CONFIG$output_dir, "99_final_incremental_value_audit")
)
invisible(lapply(OUT, safe_dir_create))

load_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing required R package: ", pkg, ". Install with install.packages('", pkg, "')")
  invisible(TRUE)
}
load_or_stop("data.table"); load_or_stop("survival"); load_or_stop("ggplot2")
suppressPackageStartupMessages(library(survival))
Surv <- survival::Surv

# v1.1 self-check: catch accidental non-exported data.table constructor typo.
if (!"data.table" %in% getNamespaceExports("data.table")) {
  stop("data.table package does not export data.table(). Check package installation.")
}

write_tsv <- function(dt, file) {
  safe_dir_create(dirname(file))
  data.table::fwrite(data.table::as.data.table(dt), file, sep = "\t", quote = FALSE, na = "NA")
}
write_text <- function(x, file) { safe_dir_create(dirname(file)); writeLines(x, file, useBytes = TRUE) }
clean_names <- function(x) { x <- trimws(as.character(x)); x <- gsub("\\s+", "_", x); gsub("[.]+", "_", x) }
read_table_auto <- function(file) {
  file <- norm_path(file)
  if (!file.exists(file)) stop("File not found: ", file)
  dt <- data.table::fread(file, sep = "auto", data.table = TRUE, showProgress = FALSE, quote = "", fill = TRUE)
  data.table::setnames(dt, clean_names(names(dt)))
  dt
}
as_numeric_safe <- function(x) { if (is.numeric(x)) return(x); suppressWarnings(as.numeric(gsub(",", "", as.character(x)))) }
safe_term <- function(x) if (grepl("^[A-Za-z.][A-Za-z0-9_.]*$", x)) x else paste0("`", gsub("`", "", x), "`")
make_formula <- function(predictors) stats::as.formula(paste0("survival::Surv(OS_time, OS_event) ~ ", paste(vapply(predictors, safe_term, character(1)), collapse = " + ")))
format_p_for_plot <- function(p, floor_value = CONFIG$p_value_floor_for_display) {
  p <- suppressWarnings(as.numeric(p)); if (!is.finite(p)) return("p = NA")
  if (p < floor_value) return(paste0("p < ", format(floor_value, scientific = TRUE, digits = 2)))
  if (p < 0.001) return(paste0("p = ", format(p, scientific = TRUE, digits = 2)))
  paste0("p = ", sprintf("%.3f", p))
}

detect_survival_time_unit <- function(time, configured = CONFIG$time_unit) {
  configured <- tolower(as.character(configured)); if (configured %in% c("days", "months")) return(configured)
  tt <- as.numeric(time); tt <- tt[is.finite(tt) & tt > 0]
  if (length(tt) == 0) return("unknown")
  if (stats::quantile(tt, 0.95, na.rm = TRUE) <= 240 && max(tt, na.rm = TRUE) <= 400) "months" else "days"
}
choose_eval_times <- function(time, unit) {
  tt <- as.numeric(time); tt <- tt[is.finite(tt) & tt > 0]
  if (length(tt) == 0) return(numeric(0))
  cand <- if (identical(unit, "months")) CONFIG$eval_times_months else CONFIG$eval_times_days
  cand[cand > min(tt, na.rm = TRUE) & cand < max(tt, na.rm = TRUE)]
}

SIDISH_COLORS <- c(red="#E64B35FF", blue="#4DBBD5FF", green="#00A087FF", navy="#3C5488FF", orange="#F39B7FFF", purple="#8491B4FF", grey="#8C8C8C", black="#000000")
pal <- function(x) unname(SIDISH_COLORS[x])
pub_theme <- function(base_size = CONFIG$plot_base_size, base_family = CONFIG$plot_base_family) {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, colour = "black", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(colour = "#333333", size = base_size),
      axis.title = ggplot2::element_text(colour = "black"), axis.text = ggplot2::element_text(colour = "black"),
      panel.grid.major = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.3), panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", linewidth = CONFIG$plot_linewidth),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = CONFIG$plot_linewidth),
      panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = CONFIG$plot_linewidth),
      legend.title = ggplot2::element_text(face = "bold", colour = "black"), legend.text = ggplot2::element_text(colour = "black"),
      plot.margin = grid::unit(c(5.5, 12, 5.5, 5.5), "pt")
    )
}
save_pub_plot <- function(p, filename_stub, width, height) {
  safe_dir_create(OUT$fig)
  pdf_file <- file.path(OUT$fig, paste0(filename_stub, ".pdf")); png_file <- file.path(OUT$fig, paste0(filename_stub, ".png"))
  pdf_status <- "not_run"; pdf_error <- NA_character_
  tryCatch({
    ggplot2::ggsave(pdf_file, plot = p, width = width, height = height, units = "in", device = grDevices::cairo_pdf, family = CONFIG$plot_base_family, bg = "white")
    pdf_status <- paste0("ok_cairo_pdf_family_", CONFIG$plot_base_family)
  }, error = function(e) {
    pdf_error <<- conditionMessage(e)
    ggplot2::ggsave(pdf_file, plot = p + ggplot2::theme(text = ggplot2::element_text(family = "sans")), width = width, height = height, units = "in", device = grDevices::cairo_pdf, family = "sans", bg = "white")
    pdf_status <<- "ok_cairo_pdf_family_sans_fallback"
  })
  ggplot2::ggsave(png_file, plot = p, width = width, height = height, units = "in", dpi = 600, bg = "white")
  log_file <- file.path(OUT$fig, "figure_export_status.tsv")
  row <- data.table::data.table(figure = filename_stub, pdf_file = pdf_file, png_file = png_file, requested_family = CONFIG$plot_base_family, pdf_status = pdf_status, pdf_error = pdf_error, width = width, height = height)
  if (file.exists(log_file)) row <- data.table::rbindlist(list(data.table::fread(log_file), row), fill = TRUE)
  write_tsv(row, log_file)
}

###############################################################################
# Data preparation
###############################################################################
msg("SIDISH Step 4.5 clinical incremental value audit started.")
if (!file.exists(CONFIG$step4_merged_file)) stop("Step 4 merged dataset not found: ", CONFIG$step4_merged_file)
dt <- read_table_auto(CONFIG$step4_merged_file)
required <- c("cohort_final", "dataset_stage", "sample_id", "OS_time", "OS_event", "locked_risk_score")
miss <- setdiff(required, names(dt)); if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))
dt[, OS_time := as_numeric_safe(OS_time)]; dt[, OS_event := as_numeric_safe(OS_event)]; dt[, locked_risk_score := as_numeric_safe(locked_risk_score)]
dt <- dt[is.finite(OS_time) & OS_time > 0 & OS_event %in% c(0,1) & is.finite(locked_risk_score)]

dev_mean <- mean(dt[dataset_stage == "development", locked_risk_score], na.rm = TRUE)
dev_sd <- stats::sd(dt[dataset_stage == "development", locked_risk_score], na.rm = TRUE)
if (!is.finite(dev_sd) || dev_sd <= 0) dev_sd <- stats::sd(dt$locked_risk_score, na.rm = TRUE)
dt[, SIDISH_risk_z_devscaled := (locked_risk_score - dev_mean) / dev_sd]
dt[, SIDISH_risk_z_within := as.numeric(scale(locked_risk_score)), by = cohort_final]
if (!"Age_z" %in% names(dt)) { dt[, Age := as_numeric_safe(Age)]; dt[, Age_z := as.numeric(scale(Age))] }
if (!"grade_binary_G4" %in% names(dt) && "WHO_grade" %in% names(dt)) dt[, grade_binary_G4 := ifelse(WHO_grade == "G4", "G4", ifelse(!is.na(WHO_grade), "G2_G3", NA_character_))]
if (!"IDH_binary" %in% names(dt) && "IDH_status" %in% names(dt)) dt[, IDH_binary := IDH_status]
if (!"1p19q_status" %in% names(dt)) dt[, `1p19q_status` := NA_character_]
if (!"Sex" %in% names(dt)) dt[, Sex := NA_character_]
for (cc in c("grade_binary_G4", "IDH_binary", "1p19q_status", "Sex")) if (cc %in% names(dt)) dt[[cc]] <- as.factor(dt[[cc]])
write_tsv(dt, file.path(OUT$data, "Step4_5_incremental_value_input_dataset.tsv"))

input_summary <- dt[, .(n=.N, events=sum(OS_event==1), Age_n=sum(!is.na(Age_z)), Grade_n=sum(!is.na(grade_binary_G4)), IDH_n=sum(!is.na(IDH_binary)), p1p19q_n=sum(!is.na(`1p19q_status`)), Sex_n=sum(!is.na(Sex)), risk_n=sum(is.finite(SIDISH_risk_z_devscaled))), by=.(dataset_stage, cohort_final)]
write_tsv(input_summary, file.path(OUT$audit, "Step4_5_input_covariate_availability.tsv"))

###############################################################################
# Modeling helpers
###############################################################################
select_clinical_covariates <- function(dd, covariates = CONFIG$clinical_covariate_priority) {
  selected <- character(0)
  for (cc in covariates) {
    if (!cc %in% names(dd)) next
    v <- dd[[cc]]; ok <- !is.na(v)
    if (sum(ok) < CONFIG$min_complete_cases) next
    if (is.numeric(v)) {
      if (stats::sd(v[ok], na.rm = TRUE) > 0) selected <- c(selected, cc)
    } else {
      tab <- table(v[ok]); if (length(tab) >= 2 && sum(tab >= CONFIG$min_factor_level_n) >= 2) selected <- c(selected, cc)
    }
  }
  selected
}
complete_model_data <- function(dd, predictors) {
  cols <- unique(c("sample_id", "cohort_final", "dataset_stage", "OS_time", "OS_event", predictors)); cols <- cols[cols %in% names(dd)]
  out <- data.table::copy(dd[, ..cols]); out <- out[is.finite(OS_time) & OS_time > 0 & OS_event %in% c(0,1)]
  out <- out[stats::complete.cases(as.data.frame(out[, predictors, with = FALSE]))]
  droplevels(as.data.frame(out))
}
fit_cox_model <- function(dd, predictors, model_name, context, evaluation_type) {
  if (length(predictors) == 0) return(NULL)
  md <- complete_model_data(dd, predictors); n <- nrow(md); events <- sum(md$OS_event == 1, na.rm = TRUE)
  if (n < CONFIG$min_complete_cases || events < CONFIG$min_events) return(NULL)
  keep <- character(0)
  for (pp in predictors) {
    vv <- md[[pp]]
    if (is.numeric(vv)) { if (stats::sd(vv, na.rm = TRUE) > 0) keep <- c(keep, pp) }
    else { vv <- droplevels(as.factor(vv)); md[[pp]] <- vv; if (nlevels(vv) >= 2) keep <- c(keep, pp) }
  }
  if (length(keep) == 0) return(NULL)
  fit <- tryCatch(survival::coxph(make_formula(keep), data = md, x = TRUE, model = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  sm <- summary(fit); lp <- as.numeric(stats::predict(fit, newdata = md, type = "lp"))
  list(fit = fit, model_data = data.table::as.data.table(md), lp = lp,
       metrics = data.table::data.table(context=context, evaluation_type=evaluation_type, model=model_name, predictors=paste(keep, collapse=";"), n=n, events=events, cindex=as.numeric(sm$concordance[1]), AIC=stats::AIC(fit), BIC=stats::BIC(fit), status="ok"))
}
marker_cindex <- function(time, event, marker) {
  dd <- data.frame(time=as.numeric(time), event=as.numeric(event), marker=as.numeric(marker))
  dd <- dd[is.finite(dd$time) & dd$time > 0 & dd$event %in% c(0,1) & is.finite(dd$marker), , drop=FALSE]
  if (nrow(dd) < CONFIG$min_complete_cases || sum(dd$event==1) < CONFIG$min_events || stats::sd(dd$marker) <= 0) return(NA_real_)
  fit <- tryCatch(survival::coxph(survival::Surv(time,event) ~ marker, data=dd), error=function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  as.numeric(summary(fit)$concordance[1])
}
lrt_compare <- function(clin_fit, combined_fit, context, evaluation_type) {
  if (is.null(clin_fit) || is.null(combined_fit)) return(data.table::data.table())
  a <- tryCatch(stats::anova(clin_fit$fit, combined_fit$fit, test="LRT"), error=function(e) NULL)
  p_lrt <- chisq <- NA_real_; if (!is.null(a) && nrow(a) >= 2) { p_lrt <- suppressWarnings(as.numeric(a[2,"P(>|Chi|)"])); chisq <- suppressWarnings(as.numeric(a[2,"Chisq"])) }
  data.table::data.table(context=context, evaluation_type=evaluation_type, clinical_model="Clinical-only", expanded_model="Clinical+SIDISH", n_clinical=clin_fit$metrics$n, n_combined=combined_fit$metrics$n, events_clinical=clin_fit$metrics$events, events_combined=combined_fit$metrics$events, cindex_clinical=clin_fit$metrics$cindex, cindex_combined=combined_fit$metrics$cindex, delta_cindex=combined_fit$metrics$cindex-clin_fit$metrics$cindex, AIC_clinical=clin_fit$metrics$AIC, AIC_combined=combined_fit$metrics$AIC, delta_AIC_combined_minus_clinical=combined_fit$metrics$AIC-clin_fit$metrics$AIC, BIC_clinical=clin_fit$metrics$BIC, BIC_combined=combined_fit$metrics$BIC, delta_BIC_combined_minus_clinical=combined_fit$metrics$BIC-clin_fit$metrics$BIC, LRT_Chisq=chisq, LRT_p=p_lrt, reporting_note=CONFIG$objective_reporting_note)
}
extract_risk_term <- function(fit, context, model, n, events) {
  sm <- summary(fit); cf <- as.data.frame(sm$coefficients); ci <- as.data.frame(sm$conf.int); terms <- rownames(cf)
  idx <- which(terms %in% c("SIDISH_risk_z_devscaled", "SIDISH_risk_z_within")); if (length(idx)==0) return(data.table::data.table())
  data.table::data.table(context=context, model=model, term=terms[idx], n=n, events=events, beta=cf[idx,"coef"], HR=ci[idx,"exp(coef)"], lower95=ci[idx,"lower .95"], upper95=ci[idx,"upper .95"], p=cf[idx,"Pr(>|z|)"])
}
bootstrap_delta_cindex <- function(dd, clinical_predictors, combined_predictors, context, B=CONFIG$bootstrap_n) {
  if (!isTRUE(CONFIG$run_bootstrap_delta_cindex)) return(data.table::data.table())
  rows <- vector("list", B); cohorts <- unique(dd$cohort_final)
  for (b in seq_len(B)) {
    idx <- integer(0); for (coh in cohorts) { ii <- which(dd$cohort_final == coh); idx <- c(idx, sample(ii, length(ii), replace = TRUE)) }
    boot <- dd[idx]; clin <- fit_cox_model(boot, clinical_predictors, "Clinical-only", context, "bootstrap"); comb <- fit_cox_model(boot, combined_predictors, "Clinical+SIDISH", context, "bootstrap")
    rows[[b]] <- data.table::data.table(context=context, bootstrap_id=b, cindex_clinical=if(!is.null(clin)) clin$metrics$cindex else NA_real_, cindex_combined=if(!is.null(comb)) comb$metrics$cindex else NA_real_, delta_cindex=if(!is.null(clin)&&!is.null(comb)) comb$metrics$cindex-clin$metrics$cindex else NA_real_)
  }
  data.table::rbindlist(rows, fill=TRUE)
}
predict_risk_at_times <- function(fit, newdata, times) {
  if (is.null(fit) || length(times)==0) return(matrix(NA_real_, nrow=nrow(newdata), ncol=length(times)))
  sf <- tryCatch(survival::survfit(fit, newdata=as.data.frame(newdata)), error=function(e) NULL)
  if (is.null(sf)) return(matrix(NA_real_, nrow=nrow(newdata), ncol=length(times)))
  ss <- tryCatch(summary(sf, times=times, extend=TRUE), error=function(e) NULL)
  if (is.null(ss)) return(matrix(NA_real_, nrow=nrow(newdata), ncol=length(times)))
  surv <- ss$surv
  if (is.null(dim(surv))) surv <- matrix(surv, nrow=nrow(newdata), ncol=length(times), byrow=FALSE)
  if (nrow(surv) != nrow(newdata) && ncol(surv) == nrow(newdata)) surv <- t(surv)
  if (ncol(surv) != length(times)) surv <- matrix(surv, nrow=nrow(newdata), ncol=length(times), byrow=FALSE)
  1 - surv
}
binary_outcome_at_time <- function(time, event, t0) { out <- rep(NA_integer_, length(time)); out[event==1 & time<=t0] <- 1L; out[time > t0] <- 0L; out }
empirical_auc <- function(marker, y) {
  keep <- is.finite(marker) & y %in% c(0L,1L); marker <- marker[keep]; y <- y[keep]
  n1 <- sum(y==1L); n0 <- sum(y==0L); if (n1<1 || n0<1) return(NA_real_)
  r <- rank(marker, ties.method="average"); (sum(r[y==1L]) - n1*(n1+1)/2)/(n1*n0)
}
compute_time_auc <- function(time, event, marker, times, unit, context, model) {
  if (length(times)==0) return(data.table::data.table())
  if (requireNamespace("timeROC", quietly=TRUE)) {
    suppressPackageStartupMessages(library(timeROC)); suppressPackageStartupMessages(library(survival)); Surv <- survival::Surv
    troc <- tryCatch(timeROC::timeROC(T=as.numeric(time), delta=as.numeric(event), marker=as.numeric(marker), cause=1, times=times, iid=FALSE), error=function(e) e)
    if (!inherits(troc, "error") && length(troc$AUC)>0) return(data.table::data.table(context=context, model=model, time_unit=unit, time=troc$times, AUC=as.numeric(troc$AUC), method="IPCW_timeROC", status="ok", reason="timeROC_completed"))
    reason <- if (inherits(troc,"error")) conditionMessage(troc) else "timeROC_no_auc"
  } else reason <- "timeROC_not_installed"
  rows <- lapply(times, function(t0) data.table::data.table(context=context, model=model, time_unit=unit, time=t0, AUC=empirical_auc(marker, binary_outcome_at_time(time,event,t0)), method="empirical_dynamic_AUC_excluding_censored_not_IPCW", status="fallback", reason=reason))
  data.table::rbindlist(rows)
}
brier_score_at_time <- function(pred, time, event, t0) {
  y <- binary_outcome_at_time(time,event,t0); keep <- y %in% c(0L,1L) & is.finite(pred)
  if (sum(keep) < CONFIG$min_complete_cases) return(list(bs=NA_real_, n=sum(keep), events=NA_integer_))
  list(bs=mean((pred[keep]-y[keep])^2), n=sum(keep), events=sum(y[keep]==1L))
}
net_benefit <- function(pred, y, thresholds) {
  rows <- lapply(thresholds, function(pt) {
    keep <- y %in% c(0L,1L) & is.finite(pred); yy <- y[keep]; pp <- pred[keep] >= pt; n <- length(yy)
    if (n==0) return(data.table::data.table(threshold=pt, net_benefit=NA_real_, n=0, TP=NA_integer_, FP=NA_integer_))
    TP <- sum(pp & yy==1L); FP <- sum(pp & yy==0L); data.table::data.table(threshold=pt, net_benefit=TP/n - FP/n*pt/(1-pt), n=n, TP=TP, FP=FP)
  })
  data.table::rbindlist(rows)
}

###############################################################################
# Within-context comparisons
###############################################################################
contexts <- list()
for (coh in CONFIG$cohort_order) if (coh %in% dt$cohort_final) contexts[[coh]] <- dt[cohort_final == coh]
if (isTRUE(CONFIG$include_development_pooled)) contexts[["DEVELOPMENT_POOLED"]] <- dt[dataset_stage == "development"]
if (isTRUE(CONFIG$include_all_pooled)) contexts[["ALL_POOLED_WITH_EXTERNAL"]] <- dt

model_metrics_all <- list(); lrt_all <- list(); risk_terms_all <- list(); boot_all <- list(); pred_all <- list(); mi <- li <- ri <- bi <- pi <- 1L
for (ctx in names(contexts)) {
  dd <- data.table::copy(contexts[[ctx]]); if (nrow(dd)<CONFIG$min_complete_cases || sum(dd$OS_event==1)<CONFIG$min_events) next
  clinical_covs <- select_clinical_covariates(dd); if (length(clinical_covs)==0) next
  sidish_cov <- if (ctx %in% CONFIG$cohort_order) "SIDISH_risk_z_within" else "SIDISH_risk_z_devscaled"
  fit_clin <- fit_cox_model(dd, clinical_covs, "Clinical-only", ctx, "within_context_fitted")
  fit_sid <- fit_cox_model(dd, sidish_cov, "SIDISH-only", ctx, "within_context_fitted")
  fit_comb <- fit_cox_model(dd, c(clinical_covs, sidish_cov), "Clinical+SIDISH", ctx, "within_context_fitted")
  for (ff in list(fit_clin, fit_sid, fit_comb)) {
    if (is.null(ff)) next
    model_metrics_all[[mi]] <- ff$metrics; mi <- mi+1L
    pred_all[[pi]] <- data.table::data.table(context=ctx, evaluation_type="within_context_fitted", model=ff$metrics$model, sample_id=ff$model_data$sample_id, cohort_final=ff$model_data$cohort_final, dataset_stage=ff$model_data$dataset_stage, OS_time=ff$model_data$OS_time, OS_event=ff$model_data$OS_event, lp=ff$lp); pi <- pi+1L
  }
  lrt <- lrt_compare(fit_clin, fit_comb, ctx, "within_context_fitted"); if (nrow(lrt)>0) { lrt_all[[li]] <- lrt; li <- li+1L }
  if (!is.null(fit_comb)) { rt <- extract_risk_term(fit_comb$fit, ctx, "Clinical+SIDISH", fit_comb$metrics$n, fit_comb$metrics$events); if (nrow(rt)>0) { risk_terms_all[[ri]] <- rt; ri <- ri+1L } }
  if (!is.null(fit_clin) && !is.null(fit_comb)) { bb <- bootstrap_delta_cindex(dd, clinical_covs, c(clinical_covs, sidish_cov), ctx); if (nrow(bb)>0) { boot_all[[bi]] <- bb; bi <- bi+1L } }
}

###############################################################################
# v1.3 schema-safe helper functions
###############################################################################

empty_lrt_schema <- function() {
  data.table::data.table(
    context = character(),
    evaluation_type = character(),
    clinical_model = character(),
    expanded_model = character(),
    n_clinical = integer(),
    n_combined = integer(),
    events_clinical = integer(),
    events_combined = integer(),
    cindex_clinical = numeric(),
    cindex_combined = numeric(),
    delta_cindex = numeric(),
    AIC_clinical = numeric(),
    AIC_combined = numeric(),
    delta_AIC_combined_minus_clinical = numeric(),
    BIC_clinical = numeric(),
    BIC_combined = numeric(),
    delta_BIC_combined_minus_clinical = numeric(),
    LRT_Chisq = numeric(),
    LRT_p = numeric(),
    reporting_note = character()
  )
}

empty_boot_schema <- function() {
  data.table::data.table(
    context = character(),
    bootstrap_n = integer(),
    delta_cindex_median = numeric(),
    delta_cindex_lower95 = numeric(),
    delta_cindex_upper95 = numeric()
  )
}

empty_auc_delta_schema <- function() {
  data.table::data.table(
    context = character(),
    evaluation_type = character(),
    time_unit = character(),
    time = numeric(),
    `Clinical-only` = numeric(),
    `Clinical+SIDISH` = numeric(),
    delta_AUC_combined_minus_clinical = numeric()
  )
}

empty_brier_delta_schema <- function() {
  data.table::data.table(
    context = character(),
    evaluation_type = character(),
    time_unit = character(),
    time = numeric(),
    `Clinical-only` = numeric(),
    `Clinical+SIDISH` = numeric(),
    delta_Brier_combined_minus_clinical = numeric()
  )
}

empty_ibs_delta_schema <- function() {
  data.table::data.table(
    context = character(),
    evaluation_type = character(),
    `Clinical-only` = numeric(),
    `Clinical+SIDISH` = numeric(),
    delta_IBS_like_combined_minus_clinical = numeric()
  )
}

safe_merge_if_possible <- function(x, y, by, all.x = TRUE) {
  x <- data.table::as.data.table(x)
  y <- data.table::as.data.table(y)
  if (nrow(x) == 0 || ncol(x) == 0) return(x)
  if (nrow(y) == 0 || ncol(y) == 0) return(x)
  if (!all(by %in% names(x)) || !all(by %in% names(y))) return(x)
  merge(x, y, by = by, all.x = all.x)
}

safe_rbindlist <- function(x, schema_fun = NULL) {
  out <- data.table::rbindlist(x, fill = TRUE)
  if ((nrow(out) == 0 || ncol(out) == 0) && !is.null(schema_fun)) out <- schema_fun()
  out
}

build_incremental_from_model_metrics <- function(model_metrics) {
  mm <- data.table::as.data.table(model_metrics)
  if (nrow(mm) == 0 || !all(c("context", "evaluation_type", "model") %in% names(mm))) {
    return(empty_lrt_schema())
  }

  clin <- mm[model == "Clinical-only"]
  comb <- mm[model == "Clinical+SIDISH"]
  if (nrow(clin) == 0 || nrow(comb) == 0) return(empty_lrt_schema())

  data.table::setnames(
    clin,
    old = intersect(names(clin), c("n", "events", "cindex", "AIC", "BIC")),
    new = paste0(intersect(names(clin), c("n", "events", "cindex", "AIC", "BIC")), "_clinical")
  )
  data.table::setnames(
    comb,
    old = intersect(names(comb), c("n", "events", "cindex", "AIC", "BIC")),
    new = paste0(intersect(names(comb), c("n", "events", "cindex", "AIC", "BIC")), "_combined")
  )

  out <- merge(
    clin[, .(context, evaluation_type, n_clinical, events_clinical, cindex_clinical, AIC_clinical, BIC_clinical)],
    comb[, .(context, evaluation_type, n_combined, events_combined, cindex_combined, AIC_combined, BIC_combined)],
    by = c("context", "evaluation_type"),
    all = FALSE
  )
  if (nrow(out) == 0) return(empty_lrt_schema())

  out[, clinical_model := "Clinical-only"]
  out[, expanded_model := "Clinical+SIDISH"]
  out[, delta_cindex := cindex_combined - cindex_clinical]
  out[, delta_AIC_combined_minus_clinical := AIC_combined - AIC_clinical]
  out[, delta_BIC_combined_minus_clinical := BIC_combined - BIC_clinical]
  out[, LRT_Chisq := NA_real_]
  out[, LRT_p := NA_real_]
  out[, reporting_note := paste0(CONFIG$objective_reporting_note, "; rebuilt_from_model_metrics_without_LRT_p")]
  data.table::setcolorder(out, names(empty_lrt_schema()))
  out[]
}


model_metrics <- data.table::rbindlist(model_metrics_all, fill=TRUE)
lrt_table <- safe_rbindlist(lrt_all, schema_fun = empty_lrt_schema)
risk_term_table <- data.table::rbindlist(risk_terms_all, fill=TRUE)
bootstrap_delta <- data.table::rbindlist(boot_all, fill=TRUE)
predictions <- data.table::rbindlist(pred_all, fill=TRUE)

if (nrow(lrt_table) == 0 || ncol(lrt_table) == 0) lrt_table <- empty_lrt_schema()

# v1.2 diagnostic: why no Clinical-only vs Clinical+SIDISH comparison was produced.
if (nrow(model_metrics) > 0 && all(c("context", "model") %in% names(model_metrics))) {
  model_presence <- unique(model_metrics[, .(context, model)])
  model_presence[, present := TRUE]
  model_presence_wide <- data.table::dcast(model_presence, context ~ model, value.var = "present", fill = FALSE)
  needed <- c("Clinical-only", "Clinical+SIDISH")
  for (nn in needed) if (!nn %in% names(model_presence_wide)) model_presence_wide[, (nn) := FALSE]
  model_presence_wide[, incremental_pair_available := get("Clinical-only") == TRUE & get("Clinical+SIDISH") == TRUE]
  model_presence_wide[, diagnostic_reason := data.table::fcase(
    incremental_pair_available == TRUE, "pair_available_but_LRT_table_empty_check_anova_or_complete_case_alignment",
    get("Clinical-only") == FALSE & get("Clinical+SIDISH") == FALSE, "both_clinical_only_and_combined_models_missing",
    get("Clinical-only") == FALSE, "clinical_only_model_missing",
    get("Clinical+SIDISH") == FALSE, "clinical_plus_SIDISH_model_missing",
    default = "unknown"
  )]
} else {
  model_presence_wide <- data.table::data.table(
    context = character(),
    `Clinical-only` = logical(),
    `Clinical+SIDISH` = logical(),
    incremental_pair_available = logical(),
    diagnostic_reason = character()
  )
}
write_tsv(model_presence_wide, file.path(OUT$audit, "Step4_5_incremental_pair_availability_diagnostic.tsv"))

# v1.4: if LRT/anova output is empty but paired model metrics exist, rebuild
# the incremental summary from model_metrics. LRT p is left NA in this fallback.
if (nrow(lrt_table) == 0 && nrow(model_metrics) > 0) {
  lrt_fallback <- build_incremental_from_model_metrics(model_metrics)
  if (nrow(lrt_fallback) > 0) {
    lrt_table <- lrt_fallback
  }
}
write_tsv(data.table::data.table(
  metric = c("lrt_table_rows_after_fallback", "lrt_fallback_used"),
  value = c(nrow(lrt_table), ifelse(any(grepl("rebuilt_from_model_metrics", lrt_table$reporting_note)), TRUE, FALSE))
), file.path(OUT$audit, "Step4_5_lrt_fallback_status.tsv"))

write_tsv(model_metrics, file.path(OUT$discrimination, "Step4_5_model_metrics_within_context.tsv"))
write_tsv(lrt_table, file.path(OUT$discrimination, "Step4_5_incremental_Cindex_AIC_BIC_LRT.tsv"))
write_tsv(risk_term_table, file.path(OUT$discrimination, "Step4_5_clinical_plus_SIDISH_risk_term.tsv"))
write_tsv(bootstrap_delta, file.path(OUT$discrimination, "Step4_5_bootstrap_delta_Cindex_long.tsv"))
write_tsv(predictions, file.path(OUT$data, "Step4_5_model_linear_predictors_long.tsv"))

if (nrow(bootstrap_delta)>0 && all(c("context","delta_cindex") %in% names(bootstrap_delta))) {
  boot_summary <- bootstrap_delta[, .(
    bootstrap_n=sum(is.finite(delta_cindex)),
    delta_cindex_median=median(delta_cindex,na.rm=TRUE),
    delta_cindex_lower95=stats::quantile(delta_cindex,0.025,na.rm=TRUE),
    delta_cindex_upper95=stats::quantile(delta_cindex,0.975,na.rm=TRUE)
  ), by=context]
} else {
  boot_summary <- empty_boot_schema()
}
write_tsv(boot_summary, file.path(OUT$discrimination, "Step4_5_bootstrap_delta_Cindex_summary.tsv"))

###############################################################################
# Development-trained external application
###############################################################################
ext_app <- data.table::data.table(); ext_pred <- data.table::data.table()
if (isTRUE(CONFIG$run_development_trained_external_application) && "GSE16011" %in% dt$cohort_final) {
  train <- data.table::copy(dt[dataset_stage == "development"]); test <- data.table::copy(dt[cohort_final == "GSE16011"])
  common_covs <- select_clinical_covariates(data.table::rbindlist(list(train, test), fill=TRUE))
  common_covs <- common_covs[vapply(common_covs, function(cc) sum(!is.na(train[[cc]]))>=CONFIG$min_complete_cases && sum(!is.na(test[[cc]]))>=CONFIG$min_complete_cases, logical(1))]
  if (length(common_covs)>0) {
    fits <- list(
      `Clinical-only`=fit_cox_model(train, common_covs, "Clinical-only", "DEVELOPMENT_TRAIN_TO_GSE16011", "development_trained_external_applied"),
      `Clinical+SIDISH`=fit_cox_model(train, c(common_covs, "SIDISH_risk_z_devscaled"), "Clinical+SIDISH", "DEVELOPMENT_TRAIN_TO_GSE16011", "development_trained_external_applied"),
      `SIDISH-only`=fit_cox_model(train, "SIDISH_risk_z_devscaled", "SIDISH-only", "DEVELOPMENT_TRAIN_TO_GSE16011", "development_trained_external_applied")
    )
    rows <- list(); k <- 1L
    for (nm in names(fits)) {
      ff <- fits[[nm]]; if (is.null(ff)) next
      pr <- tryCatch(as.numeric(stats::predict(ff$fit, newdata=as.data.frame(test), type="lp")), error=function(e) rep(NA_real_, nrow(test)))
      rows[[k]] <- data.table::data.table(context="GSE16011", evaluation_type="development_trained_external_applied", model=nm, n=sum(is.finite(pr)), events=sum(test$OS_event[is.finite(pr)]==1, na.rm=TRUE), cindex=marker_cindex(test$OS_time, test$OS_event, pr), predictors=ff$metrics$predictors, training_context="development_pooled"); k <- k+1L
      ext_pred <- data.table::rbindlist(list(ext_pred, data.table::data.table(context="GSE16011", evaluation_type="development_trained_external_applied", model=nm, sample_id=test$sample_id, cohort_final=test$cohort_final, dataset_stage=test$dataset_stage, OS_time=test$OS_time, OS_event=test$OS_event, lp=pr)), fill=TRUE)
    }
    ext_app <- data.table::rbindlist(rows, fill=TRUE)
    if (all(c("Clinical-only", "Clinical+SIDISH") %in% ext_app$model)) ext_app[, delta_cindex_vs_clinical_only := ifelse(model=="Clinical+SIDISH", cindex - ext_app[model=="Clinical-only", cindex][1], NA_real_)]
  }
}
write_tsv(ext_app, file.path(OUT$discrimination, "Step4_5_development_trained_external_application_GSE16011.tsv"))
write_tsv(ext_pred, file.path(OUT$data, "Step4_5_development_trained_external_predictions_GSE16011.tsv"))
predictions_all <- data.table::rbindlist(list(predictions, ext_pred), fill=TRUE)
write_tsv(predictions_all, file.path(OUT$data, "Step4_5_all_model_linear_predictors_long.tsv"))

###############################################################################
# Time AUC, Brier, calibration, reclassification, DCA
###############################################################################
auc_rows <- list(); brier_rows <- list(); cal_rows <- list(); reclass_rows <- list(); inter_rows <- list(); dca_rows <- list(); impact_rows <- list()
ai <- br <- ca <- rc <- ir <- dc <- ip <- 1L
for (ctx in unique(predictions_all$context)) {
  dd_ctx <- predictions_all[context == ctx]; unit <- detect_survival_time_unit(dd_ctx$OS_time); times <- choose_eval_times(dd_ctx$OS_time, unit); if (length(times)==0) next
  for (evtype in unique(dd_ctx$evaluation_type)) {
    dd_ev <- dd_ctx[evaluation_type == evtype]
    pred_risk_store <- list(); base_store <- list()
    for (mod in unique(dd_ev$model)) {
      pm <- dd_ev[model == mod]; if (nrow(pm)<CONFIG$min_complete_cases || sum(pm$OS_event==1)<CONFIG$min_events) next
      auc <- compute_time_auc(pm$OS_time, pm$OS_event, pm$lp, times, unit, ctx, mod); if (nrow(auc)>0) { auc[, evaluation_type := evtype]; auc_rows[[ai]] <- auc; ai <- ai+1L }
      cal_fit <- tryCatch(survival::coxph(survival::Surv(OS_time, OS_event) ~ lp, data=as.data.frame(pm), x=TRUE), error=function(e) NULL); if (is.null(cal_fit)) next
      pred_mat <- predict_risk_at_times(cal_fit, as.data.frame(pm), times); pred_risk_store[[mod]] <- pred_mat; base_store[[mod]] <- pm
      for (ti in seq_along(times)) {
        t0 <- times[ti]; pr <- as.numeric(pred_mat[, ti]); bs <- brier_score_at_time(pr, pm$OS_time, pm$OS_event, t0)
        brier_rows[[br]] <- data.table::data.table(context=ctx, evaluation_type=evtype, model=mod, time_unit=unit, time=t0, n_evaluable=bs$n, events_by_time=bs$events, Brier=bs$bs, method="binary_timepoint_Brier_excluding_censored_before_or_at_time"); br <- br+1L
        y <- binary_outcome_at_time(pm$OS_time, pm$OS_event, t0); cdt <- data.table::data.table(pred_risk=pr, y=y); cdt <- cdt[y %in% c(0L,1L) & is.finite(pred_risk)]
        if (nrow(cdt)>=40 && data.table::uniqueN(cdt$y)>=2) { cdt[, rank := data.table::frank(pred_risk)]; breaks <- unique(stats::quantile(cdt$rank, probs=seq(0,1,length.out=5), na.rm=TRUE)); if (length(breaks)>=3) { cdt[, bin := cut(rank, breaks=breaks, include.lowest=TRUE)]; cal <- cdt[, .(n=.N, mean_predicted_risk=mean(pred_risk), observed_event_rate=mean(y==1L)), by=bin]; cal[, `:=`(context=ctx, evaluation_type=evtype, model=mod, time_unit=unit, time=t0)]; cal_rows[[ca]] <- cal; ca <- ca+1L } }
        nb <- net_benefit(pr, y, CONFIG$dca_thresholds); nb[, `:=`(context=ctx, evaluation_type=evtype, model=mod, time_unit=unit, time=t0, method_note="DCA_style_binary_timepoint_excluding_censored_before_or_at_time")]; dca_rows[[dc]] <- nb; dc <- dc+1L
        valid <- y %in% c(0L,1L) & is.finite(pr)
        for (pt in CONFIG$dca_thresholds) { ph <- pr[valid] >= pt; impact_rows[[ip]] <- data.table::data.table(context=ctx, evaluation_type=evtype, model=mod, time_unit=unit, time=t0, threshold=pt, n_evaluable=sum(valid), predicted_high_per_100=100*mean(ph), true_events_among_predicted_high_per_100=100*mean(ph & y[valid]==1L), observed_event_rate_by_time=mean(y[valid]==1L)); ip <- ip+1L }
      }
    }
    if (all(c("Clinical-only", "Clinical+SIDISH") %in% names(pred_risk_store))) {
      bd <- base_store[["Clinical-only"]]
      for (ti in seq_along(times)) {
        t0 <- times[ti]; cr <- pred_risk_store[["Clinical-only"]][,ti]; cm <- pred_risk_store[["Clinical+SIDISH"]][,ti]
        y <- binary_outcome_at_time(bd$OS_time, bd$OS_event, t0); keep <- y %in% c(0L,1L) & is.finite(cr) & is.finite(cm); if (sum(keep)<CONFIG$min_complete_cases) next
        delta <- cm[keep] - cr[keep]; yk <- y[keep]
        nri_event <- mean(delta[yk==1L] > 0) - mean(delta[yk==1L] < 0); nri_nonevent <- mean(delta[yk==0L] < 0) - mean(delta[yk==0L] > 0)
        idi <- (mean(cm[keep][yk==1L]) - mean(cm[keep][yk==0L])) - (mean(cr[keep][yk==1L]) - mean(cr[keep][yk==0L]))
        reclass_rows[[rc]] <- data.table::data.table(context=ctx, evaluation_type=evtype, time_unit=unit, time=t0, n_evaluable=sum(keep), n_events_by_time=sum(yk==1L), n_nonevents_by_time=sum(yk==0L), continuous_NRI=nri_event+nri_nonevent, NRI_event_component=nri_event, NRI_nonevent_component=nri_nonevent, IDI=idi, method_note="continuous_NRI_IDI_binary_timepoint_excluding_censored_before_or_at_time"); rc <- rc+1L
        qlo <- as.numeric(stats::quantile(cr[keep], CONFIG$intermediate_lower_quantile, na.rm=TRUE)); qhi <- as.numeric(stats::quantile(cr[keep], CONFIG$intermediate_upper_quantile, na.rm=TRUE)); int_idx <- keep & cr >= qlo & cr <= qhi
        int_dt <- data.table::copy(bd[int_idx]); if (nrow(int_dt)>=CONFIG$min_complete_cases && sum(int_dt$OS_event==1)>=CONFIG$min_events && "locked_risk_group" %in% names(dt)) {
          rg <- unique(dt[, .(sample_id, cohort_final, locked_risk_group)]); int_dt <- merge(int_dt, rg, by=c("sample_id","cohort_final"), all.x=TRUE)
          if (data.table::uniqueN(int_dt$locked_risk_group[!is.na(int_dt$locked_risk_group)])>=2) { int_dt[, group01 := ifelse(locked_risk_group=="High",1,0)]; fit_int <- tryCatch(survival::coxph(survival::Surv(OS_time,OS_event) ~ group01, data=as.data.frame(int_dt)), error=function(e) NULL); if (!is.null(fit_int)) { sm <- summary(fit_int); inter_rows[[ir]] <- data.table::data.table(context=ctx, evaluation_type=evtype, time_unit=unit, time=t0, clinical_intermediate_definition=paste0("clinical_only_middle_tertile_", round(qlo,4), "_to_", round(qhi,4)), n_intermediate=nrow(int_dt), events_intermediate=sum(int_dt$OS_event==1), high_n=sum(int_dt$locked_risk_group=="High"), low_n=sum(int_dt$locked_risk_group=="Low"), HR_SIDISH_high_vs_low=exp(as.numeric(stats::coef(fit_int)[1])), p=sm$coefficients[1,"Pr(>|z|)"], cindex=as.numeric(sm$concordance[1])); ir <- ir+1L } }
        }
      }
    }
  }
}
auc_table <- data.table::rbindlist(auc_rows, fill=TRUE); brier_table <- data.table::rbindlist(brier_rows, fill=TRUE); calibration_table <- data.table::rbindlist(cal_rows, fill=TRUE); reclassification_table <- data.table::rbindlist(reclass_rows, fill=TRUE); intermediate_table <- data.table::rbindlist(inter_rows, fill=TRUE); dca_table <- data.table::rbindlist(dca_rows, fill=TRUE); impact_table <- data.table::rbindlist(impact_rows, fill=TRUE)
write_tsv(auc_table, file.path(OUT$discrimination, "Step4_5_time_dependent_AUC_by_model.tsv"))
write_tsv(brier_table, file.path(OUT$prediction_error, "Step4_5_Brier_score_by_model_time.tsv"))
write_tsv(calibration_table, file.path(OUT$prediction_error, "Step4_5_calibration_by_model_quantile.tsv"))
write_tsv(reclassification_table, file.path(OUT$reclassification, "Step4_5_continuous_NRI_IDI_reclassification.tsv"))
write_tsv(intermediate_table, file.path(OUT$reclassification, "Step4_5_clinical_intermediate_risk_rescue_analysis.tsv"))
write_tsv(dca_table, file.path(OUT$utility, "Step4_5_DCA_style_net_benefit_by_model.tsv"))
write_tsv(impact_table, file.path(OUT$utility, "Step4_5_clinical_impact_curve_data.tsv"))
auc_w <- if (nrow(auc_table)>0) data.table::dcast(auc_table, context + evaluation_type + time_unit + time ~ model, value.var="AUC") else data.table::data.table(); if (nrow(auc_w)>0 && all(c("Clinical-only","Clinical+SIDISH") %in% names(auc_w))) auc_w[, delta_AUC_combined_minus_clinical := `Clinical+SIDISH` - `Clinical-only`]
write_tsv(auc_w, file.path(OUT$discrimination, "Step4_5_delta_time_dependent_AUC.tsv"))
brier_w <- if (nrow(brier_table)>0) data.table::dcast(brier_table, context + evaluation_type + time_unit + time ~ model, value.var="Brier") else data.table::data.table(); if (nrow(brier_w)>0 && all(c("Clinical-only","Clinical+SIDISH") %in% names(brier_w))) brier_w[, delta_Brier_combined_minus_clinical := `Clinical+SIDISH` - `Clinical-only`]
write_tsv(brier_w, file.path(OUT$prediction_error, "Step4_5_delta_Brier_score.tsv"))
ibs <- if (nrow(brier_table)>0) brier_table[, .(IBS_like_mean_Brier=mean(Brier, na.rm=TRUE)), by=.(context,evaluation_type,model)] else data.table::data.table(); ibs_w <- if (nrow(ibs)>0) data.table::dcast(ibs, context + evaluation_type ~ model, value.var="IBS_like_mean_Brier") else data.table::data.table(); if (nrow(ibs_w)>0 && all(c("Clinical-only","Clinical+SIDISH") %in% names(ibs_w))) ibs_w[, delta_IBS_like_combined_minus_clinical := `Clinical+SIDISH` - `Clinical-only`]
write_tsv(ibs, file.path(OUT$prediction_error, "Step4_5_IBS_like_mean_Brier_by_model.tsv")); write_tsv(ibs_w, file.path(OUT$prediction_error, "Step4_5_delta_IBS_like_mean_Brier.tsv"))

###############################################################################
# Figures
###############################################################################
if (isTRUE(CONFIG$make_publication_figures)) {
  msg("Generating Step 4.5 figures...")
  fig1 <- data.table::copy(lrt_table); if (nrow(fig1)>0 && "context" %in% names(fig1)) { fig1 <- safe_merge_if_possible(fig1, boot_summary, by="context", all.x=TRUE); fig1[, context := factor(context, levels=rev(unique(context)))]; p1 <- ggplot2::ggplot(fig1, ggplot2::aes(x=delta_cindex, y=context)) + ggplot2::geom_vline(xintercept=0, linetype="dashed", linewidth=CONFIG$plot_linewidth) + ggplot2::geom_errorbarh(ggplot2::aes(xmin=delta_cindex_lower95, xmax=delta_cindex_upper95), height=0.18, linewidth=CONFIG$plot_linewidth, na.rm=TRUE, color=pal("grey")) + ggplot2::geom_point(size=CONFIG$plot_point_size+0.4, color=pal("red")) + ggplot2::labs(title="Step 4.5 | Incremental C-index", subtitle="Clinical+SIDISH minus clinical-only", x="ΔC-index", y=NULL) + pub_theme(); save_pub_plot(p1, "Figure_1_incremental_delta_Cindex", 5.2, max(3.3, 0.25*nrow(fig1)+1.4)) }
  fig2 <- data.table::copy(auc_w); if (nrow(fig2)>0 && "delta_AUC_combined_minus_clinical" %in% names(fig2)) { p2 <- ggplot2::ggplot(fig2, ggplot2::aes(x=paste0(time," ",time_unit), y=context, fill=delta_AUC_combined_minus_clinical)) + ggplot2::geom_tile(color="white", linewidth=0.25) + ggplot2::geom_text(ggplot2::aes(label=sprintf("%.3f", delta_AUC_combined_minus_clinical)), size=2.7, family=CONFIG$plot_base_family) + ggplot2::scale_fill_gradient2(low=pal("blue"), mid="white", high=pal("red"), midpoint=0, na.value="#F2F2F2") + ggplot2::labs(title="Step 4.5 | Incremental time-dependent AUC", subtitle="Clinical+SIDISH minus clinical-only", x="Time", y=NULL, fill="ΔAUC") + pub_theme(); save_pub_plot(p2, "Figure_2_delta_time_dependent_AUC_heatmap", 5.8, max(3.4, 0.25*data.table::uniqueN(fig2$context)+1.6)) }
  fig3 <- data.table::copy(brier_w); if (nrow(fig3)>0 && "delta_Brier_combined_minus_clinical" %in% names(fig3)) { p3 <- ggplot2::ggplot(fig3, ggplot2::aes(x=paste0(time," ",time_unit), y=delta_Brier_combined_minus_clinical, group=context, color=context)) + ggplot2::geom_hline(yintercept=0, linetype="dashed", linewidth=CONFIG$plot_linewidth) + ggplot2::geom_line(linewidth=CONFIG$plot_linewidth) + ggplot2::geom_point(size=CONFIG$plot_point_size) + ggplot2::labs(title="Step 4.5 | Prediction-error change", subtitle="ΔBrier = Clinical+SIDISH minus clinical-only; <0 indicates lower error", x="Time", y="ΔBrier", color="Context") + pub_theme(); save_pub_plot(p3, "Figure_3_delta_Brier_prediction_error", 5.8, 3.6) }
  fig4 <- data.table::copy(intermediate_table); if (nrow(fig4)>0) { fig4[, label := paste(context, paste0(time," ",time_unit), sep=" | ")]; fig4[, label := make.unique(as.character(label), sep=" #")]; fig4[, label := factor(label, levels=rev(unique(label)))]; p4 <- ggplot2::ggplot(fig4, ggplot2::aes(x=HR_SIDISH_high_vs_low, y=label, color=evaluation_type)) + ggplot2::geom_vline(xintercept=1, linetype="dashed", linewidth=CONFIG$plot_linewidth) + ggplot2::geom_point(size=CONFIG$plot_point_size+0.4) + ggplot2::scale_x_log10() + ggplot2::labs(title="Step 4.5 | Intermediate-risk rescue", subtitle="SIDISH high vs low within clinical-only intermediate-risk patients", x="Hazard ratio, log scale", y=NULL, color="Evaluation") + pub_theme(); save_pub_plot(p4, "Figure_4_intermediate_risk_rescue_HR", 6.0, max(3.4, 0.22*nrow(fig4)+1.4)) }
  fig5 <- dca_table[model %in% c("Clinical-only","Clinical+SIDISH","SIDISH-only") & context %in% c("TCGA","CGGA325","CGGA693","GSE16011","DEVELOPMENT_POOLED")]; if (nrow(fig5)>0) { p5 <- ggplot2::ggplot(fig5, ggplot2::aes(x=threshold, y=net_benefit, color=model, linetype=model)) + ggplot2::geom_line(linewidth=CONFIG$plot_linewidth) + ggplot2::facet_grid(context ~ paste0(time," ",time_unit), scales="free_y") + ggplot2::scale_color_manual(values=c("Clinical-only"=pal("navy"), "Clinical+SIDISH"=pal("red"), "SIDISH-only"=pal("orange"))) + ggplot2::labs(title="Step 4.5 | DCA-style clinical utility", subtitle="Binary timepoint approximation; censored before/at time excluded", x="Threshold probability", y="Net benefit", color="Model", linetype="Model") + pub_theme(); save_pub_plot(p5, "Figure_5_DCA_style_clinical_utility", 7.6, 6.4) }
  fig6 <- impact_table[context=="GSE16011" & model %in% c("Clinical-only","Clinical+SIDISH")]; if (nrow(fig6)>0) { mt <- sort(unique(fig6$time))[ceiling(length(unique(fig6$time))/2)]; fig6 <- fig6[time==mt]; p6 <- ggplot2::ggplot(fig6, ggplot2::aes(x=threshold, color=model)) + ggplot2::geom_line(ggplot2::aes(y=predicted_high_per_100), linewidth=CONFIG$plot_linewidth) + ggplot2::geom_line(ggplot2::aes(y=true_events_among_predicted_high_per_100), linewidth=CONFIG$plot_linewidth, linetype="dashed") + ggplot2::scale_color_manual(values=c("Clinical-only"=pal("navy"), "Clinical+SIDISH"=pal("red"))) + ggplot2::labs(title="Step 4.5 | Clinical impact curve in GSE16011", subtitle=paste0("Solid: predicted high-risk per 100; dashed: true events among predicted high-risk | time=", mt, " ", unique(fig6$time_unit)[1]), x="Threshold probability", y="Patients per 100", color="Model") + pub_theme(); save_pub_plot(p6, "Figure_6_GSE16011_clinical_impact_curve", 5.8, 3.6) }
  fig7 <- calibration_table[model %in% c("Clinical-only","Clinical+SIDISH") & context %in% c("GSE16011","DEVELOPMENT_POOLED","TCGA","CGGA325","CGGA693")]; if (nrow(fig7)>0) { p7 <- ggplot2::ggplot(fig7, ggplot2::aes(x=mean_predicted_risk, y=observed_event_rate, color=model)) + ggplot2::geom_abline(intercept=0, slope=1, linetype="dashed", linewidth=CONFIG$plot_linewidth) + ggplot2::geom_point(size=CONFIG$plot_point_size) + ggplot2::geom_line(linewidth=CONFIG$plot_linewidth) + ggplot2::facet_grid(context ~ paste0(time," ",time_unit), scales="free") + ggplot2::scale_color_manual(values=c("Clinical-only"=pal("navy"), "Clinical+SIDISH"=pal("red"))) + ggplot2::labs(title="Step 4.5 | Calibration comparison", subtitle="Observed vs predicted event risk by predicted-risk quantile", x="Mean predicted risk", y="Observed event rate", color="Model") + pub_theme(); save_pub_plot(p7, "Figure_7_calibration_comparison", 7.2, 6.2) }
}

###############################################################################
# Final summaries
###############################################################################
final_summary <- data.table::data.table(metric=c("n_input_samples","n_development_samples","n_external_samples","n_incremental_comparison_rows","n_auc_rows","n_brier_rows","n_reclassification_rows","n_intermediate_rescue_rows","n_dca_rows","n_clinical_impact_rows","bootstrap_delta_cindex_n","reporting_policy"), value=c(nrow(dt), nrow(dt[dataset_stage=="development"]), nrow(dt[dataset_stage=="external"]), nrow(lrt_table), nrow(auc_table), nrow(brier_table), nrow(reclassification_table), nrow(intermediate_table), nrow(dca_table), nrow(impact_table), CONFIG$bootstrap_n, CONFIG$objective_reporting_note))
write_tsv(final_summary, file.path(OUT$final, "SIDISH_Step4_5_incremental_value_final_summary.tsv"))
main_table <- data.table::copy(lrt_table)
if (nrow(main_table) > 0 && "context" %in% names(main_table)) {
  main_table <- safe_merge_if_possible(main_table, boot_summary, by="context", all.x=TRUE)

  if (nrow(auc_w)>0 && all(c("context","delta_AUC_combined_minus_clinical") %in% names(auc_w))) {
    auc_delta <- auc_w[, .(
      mean_delta_AUC=mean(delta_AUC_combined_minus_clinical, na.rm=TRUE),
      min_delta_AUC=min(delta_AUC_combined_minus_clinical, na.rm=TRUE),
      max_delta_AUC=max(delta_AUC_combined_minus_clinical, na.rm=TRUE)
    ), by=context]
    main_table <- safe_merge_if_possible(main_table, auc_delta, by="context", all.x=TRUE)
  }

  if (nrow(ibs_w)>0 && all(c("context","evaluation_type","delta_IBS_like_combined_minus_clinical") %in% names(ibs_w))) {
    ibs_delta <- ibs_w[, .(context, evaluation_type, delta_IBS_like=delta_IBS_like_combined_minus_clinical)]
    main_table <- safe_merge_if_possible(main_table, ibs_delta, by=c("context","evaluation_type"), all.x=TRUE)
  }
} else {
  main_table <- data.table::data.table(
    context = character(),
    evaluation_type = character(),
    note = character()
  )
}

write_tsv(main_table, file.path(OUT$final, "SIDISH_Step4_5_main_incremental_value_summary.tsv"))

if (nrow(main_table) == 0) {
  write_tsv(data.table::data.table(
    issue = "main_incremental_value_summary_empty",
    likely_reason = "No valid Clinical-only and Clinical+SIDISH paired Cox comparison was generated in lrt_table.",
    inspect_first = file.path(OUT$audit, "Step4_5_incremental_pair_availability_diagnostic.tsv"),
    inspect_second = file.path(OUT$discrimination, "Step4_5_model_metrics_within_context.tsv"),
    inspect_third = file.path(OUT$discrimination, "Step4_5_incremental_Cindex_AIC_BIC_LRT.tsv")
  ), file.path(OUT$final, "SIDISH_Step4_5_empty_incremental_summary_diagnostic.tsv"))
}
readme <- c("SIDISH Step 4.5 Clinical incremental value and translational utility audit v1.4", "", "This step does not alter the locked SIDISH gene model.", "It quantifies whether SIDISH adds information beyond a clinical-only model.", "", "Core comparison: Clinical-only vs Clinical+SIDISH", "", "Key outputs:", file.path(OUT$discrimination, "Step4_5_incremental_Cindex_AIC_BIC_LRT.tsv"), file.path(OUT$discrimination, "Step4_5_bootstrap_delta_Cindex_summary.tsv"), file.path(OUT$discrimination, "Step4_5_delta_time_dependent_AUC.tsv"), file.path(OUT$prediction_error, "Step4_5_delta_Brier_score.tsv"), file.path(OUT$prediction_error, "Step4_5_delta_IBS_like_mean_Brier.tsv"), file.path(OUT$reclassification, "Step4_5_continuous_NRI_IDI_reclassification.tsv"), file.path(OUT$reclassification, "Step4_5_clinical_intermediate_risk_rescue_analysis.tsv"), file.path(OUT$utility, "Step4_5_DCA_style_net_benefit_by_model.tsv"), file.path(OUT$utility, "Step4_5_clinical_impact_curve_data.tsv"), file.path(OUT$final, "SIDISH_Step4_5_main_incremental_value_summary.tsv"), "", "DCA/Brier/reclassification are binary timepoint audits excluding subjects censored before/at each time point. Report raw delta metrics objectively.
Version v1.4 rebuilds the Clinical-only vs Clinical+SIDISH incremental summary from model_metrics when anova/LRT output is empty; LRT p is left NA in that fallback, while ΔC-index, ΔAIC and ΔBIC are preserved.")
write_text(readme, file.path(OUT$root, "README_SIDISH_Step4_5_Clinical_Incremental_Value_Audit_v1_4.txt"))
write_text(capture.output(sessionInfo()), file.path(OUT$audit, "sessionInfo.txt"))
msg("SIDISH Step 4.5 clinical incremental value audit complete.")
msg("Output directory: ", OUT$root)
