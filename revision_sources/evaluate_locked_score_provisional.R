# Revision-stage validation code, written 15 October 2026.
# This is NOT the historical code that generated the submitted manuscript.
# Results are PROVISIONAL until the source-provenance audit is complete.
# The nine coefficients and development cutoff are never re-estimated.
invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub('^--file=', '', args[grepl('^--file=', args)])
task_dir <- dirname(normalizePath(script_path, winslash = '/', mustWork = TRUE))
suppressPackageStartupMessages(library(survival))
library(jsonlite)
set.seed(20261015)
cutoff <- -0.0269027400037899
coef_tab <- read.csv(file.path(task_dir, 'locked_nine_gene_coefficients.csv'), stringsAsFactors = FALSE)
genes <- coef_tab$gene
beta <- setNames(coef_tab$coefficient, genes)
stopifnot(length(beta) == 9L, all(is.finite(beta)))
out_dir <- file.path(task_dir, 'provisional_results')
dir.create(out_dir, showWarnings = FALSE)
performance <- list(); associations <- list(); attrs <- list(); counter <- 0L

assoc <- function(d, rhs, label) {
  f <- as.formula(paste('Surv(time_months, event) ~', rhs))
  fit <- coxph(f, data = d, ties = 'efron', x = TRUE, model = TRUE)
  ss <- summary(fit)
  if (!'score_sd' %in% rownames(ss$coefficients)) stop('Score term missing')
  v <- ss$coefficients['score_sd', ]
  ci <- ss$conf.int['score_sd', ]
  z <- tryCatch(cox.zph(fit)$table['score_sd', 'p'], error = function(e) NA_real_)
  data.frame(context = label, formula = paste(deparse(f), collapse = ' '), n = fit$n, events = fit$nevent,
             HR_per_score_SD = ci['exp(coef)'], lower95 = ci['lower .95'],
             upper95 = ci['upper .95'], p_Wald = v['Pr(>|z|)'],
             score_PH_test_p = z, row.names = NULL,
             interpretation = 'exploratory association; not validation of a clinical prediction model')
}

evaluate <- function(d, context, covariates) {
  if (nrow(d) < 15 || sum(d$event) < 10) stop('Insufficient evaluable observations')
  d$score_sd <- as.numeric(scale(d$score))
  cc <- concordance(Surv(time_months, event) ~ score, data = d, reverse = TRUE)
  boot <- replicate(1000, {
    ind <- sample.int(nrow(d), replace = TRUE)
    dd <- d[ind, ]
    tryCatch(concordance(Surv(time_months, event) ~ score, data = dd, reverse = TRUE)$concordance,
             error = function(e) NA_real_)
  })
  q <- quantile(boot[is.finite(boot)], c(.025, .975), names = FALSE)
  d$risk_group <- factor(ifelse(d$score > cutoff, 'High', 'Low'), levels = c('Low', 'High'))
  split_hr <- split_lo <- split_hi <- logrank <- NA_real_
  if (all(table(d$risk_group) > 0)) {
    fg <- coxph(Surv(time_months, event) ~ risk_group, data = d)
    ci <- summary(fg)$conf.int[1, ]
    split_hr <- unname(ci['exp(coef)']); split_lo <- unname(ci['lower .95']); split_hi <- unname(ci['upper .95'])
    lr <- survdiff(Surv(time_months, event) ~ risk_group, data = d)
    logrank <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
  }
  per <- data.frame(context = context, n = nrow(d), events = sum(d$event),
                    cindex_locked_score = cc$concordance, cindex_bootstrap_lower95 = q[1],
                    cindex_bootstrap_upper95 = q[2], bootstrap_valid = sum(is.finite(boot)),
                    fixed_cutoff = cutoff, n_high = sum(d$risk_group == 'High'), n_low = sum(d$risk_group == 'Low'),
                    HR_high_vs_low = split_hr, HR_high_vs_low_lower95 = split_lo,
                    HR_high_vs_low_upper95 = split_hi, logrank_p = logrank,
                    status = 'PROVISIONAL; source-provenance audit not yet closed')
  aa <- assoc(d, 'score_sd', paste0(context, ': unadjusted'))
  if (length(covariates)) {
    aa <- rbind(aa, assoc(d, paste(c('score_sd', covariates), collapse = ' + '), paste0(context, ': covariate-adjusted')))
  }
  list(performance = per, associations = aa, data = d)
}

for (acc in c('GSE43378', 'GSE7696')) {
  raw <- read.csv(file.path(task_dir, paste0(acc, '_prepared.csv')), stringsAsFactors = FALSE, na.strings = c('', 'NA', '?'))
  d <- raw[raw$candidate_eligible == 1, ]
  stopifnot(!anyDuplicated(d$patient_id), all(d$event %in% 0:1), all(d$time_months > 0))
  x <- as.matrix(d[, genes]);storage.mode(x) <- 'double'
  stopifnot(all(is.finite(x)), all(apply(x, 2, sd) > 0))
  z <- scale(x)  # sample SD, separately within each candidate cohort
  d$score <- as.numeric(z %*% beta)
  stopifnot(max(abs(d$score - rowSums(sweep(z, 2, beta, '*')))) < 1e-12)
  d$age10 <- as.numeric(d$age) / 10
  d$sex <- factor(tolower(d$sex))
  d$kps10 <- as.numeric(d$kps) / 10
  d$grade <- factor(d$grade)
  d$treatment <- factor(d$treatment)
  d$mgmt <- factor(d$mgmt, levels = c('U', 'M'))
  refs <- data.frame(gene = genes, expression_mean = attr(z, 'scaled:center'), expression_sample_SD = attr(z, 'scaled:scale'), row.names = NULL)
  write.csv(refs, file.path(out_dir, paste0(acc, '_cohort_reference.csv')), row.names = FALSE)
  write.csv(d, file.path(out_dir, paste0(acc, '_scored_patients.csv')), row.names = FALSE)
  attrs[[acc]] <- data.frame(accession = acc, input_n = nrow(raw), candidate_eligible_n = nrow(d), excluded_n = nrow(raw) - nrow(d),
                             event_coding = if (acc == 'GSE43378') 'DEAD=1; ALIVE=0' else 'GEO status assumed 1=death; verify before final use',
                             standardization = 'all candidate-eligible patients within accession; reference frozen for subset analyses')
  if (acc == 'GSE43378') {
    # Full series is a case-mix analysis. Primary-GBM subset has explicit baseline support in the source paper.
    tasks <- list(list(d = d, label = 'GSE43378_all_50_provisional_case_mix', covars = c('age10','sex','grade','kps10')),
                  list(d = d[d$histology == 'GBM' & d$grade == '4', ], label = 'GSE43378_primary_GBM_32', covars = c('age10','sex','kps10')))
  } else {
    tasks <- list(list(d = d, label = 'GSE7696_primary_GBM_70', covars = c('age10','sex','treatment')),
                  list(d = d[complete.cases(d[, c('age10','sex','treatment','mgmt')]), ], label = 'GSE7696_MGMT_complete_cases', covars = c('age10','sex','treatment','mgmt')))
  }
  for (task in tasks) {
    a <- evaluate(task$d, task$label, task$covars);counter <- counter + 1L
    performance[[counter]] <- a$performance;associations[[counter]] <- a$associations
    write.csv(a$data, file.path(out_dir, paste0(task$label, '_evaluated.csv')), row.names = FALSE)
  }
}
write.csv(do.call(rbind, performance), file.path(out_dir, 'locked_score_performance.csv'), row.names = FALSE)
write.csv(do.call(rbind, associations), file.path(out_dir, 'score_survival_associations.csv'), row.names = FALSE)
write.csv(do.call(rbind, attrs), file.path(out_dir, 'sample_attrition.csv'), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, 'R_sessionInfo.txt'))
write_json(list(status = 'PROVISIONAL; not inserted into manuscript', seed = 20261015,
                coefficient_rule = 'nine submitted weights unchanged; no feature selection or coefficient refit',
                cutoff_rule = 'development cutoff unchanged; high risk is score > cutoff',
                bootstrap_rule = '1000 patient bootstrap resamples; fixed observed cohort expression reference; percentile intervals',
                subgroup_reference_rule = 'all-eligible accession reference retained for every subset; no subgroup rescaling of gene expression',
                uncertainty = c('Source independence and baseline sampling audit not closed',
                                'GSE7696 status coding must be confirmed before final use',
                                'Small subgroups, covariate missingness, and multiple exploratory analyses limit interpretation',
                                'Adjusted Cox associations are not validated clinical utility estimates')),
           file.path(out_dir, 'analysis_status.json'), pretty = TRUE, auto_unbox = TRUE)
print(do.call(rbind, performance))
print(do.call(rbind, associations))
