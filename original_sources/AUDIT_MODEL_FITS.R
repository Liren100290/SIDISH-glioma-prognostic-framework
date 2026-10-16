.libPaths(c(file.path(R.home(),"library"),.libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(survival))
args <- commandArgs(trailingOnly=TRUE)
base <- if(length(args)) args[1] else getwd()
root <- file.path(base,"original_project")
out <- file.path(base,"verification_results")
stopifnot(dir.exists(out))
step2 <- file.path(root,"SIDISH_Step2_ML175_DirectionStable_Model_Development_v2_6_strict_all3")
source <- paste0(step2,".R")
parsed <- parse(source,encoding="UTF-8",keep.source=TRUE)
allowed <- c("eval_risk","safe_univ_gene_stats","cohort_consensus_stats","rank_genes_for_training",
             "parse_learner_config","safe_gene_terms","fit_predict_candidate")
audit_env <- new.env(parent=globalenv())
loaded <- character()
for (expr in parsed) {
  if (is.call(expr) && length(expr)==3 && identical(expr[[1]],as.name("<-"))) {
    if (!is.symbol(expr[[2]])) next
    name <- as.character(expr[[2]])
    if (identical(name,"CONFIG")) {
      stopifnot(is.call(expr[[3]]),identical(expr[[3]][[1]],as.name("list")))
      eval(expr,audit_env)
    } else if (name %in% allowed) {
      stopifnot(is.call(expr[[3]]),identical(expr[[3]][[1]],as.name("function")))
      eval(expr,audit_env);loaded <- c(loaded,name)
    }
  }
}
stopifnot(setequal(loaded,allowed),exists("CONFIG",audit_env,inherits=FALSE))
dev <- fread(file.path(step2,"01_development_data","Step2_development_matrix_zscored.tsv"))
locked <- fread(file.path(step2,"99_locked_model","SIDISH_Step2_locked_model_coefficients.tsv"))
stopifnot(nrow(dev)==1580,nrow(locked)==9,!anyDuplicated(paste(dev$cohort,dev$sample_id)))
form <- as.formula(paste("Surv(OS_time,OS_event) ~",paste(locked$gene,collapse=" + ")))
fit <- coxph(form,data=dev)
check <- data.table(gene=names(coef(fit)),recomputed_coefficient=as.numeric(coef(fit)))
check <- merge(check,locked[,.(gene,locked_coefficient)],by="gene")
check[,difference:=recomputed_coefficient-locked_coefficient]
fwrite(check,file.path(out,"independent_nine_gene_Cox_refit.csv"))
stopifnot(max(abs(check$difference))<1e-8)

classa <- fread(file.path(root,"SIDISH_Step1_Core_marker_refinement_v3_5_flat_inputs","99_final_outputs","SIDISH_Step1_ClassA_signed_core_marker_pool_for_Step2_ML.tsv"))
gate <- fread(file.path(step2,"02_model_development","Step2_v2_6_gene_direction_stability_summary.tsv"))
genes <- gate[direction_stable_for_modeling==TRUE,gene]
stopifnot(length(genes)==24)
audit_env$dev_dt <- dev
audit_env$x <- as.matrix(dev[,..genes])
audit_env$y <- Surv(dev$OS_time,dev$OS_event)
audit_env$classA <- classa
selected <- fread(file.path(step2,"02_model_development","Step2_selected_model_configuration.tsv"))
stopifnot(selected$candidate_id=="ML0133",selected$selector=="positive_beta",selected$feature_size==15)
prior_loco <- fread(file.path(step2,"02_model_development","Step2_multiconfig_LOCO_results.tsv"))[candidate_id=="ML0133"]
results <- list()
for (cohort in audit_env$CONFIG$expected_development_cohorts) {
  tr <- which(dev$cohort!=cohort);te <- which(dev$cohort==cohort)
  ranked <- audit_env$rank_genes_for_training(tr,selected$selector,genes,selected$feature_size)
  fitted <- audit_env$fit_predict_candidate(tr,te,ranked$genes,selected$learner_config)
  stopifnot(!is.null(fitted))
  metric <- audit_env$eval_risk(dev$OS_time[te],dev$OS_event[te],fitted$risk)
  fixed_c <- concordance(Surv(OS_time,OS_event)~risk,data=data.frame(dev[te],risk=fitted$risk),reverse=TRUE)$concordance
  previous <- prior_loco[holdout_cohort==cohort]
  stopifnot(nrow(previous)==1,identical(paste(ranked$genes,collapse=";"),previous$selected_gene_set))
  results[[cohort]] <- data.table(cohort=cohort,n=length(te),deaths=sum(dev$OS_event[te]),
    source_cindex=previous$cindex,recomputed_original_eval_cindex=as.numeric(metric["cindex"]),
    independently_oriented_concordance=fixed_c,
    maximum_absolute_difference=abs(as.numeric(metric["cindex"])-previous$cindex),
    selected_gene_order_identical=TRUE,nonzero_genes=fitted$n_nonzero)
}
loco <- rbindlist(results);stopifnot(max(loco$maximum_absolute_difference)<1e-8)
fwrite(loco,file.path(out,"selected_ML0133_LOCO_replay.csv"))
report <- list(R_version=R.version.string,survival_version=as.character(packageVersion("survival")),
  nine_gene_Cox_refit_maximum_coefficient_difference=max(abs(check$difference)),
  selected_ML0133_LOCO_maximum_Cindex_difference=max(loco$maximum_absolute_difference),
  selected_ML0133_all_fold_gene_rankings_identical=TRUE,
  original_245_configuration_search_rerun=FALSE,original_100_bootstrap_runs_rerun=FALSE,
  raw_expression_preprocessing_rerun=FALSE,full_original_script_executed=FALSE,
  exact_source_functions_evaluated=allowed,
  development_validation_limit="The 24-gene direction gate was computed using all three development cohorts before LOCO. This replay is of the archived selection workflow and is not an unbiased fully nested validation estimate.")
writeLines(jsonlite::toJSON(report,auto_unbox=TRUE,pretty=TRUE),file.path(out,"independent_R_model_audit.json"))
writeLines(capture.output(sessionInfo()),file.path(out,"independent_R_model_audit_sessionInfo.txt"))
cat(jsonlite::toJSON(report,auto_unbox=TRUE),"\n")
