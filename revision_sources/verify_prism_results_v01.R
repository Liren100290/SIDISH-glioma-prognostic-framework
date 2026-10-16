# Independent R arithmetic/rank/BH verification of the Python PRISM analysis.
# Numerical implementation checking is separate from biological validation.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
av<-commandArgs(trailingOnly=FALSE)
here<-dirname(normalizePath(sub('^--file=','',av[grepl('^--file=',av)]),winslash='/',mustWork=TRUE))
out<-file.path(here,'results_v01')
library(jsonlite)
d<-read.csv(file.path(out,'glioma_parent_expression_and_scores.csv'),check.names=FALSE)
for(v in c('documented_GBM','documented_adult','documented_BRAF_model_excluded_in_sensitivity'))d[[v]]<-as.logical(d[[v]])
w<-read.csv(file.path(out,'frozen_expression_reference_and_weights.csv'),check.names=FALSE)
x<-as.matrix(d[,w$gene]);storage.mode(x)<-'double'
stopifnot(nrow(d)==39,length(w$gene)==9,anyDuplicated(d$depmap_id)==0)
z<-scale(x);scores<-as.numeric(z%*%w$locked_coefficient)
score_gap<-max(abs(scores-d$score))
reference_gap<-max(abs(colMeans(x)-w$parent_mean_log2_TPM1),abs(apply(x,2,sd)-w$parent_sample_SD))
stopifnot(score_gap<1e-12,reference_gap<1e-12)
r<-read.csv(file.path(out,'model_level_observed_responses_and_scores.csv'),check.names=FALSE)
t<-read.csv(file.path(out,'all_context_screen_tests.csv'),check.names=FALSE)
t$inferentially_eligible<-as.logical(t$inferentially_eligible)
p<-read.csv(file.path(out,'primary_16_candidate_results.csv'),check.names=FALSE)
stopifnot(nrow(t)==96,nrow(p)==16,anyDuplicated(r[,c('broad_id','screen_id','depmap_id')])==0)
gaps<-numeric(nrow(t))
for(i in seq_len(nrow(t))){
 ids<-switch(t$context[i],all_glioma=d$depmap_id,
  documented_GBM=d$depmap_id[d$documented_GBM],
  documented_adult=d$depmap_id[d$documented_adult],
  excluding_two_documented_BRAF_models=d$depmap_id[!d$documented_BRAF_model_excluded_in_sensitivity])
 dd<-r[r$broad_id==t$broad_id[i]&r$screen_id==t$screen_id[i]&r$depmap_id%in%ids,,drop=FALSE]
 stopifnot(nrow(dd)==t$n_models[i])
 if(t$inferentially_eligible[i])gaps[i]<-abs(cor(dd$score,dd$auc,method='spearman')-t$Spearman_rho[i])
}
qgap<-max(abs(p.adjust(t$permutation_p_two_sided,method='BH')-t$BH_q_sensitivity_family),
          abs(p.adjust(p$permutation_p_two_sided,method='BH')-p$BH_q_primary_16))
stopifnot(max(gaps)<1e-12,qgap<1e-12,sum(p$BH_q_primary_16<.05)==0)
source_check<-list(performed=FALSE)
raw_response<-file.path(here,'primary_sources','PRISM_19Q4_secondary_dose_response.csv')
raw_expression<-file.path(here,'primary_sources','CCLE_19Q4_expression.csv')
if(file.exists(raw_response)&&file.exists(raw_expression)){
 suppressPackageStartupMessages(library(data.table))
 cols<-names(fread(raw_expression,nrows=0))
 er<-fread(raw_expression,select=c(1L,match(w$source_column,cols)),check.names=FALSE)
 er<-er[match(d$depmap_id,er[[1]]),]
 xr<-as.matrix(er[,-1]);storage.mode(xr)<-'double'
 raw_scores<-as.numeric(scale(xr)%*%w$locked_coefficient)
 sr<-fread(raw_response,select=c('broad_id','depmap_id','screen_id','auc'))
 sr<-merge(as.data.table(r[,c('broad_id','depmap_id','screen_id','auc')]),sr,
           by=c('broad_id','depmap_id','screen_id'),suffixes=c('_prepared','_source'))
 stopifnot(nrow(sr)==687)
 raw_gap<-max(abs(sr$auc_prepared-sr$auc_source))
 raw_score_gap<-max(abs(raw_scores-d$score))
 ref_rhos<-numeric(nrow(p))
 for(i in seq_len(nrow(p))){
  xx<-sr[broad_id==p$broad_id[i]&screen_id==p$screen_id[i]]
  xx$raw_score<-raw_scores[match(xx$depmap_id,d$depmap_id)]
  ref_rhos[i]<-cor(xx$raw_score,xx$auc_source,method='spearman')
 }
 source_check<-list(performed=TRUE,response_rows=nrow(sr),source_AUC_max_gap=raw_gap,
  source_score_max_gap=raw_score_gap,source_primary_rho_max_gap=max(abs(ref_rhos-p$Spearman_rho)),
  source_primary_rhos=as.list(setNames(ref_rhos,p$original_candidate)))
 stopifnot(raw_gap<1e-12,raw_score_gap<1e-12)
}
# A separate R Monte Carlo check of the prespecified main experimental drug.
td<-r[r$original_candidate=='Trametinib'&r$screen_id=='MTS010',]
stopifnot(nrow(td)==29)
set.seed(20261015)
a<-rank(td$score);a<-a-mean(a);b<-rank(td$auc);b<-b-mean(b)
rho<-sum(a*b)/sqrt(sum(a*a)*sum(b*b))
perm<-replicate(20000,sum(a*sample(b))/sqrt(sum(a*a)*sum(b*b)))
rp<-(1+sum(abs(perm)>=abs(rho)-1e-14))/20001
py<-p$permutation_p_two_sided[p$original_candidate=='Trametinib']
stopifnot(abs(rp-py)<.025)
write_json(list(independent_R_score_checks_passed=TRUE,score_max_gap=score_gap,
 reference_max_gap=reference_gap,all96_rank_correlations_checked=TRUE,
 maximum_rho_gap=max(gaps),BH_max_gap=qgap,primary_significant_candidates=0,
 primary_source_independent_read=source_check,
 trametinib_independent_R_permutations=20000,trametinib_R_permutation_p=rp,
 trametinib_Python_permutation_p=py,scope='Implementation verification only; no independent biological-validation claim.'),
 file.path(out,'Independent_R_verification_v01.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat('Independent R verification passed for39 parent scores,96 rank correlations and both BH families.\n')
