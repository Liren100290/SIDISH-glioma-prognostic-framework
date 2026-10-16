# Publish the complete coefficients of the newly fitted within-cohort models.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
av<-commandArgs(trailingOnly=FALSE);here<-dirname(normalizePath(sub('^--file=','',av[grepl('^--file=',av)]),winslash='/',mustWork=TRUE))
suppressPackageStartupMessages(library(survival));library(jsonlite)
out<-file.path(here,'clinical_comparison_v04')
d<-read.csv(file.path(here,'source_results_v04','GSE16011_revision_source_all_scored.csv'),stringsAsFactors=FALSE,na.strings=c('','NA'))
for(v in c('sex','grade','idh1','codeletion','rt_record','chemo_record','resection_complete'))d[[v]]<-factor(d[[v]])
audit<-fromJSON(file.path(out,'Clinical_comparison_audit_v04.json'),simplifyVector=FALSE)
rows<-list();baseline<-list()
for(a in audit$contexts){
 dd<-d[as.character(d$patient_id)%in%unlist(a$patient_ids),]
 stopifnot(nrow(dd)==a$n)
 for(type in c('clinical','combined')){
  formula<-as.formula(a[[paste0(type,'_formula')]])
  fit<-coxph(formula,data=dd,ties='efron',x=TRUE,model=TRUE,singular.ok=FALSE)
  s<-summary(fit)
  for(term in rownames(s$coefficients))rows[[length(rows)+1L]]<-data.frame(context=a$context,model=type,n=nrow(dd),events=sum(dd$event),
    term=term,beta=s$coefficients[term,'coef'],SE=s$coefficients[term,'se(coef)'],HR=s$conf.int[term,'exp(coef)'],
    lower95=s$conf.int[term,'lower .95'],upper95=s$conf.int[term,'upper .95'],p=s$coefficients[term,'Pr(>|z|)'])
  bh<-basehaz(fit,centered=FALSE);baseline[[length(baseline)+1L]]<-data.frame(context=a$context,model=type,bh)
 }
}
write.csv(do.call(rbind,rows),file.path(out,'clinical_model_coefficients.csv'),row.names=FALSE)
write.csv(do.call(rbind,baseline),file.path(out,'clinical_model_baseline_hazards.csv'),row.names=FALSE)
cat('Exported',length(rows),'coefficients from8 newly fitted within-cohort models\n')
