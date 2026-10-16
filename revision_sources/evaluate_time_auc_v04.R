# Descriptive time-AUC extension for all source-audited parent populations.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
a<-commandArgs(trailingOnly=FALSE);here<-dirname(normalizePath(sub('^--file=','',a[grepl('^--file=',a)]),winslash='/',mustWork=TRUE))
suppressPackageStartupMessages(library(survival));suppressPackageStartupMessages(library(timeROC));library(jsonlite)
out<-file.path(here,'source_results_v04')
protocol<-file.path(here,'time_auc_extension_protocol_v04.json')
parents<-c('GSE16011_revision_source','GSE108474_timing_uncertain','GSE43378_case_mix')
if(!file.exists(protocol))write_json(list(created_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
 stage='Before any revised time-AUC calculation; after concordance and other source analyses',
 contexts=parents,times_months=c(12,36,60),method='timeROC, cumulative/dynamic IPCW AUC, marginal censoring weights',
 report='All nine parent/time combinations, descriptive point estimates; no favorable-time selection, model update or cutoff optimization'),protocol,pretty=TRUE,auto_unbox=TRUE)
rows<-list()
for(parent in parents){
 d<-read.csv(file.path(out,paste0(parent,'_all_scored.csv')))
 roc<-timeROC::timeROC(T=d$time_months,delta=d$event,marker=d$score,cause=1,weighting='marginal',times=c(12,36,60),iid=FALSE)
 for(i in seq_along(roc$times))rows[[length(rows)+1L]]<-data.frame(context=parent,n=nrow(d),month=roc$times[i],AUC=roc$AUC[i],
  cases=roc$Stats[i,'Cases'],survivors=roc$Stats[i,'survivor at t'],censored=roc$Stats[i,'Censored at t'],
  interpretation='Descriptive cumulative/dynamic AUC point estimate; see source/timing limitations')
}
write.csv(do.call(rbind,rows),file.path(out,'descriptive_time_dependent_AUC.csv'),row.names=FALSE)
cat('All',length(rows),'parent/time AUC estimates saved. timeROC version',as.character(packageVersion('timeROC')),'\n')
