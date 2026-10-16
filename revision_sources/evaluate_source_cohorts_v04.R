# New revision-stage evaluation of the original locked gene formula.
# Original development code is not claimed to have been recovered.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
av<-commandArgs(trailingOnly=FALSE)
script<-sub('^--file=','',av[grepl('^--file=',av)])
here<-dirname(normalizePath(script,winslash='/',mustWork=TRUE))
input<-file.path(here,'primary_prepared_v04')
out<-file.path(here,'source_results_v04');dir.create(out,showWarnings=FALSE)
suppressPackageStartupMessages(library(survival));library(jsonlite)
stopifnot(file.exists(file.path(here,'source_scoring_plan_v04.json')))
set.seed(20261015)
weights<-read.csv(file.path(here,'locked_nine_gene_coefficients.csv'),stringsAsFactors=FALSE)
genes<-weights$gene;beta<-setNames(weights$coefficient,genes);cutoff<--0.0269027400037899
performance<-list();associations<-list();fit_audit<-list();sg_counts<-list();k<-0L;a<-0L;s<-0L

association<-function(d,terms,context,parent,label){
 variables<-c('time_months','event','score_sd',terms)
 dd<-d[complete.cases(d[,variables,drop=FALSE]),,drop=FALSE]
 constants<-terms[vapply(terms,function(x)length(unique(dd[[x]]))<2L,logical(1))]
 rhs<-setdiff(terms,constants)
 fit_audit[[length(fit_audit)+1L]]<<-list(parent=parent,context=context,model=label,parent_n=nrow(d),complete_n=nrow(dd),events=sum(dd$event),omitted_invariant=constants)
 if(nrow(dd)<30L||sum(dd$event)<20L)return(NULL)
 f<-reformulate(c('score_sd',rhs),response='Surv(time_months,event)')
 fit<-coxph(f,data=dd,ties='efron',x=TRUE,model=TRUE)
 ss<-summary(fit);ci<-ss$conf.int['score_sd',]
 ph<-tryCatch(cox.zph(fit)$table['score_sd','p'],error=function(e)NA_real_)
 data.frame(parent=parent,context=context,model=label,n=nrow(dd),events=sum(dd$event),
            HR_per_parent_score_SD=ci['exp(coef)'],lower95=ci['lower .95'],upper95=ci['upper .95'],
            p_Wald=ss$coefficients['score_sd','Pr(>|z|)'],score_PH_p=ph,
            formula=paste(deparse(f),collapse=' '),omitted_invariant=paste(constants,collapse=';'),
            interpretation='Exploratory within-cohort association with fixed nine internal weights; not externally validated clinical utility')
}

evaluate<-function(d,parent,context){
 k<<-k+1L
 cc<-concordance(Surv(time_months,event)~score,data=d,reverse=TRUE)$concordance
 bs<-replicate(1000,{
  ix<-sample.int(nrow(d),replace=TRUE);dd<-d[ix,]
  tryCatch(concordance(Surv(time_months,event)~score,data=dd,reverse=TRUE)$concordance,error=function(e)NA_real_)
 })
 ci<-quantile(bs[is.finite(bs)],c(.025,.975),names=FALSE)
 d$risk_group<-factor(ifelse(d$score>cutoff,'High','Low'),levels=c('Low','High'))
 hr<-lo<-hi<-lp<-NA_real_
 if(all(table(d$risk_group)>0)){
  ff<-coxph(Surv(time_months,event)~risk_group,data=d)
  conf<-summary(ff)$conf.int[1,];hr<-conf['exp(coef)'];lo<-conf['lower .95'];hi<-conf['upper .95']
  lr<-survdiff(Surv(time_months,event)~risk_group,data=d);lp<-pchisq(lr$chisq,1,lower.tail=FALSE)
 }
 performance[[k]]<<-data.frame(parent=parent,context=context,n=nrow(d),events=sum(d$event),
  Cindex=cc,lower95=ci[1],upper95=ci[2],bootstrap_valid=sum(is.finite(bs)),cutoff=cutoff,
  high=sum(d$risk_group=='High'),low=sum(d$risk_group=='Low'),HR_high_low=hr,HR_lower95=lo,HR_upper95=hi,logrank_p=lp)
 contexts<-list(unadjusted=character(0),age_sex_grade=c('age10','sex','grade'))
 if('kps10'%in%names(d))contexts$age_sex_grade_KPS<-c('age10','sex','grade','kps10')
 if(all(c('idh1','codeletion')%in%names(d)))contexts$molecular<-c('age10','sex','grade','idh1','codeletion')
 if(all(c('rt_record','chemo_record','resection_complete','kps10')%in%names(d)))
  contexts$treatment_records<-c('age10','sex','grade','kps10','rt_record','chemo_record','resection_complete')
 for(label in names(contexts)){
  res<-association(d,contexts[[label]],context,parent,label)
  if(!is.null(res)){a<<-a+1L;associations[[a]]<<-res}
 }
 write.csv(d,file.path(out,paste0(context,'_evaluated.csv')),row.names=FALSE)
 cat('Completed',context,'n=',nrow(d),'C=',round(cc,6),'\n');flush.console()
}

files<-c(GSE16011_revision_source='GSE16011_eligible_revision_audit.csv',
         GSE108474_timing_uncertain='GSE108474_eligible_timing_uncertain.csv',
         GSE43378_case_mix='GSE43378_prepared.csv')
for(parent in names(files)){
 path<-if(grepl('GSE43378',parent))file.path(here,files[[parent]]) else file.path(input,files[[parent]])
 d<-read.csv(path,stringsAsFactors=FALSE,na.strings=c('','NA','?'))
 if('candidate_eligible'%in%names(d))d<-d[d$candidate_eligible==1,]
 stopifnot(!anyDuplicated(d$patient_id),all(d$event%in%0:1),all(d$time_months>0))
 x<-as.matrix(d[,genes]);storage.mode(x)<-'double'
 stopifnot(all(is.finite(x)),all(apply(x,2,sd)>0))
 z<-scale(x);d$score<-as.numeric(z%*%beta)
 d$score_sd<-(d$score-mean(d$score))/sd(d$score)
 d$age10<-as.numeric(d$age)/10;d$sex<-factor(tolower(d$sex));d$grade<-factor(as.integer(d$grade),levels=2:4)
 if('kps'%in%names(d))d$kps10<-as.numeric(d$kps)/10
 for(v in intersect(c('idh1','codeletion','rt_record','chemo_record','resection_complete'),names(d)))d[[v]]<-factor(d[[v]])
 write.csv(data.frame(gene=genes,mean=attr(z,'scaled:center'),sample_SD=attr(z,'scaled:scale')),
           file.path(out,paste0(parent,'_reference.csv')),row.names=FALSE)
 write.csv(d,file.path(out,paste0(parent,'_all_scored.csv')),row.names=FALSE)
 evaluate(d,parent,parent)
 for(g in levels(d$grade)){
  dd<-d[!is.na(d$grade)&d$grade==g,,drop=FALSE]
  s<-s+1L;sg_counts[[s]]<-data.frame(parent=parent,variable='grade',level=g,n=nrow(dd),events=sum(dd$event),evaluated=nrow(dd)>=20&&sum(dd$event)>=10)
  if(nrow(dd)>=20L&&sum(dd$event)>=10L)evaluate(dd,parent,paste0(parent,'_grade',g))
 }
 if('idh1'%in%names(d))for(g in levels(d$idh1)){
  dd<-d[!is.na(d$idh1)&d$idh1==g,,drop=FALSE]
  s<-s+1L;sg_counts[[s]]<-data.frame(parent=parent,variable='IDH1_R132',level=g,n=nrow(dd),events=sum(dd$event),evaluated=nrow(dd)>=20&&sum(dd$event)>=10)
  if(nrow(dd)>=20L&&sum(dd$event)>=10L)evaluate(dd,parent,paste0(parent,'_IDH1_',g))
 }
}
write.csv(do.call(rbind,performance),file.path(out,'locked_score_performance.csv'),row.names=FALSE)
write.csv(do.call(rbind,associations),file.path(out,'score_survival_associations.csv'),row.names=FALSE)
write.csv(do.call(rbind,sg_counts),file.path(out,'subgroup_eligibility.csv'),row.names=FALSE)
write_json(fit_audit,file.path(out,'association_completeness_audit.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
write_json(list(created_utc=format(Sys.time(),tz='UTC',usetz=TRUE),parents=names(files),seed=20261015,bootstrap=1000,
                nine_weights_refit=FALSE,cutoff_optimized=FALSE,subgroup_restandardization=FALSE,
                status='Revision-source evaluations; consult source/timing/independence audits before interpreting as final external validation',
                score_HR_unit='One SD of the full eligible parent score, retained in subgroups',
                GSE7696='Earlier provisional diagnostics remain retained separately; not promoted to verified validation because numeric survival-status meaning remains unresolved'),
           file.path(out,'Evaluation_status_v04.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
