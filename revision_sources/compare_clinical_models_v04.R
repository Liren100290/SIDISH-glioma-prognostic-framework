# Revision-stage within-cohort clinical comparison; see registered protocol.
# This is not external validation of a development-trained combined model.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
av<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',av[grepl('^--file=',av)])
here<-dirname(normalizePath(script,winslash='/',mustWork=TRUE))
out<-file.path(here,'clinical_comparison_v04');dir.create(out,showWarnings=FALSE)
suppressPackageStartupMessages(library(survival));library(jsonlite)
stopifnot(file.exists(file.path(here,'clinical_comparison_protocol_v04.json')))
all<-read.csv(file.path(here,'source_results_v04','GSE16011_revision_source_all_scored.csv'),
              stringsAsFactors=FALSE,na.strings=c('','NA'))
for(v in c('sex','grade','idh1','codeletion','rt_record','chemo_record','resection_complete'))all[[v]]<-factor(all[[v]])
contexts<-list(age_sex_grade=c('age10','sex','grade'),
               age_sex_grade_KPS=c('age10','sex','grade','kps10'),
               molecular=c('age10','sex','grade','idh1','codeletion'),
               treatment_records=c('age10','sex','grade','kps10','rt_record','chemo_record','resection_complete'))
performance<-list();briers<-list();netbenefit<-list();calibration<-list();audits<-list();k<-0
times<-c(12,36,60)
warnings_seen<-character(0)
quiet_fit<-function(f,d){
 withCallingHandlers(coxph(f,data=d,ties='efron',x=TRUE,model=TRUE,singular.ok=FALSE),
   warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')})
}
cc<-function(time,event,risk){
 if(!all(is.finite(risk)))return(NA_real_)
 concordance(Surv(time,event)~risk,reverse=TRUE)$concordance
}
survival_predictions<-function(fit,new,times){
 lp<-as.numeric(predict(fit,newdata=new,type='lp',reference='zero'))
 bh<-basehaz(fit,centered=FALSE)
 h0<-stepfun(bh$time,c(0,bh$hazard))(times)
 lch<-outer(lp,log(h0),'+')
 p<- -expm1(-exp(lch))
 list(prob=p,lch=lch)
}
# IPCW implementation: censored-before-horizon cases retain zero weight and the
# denominator remains all n, rather than dropping them as known negative cases.
ipcw<-function(time,event,t){
 fit<-survfit(Surv(time,1-event)~1)
 G<-function(x,left=FALSE){vapply(x,function(z){
  ii<-which(if(left)fit$time<z else fit$time<=z)
  if(length(ii))fit$surv[max(ii)] else 1
 },numeric(1))}
 gt<-G(t);gleft<-G(time,left=TRUE)
 cases<-event==1 & time<=t;noncases<-time>t
 w1<-ifelse(cases,1/pmax(gleft,1e-12),0)
 w0<-ifelse(noncases,1/pmax(gt,1e-12),0)
 list(cases=cases,noncases=noncases,wcase=w1,wnoncase=w0,G=gt)
}
# A mathematical boundary check: without censoring, the IPCW Brier score must
# equal the ordinary fully observed binary Brier score.
toy_t<-c(1,2,4,5);toy_e<-rep(1,4);toy_p<-c(.8,.7,.2,.1);toy<-ipcw(toy_t,toy_e,3)
stopifnot(abs(mean(toy$wcase*(1-toy_p)^2+toy$wnoncase*toy_p^2)-mean((as.numeric(toy_t<=3)-toy_p)^2))<1e-12)

for(context in names(contexts)){
 k<-k+1;terms<-contexts[[context]]
 d<-all[complete.cases(all[,c('time_months','event','score_sd',terms)]),,drop=FALSE]
 omitted<-terms[vapply(terms,function(v)length(unique(d[[v]]))<2,logical(1))]
 terms<-setdiff(terms,omitted)
 stopifnot(nrow(d)>=30,sum(d$event)>=20,!anyDuplicated(d$patient_id))
 f0<-reformulate(terms,response='Surv(time_months,event)')
 f1<-reformulate(c(terms,'score_sd'),response='Surv(time_months,event)')
 warnings_before<-length(warnings_seen)
 fit0<-quiet_fit(f0,d);fit1<-quiet_fit(f1,d)
 lp0<-predict(fit0,type='lp');lp1<-predict(fit1,type='lp')
 c0<-cc(d$time_months,d$event,lp0);c1<-cc(d$time_months,d$event,lp1)
 set.seed(20261015+k)
 boots<-replicate(1000,{
  ix<-sample.int(nrow(d),replace=TRUE);dd<-d[ix,,drop=FALSE]
  tryCatch({
   b0<-quiet_fit(f0,dd);b1<-quiet_fit(f1,dd)
   cc(dd$time_months,dd$event,predict(b1,type='lp'))-cc(dd$time_months,dd$event,predict(b0,type='lp'))
  },error=function(e)NA_real_)
 })
 ci<-quantile(boots[is.finite(boots)],c(.025,.975),names=FALSE)
 write.csv(data.frame(replicate=seq_along(boots),apparent_delta_C=boots),file.path(out,paste0(context,'_paired_bootstrap.csv')),row.names=FALSE)
 cv<-matrix(NA_real_,10,3,dimnames=list(NULL,c('clinical','combined','delta')))
 ps0<-array(NA_real_,c(nrow(d),3,10));ps1<-ps0
 fold_records<-list()
 for(rep in 1:10){
  set.seed(20261015+rep)
  folds<-integer(nrow(d))
  strata<-interaction(d$event,d$grade,drop=TRUE)
  for(g in levels(strata)){
   ii<-which(strata==g)
   folds[ii]<-sample(rep(1:5,length.out=length(ii)))
  }
  lch0<-lch1<-rep(NA_real_,nrow(d))
  for(fold in 1:5){
   test<-which(folds==fold);train<-which(folds!=fold)
   stopifnot(length(intersect(d$patient_id[test],d$patient_id[train]))==0)
   fit0cv<-quiet_fit(f0,d[train,,drop=FALSE]);fit1cv<-quiet_fit(f1,d[train,,drop=FALSE])
   p0<-survival_predictions(fit0cv,d[test,,drop=FALSE],times)
   p1<-survival_predictions(fit1cv,d[test,,drop=FALSE],times)
   ps0[test,,rep]<-p0$prob;ps1[test,,rep]<-p1$prob
   lch0[test]<-p0$lch[,2];lch1[test]<-p1$lch[,2]
  }
  stopifnot(all(is.finite(lch0)),all(is.finite(lch1)))
  cv[rep,1]<-cc(d$time_months,d$event,lch0);cv[rep,2]<-cc(d$time_months,d$event,lch1);cv[rep,3]<-cv[rep,2]-cv[rep,1]
  fold_records[[rep]]<-data.frame(patient_id=d$patient_id,repeat_id=rep,fold=folds)
 }
 write.csv(do.call(rbind,fold_records),file.path(out,paste0(context,'_crossvalidation_folds.csv')),row.names=FALSE)
 write.csv(data.frame(repeat_id=1:10,cv),file.path(out,paste0(context,'_crossvalidation_concordance.csv')),row.names=FALSE)
 p0<-apply(ps0,c(1,2),mean);p1<-apply(ps1,c(1,2),mean)
 pred<-data.frame(patient_id=d$patient_id,time_months=d$time_months,event=d$event)
 for(j in seq_along(times)){
  t<-times[j];ip<-ipcw(d$time_months,d$event,t)
  pred[[paste0('clinical_p',t)]]<-p0[,j];pred[[paste0('combined_p',t)]]<-p1[,j]
  if(ip$G<.1)next
  bs0<-mean(ip$wcase*(1-p0[,j])^2+ip$wnoncase*p0[,j]^2)
  bs1<-mean(ip$wcase*(1-p1[,j])^2+ip$wnoncase*p1[,j]^2)
  briers[[length(briers)+1]]<-data.frame(context=context,n=nrow(d),month=t,censoring_survival=ip$G,
    deaths_by_t=sum(ip$cases),observed_beyond_t=sum(ip$noncases),early_censored=sum(d$event==0&d$time_months<=t),
    clinical_IPCW_Brier=bs0,combined_IPCW_Brier=bs1,combined_minus_clinical=bs1-bs0)
  for(pt in seq(.05,.80,by=.05)){
   nb<-function(p){mean(ip$wcase*(p>=pt))-mean(ip$wnoncase*(p>=pt))*pt/(1-pt)}
   netbenefit[[length(netbenefit)+1]]<-data.frame(context=context,month=t,threshold=pt,clinical=nb(p0[,j]),combined=nb(p1[,j]),
     treat_all=mean(ip$wcase)-mean(ip$wnoncase)*pt/(1-pt),treat_none=0)
  }
  for(model in c('clinical','combined')){
   prob<-if(model=='clinical')p0[,j] else p1[,j]
   group<-integer(length(prob));group[order(prob,d$patient_id)]<-pmin(5L,ceiling(seq_along(prob)*5/length(prob)))
   for(g in 1:5){
    ii<-which(group==g);sf<-survfit(Surv(time_months,event)~1,data=d[ii,])
    st<-summary(sf,times=t,extend=TRUE)
    calibration[[length(calibration)+1]]<-data.frame(context=context,model=model,month=t,group=g,n=length(ii),
      mean_predicted=mean(prob[ii]),KM_observed=1-st$surv[1],observed_lower95=1-st$upper[1],observed_upper95=1-st$lower[1],at_risk=st$n.risk[1])
   }
  }
 }
 write.csv(pred,file.path(out,paste0(context,'_mean_out_of_fold_predictions.csv')),row.names=FALSE)
 un<-coxph(Surv(time_months,event)~score_sd,data=d)
 uci<-summary(un)$conf.int['score_sd',];aci<-summary(fit1)$conf.int['score_sd',]
 performance[[k]]<-data.frame(context=context,n=nrow(d),events=sum(d$event),
  apparent_C_clinical=c0,apparent_C_combined=c1,apparent_delta_C=c1-c0,
  apparent_delta_boot_lower95=ci[1],apparent_delta_boot_upper95=ci[2],bootstrap_valid=sum(is.finite(boots)),
  crossvalidated_C_clinical_mean=mean(cv[,1]),crossvalidated_C_combined_mean=mean(cv[,2]),
  crossvalidated_delta_mean=mean(cv[,3]),crossvalidated_delta_min=min(cv[,3]),crossvalidated_delta_max=max(cv[,3]),
  matched_unadjusted_HR=uci['exp(coef)'],matched_unadjusted_lower95=uci['lower .95'],matched_unadjusted_upper95=uci['upper .95'],
  adjusted_HR=aci['exp(coef)'],adjusted_lower95=aci['lower .95'],adjusted_upper95=aci['upper .95'],
  omitted_invariant=paste(omitted,collapse=';'),warnings=length(warnings_seen)-warnings_before)
 audits[[k]]<-list(context=context,n=nrow(d),patient_ids=as.character(d$patient_id),clinical_formula=paste(deparse(f0),collapse=' '),combined_formula=paste(deparse(f1),collapse=' '),omitted_invariant=omitted)
 cat('Completed',context,'n=',nrow(d),'apparent delta=',c1-c0,'mean CV delta=',mean(cv[,3]),'\n');flush.console()
}
write.csv(do.call(rbind,performance),file.path(out,'paired_clinical_performance.csv'),row.names=FALSE)
write.csv(do.call(rbind,briers),file.path(out,'out_of_fold_IPCW_Brier.csv'),row.names=FALSE)
write.csv(do.call(rbind,netbenefit),file.path(out,'out_of_fold_IPCW_net_benefit.csv'),row.names=FALSE)
write.csv(do.call(rbind,calibration),file.path(out,'out_of_fold_KM_calibration.csv'),row.names=FALSE)
saveRDS(list(contexts=audits,warnings=warnings_seen),file.path(out,'clinical_fit_audit_checkpoint.rds'))
write_json(list(created_utc=format(Sys.time(),tz='UTC',usetz=TRUE),population='Source-reconstructed GSE16011, not asserted identical to original240',
 contexts=audits,repeats=10,folds=5,paired_bootstrap=1000,score_internal_weights_changed=FALSE,
 model_type='New within-cohort clinical models, apparent assessment plus repeated cross-validation',
 warnings=as.list(table(warnings_seen)),IPCW_fully_observed_boundary_check=TRUE,
 limitations='Cross-validation does not replace development-trained external clinical prediction validation. Treatment timing, complete cases and historical molecular data remain limitations.'),
 file.path(out,'Clinical_comparison_audit_v04.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
