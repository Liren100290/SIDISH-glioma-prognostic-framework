# R-only scientific figures from retained revision outputs; see figure contracts.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
# Prefer this runtime's coherent library. The inherited user library mixes a
# newer S7 ggplot build with older runtimes; do not alter global settings.
.libPaths(c(file.path(R.home(),'library'),.libPaths()))
av<-commandArgs(trailingOnly=FALSE);script<-sub('^--file=','',av[grepl('^--file=',av)])
out<-dirname(normalizePath(script,winslash='/',mustWork=TRUE));base<-dirname(out)
suppressPackageStartupMessages(library(ggplot2));suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(survival));library(grid)
src<-file.path(base,'05_Analyses','external_validation','source_results_v04')
clin<-file.path(base,'05_Analyses','external_validation','clinical_comparison_v04')
fun<-file.path(base,'05_Analyses','functional_enrichment')
blue<-'#365F7F';grey<-'#7F8C97';light<-'#E9F0F5';ink<-'#263746'
theme_nature_contract<-function(){theme_classic(base_size=8,base_family='Arial')+
 theme(axis.line=element_line(linewidth=.35,colour=ink),axis.ticks=element_line(linewidth=.3),
       legend.title=element_blank(),legend.text=element_text(size=7),legend.position='bottom',
       plot.title=element_text(size=9,face='bold',colour=ink),plot.subtitle=element_text(size=7.5),
       axis.title=element_text(size=8),axis.text=element_text(size=7,colour=ink),
       plot.margin=margin(5,7,5,5),strip.text=element_text(size=8,face='bold'),
       strip.background=element_rect(fill=light,colour=NA),panel.grid=element_blank())}
theme_set(theme_nature_contract())
save_pub_r<-function(plot,name,width_mm,height_mm){
 selected<-commandArgs(trailingOnly=TRUE)
 if(length(selected)>0L && !name%in%selected)return(invisible(NULL))
 w<-width_mm/25.4;h<-height_mm/25.4
 svglite::svglite(file.path(out,paste0(name,'.svg')),width=w,height=h,bg='white');print(plot);dev.off()
 grDevices::cairo_pdf(file.path(out,paste0(name,'.pdf')),width=w,height=h,family='Arial',bg='white');print(plot);dev.off()
 ragg::agg_tiff(file.path(out,paste0(name,'.tiff')),width=w,height=h,units='in',res=600,compression='lzw',background='white');print(plot);dev.off()
 ragg::agg_png(file.path(out,paste0(name,'.png')),width=w,height=h,units='in',res=180,background='white');print(plot);dev.off()
 cat('Exported',name,'\n');flush.console()
}

# Figure1: each quantitative anchor is source-verified; no asserted clinical gain.
panels<-data.frame(letter=LETTERS[1:6],title=c('Risk-state discovery','Patient evidence','Locked candidate model','External evaluation','Clinical assessment','Biological hypotheses'),
 x=rep(c(.015,.345,.675),2),y=rep(c(.56,.055),each=3),
 body=c('Published SIDISH implementation\nGBmap + TCGA training inputs\n847,201 analyzed cells\n18,620 saved high-risk cells\n245 fixed candidate markers',
        '238 source donors; 100 eligible\n4 strict + 29 supportive markers\nMatched checks: 81 donors\n159 donor/type/assay strata\nTechnical effects remain possible',
        'TCGA + CGGA325 + CGGA693\n33 Class A -> 24 direction-stable\nNine-gene formula locked\nOriginal internal weights retained\nFixed cutoff: -0.02690274',
        'GSE16011 source reconstruction\nAdditional public-cohort evaluation\nGSE43378 and REMBRANDT\nFull cohorts + within-grade checks\nNo gene-weight or cutoff tuning',
        'Matched clinical comparisons\nMissingness + therapy-record audit\nBootstrap + cross-validation\nCensoring-aware error/calibration\nAdded clinical value is uncertain',
        'Whole-transcriptome gene-set test\n81 donors; 50 Hallmark sets\nPan-cancer drug predictions\nExploratory associations\nCNS-specific validation needed'))
ar<-data.frame(x=c(.325,.655,.825,.825,.165,.325,.655),xend=c(.342,.672,.825,.165,.165,.342,.672),
               y=c(.77,.77,.56,.525,.525,.265,.265),yend=c(.77,.77,.525,.525,.478,.265,.265))
p1<-ggplot()+geom_rect(data=panels,aes(xmin=x,xmax=x+.31,ymin=y,ymax=y+.42),fill='#F7FAFC',colour='#8EA7BA',linewidth=.45)+
 geom_rect(data=panels,aes(xmin=x,xmax=x+.31,ymin=y+.34,ymax=y+.42),fill=light,colour=NA)+
 geom_segment(data=ar,aes(x=x,y=y,xend=xend,yend=yend),colour=grey,linewidth=.45,
  arrow=arrow(length=unit(1.5,'mm'),type='closed'))+
 geom_text(data=panels,aes(x=x+.015,y=y+.38,label=letter),hjust=0,fontface='bold',size=3.3,colour=blue,family='Arial')+
 geom_text(data=panels,aes(x=x+.045,y=y+.38,label=title),hjust=0,fontface='bold',size=2.85,colour=ink,family='Arial')+
 geom_text(data=panels,aes(x=x+.016,y=y+.295,label=body),hjust=0,vjust=1,size=2.65,lineheight=1.55,colour=ink,family='Arial')+
 coord_cartesian(xlim=c(0,1),ylim=c(.025,1),expand=FALSE,clip='off')+theme_void()+theme(plot.margin=margin(5,5,5,5))
save_pub_r(p1,'Figure1_revision_v04',183,132)
write.csv(panels,file.path(out,'Figure1_source_text.csv'),row.names=FALSE)

perf<-read.csv(file.path(src,'locked_score_performance.csv'),stringsAsFactors=FALSE)
assoc<-read.csv(file.path(src,'score_survival_associations.csv'),stringsAsFactors=FALSE)
parents<-c('GSE16011_revision_source','GSE108474_timing_uncertain','GSE43378_case_mix')
parent_labels<-c('GSE16011 source audit','REMBRANDT','GSE43378')
colours<-setNames(c(blue,'#6F94AC','#9CA9B2'),parent_labels)
coeff<-read.csv(file.path(base,'05_Analyses','external_validation','locked_nine_gene_coefficients.csv'))
coeff$gene<-factor(coeff$gene,levels=coeff$gene[order(coeff$coefficient)])
pa<-ggplot(coeff,aes(x=coefficient,y=gene))+geom_vline(xintercept=0,colour='#B7C1C8',linewidth=.35)+
 geom_segment(aes(x=0,xend=coefficient,yend=gene),linewidth=1,colour=blue)+geom_point(size=1.8,colour=blue)+
 labs(title='A  Locked nine-gene weights',x='Original Cox coefficient',y=NULL)
mix<-do.call(rbind,lapply(seq_along(parents),function(i){
 d<-read.csv(file.path(src,paste0(parents[i],'_all_scored.csv')))
 tab<-as.data.frame(table(factor(d$grade,levels=2:4)))
 data.frame(cohort=parent_labels[i],grade=tab$Var1,n=tab$Freq)
}))
mix$cohort<-factor(mix$cohort,levels=rev(parent_labels))
pb<-ggplot(mix,aes(x=n,y=cohort,fill=grade))+geom_col(position='fill',width=.55)+
 scale_fill_manual(values=c('2'='#DCE6ED','3'='#91AFC4','4'=blue),labels=c('II','III','IV'))+
 scale_x_continuous(labels=function(x)paste0(round(x*100),'%'))+
 labs(title='B  Historical grade composition',x='Patients within each parent cohort',y=NULL)
d<-read.csv(file.path(src,'GSE16011_revision_source_all_scored.csv'))
d$risk<-factor(ifelse(d$score>-.0269027400037899,'High','Low'),levels=c('Low','High'))
sf<-survfit(Surv(time_months,event)~risk,data=d)
km<-data.frame(time=sf$time,surv=sf$surv,lower=sf$lower,upper=sf$upper,
 risk=rep(sub('risk=','',names(sf$strata)),sf$strata))
km<-rbind(data.frame(time=0,surv=1,lower=1,upper=1,risk=c('Low','High')),km)
pc<-ggplot(km,aes(time,surv,colour=risk))+geom_step(linewidth=.65)+
 scale_colour_manual(values=c(Low=grey,High=blue))+
 scale_x_continuous(limits=c(0,ceiling(max(d$time_months)/60)*60),breaks=seq(0,ceiling(max(d$time_months)/60)*60,60))+scale_y_continuous(limits=c(0,1))+
 labs(title='C  Fixed-cutoff survival separation',subtitle='254 adults; 230 deaths; source reconstruction',x='Months',y='Overall survival')
forest_c<-function(dd,title){
 dd$cohort<-factor(dd$cohort,levels=rev(parent_labels))
 ggplot(dd,aes(x=Cindex,y=cohort,colour=cohort))+
  geom_vline(xintercept=.5,linetype=2,colour=grey,linewidth=.35)+
  geom_segment(aes(x=lower95,xend=upper95,yend=cohort),linewidth=.65)+geom_point(size=2)+
  scale_colour_manual(values=colours,guide='none')+scale_x_continuous(limits=c(.35,.80),breaks=c(.4,.5,.6,.7,.8))+
  labs(title=title,x='C-index (95% bootstrap interval)',y=NULL)
}
full<-perf[match(parents,perf$context),];full$cohort<-parent_labels
pd<-forest_c(full,'D  Parent-cohort discrimination')
gbm<-perf[match(paste0(parents,'_grade4'),perf$context),];gbm$cohort<-parent_labels
pe<-forest_c(gbm,'E  Grade-IV subset discrimination')
adjusted<-assoc[assoc$context%in%parents & assoc$model=='age_sex_grade',]
adjusted$cohort<-factor(parent_labels[match(adjusted$context,parents)],levels=rev(parent_labels))
pf<-ggplot(adjusted,aes(HR_per_parent_score_SD,cohort,colour=cohort))+
 geom_vline(xintercept=1,linetype=2,colour=grey,linewidth=.35)+
 geom_segment(aes(x=lower95,xend=upper95,yend=cohort),linewidth=.65)+geom_point(size=2)+
 scale_colour_manual(values=colours,guide='none')+scale_x_log10(limits=c(.6,2),breaks=c(.75,1,1.5,2))+
 labs(title='F  Age/sex/grade-adjusted association',subtitle='Within-cohort Cox fits; complete cases',x='HR per parent-score SD (95% CI)',y=NULL)
p4<-(pa+pb)/(pc+pd)/(pe+pf)+plot_layout(guides='collect')+
 plot_annotation(caption=paste(strwrap('Internal weights and cutoff remain fixed. Sampling timing is uncertain for REMBRANDT and part of GSE43378. Adjusted associations do not establish clinical utility.',width=155),collapse='\n'),
  theme=theme(plot.caption=element_text(size=6.8,hjust=0),legend.position='bottom'))
save_pub_r(p4,'Figure4_revision_v04',183,176)

cp<-read.csv(file.path(clin,'paired_clinical_performance.csv'),stringsAsFactors=FALSE)
keys<-c('age_sex_grade','age_sex_grade_KPS','molecular','treatment_records')
labels<-c('Age / sex / grade','+ KPS','+ IDH1 / 1p19q','KPS / therapy / surgery')
cp$label<-factor(labels[match(cp$context,keys)],levels=rev(labels))
qa<-ggplot(cp,aes(n,label))+geom_segment(aes(x=0,xend=n,yend=label),linewidth=1,colour='#9EB7C9')+
 geom_point(size=2.2,colour=blue)+geom_text(aes(label=paste0(n,' / ',events)),hjust=-.15,size=2.5,colour=ink)+
 scale_x_continuous(limits=c(0,310),breaks=c(0,100,200,300))+
 labs(title='A  Complete-case clinical contexts',subtitle='Labels show patients / deaths',x='Patients',y=NULL)
hrs<-rbind(data.frame(cp[c('label','context')],model='Unadjusted',hr=cp$matched_unadjusted_HR,lo=cp$matched_unadjusted_lower95,hi=cp$matched_unadjusted_upper95),
 data.frame(cp[c('label','context')],model='Adjusted',hr=cp$adjusted_HR,lo=cp$adjusted_lower95,hi=cp$adjusted_upper95))
dodge<-position_dodge(width=.45)
qb<-ggplot(hrs,aes(hr,label,colour=model))+geom_vline(xintercept=1,linetype=2,colour=grey,linewidth=.35)+
 geom_errorbar(aes(xmin=lo,xmax=hi),orientation='y',width=.15,position=dodge,linewidth=.5)+
 geom_point(position=dodge,size=1.7)+scale_colour_manual(values=c(Unadjusted=grey,Adjusted=blue))+
 scale_x_log10(limits=c(.8,2.1),breaks=c(1,1.5,2))+
 labs(title='B  Same-patient score association',x='HR per parent-score SD (95% CI)',y=NULL)
qc<-ggplot(cp,aes(apparent_delta_C,label))+geom_vline(xintercept=0,linetype=2,colour=grey,linewidth=.35)+
 geom_segment(aes(x=apparent_delta_boot_lower95,xend=apparent_delta_boot_upper95,yend=label),linewidth=.6,colour=blue)+
 geom_point(size=2,colour=blue)+scale_x_continuous(limits=c(-.012,.05),breaks=c(-.01,0,.02,.04))+
 labs(title='C  Apparent discrimination increment',subtitle='Paired bootstrap; 1,000 resamples',x='C-index difference (95% interval)',y=NULL)
cv<-do.call(rbind,lapply(keys,function(key){
 z<-read.csv(file.path(clin,paste0(key,'_crossvalidation_concordance.csv')));z$context<-key;z
}))
cv$label<-factor(labels[match(cv$context,keys)],levels=rev(labels))
set.seed(20261015)
qd<-ggplot(cv,aes(delta,label))+geom_vline(xintercept=0,linetype=2,colour=grey,linewidth=.35)+
 geom_point(position=position_jitter(width=0,height=.11),size=1.4,colour=blue,alpha=.75)+
 scale_x_continuous(limits=c(-.012,.05),breaks=c(-.01,0,.02,.04))+
 labs(title='D  Repeated cross-validation',subtitle='10 repeats of 5 folds; each point is one repeat',x='Out-of-fold C-index difference',y=NULL)
br<-read.csv(file.path(clin,'out_of_fold_IPCW_Brier.csv'));br$label<-factor(labels[match(br$context,keys)],levels=labels)
qe<-ggplot(br,aes(month,combined_minus_clinical,colour=label,shape=label,linetype=label))+
 geom_hline(yintercept=0,linetype=2,colour=grey,linewidth=.35)+geom_line(linewidth=.5)+geom_point(size=1.7)+
 scale_colour_manual(values=c(blue,'#6C92AA','#8A6D8D','#6F8B78'))+
 scale_x_continuous(breaks=c(12,36,60))+labs(title='E  Out-of-fold prediction error',subtitle='Negative difference favors clinical + score',x='Months',y='IPCW Brier-score difference')+
 theme(legend.text=element_text(size=6.1),legend.key.width=unit(7,'mm'))+guides(colour=guide_legend(nrow=2))
cal<-read.csv(file.path(clin,'out_of_fold_KM_calibration.csv'));cal<-cal[cal$context=='molecular',]
cal$model<-factor(cal$model,levels=c('clinical','combined'),labels=c('Clinical','Clinical + score'))
cal_plot<-function(z,title){ggplot(z,aes(mean_predicted,KM_observed,colour=model))+
 geom_abline(slope=1,intercept=0,colour='#A6B1B9',linetype=2,linewidth=.35)+
 geom_errorbar(aes(ymin=observed_lower95,ymax=observed_upper95),width=.015,linewidth=.4)+
 geom_line(linewidth=.45)+geom_point(size=1.6)+
 scale_colour_manual(values=c('Clinical'=grey,'Clinical + score'=blue))+
 scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+
 labs(title=title,x='Predicted death probability',y='Observed death probability')}
qf<-cal_plot(cal[cal$month==36,],'F  Three-year calibration')+labs(subtitle='Molecular complete cases; out-of-fold predictions')
p5<-(qa+qb)/(qc+qd)/(qe+qf)+plot_annotation(caption=paste(strwrap('All clinical models were newly fitted within source-reconstructed GSE16011. Cross-validation points are resampling repeats, not confidence intervals. Therapy timing remains incomplete.',width=155),collapse='\n'),
 theme=theme(plot.caption=element_text(size=6.8,hjust=0)))
save_pub_r(p5,'Figure5_revision_v04',183,182)

# Supplementary Figure2: actual survival calibration and incremental IPCW NB.
nb<-read.csv(file.path(clin,'out_of_fold_IPCW_net_benefit.csv'));nb<-nb[nb$context=='molecular',]
nb$increment<-nb$combined-nb$clinical
tops<-lapply(c(12,36,60),function(t)cal_plot(cal[cal$month==t,],paste0(t,' months')))
bottoms<-lapply(c(12,36,60),function(t){ggplot(nb[nb$month==t,],aes(threshold,increment))+
 geom_hline(yintercept=0,linetype=2,colour=grey,linewidth=.35)+geom_line(colour=blue,linewidth=.6)+
 geom_point(colour=blue,size=1.1)+scale_x_continuous(breaks=c(.1,.4,.8))+
 labs(title=paste0(t,' months'),x='Hypothetical risk threshold',y='Incremental IPCW net benefit')})
s2<-wrap_plots(c(tops,bottoms),ncol=3,guides='collect')+
 plot_annotation(title='Molecular context: 127 patients, 108 deaths',
  caption='Upper row: out-of-fold calibration with KM/Greenwood intervals. Lower row: clinical + score minus clinical net benefit. No decision threshold or patient benefit is established.',
  theme=theme(plot.title=element_text(size=10,face='bold'),plot.caption=element_text(size=6.7,hjust=0),legend.position='bottom'))
save_pub_r(s2,'FigureS2_revision_v04',183,136)

# Supplementary Figure3 shows all tested sets, rather than only the six positives.
gs<-read.csv(file.path(fun,'hallmark_camera_all_results.csv'),stringsAsFactors=FALSE)
gs$label<-gsub('_',' ',sub('^HALLMARK_','',gs$pathway))
gs$signed_logq<- -log10(gs$BH_q)*ifelse(gs$Direction=='Up',1,-1)
gs$label<-factor(gs$label,levels=gs$label[order(gs$signed_logq)])
gs$significant<-gs$BH_q<.05
s3<-ggplot(gs,aes(signed_logq,label,colour=significant,shape=Direction))+
 geom_vline(xintercept=0,colour='#A6B1B9',linewidth=.35)+
 geom_vline(xintercept=c(-1,1)*(-log10(.05)),linetype=2,colour='#A6B1B9',linewidth=.35)+
 geom_segment(aes(x=0,xend=signed_logq,yend=label),linewidth=.35,show.legend=FALSE)+geom_point(size=1.8)+
 scale_colour_manual(values=c('FALSE'='#A8B1B9','TRUE'=blue),labels=c('BH q >= 0.05','BH q < 0.05'))+
 labs(title='Exploratory Hallmark associations after matching',subtitle='81 donors | 159 donor / cell-type / assay strata | 21,667 background genes',
      x='Signed -log10(BH q): Down < 0; Up > 0',y=NULL,
      caption=paste(strwrap('CAMERA estimates within-set gene correlation; BH adjustment covers all 50 sets. Direction is relative to the background transcriptome and does not demonstrate pathway activation.',width=150),collapse='\n'))+
 theme(axis.text.y=element_text(size=7),plot.caption=element_text(size=6.6,hjust=0),plot.caption.position='plot',legend.position='bottom')
save_pub_r(s3,'FigureS3_revision_v04',183,208)
write.csv(mix,file.path(out,'Figure4_grade_composition_source.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
