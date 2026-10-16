# R-only figure: retain all sixteen frozen candidates and unfavorable findings.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
.libPaths(c(file.path(R.home(),'library'),.libPaths()))
av<-commandArgs(trailingOnly=FALSE)
here<-dirname(normalizePath(sub('^--file=','',av[grepl('^--file=',av)]),winslash='/',mustWork=TRUE))
base<-dirname(here)
suppressPackageStartupMessages(library(ggplot2));suppressPackageStartupMessages(library(patchwork))
library(grid)
src<-file.path(base,'05_Analyses','glioma_drug_audit','results_v01')
p<-read.csv(file.path(src,'primary_16_candidate_results.csv'),check.names=FALSE)
r<-read.csv(file.path(src,'model_level_observed_responses_and_scores.csv'),check.names=FALSE)
d<-read.csv(file.path(src,'glioma_parent_expression_and_scores.csv'),check.names=FALSE)
stopifnot(nrow(p)==16,nrow(d)==39,all(p$BH_q_primary_16>=.05))
blue<-'#365F7F';ink<-'#263746';grey<-'#8A949D'
theme_set(theme_classic(base_size=8,base_family='Arial')+
 theme(axis.line=element_line(linewidth=.35,colour=ink),axis.ticks=element_line(linewidth=.3),
       axis.text=element_text(size=7,colour=ink),axis.title=element_text(size=8),
       plot.title=element_text(size=9,face='bold',colour=ink),plot.subtitle=element_text(size=7.5),
       plot.margin=margin(5,7,5,5),legend.title=element_blank(),legend.text=element_text(size=7),
       legend.position='bottom'))
p$original_candidate<-factor(p$original_candidate,levels=rev(p$original_candidate))
p$label<-sprintf('n=%d; q=%.3f',p$n_models,p$BH_q_primary_16)
forest<-ggplot(p,aes(x=Spearman_rho,y=original_candidate))+
 geom_vline(xintercept=0,linetype=2,colour=grey,linewidth=.35)+
 geom_errorbarh(aes(xmin=bootstrap_lower95,xmax=bootstrap_upper95),height=.15,linewidth=.5,colour=blue)+
 geom_point(size=1.8,colour=blue)+
 geom_text(aes(x=1.02,label=label),hjust=0,size=2.45,family='Arial',colour=ink)+
 scale_x_continuous(breaks=seq(-1,1,.5),limits=c(-1,1.46))+
 labs(title='A  All 16 submitted curated compounds',
      subtitle='Source-recommended preferred screens; no association met BH q < 0.05',
      x='Score versus measured viability AUC: Spearman rho (95% bootstrap interval)',y=NULL)+
 theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(6,10,5,5))

plots<-list()
for(i in seq_along(c('Selumetinib','Trametinib','Dabrafenib'))){
 drug<-c('Selumetinib','Trametinib','Dabrafenib')[i]
 pr<-p[as.character(p$original_candidate)==drug,]
 x<-r[r$original_candidate==drug&r$screen_id==pr$screen_id&r$broad_id==pr$broad_id,]
 x$BRAF_documented<-as.logical(x$documented_BRAF_model_excluded_in_sensitivity)
 plots[[i]]<-ggplot(x,aes(x=score,y=auc))+
  geom_point(aes(shape=BRAF_documented),size=1.6,colour=blue,alpha=.8)+
  scale_shape_manual(values=c('FALSE'=16,'TRUE'=17),guide='none')+
  labs(title=paste0(LETTERS[i+1],'  ',drug),
       subtitle=sprintf('n=%d; rho=%.2f; q=%.3f',pr$n_models,pr$Spearman_rho,pr$BH_q_primary_16),
       x='Fixed nine-gene score',y=if(i==1)'Measured viability AUC' else NULL)+
  theme(plot.title=element_text(size=8.5),plot.subtitle=element_text(size=7),plot.margin=margin(5,5,5,5))
}
combo<-forest/(plots[[1]]|plots[[2]]|plots[[3]])+
 plot_layout(heights=c(2.35,1))+
 plot_annotation(caption=paste(strwrap(paste(
  '39 identity-QC-passed glioma models with all nine expression measurements define the frozen reference; each preferred drug screen has 24-31 models.',
  'Intervals: 1,000 model bootstraps. Two-sided permutation tests: 10,000 permutations; BH adjustment across all 16 primary candidates.',
  'Negative rho is the originally hypothesized sensitivity direction. Triangles identify two models documented as BRAF V600E in a prior study; current genotype was not reassayed.',
  'AUC and source curve-fit quality limit interpretation. These exploratory data did not confirm the proposed association between score and sensitivity or establish clinical treatment selection.',sep=' '),width=132),collapse='\n'),
  theme=theme(plot.caption=element_text(size=7,hjust=0,lineheight=1.15),
              plot.caption.position='plot',plot.margin=margin(5,5,5,5)))
name<-'FigureS4_revision_v05';w<-183/25.4;h<-235/25.4
svglite::svglite(file.path(here,paste0(name,'.svg')),width=w,height=h,bg='white');print(combo);dev.off()
grDevices::cairo_pdf(file.path(here,paste0(name,'.pdf')),width=w,height=h,family='Arial',bg='white');print(combo);dev.off()
ragg::agg_tiff(file.path(here,paste0(name,'.tiff')),width=w,height=h,units='in',res=600,compression='lzw',background='white');print(combo);dev.off()
ragg::agg_png(file.path(here,paste0(name,'.png')),width=w,height=h,units='in',res=180,background='white');print(combo);dev.off()
capture.output(sessionInfo(),file=file.path(here,'PRISM_Figure_R_sessionInfo_v05.txt'))
cat('Exported all16-candidate measured-response figure in four formats.\n')

# Main Figure6: source-reported predictions and subjective curation are separated
# from the new measured public-screen findings. Unsupported lower-tier tallies
# and gene-program heatmaps without retained numerical inputs are not recreated.
curated<-read.csv(file.path(src,'Source_curated_compounds_from_submitted_ST8.csv'),check.names=FALSE)
stopifnot(nrow(curated)==16,identical(curated$curated_drug_name,rev(levels(p$original_candidate))))
curated$drug<-factor(curated$curated_drug_name,levels=levels(p$original_candidate))
pa<-ggplot(curated,aes(x=mean_standardized_delta,y=drug))+
 geom_vline(xintercept=0,colour=grey,linewidth=.35)+
 geom_segment(aes(x=0,xend=mean_standardized_delta,yend=drug),colour=blue,linewidth=.55)+
 geom_point(colour=blue,size=1.7)+
 scale_x_continuous(breaks=c(-1,-.5,0),limits=c(-1.15,.05))+
 labs(title='A  Submitted predictions',subtitle='Reported summaries\nPan-cancer resource',
      x='Standardized predicted difference',y=NULL)+
 theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
fields<-c('CNS_score','glioma_score','readiness_score','feasibility_score','mechanism_score')
mm<-do.call(rbind,lapply(fields,function(f)data.frame(drug=curated$drug,criterion=f,value=curated[[f]])))
mm$criterion<-factor(mm$criterion,levels=fields,labels=c('CNS','Glioma','Readiness','Feasibility','Mechanism'))
pb<-ggplot(mm,aes(x=criterion,y=drug,fill=value))+
 geom_tile(colour='white',linewidth=.4)+
 geom_text(aes(label=value,colour=value>=2.5),size=2.4,family='Arial')+
 scale_fill_gradient(low='#EDF2F6',high=blue,limits=c(0,3),guide='none')+
 scale_colour_manual(values=c('FALSE'=ink,'TRUE'='white'),guide='none')+
 labs(title='B  Original curation',subtitle='Ordinal judgements\nNo measured CNS exposure',x=NULL,y=NULL)+
 theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line=element_blank(),
       axis.ticks.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1))
pc<-ggplot(p,aes(x=Spearman_rho,y=original_candidate))+
 geom_vline(xintercept=0,linetype=2,colour=grey,linewidth=.35)+
 geom_errorbarh(aes(xmin=bootstrap_lower95,xmax=bootstrap_upper95),height=.14,colour=blue,linewidth=.45)+
 geom_point(colour=blue,size=1.7)+
 scale_x_continuous(breaks=c(-1,-.5,0,.5,1),limits=c(-1,1))+
 labs(title='C  Measured glioma response',subtitle='PRISM preferred screens\nAll 16 primary BH q >= 0.315',
      x='Score versus AUC: rho (95% interval)',y=NULL)+
 theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())
main6<-pa+pb+pc+plot_layout(widths=c(1.15,1,1.1))+
 plot_annotation(caption=paste(strwrap(paste(
  'All 16 named compounds in the submitted final-curation table are retained. Panel A reports the original computational summaries; the historical fast_pcr implementation remains source-unmatched.',
  'Panel B reproduces author-assigned prioritization scores and does not establish blood-brain-barrier exposure or glioma efficacy.',
  'Panel C uses source-recommended PRISM screens and unchanged nine-gene weights; 24-31 models per compound, 1,000 model-bootstrap intervals, and BH across 16 tests.',
  'The measured screen does not confirm the proposed association between score and drug sensitivity. All screens and contexts: Supplementary Table 19.',sep=' '),width=134),collapse='\n'),
  theme=theme(plot.caption=element_text(size=7,hjust=0,lineheight=1.15),
              plot.caption.position='plot',plot.margin=margin(5,5,5,5)))
name<-'Figure6_revision_v05';w<-183/25.4;h<-169/25.4
svglite::svglite(file.path(here,paste0(name,'.svg')),width=w,height=h,bg='white');print(main6);dev.off()
grDevices::cairo_pdf(file.path(here,paste0(name,'.pdf')),width=w,height=h,family='Arial',bg='white');print(main6);dev.off()
ragg::agg_tiff(file.path(here,paste0(name,'.tiff')),width=w,height=h,units='in',res=600,compression='lzw',background='white');print(main6);dev.off()
ragg::agg_png(file.path(here,paste0(name,'.png')),width=w,height=h,units='in',res=180,background='white');print(main6);dev.off()
cat('Exported corrected source-supported main Figure6 in four formats.\n')
