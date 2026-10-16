# Reproduce the functional table from the retained full-transcriptome donor inputs.
# This checks numerical reproducibility, not independent biological validation.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
a<-commandArgs(trailingOnly=FALSE);here<-dirname(normalizePath(sub('^--file=','',a[grepl('^--file=',a)]),winslash='/',mustWork=TRUE))
suppressPackageStartupMessages(library(limma));library(jsonlite)
x<-readRDS(file.path(here,'whole_transcriptome_matched_effects.rds'))
stopifnot(ncol(x$donor_delta)==81,nrow(x$pairs)==159,nrow(x$donor_delta)==25987)
sdv<-apply(x$donor_delta,1,sd)
keep<-rowSums(x$donor_mean_expression>0)>=10 & is.finite(sdv)&sdv>1e-12
y<-x$donor_delta[keep,,drop=FALSE]
stopifnot(nrow(y)==21667)
m<-read.csv(file.path(here,'hallmark_membership_and_mapping.csv'),stringsAsFactors=FALSE)
indices<-lapply(split(m$gene_symbol,m$gs_name),function(g)which(rownames(y)%in%unique(g)))
indices<-indices[lengths(indices)>=15&lengths(indices)<=500]
z<-camera(y,indices,design=matrix(1,ncol(y),1),contrast=1,use.ranks=FALSE,allow.neg.cor=FALSE,
          inter.gene.cor=NA,trend.var=FALSE,sort=FALSE)
old<-read.csv(file.path(here,'hallmark_camera_all_results.csv'),stringsAsFactors=FALSE)
old<-old[match(rownames(z),old$pathway),]
gap<-max(abs(z$PValue-old$PValue),abs(z$FDR-old$BH_q),abs(z$Correlation-old$Correlation))
stopifnot(gap<1e-12,all(z$Direction==old$Direction),nrow(z)==50)
write_json(list(aggregate_reproduction_passed=TRUE,sets=50,background_genes=nrow(y),donors=ncol(y),max_absolute_discrepancy=gap,
               numerical_check_not_independent_validation=TRUE,limma_version=as.character(packageVersion('limma'))),
           file.path(here,'Aggregate_reproduction_check_v04.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat('All50 gene-set results reproduced; max absolute discrepancy',gap,'\n')
