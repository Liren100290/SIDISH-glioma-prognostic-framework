# Revision-stage analysis, authored 15 October 2026; not historical study code.
# All inputs are read-only. See analysis_protocol_v04.json written before effects.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
argv <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=','',argv[grepl('^--file=',argv)])
out <- dirname(normalizePath(script,winslash='/',mustWork=TRUE))
analysis <- dirname(out); revision <- dirname(analysis); workspace <- dirname(revision)
robust <- file.path(analysis,'patient_robustness')
input <- file.path(workspace,'SIDISH','SIDISH_\u5168\u91cf\u4ea4\u63a5\u5305_20261015',
                   '01_\u9879\u76ee\u539f\u4ef6','01_SIDISH','GBmap_Extend_scRNA_TCGAbulk','02_Filtered_Expression_Matrix.rds')
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(jsonlite))
stopifnot(file.exists(file.path(out,'analysis_protocol_v04.json')))
start <- proc.time()[['elapsed']]
groups <- as.data.frame(fread(file.path(robust,'expression_assay_group_metadata.csv')))
groups$rowcode <- seq_len(nrow(groups))
stopifnot(setequal(unique(groups$state),c('h','b')))
pairs <- dcast(as.data.table(groups[groups$n_cells>=10,]),patient_id+celltype+assay~state,value.var='rowcode')
pairs <- as.data.frame(pairs[!is.na(h)&!is.na(b)])
reference <- as.data.frame(fread(file.path(robust,'assay_matched_donor_delta.csv')))
donors <- as.character(reference$patient_id)
stopifnot(nrow(pairs)==159L,length(unique(pairs$patient_id))==81L,setequal(donors,pairs$patient_id))
gids <- sort(unique(c(pairs$h,pairs$b)))
stopifnot(length(gids)==318L,!anyDuplicated(pairs[c('patient_id','celltype','assay')]))
allcells <- fread(file.path(robust,'trained_cell_assay_group_keys.csv'))
cells <- allcells[(group_code+1L)%in%gids,]
gcodes <- match(cells$group_code+1L,gids)
stopifnot(identical(as.integer(tabulate(gcodes,nbins=length(gids))),as.integer(groups$n_cells[gids])))
write.csv(pairs,file.path(out,'eligible_matched_strata.csv'),row.names=FALSE)
cat('Selected',nrow(cells),'cells,',nrow(pairs),'strata and',length(donors),'donors\n');flush.console()
cat('Reading cleaned expression source\n');flush.console()
mat <- readRDS(input)
stopifnot(inherits(mat,'sparseMatrix'),nrow(mat)==25987L,ncol(mat)==924234L,!anyDuplicated(rownames(mat)))
ix <- match(cells$trained_cell_id,colnames(mat))
stopifnot(!anyNA(ix),!anyDuplicated(ix))
genes <- rownames(mat)
lib <- Matrix::colSums(mat)[ix]
stopifnot(all(is.finite(lib)),all(lib>0))
design <- sparseMatrix(i=seq_along(gcodes),j=gcodes,x=1,dims=c(length(gcodes),length(gids)))
gm <- matrix(NA_real_,nrow=nrow(mat),ncol=length(gids),dimnames=list(genes,as.character(gids)))
for (a in seq.int(1L,nrow(mat),by=750L)) {
  ii <- a:min(a+749L,nrow(mat))
  block <- as(mat[ii,ix,drop=FALSE],'dgCMatrix')
  block@x <- log1p(block@x*rep.int(10000/lib,diff(block@p)))
  gm[ii,] <- sweep(as.matrix(block %*% design),2,groups$n_cells[gids],'/')
  cat('Aggregated genes through',max(ii),'\n');flush.console()
  block <- NULL;invisible(gc())
}
mat <- NULL;allcells <- NULL;cells <- NULL;design <- NULL;invisible(gc())
h <- match(pairs$h,gids);b <- match(pairs$b,gids)
stratum_delta <- gm[,h,drop=FALSE]-gm[,b,drop=FALSE]
stratum_average <- (gm[,h,drop=FALSE]+gm[,b,drop=FALSE])/2
delta <- vapply(donors,function(d)rowMeans(stratum_delta[,pairs$patient_id==d,drop=FALSE]),numeric(length(genes)))
avg <- vapply(donors,function(d)rowMeans(stratum_average[,pairs$patient_id==d,drop=FALSE]),numeric(length(genes)))
rownames(delta)<-rownames(avg)<-genes;colnames(delta)<-colnames(avg)<-donors
# Directly reconcile every prior marker contrast and every matched group mean.
refgenes <- setdiff(names(reference),'patient_id')
stopifnot(all(refgenes%in%genes))
marker_discrepancy <- max(abs(t(delta[refgenes,donors,drop=FALSE])-as.matrix(reference[,refgenes,drop=FALSE])))
oldgm <- as.data.frame(fread(file.path(robust,'cleaned_assay_group_expression_means.csv.gz')))
group_discrepancy <- max(abs(t(gm[refgenes,,drop=FALSE])-as.matrix(oldgm[gids,refgenes,drop=FALSE])))
stopifnot(marker_discrepancy<1e-10,group_discrepancy<1e-10)
saveRDS(list(donor_delta=delta,donor_mean_expression=avg,matched_group_mean=gm,pairs=pairs),
        file.path(out,'whole_transcriptome_matched_effects.rds'),compress=TRUE)
positive <- rowSums(avg>0)
sds <- apply(delta,1,sd)
eligible <- positive>=10L & is.finite(sds) & sds>1e-12
gene_results <- data.frame(gene=genes,n_donors=length(donors),donors_expressing=positive,
                           mean_delta=rowMeans(delta),median_delta=apply(delta,1,median),
                           sd_delta=sds,eligible_background=eligible,
                           exclusion_reason=ifelse(eligible,'',ifelse(positive<10L,'expression in fewer than 10 donors','zero or numerically negligible contrast SD')))
write.csv(gene_results,file.path(out,'whole_transcriptome_gene_effects.csv'),row.names=FALSE)
cat('Loading prespecified Hallmark sets\n');flush.console()
sets <- as.data.frame(msigdbr::msigdbr(db_species='HS',species='Homo sapiens',collection='H'))
stopifnot(length(unique(sets$gs_name))==50L,identical(unique(sets$db_version),'2025.1.Hs'))
sets <- unique(sets[c('gs_name','gene_symbol','db_version','gs_description','gs_pmid','gs_url')])
sets$eligible_background <- sets$gene_symbol%in%genes[eligible]
write.csv(sets,file.path(out,'hallmark_membership_and_mapping.csv'),row.names=FALSE)
y <- delta[eligible,,drop=FALSE]
indices <- lapply(split(sets$gene_symbol,sets$gs_name),function(g)which(rownames(y)%in%unique(g)))
set_audit <- data.frame(pathway=names(indices),mapped_genes=lengths(indices),stringsAsFactors=FALSE)
set_audit$tested <- set_audit$mapped_genes>=15L & set_audit$mapped_genes<=500L
write.csv(set_audit,file.path(out,'hallmark_set_eligibility.csv'),row.names=FALSE)
idx <- indices[set_audit$tested]
camera <- limma::camera(y,index=idx,design=matrix(1,ncol(y),1,dimnames=list(colnames(y),'High_minus_background')),
                       contrast=1,use.ranks=FALSE,allow.neg.cor=FALSE,inter.gene.cor=NA,trend.var=FALSE,sort=FALSE)
camera$pathway <- rownames(camera)
camera$BH_q <- p.adjust(camera$PValue,method='BH')
stopifnot(max(abs(camera$BH_q-camera$FDR))<1e-12,nrow(camera)==length(idx))
camera <- camera[,c('pathway',setdiff(names(camera),'pathway'))]
write.csv(camera,file.path(out,'hallmark_camera_all_results.csv'),row.names=FALSE)
summary <- list(created_utc=format(Sys.time(),tz='UTC',usetz=TRUE),source=input,source_modified=FALSE,
                donors=length(donors),matched_strata=nrow(pairs),selected_cells=length(ix),source_genes=length(genes),
                eligible_background_genes=sum(eligible),gene_sets_available=length(indices),gene_sets_tested=nrow(camera),
                significant_BH_005=sum(camera$BH_q<0.05),
                significant_pathways=camera[camera$BH_q<0.05,,drop=FALSE],
                prior_245_donor_delta_max_abs_discrepancy=marker_discrepancy,prior_matched_group_mean_max_abs_discrepancy=group_discrepancy,
                package_versions=as.list(sapply(c('limma','msigdbr','Matrix','data.table'),function(p)as.character(packageVersion(p)))),
                msigdb_version=unique(sets$db_version),elapsed_seconds=proc.time()[['elapsed']]-start,
                interpretation='Exploratory competitive gene-set association after donor/type/assay matching; not independent validation, absolute pathway activation, mechanism or drug-response proof.')
write_json(summary,file.path(out,'Functional_enrichment_results_v04.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=16)
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
cat('Completed',nrow(camera),'sets;',sum(camera$BH_q<0.05),'BH q<0.05; seconds',proc.time()[['elapsed']]-start,'\n')
