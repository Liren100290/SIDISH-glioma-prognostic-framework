# Independent implementation check, not an independent biological validation.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
args <- commandArgs(trailingOnly=FALSE)
f <- sub('^--file=','',args[grepl('^--file=',args)])
o <- file.path(dirname(normalizePath(f,winslash='/')),'patient_robustness')
suppressPackageStartupMessages(library(data.table));library(jsonlite)
d <- fread(file.path(o,'cleaned_assay_group_expression_means.csv.gz'))
genes <- readLines(file.path(o,'submitted_marker_genes.txt'))
parent <- fread(file.path(o,'patient_eligibility.csv'))[eligible_ge10==TRUE,analysis_patient_id]
h <- d[state=='h' & n_cells>=10 & patient_id %in% parent]
b <- d[state=='b' & n_cells>=10 & patient_id %in% parent]
keys <- c('patient_id','celltype','assay')
joined <- merge(h,b,by=keys,suffixes=c('_h','_b'))
delta <- as.matrix(joined[,paste0(genes,'_h'),with=FALSE])-as.matrix(joined[,paste0(genes,'_b'),with=FALSE])
colnames(delta)<-genes
ids <- joined$patient_id
result <- rowsum(delta,ids)/as.numeric(table(factor(ids,levels=sort(unique(ids)))))
py <- fread(file.path(o,'assay_matched_donor_delta.csv'))
err <- max(abs(result[py$patient_id,genes]-as.matrix(py[,genes,with=FALSE])))
stopifnot(nrow(joined)==159,nrow(result)==81,err<1e-12)
z <- result;z[abs(z)<=1e-12]<-0
stat <- fread(file.path(o,'assay_matched_statistics.csv'))
med <- apply(z,2,median);frac <- colMeans(z>0)
probs <- apply(z,2,function(x) if(all(x==0)) 1 else wilcox.test(x,exact=FALSE,correct=FALSE)$p.value)
mederr <- max(abs(med[stat$gene]-stat$median_delta))
fracerr <- max(abs(frac[stat$gene]-stat$positive_fraction))
perr <- max(abs(probs[stat$gene]-stat$wilcoxon_p_exploratory))
stopifnot(mederr<1e-12,fracerr<1e-12,perr<1e-10)
qerr <- max(abs(p.adjust(probs,method='BH')[stat$gene]-stat$bh_q_exploratory))
stopifnot(qerr<1e-10)
write_json(list(donor_type_assay_strata=nrow(joined),donors=nrow(result),genes=length(genes),
  maximum_donor_delta_discrepancy=err,maximum_median_discrepancy=mederr,
  maximum_positive_fraction_discrepancy=fracerr,maximum_wilcoxon_p_discrepancy=perr,
  maximum_BH_q_discrepancy=qerr,independent_implementation='R data.table/base statistics vs Python pandas/scipy',
  original_expression_matrices_modified=FALSE),file.path(o,'Independent_R_method_check.json'),pretty=TRUE,auto_unbox=TRUE)
cat('Independent R check passed for 245 genes, 81 donors and 159 assay-matched strata.\n')
