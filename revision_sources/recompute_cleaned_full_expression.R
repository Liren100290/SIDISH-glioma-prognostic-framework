# Reconstruct the patient-expression gate from the author's full post-QC RDS.
# This is new revision audit code, not a claimed historical analysis script.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
args <- commandArgs(trailingOnly=FALSE)
this_file <- sub('^--file=','',args[grepl('^--file=',args)])
analysis_dir <- dirname(normalizePath(this_file,winslash='/',mustWork=TRUE))
revision_dir <- dirname(analysis_dir)
workspace_dir <- dirname(revision_dir)
input_path <- file.path(workspace_dir,'SIDISH','SIDISH_\u5168\u91cf\u4ea4\u63a5\u5305_20261015',
                        '01_\u9879\u76ee\u539f\u4ef6','01_SIDISH','GBmap_Extend_scRNA_TCGAbulk','02_Filtered_Expression_Matrix.rds')
out_dir <- file.path(analysis_dir,'patient_robustness')
assay_mode <- '--assay' %in% commandArgs(trailingOnly=TRUE)
key_file <- if (assay_mode) 'trained_cell_assay_group_keys.csv' else 'trained_cell_group_keys.csv'
group_file <- if (assay_mode) 'expression_assay_group_metadata.csv' else 'expression_group_metadata.csv'
prefix <- if (assay_mode) 'cleaned_assay' else 'cleaned_full'
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(data.table))
library(jsonlite)
cells <- fread(file.path(out_dir,key_file))
groups <- fread(file.path(out_dir,group_file))
genes <- readLines(file.path(out_dir,'submitted_marker_genes.txt'),warn=FALSE)
started <- proc.time()[['elapsed']]
cat('Reading original full expression RDS:',input_path,'\n');flush.console()
mat <- readRDS(input_path)
cat('RDS class:',paste(class(mat),collapse=','),'dimensions:',paste(dim(mat),collapse=' x '),'\n');flush.console()
stopifnot(inherits(mat,'sparseMatrix'))
idx <- match(cells$trained_cell_id,colnames(mat));gi <- match(genes,rownames(mat))
stopifnot(!anyNA(idx),!anyDuplicated(idx),!anyNA(gi),length(idx)==847201,length(gi)==245)
gcodes <- cells$group_code+1L
stopifnot(identical(as.integer(tabulate(gcodes,nbins=nrow(groups))),as.integer(groups$n_cells)))
cat('Computing full-transcriptome library sizes and extracting the fixed 245 genes\n');flush.console()
lib_totals <- Matrix::colSums(mat)[idx]
small <- as(mat[gi,idx,drop=FALSE],'dgCMatrix')
stopifnot(all(lib_totals>0),all(is.finite(small@x)),all(small@x>=0))
mat <- NULL;invisible(gc())
small@x <- log1p(small@x * rep.int(10000/lib_totals,diff(small@p)))
design <- sparseMatrix(i=seq_along(gcodes),j=gcodes,x=1,dims=c(length(gcodes),nrow(groups)))
group_sums <- small %*% design
means <- sweep(as.matrix(t(group_sums)),1,groups$n_cells,'/')
colnames(means) <- genes
res <- cbind(as.data.frame(groups),as.data.frame(means,check.names=FALSE))
write.csv(res,gzfile(file.path(out_dir,paste0(prefix,'_group_expression_means.csv.gz'))),row.names=FALSE)
if (!assay_mode) write.csv(data.frame(trained_cell_id=cells$trained_cell_id,full_library_total=lib_totals),
          gzfile(file.path(out_dir,'cleaned_full_library_totals.csv.gz')),row.names=FALSE)
write_json(list(source=input_path,mode='cleaned_full',cells=length(idx),genes=length(genes),
                normalization='Full post-QC/decontX transcriptome library sum, target 10000, natural log1p; aggregate fixed marker genes by true donor/type/state',
                elapsed_seconds=proc.time()[['elapsed']]-started,
                source_modified=FALSE,coefficient_or_cutoff_refit=FALSE),
           file.path(out_dir,paste0(prefix,'_expression_computation.json')),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out_dir,paste0(prefix,'_R_sessionInfo.txt')))
cat('Full-matrix reconstruction complete. Seconds:',proc.time()[['elapsed']]-started,'\n')
