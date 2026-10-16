# Data sources and access

The release retains source-derived inputs required by the supplied checks and revision replay.
It does not claim to redistribute every original raw dataset or historical training input.

| Resource | Access and role |
|---|---|
| GBmap | Source atlas linked in the retained provenance; cell/donor/type/assay summaries support fixed-label sensitivity analysis. |
| TCGA and CGGA | Development cohort inputs and retained standardized matrices are identified in the original-workflow archive. Obtain full original data from their source portals. |
| GSE16011 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16011 ; source-reconstructed adult diffuse-glioma survival/clinical evaluation. |
| REMBRANDT / GSE108474 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108474 ; official clinical records and specimen mapping are documented in the source ledger. |
| GSE43378 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43378 ; parent and documented untreated grade-IV subsets are distinguished. |
| PRISM / DepMap | https://depmap.org/portal/data_page/ ; source-version metadata, expression/reference subsets and observed glioma drug response are retained with file-level provenance. |
| MSigDB Hallmark 2025.1.Hs | https://www.gsea-msigdb.org/gsea/msigdb/ ; all 50 Hallmark sets, source membership and matched donor-level contrasts are documented in the functional-analysis package. |
| Historical GDSC / CLUE | Original workflow sources and small retained outputs are provided; full historical training and query resources are not claimed as replayed. |

Refer to each package's manifests, retrieval ledgers and data dictionaries for exact filenames,
versions, units, exclusions and hashes. Expired S3 signing parameters are omitted from the
public ledgers; stable source addresses and payload checksums remain available. Source
dataset citations and usage terms continue to apply to reused or derived material.
