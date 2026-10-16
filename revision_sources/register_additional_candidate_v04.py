from pathlib import Path
from datetime import datetime,timezone
import json
out=Path(__file__).resolve().parent
record={
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'stage':'After inspection of the public GEO series description, before downloading REMBRANDT patient clinical or expression data and before any scoring',
 'candidate':'GSE108474 (REMBRANDT)',
 'reason':'The original prediction target is adult diffuse glioma; this accession offers a broader glioma case mix and primary-source clinical data, biospecimen mapping and a data dictionary. The earlier GBM-focused candidate results do not determine inclusion or exclusion of this cohort.',
 'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108474',
 'eligibility':['Age >=18 years for the adult-target evaluation; missing or younger ages excluded and counted','Recorded diffuse glioma histology and grade II-IV where supported; do not infer unknown grade from an unspecified astrocytoma/oligodendroglioma label','Positive overall survival and explicitly verified event coding from primary data dictionary','One independent observation per patient, linked via the official biospecimen mapping','Prefer a documented initial untreated/primary specimen. Exclude documented recurrence; if baseline timing is unavailable, label that population as sampling-timing uncertain, not verified pretreatment validation','All nine gene values must be available; no gene substitution'],
 'duplicates':'If multiple documented baseline specimens per patient remain without a date-based ordering, average the gene expression measurements within patient before cohort standardization; report the count. Never select a specimen by association strength.',
 'scoring':'Use the original nine coefficients and cutoff -0.0269027400037899; outcome-independent exact single-symbol GPL570 probe mapping, average probes; follow deposited preprocessing. Standardize each gene across the eligible accession patient reference, then freeze that reference for every subgroup.',
 'evaluation':['Report every eligible cohort/subgroup regardless of direction or significance','C-index and 1000 patient-bootstrap percentile intervals, seed 20261015','Continuous score HR/SD and fixed-cutoff high/low contrast','Available-covariate adjustment and treatment-aware association when fields are interpretable; complete-case attrition explicit','Grade/diagnosis-specific checks are exploratory, without restandardization or coefficient fitting'],
 'independence':'Compare patient/sample identifiers and recruitment provenance with the development and original external resources. GSE68848 is a related REMBRANDT accession and must not be counted as an additional independent cohort.',
 'boundaries':['This is a revision-stage extension, not historical prespecification','No favorable-result-based candidate selection, coefficient tuning or cutoff optimization','If eligibility/outcome/independence cannot be established, retain the audit and do not call it a successful independent validation']
}
p=out/'additional_candidate_protocol_v04.json'
if p.exists():raise SystemExit('Preserve existing protocol; write a separate amendment if needed')
p.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
print(p)
