from pathlib import Path
from datetime import datetime,timezone
import json
p=Path(__file__).resolve().parent/'gse16011_reconstruction_protocol_v04.json'
if p.exists():raise SystemExit('Preserve the existing protocol')
d={
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'stage':'After source-table schema inspection, before extracting the complete clinical table, cohort selection or calculating any new GSE16011 score',
 'purpose':'Revision-stage reconstruction from original GEO expression and publisher-hosted clinical/molecular/repeat-sample tables. This is not represented as recovered historical input or code.',
 'primary_source':'https://doi.org/10.1158/0008-5472.22383048.v1',
 'eligibility':['Map GEO sample titles to the database numbers in Supplementary Table 1','Exclude controls, grade I/non-diffuse diagnoses and age below 18 or unresolvable age','Exclude the repeat specimens explicitly listed in Supplementary Table 5, keeping the documented primary specimen only','Retain positive follow-up; Dead=1, Alive=0, Lost to follow-up=0 at documented last follow-up; unknown status/time excluded','Use all nine genes from the official GPL8542 custom Entrez probe annotation, not the physical GPL570 array annotation'],
 'scoring':'Deposited log2 gene expression, original nine coefficients and fixed cutoff; standardize across all eligible patients and freeze that reference for subgroup and complete-case analyses',
 'cohort_identity':'Do not alter inclusion rules to force n=240, events=213 or any submitted performance statistic. Compare identities/numerical summaries if original identifiers become available; otherwise state the reconstruction is a distinct revision-stage audit of the same public cohort.',
 'clinical_contexts':['Original clinical-variable set: age, sex, grade, IDH1 and 1p/19q where measured; complete cases reported','Treatment-record sensitivity: age, sex, grade, radiotherapy, chemotherapy, KPS and surgical extent where sufficiently observed and variable; original score formula fixed'],
 'treatment_rules':['Missing treatment fields remain missing, never no treatment','Any recorded chemotherapy other than explicit no is recorded yes, with separate TMZ/PCV/other descriptions retained','Complete resection versus other recorded surgical procedures; unknown stays missing','No treatment start dates are supplied: analyses using recorded therapy are exploratory conditional associations, vulnerable to treatment timing/selection bias. Do not claim causal effects or time-dependent treatment adjustment.'],
 'evaluation':['Full attrition and covariate completeness','Fixed-score concordance with patient-bootstrap uncertainty','Continuous score HR/SD, with available-covariate and treatment-record adjustment','Grade/IDH subgroup checks keep the full eligible reference unchanged','Clinical model comparisons are fitted within this source cohort and must not be called no-refit external validation of combined models'],
 'reporting':'Retain all eligible results, favorable or unfavorable; no gene, coefficient, cutoff or score-direction optimization.'
}
p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
print(p)
