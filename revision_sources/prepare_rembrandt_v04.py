"""Prepare the registered REMBRANDT candidate without score-survival calculations."""
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
import csv,gzip,json,hashlib,statistics,re
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent;SRC=HERE/'primary_sources_v04';OUT=HERE/'primary_prepared_v04';OUT.mkdir(exist_ok=True)
assert (HERE/'additional_candidate_protocol_v04.json').exists() and (HERE/'additional_candidate_schema_amendment_v04.json').exists()
GENES=['GFPT2','KAZN','NTN1','PTPRO','SLC22A23','STK32A','STK32B','TCF7L2','TMEM132B']
mapping=pd.read_csv(HERE/'nine_gene_probe_mapping.csv',dtype=str)
probe_gene=dict(zip(mapping.probe,mapping.gene))
profiles={};nrows=0
with gzip.open(SRC/'GSE108474_GeneExpression.txt.gz','rt',encoding='utf-8') as f:
 r=csv.reader(f,delimiter='\t');specimens=next(r)
 assert len(specimens)==len(set(specimens))
 for row in r:
  nrows+=1
  assert len(row)==len(specimens)+1,(nrows,len(row))
  if row[0] in probe_gene:profiles[row[0]]=np.array([float(v) for v in row[1:]])
assert set(profiles)==set(probe_gene)
expr=pd.DataFrame({'expression_specimen_id':specimens})
for g in GENES:expr[g]=np.mean([profiles[p] for p in mapping.loc[mapping.gene.eq(g),'probe']],axis=0)
assert np.isfinite(expr[GENES].to_numpy()).all()
c=pd.read_csv(SRC/'GSE108474_clinical.data.txt.gz',sep='\t',dtype=str,keep_default_na=False)
b=pd.read_csv(SRC/'GSE108474_biospecimen_mapping_GEO.txt.gz',sep='\t',dtype=str,keep_default_na=False)
assert c.SUBJECT_ID.is_unique and b.SUBJECT_ID.is_unique and b.BIOSPECIMEN_ID.is_unique
b['expression_specimen_id']=b.BIOSPECIMEN_ID.map(lambda s:s[1:] if s.startswith('X') and s[1:] in set(specimens) else s)
assert set(b.expression_specimen_id).issubset(specimens)
unmapped_specimens=sorted(set(specimens)-set(b.expression_specimen_id))
expr=expr[expr.expression_specimen_id.isin(b.expression_specimen_id)]
d=c.merge(b,on='SUBJECT_ID',how='left',validate='one_to_one')
assert d.loc[d.expression_specimen_id.notna(),'expression_specimen_id'].is_unique
# Multiple patients legitimately have no mapped expression; retain those NA rows
# for the attrition denominator without treating them as duplicated specimens.
d=d.merge(expr,on='expression_specimen_id',how='left',validate='many_to_one')
assert len(d)==671
d['patient_id']=d.SUBJECT_ID;d['sample_id']=d.expression_specimen_id
bands=d.AGE_RANGE.str.extract(r'^(\d+)-(\d+)$').astype(float)
d['age_band_lower']=bands[0];d['age_band_upper']=bands[1];d['age']=(bands[0]+bands[1])/2
d['sex']=d.GENDER.str.strip().str.upper().map({'MALE':'male','M':'male','FEMALE':'female','F':'female'})
d['grade']=d.WHO_GRADE.map({'II':2,'III':3,'IV':4})
d['kps']=pd.to_numeric(d.KARNOFSKY.str.rstrip('%'),errors='coerce')
d['event']=d.EVENT_OS.map({'1':1,'0':0})
d['time_months']=pd.to_numeric(d.OVERALL_SURVIVAL_MONTHS,errors='coerce')
d['histology']=d.DISEASE_TYPE
d['specimen_timing']='not documented in the available patient/specimen map'
d['recorded_prior_radiotherapy']=d.PRIOR_THERAPY_RADIATION_TYPE.str.strip().replace('',np.nan).notna()
d['recorded_prior_chemotherapy']=d.PRIOR_THERAPY_CHEMO_AGENT_NAME.str.strip().replace('',np.nan).notna()
d['exclusion_reason']=''
rules=[(d.sample_id.isna(),'no expression biospecimen'),(~d.histology.isin(['GBM','ASTROCYTOMA','OLIGODENDROGLIOMA','MIXED'])|d.grade.isna(),'histology or grade outside documented diffuse II-IV'),(~d.age_band_lower.ge(18),'adult age cannot be established'),(~d.time_months.gt(0)|d.event.isna(),'unresolved or nonpositive overall survival')]
for mask,reason in rules:d.loc[mask&d.exclusion_reason.eq(''),'exclusion_reason']=reason
d['eligible']=d.exclusion_reason.eq('')
assert not d.loc[d.eligible,'patient_id'].duplicated().any()
assert np.isfinite(d.loc[d.eligible,GENES].to_numpy()).all()
d.to_csv(OUT/'GSE108474_audited_all_patients.csv',index=False)
d[d.eligible].to_csv(OUT/'GSE108474_eligible_timing_uncertain.csv',index=False)
audit={'created_utc':datetime.now(timezone.utc).isoformat(),'stage':'Input preparation; no REMBRANDT score-survival association calculated',
 'clinical_patients':len(c),'expression_columns':len(specimens),'mapped_expression_patients':len(b),'unmapped_expression_specimens_excluded':unmapped_specimens,'expression_probe_rows':nrows,'nine_gene_probes':len(profiles),
 'eligible_patients':int(d.eligible.sum()),'eligible_events':int(d.loc[d.eligible,'event'].sum()),
 'disjoint_exclusions':d.loc[~d.eligible,'exclusion_reason'].value_counts().to_dict(),
 'grade_counts':d.loc[d.eligible,'grade'].value_counts().to_dict(),
 'sex_completeness':int(d.loc[d.eligible,'sex'].notna().sum()),'kps_completeness':int(d.loc[d.eligible,'kps'].notna().sum()),
 'age_note':'Five-year band midpoints; ambiguous 15-19 bands excluded, not assigned an invented exact age',
 'event_basis':'Primary EVENT_OS event-indicator field; dictionary defines event as death and OS as time to death/last follow-up; 1=event, 0=no recorded event',
 'timing_note':'No patient/specimen timing is provided. This is an external evaluation with uncertain sampling timing, not confirmed pretreatment validation.',
 'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [SRC/'GSE108474_clinical.data.txt.gz',SRC/'GSE108474_biospecimen_mapping_GEO.txt.gz',SRC/'GSE108474_GeneExpression.txt.gz',SRC/'GSE108474_data.dictionary.xlsx']}}
(OUT/'GSE108474_input_preparation_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
