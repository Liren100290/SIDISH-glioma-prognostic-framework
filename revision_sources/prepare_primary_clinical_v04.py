"""Extract primary-source clinical and expression inputs without scoring outcomes.

Source PDFs, GEO matrices, annotations and original coefficients remain read-only.
All exclusions and any differences from the submitted cohort are explicit.
"""
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
import csv,gzip,re,json,hashlib,statistics
import numpy as np
import pandas as pd
import pdfplumber

HERE=Path(__file__).resolve().parent
SRC=HERE/'primary_sources_v04'
OUT=HERE/'primary_prepared_v04';OUT.mkdir(exist_ok=True)
GENES=['GFPT2','KAZN','NTN1','PTPRO','SLC22A23','STK32A','STK32B','TCF7L2','TMEM132B']
PDF=SRC/'GSE16011_00085472can092307-sup-stabs_1-6.pdf'
def clean(s):
 s='' if s is None else str(s)
 s=re.sub(r'\(cid:(?:2|3)\)',' ',s)
 s=s.replace('(cid:882)','-')
 return ' '.join(s.split())
def num(s):
 try:return float(clean(s).replace(',','.').rstrip('%'))
 except ValueError:return np.nan
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def platform_rows(p,compressed=False):
 op=gzip.open if compressed else open
 with op(p,'rt',encoding='utf-8') as f:
  for s in f:
   if s.startswith('!platform_table_begin'):break
  for row in csv.DictReader(f,delimiter='\t'):
   if row and not str(next(iter(row.values()))).startswith('!platform_table_end'):yield row

assert (HERE/'gse16011_reconstruction_protocol_v04.json').exists()
clinical=[];molecular=[];repeats=[]
with pdfplumber.open(PDF) as doc:
 for i,page in enumerate(doc.pages):
  if i not in list(range(10))+[18]+list(range(19,25)):continue
  page=page.dedupe_chars(tolerance=1,extra_attrs=())
  for table in page.extract_tables():
   for j,row in enumerate(table):
    a=[clean(s) for s in row]
    if not a or not a[0].isdigit():continue
    if i<10:
     assert len(a)==10,(i,j,a)
     clinical.append(dict(zip(['db_number','sex','histology','age_raw','kps_raw','surgery','radiotherapy','chemotherapy','vital_status','survival_years_raw'],a),source_page=i+1,source_row=j+1))
    elif i==18:
     assert len(a)==6
     repeats.append(dict(primary_db_number=a[0],primary_histology=a[1],repeat_db_number=a[3],repeat_histology=a[4],source_page=i+1))
    else:
     assert len(a)==5
     molecular.append(dict(zip(['db_number','chr1p','chr19q','idh1_raw','egfr_raw'],a),source_page=i+1))
c=pd.DataFrame(clinical);m=pd.DataFrame(molecular);r=pd.DataFrame(repeats)
assert len(c)==284 and c.db_number.is_unique,(len(c),c[c.db_number.duplicated()].to_dict('records'))
assert len(m)==284 and m.db_number.is_unique and len(r)==7
c.to_csv(OUT/'GSE16011_source_clinical_table1.csv',index=False)
m.to_csv(OUT/'GSE16011_source_molecular_table6.csv',index=False)
r.to_csv(OUT/'GSE16011_documented_repeat_pairs.csv',index=False)

meta=defaultdict(list);expression={}
entrez=defaultdict(set)
for row in platform_rows(HERE/'GPL570.annot.gz',True):
 if row.get('Gene symbol','') in GENES:
  eid=row.get('Gene ID',row.get('Entrez Gene',''))
  if eid and re.fullmatch(r'\d+',eid):entrez[row['Gene symbol']].add(eid)
assert set(entrez)==set(GENES) and all(len(x)==1 for x in entrez.values()),dict(entrez)
eid_gene={next(iter(eids)):g for g,eids in entrez.items()}
g8542={};mapping=[]
for row in platform_rows(SRC/'GPL8542_platform.soft.txt'):
 eid=row['ID'].removesuffix('_at')
 if eid in eid_gene:
  g8542[row['ID']]=eid_gene[eid]
  mapping.append(dict(gene=eid_gene[eid],entrez_id=eid,probe_id=row['ID'],original_annotation=row))
assert set(g8542.values())==set(GENES) and len(g8542)==9
with gzip.open(SRC/'GSE16011_series_matrix.txt.gz','rt',encoding='utf-8') as f:
 for s in f:
  if s.startswith('!series_matrix_table_begin'):break
  a=next(csv.reader([s],delimiter='\t'))
  if a:meta[a[0]].append(a[1:])
 reader=csv.reader(f,delimiter='\t');ids=next(reader)[1:]
 for row in reader:
  if row and row[0] in g8542:expression[g8542[row[0]]]=[float(v) for v in row[1:]]
assert set(expression)==set(GENES)
rows=[]
for j,gsm in enumerate(ids):
 title=meta['!Sample_title'][0][j]
 match=re.fullmatch(r'(?:glioma|control)\s+(\d+)',title,re.I)
 assert match,title
 a={'sample_id':gsm,'db_number':match[1],'source_title':title}
 for ar in meta['!Sample_characteristics_ch1']:
  if ':' in ar[j]:
   key,val=ar[j].split(':',1);a['geo_'+key.strip()]=val.strip()
 a.update({g:expression[g][j] for g in GENES});rows.append(a)
d=pd.DataFrame(rows).merge(c,on='db_number',how='left',validate='one_to_one').merge(m.drop(columns='source_page'),on='db_number',how='left',validate='one_to_one')
assert len(d)==284 and d.source_page.notna().all()
d['age']=d.age_raw.map(num);d['kps']=d.kps_raw.map(num)
d['time_months']=d.survival_years_raw.map(num)*12
d['event']=d.vital_status.map({'Dead':1,'Alive':0,'Lost to follow-up':0})
d['grade']=d.histology.str.extract(r'\(grade\s+(IV|III|II|I)\)')[0].map({'II':2,'III':3,'IV':4})
d['idh1']=d.idh1_raw.map({'mutation':'mutant','no mutation':'wildtype'})
d['codeletion']=np.where((d.chr1p=='LOH')&(d.chr19q=='LOH'),'codeleted',np.where((d.chr1p=='no LOH')|(d.chr19q=='no LOH'),'not_codeleted',''))
d['rt_record']=d.radiotherapy.str.lower().map({'yes':'yes','no':'no'})
d['chemo_record']=d.chemotherapy.map(lambda s: '' if not s else ('no' if s.lower()=='no' else 'yes'))
d['chemo_class']=d.chemotherapy.map(lambda s: '' if not s else ('none' if s.lower()=='no' else ('TMZ' if 'temozolomide' in s.lower() else ('PCV' if 'pcv' in s.lower() else 'other'))))
d['resection_complete']=d.surgery.map(lambda s:'' if not s else ('yes' if s=='Complete resection' else 'no'))
d['documented_repeat']=d.db_number.isin(r.repeat_db_number)
d['patient_id']=d.db_number
for _,pair in r.iterrows():d.loc[d.db_number==pair.repeat_db_number,'patient_id']=pair.primary_db_number
d['exclusion_reason']=''
rules=[(d.histology.str.lower().eq('control'),'control'),(~d.histology.str.match(r'^(?:A|OD|OA|GBM) \(grade (?:II|III|IV)\)$'),'not documented diffuse grade II-IV'),(d.documented_repeat,'documented repeat specimen'),(~d.age.ge(18),'age below 18 or unavailable'),(~d.time_months.gt(0)|d.event.isna(),'unresolved follow-up')]
for mask,reason in rules:d.loc[mask & d.exclusion_reason.eq(''),'exclusion_reason']=reason
d['eligible']=d.exclusion_reason.eq('')
age_checks=[];sex_checks=[];hist_checks=[]
for _,a in d.iterrows():
 if a.histology.lower()=='control':continue
 age_checks.append(abs(num(a.get('geo_age at diagnosis',''))-a.age))
 sex_checks.append(str(a.get('geo_gender','')).lower()==a.sex.lower())
 hist_checks.append(str(a.get('geo_histology',''))==a.histology)
assert max(age_checks)<0.011 and all(sex_checks) and all(hist_checks),(max(age_checks),sum(sex_checks),sum(hist_checks))
assert d.loc[d.eligible,'patient_id'].is_unique
d.to_csv(OUT/'GSE16011_audited_all_samples.csv',index=False)
d[d.eligible].to_csv(OUT/'GSE16011_eligible_revision_audit.csv',index=False)
(OUT/'GSE16011_gene_annotation_mapping.json').write_text(json.dumps(mapping,ensure_ascii=False,indent=2),encoding='utf-8')
audit={'created_utc':datetime.now(timezone.utc).isoformat(),'stage':'Input preparation; no new score-survival association calculated',
       'clinical_table_rows':len(c),'molecular_table_rows':len(m),'documented_repeat_pairs':len(r),
       'geo_samples':len(d),'eligible_revision_audit':int(d.eligible.sum()),'eligible_events':int(d.loc[d.eligible,'event'].sum()),
       'disjoint_exclusions':d.loc[~d.eligible,'exclusion_reason'].value_counts().to_dict(),
       'grade_counts':d.loc[d.eligible,'grade'].value_counts().to_dict(),'status_counts':d.loc[d.eligible,'vital_status'].value_counts().to_dict(),
       'source_gender_histology_agree_with_geo':True,'age_max_abs_difference':max(age_checks),
       'covariate_completeness':{x:int(d.loc[d.eligible,x].replace('',np.nan).notna().sum()) for x in ['age','sex','grade','kps','idh1','codeletion','rt_record','chemo_record','resection_complete']},
       'source_pdf_sha256':sha(PDF),'matrix_sha256':sha(SRC/'GSE16011_series_matrix.txt.gz'),
       'historical_input_identity':'Not established; do not force the source reconstruction to match the submitted n=240 or describe this as recovered historical code.',
       'clinical_timing_limits':'Post-diagnosis therapy records lack start dates; recorded-treatment sensitivity cannot establish treatment independence or causal effects.'}
(OUT/'GSE16011_input_preparation_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
