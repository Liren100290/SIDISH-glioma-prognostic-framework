"""Revision-stage data preparation. This is not the original study code.

Input files are official GEO matrices and GPL570 annotation downloaded locally.
Probe mapping is outcome-independent: exact single-symbol matches, mean collapse.
All downloaded inputs and the submitted coefficient source are hash recorded.
No score-survival association is computed here.
"""
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone
import csv, gzip, hashlib, json, math, statistics
import openpyxl

HERE=Path(__file__).resolve().parent
REVISION=HERE.parent.parent
SOURCE=REVISION/'00_Submitted'/'Supplementary_Data.xlsx'
GENES=['GFPT2','KAZN','NTN1','PTPRO','SLC22A23','STK32A','STK32B','TCF7L2','TMEM132B']

def write_csv(path,rows,fields):
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

book=openpyxl.load_workbook(SOURCE,read_only=True,data_only=True)
it=iter(book['ST6_9GeneModel'].values);header=next(it)
coefs=[dict(zip(header,row)) for row in it]
assert {x['gene'] for x in coefs}==set(GENES)
write_csv(HERE/'locked_nine_gene_coefficients.csv',
          [{'gene':g,'coefficient':next(x['locked_coefficient'] for x in coefs if x['gene']==g)} for g in GENES],
          ['gene','coefficient'])
book.close()

mapping=defaultdict(list)
with gzip.open(HERE/'GPL570.annot.gz','rt',encoding='utf-8',errors='strict') as f:
    for line in f:
        if line.startswith('!platform_table_begin'):break
    for row in csv.DictReader(f,delimiter='\t'):
        g=row.get('Gene symbol','')
        if g in GENES:mapping[g].append(row['ID'])
assert set(mapping)==set(GENES)
probe_gene={p:g for g,probes in mapping.items() for p in probes}
write_csv(HERE/'nine_gene_probe_mapping.csv',
          [{'gene':g,'probe':p,'resource':'official GPL570 annotation; exact single-symbol assignment'} for g in GENES for p in mapping[g]],
          ['gene','probe','resource'])

audit={'created_utc':datetime.now(timezone.utc).isoformat(),
       'analysis_phase':'Input preparation; no score-survival association computed',
       'source_supplement_sha256':digest(SOURCE),
       'annotation_sha256':digest(HERE/'GPL570.annot.gz'),
       'annotation_url':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz',
       'probe_rule':'Exact single-symbol annotation only; average all assigned probes without outcome filtering',
       'additional_log_transform':False,
       'normalization_basis':'GSE43378 deposited data are log2 MAS5/global-scaled; GSE7696 deposited data are RMA normalized',
       'fixed_development_cutoff':-0.0269027400037899,'datasets':{}}

for acc in ['GSE43378','GSE7696']:
    meta=defaultdict(list);expr={}
    with gzip.open(HERE/(acc+'_series_matrix.txt.gz'),'rt',encoding='utf-8',errors='strict') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):break
            row=next(csv.reader([line],delimiter='\t'))
            if row:meta[row[0]].append(row[1:])
        reader=csv.reader(f,delimiter='\t');ids=next(reader)[1:]
        for row in reader:
            if not row or row[0].startswith('!series_matrix_table_end'):break
            if row[0] in probe_gene:
                expr[row[0]]=[float(v) for v in row[1:]]
    assert ids==meta['!Sample_geo_accession'][0]
    assert len(ids)==len(set(ids))
    assert set(expr)==set(probe_gene)
    assert all(len(v)==len(ids) and all(math.isfinite(x) for x in v) for v in expr.values())
    genes={g:[statistics.mean(expr[p][j] for p in mapping[g]) for j in range(len(ids))] for g in GENES}
    rows=[]
    for j,gsm in enumerate(ids):
        d={'sample_id':gsm,'title':meta['!Sample_title'][0][j]}
        for values in meta['!Sample_characteristics_ch1']:
            val=values[j]
            if ':' in val:
                key,value=val.split(':',1);d[key.strip().lower()]=value.strip()
        if acc=='GSE43378':
            event={'DEAD':1,'ALIVE':0}.get(d['outcome'])
            row={'sample_id':gsm,'patient_id':gsm,'title':d['title'],'age':d['age (years)'],
                 'sex':d['gender'],'grade':d['who grade'],'kps':d['kps'],
                 'histology':d['pathological diagnosis'],'time_months':float(d['survival time (days)'])/30.4375,
                 'event':event,'event_source':d['outcome'],'treatment':'not provided per patient',
                 'mgmt':'NA','disease_status':'primary GBM supported for grade-IV subset by source publication; other sampling status requires verification'}
            eligible=event is not None and float(row['time_months'])>0 and d['who grade'] in ['2','3','4']
        else:
            val=d.get('survival status','NA')
            event=int(val) if val in ['0','1'] else None
            row={'sample_id':gsm,'patient_id':d.get('patient',''),'title':d['title'],'age':d.get('age','NA'),
                 'sex':d.get('gender','NA'),'grade':'4' if d.get('disease status')=='GBM' else 'NA','kps':'NA',
                 'histology':d.get('disease status',''),'time_months':d.get('survival time in months','NA'),
                 'event':event,'event_source':'GEO survival status='+val+'; confirm 1=death coding before final use',
                 'treatment':d.get('treatment','NA'),'mgmt':d.get('mgmt status','NA'),
                 'disease_status':d.get('disease status','')}
            eligible=d.get('disease status')=='GBM' and event is not None and float(row['time_months'])>0
        row['candidate_eligible']=int(eligible)
        row.update({g:genes[g][j] for g in GENES});rows.append(row)
    write_csv(HERE/(acc+'_prepared.csv'),rows,list(rows[0]))
    keep=[x for x in rows if x['candidate_eligible']]
    ids_keep=[x['patient_id'] for x in keep]
    assert len(ids_keep)==len(set(ids_keep)),(acc,'duplicate eligible patient')
    audit['datasets'][acc]={'deposited_samples':len(rows),'candidate_eligible_n':len(keep),
       'events_if_coding_confirmed':sum(int(x['event']) for x in keep),
       'grades':dict((k,sum(x['grade']==k for x in keep)) for k in sorted({x['grade'] for x in keep})),
       'treatments':dict((k,sum(x['treatment']==k for x in keep)) for k in sorted({x['treatment'] for x in keep})),
       'nine_genes_measured':True,'ambiguous_probe_assignments_excluded':True,
       'unique_eligible_patient_ids':True,'matrix_sha256':digest(HERE/(acc+'_series_matrix.txt.gz')),
       'prepared_sha256':digest(HERE/(acc+'_prepared.csv')),
       'readiness':'provisional; source sampling/independence audit and GSE7696 event coding must be closed before manuscript use'}
(HERE/'input_preparation_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit['datasets'],ensure_ascii=False,indent=2))
