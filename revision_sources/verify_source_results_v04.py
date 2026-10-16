"""Independent locked-score arithmetic and pairwise concordance verification."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib
import numpy as np
import pandas as pd
import openpyxl

HERE=Path(__file__).resolve().parent;OUT=HERE/'source_results_v04'
w=openpyxl.load_workbook(HERE.parent.parent/'00_Submitted/Supplementary_Data.xlsx',read_only=True,data_only=True)
it=iter(w['ST6_9GeneModel'].values);header=next(it);raw=[dict(zip(header,row)) for row in it];w.close()
beta={r['gene']:float(r['locked_coefficient']) for r in raw};genes=list(beta)
metrics=pd.read_csv(OUT/'locked_score_performance.csv')
parents={};checks=[]
for parent in metrics.parent.unique():
 d=pd.read_csv(OUT/(parent+'_all_scored.csv'))
 x=d[genes].to_numpy(float);score=((x-x.mean(axis=0))/x.std(axis=0,ddof=1))@np.array([beta[g] for g in genes])
 discrepancy=float(np.max(np.abs(score-d.score.to_numpy())))
 assert discrepancy<1e-12,(parent,discrepancy)
 parents[parent]={'score_max_absolute_discrepancy':discrepancy,'n':len(d),'events':int(d.event.sum()),'parent_score_sd':float(score.std(ddof=1))}
for _,m in metrics.iterrows():
 d=pd.read_csv(OUT/(m.context+'_evaluated.csv'))
 t=d.time_months.to_numpy();e=d.event.to_numpy();s=d.score.to_numpy()
 credit=0.;pairs=0
 for i in range(len(d)):
  for j in range(i+1,len(d)):
   if t[i]<t[j] and e[i]==1:a,b=i,j
   elif t[j]<t[i] and e[j]==1:a,b=j,i
   elif t[i]==t[j] and e[i]!=e[j]:a,b=(i,j) if e[i]==1 else (j,i)
   else:continue
   pairs+=1;credit+=1 if s[a]>s[b] else (.5 if s[a]==s[b] else 0)
 c=credit/pairs
 assert abs(c-m.Cindex)<1e-12,(m.context,c,m.Cindex)
 assert int((s>m.cutoff).sum())==m.high
 parent=parents[m.parent]
 assert abs(d.score_sd.std(ddof=1)-d.score.std(ddof=1)/parent['parent_score_sd'])<1e-12
 checks.append({'context':m.context,'R_Cindex':m.Cindex,'independent_pairwise_Cindex':c,'comparable_pairs':pairs})
report={'created_utc':datetime.now(timezone.utc).isoformat(),'coefficient_source':'Submitted ST6_9GeneModel; all nine weights independently reread',
        'all_score_checks_passed':True,'all_concordance_checks_passed':True,'subgroup_reference_frozen':True,'parents':parents,'contexts':checks,
        'scope':'Numerical implementation verification; source/timing/independence limits remain separate',
        'evaluation_script_sha256':hashlib.sha256((HERE/'evaluate_source_cohorts_v04.R').read_bytes()).hexdigest()}
(OUT/'Independent_numeric_verification_v04.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:report[k] for k in ['all_score_checks_passed','all_concordance_checks_passed','subgroup_reference_frozen','parents']},ensure_ascii=False,indent=2))
