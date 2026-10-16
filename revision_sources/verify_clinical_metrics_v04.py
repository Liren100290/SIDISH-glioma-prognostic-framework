"""Independently verify censoring weights, Brier scores, net benefit and CV splits."""
from pathlib import Path
from datetime import datetime,timezone
import json
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent/'clinical_comparison_v04'
b=pd.read_csv(OUT/'out_of_fold_IPCW_Brier.csv')
n=pd.read_csv(OUT/'out_of_fold_IPCW_net_benefit.csv')
errors=[];checks=[]
for context in b.context.unique():
 d=pd.read_csv(OUT/(context+'_mean_out_of_fold_predictions.csv'))
 t=d.time_months.to_numpy();e=d.event.to_numpy();uniq=np.unique(t)
 left={};right={};g=1.
 for z in uniq:
  left[z]=g;g*=1-((t==z)&(e==0)).sum()/(t>=z).sum();right[z]=g
 folds=pd.read_csv(OUT/(context+'_crossvalidation_folds.csv'))
 assert folds.groupby(['repeat_id','patient_id']).size().eq(1).all()
 assert folds.groupby('repeat_id').size().eq(len(d)).all()
 assert folds.repeat_id.nunique()==10 and set(folds.fold)==set(range(1,6))
 for _,row in b[b.context.eq(context)].iterrows():
  horizon=row.month;gt=right[max(uniq[uniq<=horizon])] if np.any(uniq<=horizon) else 1.
  wc=np.array([1/left[x] if ev==1 and x<=horizon else 0 for x,ev in zip(t,e)])
  wn=(t>horizon).astype(float)/gt
  vals={}
  for model in ['clinical','combined']:
   prob=d[f'{model}_p{int(horizon)}'].to_numpy()
   val=float(np.mean(wc*(1-prob)**2+wn*prob**2));vals[model]=val
   errors.append(abs(val-row[f'{model}_IPCW_Brier']))
   for _,nb in n[n.context.eq(context)&n.month.eq(horizon)].iterrows():
    treat=prob>=nb.threshold
    v=float(np.mean(wc*treat)-np.mean(wn*treat)*nb.threshold/(1-nb.threshold))
    errors.append(abs(v-nb[model]))
  assert int(((e==0)&(t<=horizon)).sum())==row.early_censored
  checks.append(dict(context=context,month=int(horizon),n=len(d),clinical_Brier=vals['clinical'],combined_Brier=vals['combined']))
report={'created_utc':datetime.now(timezone.utc).isoformat(),'independent_IPCW_checks_passed':max(errors)<1e-12,
        'max_absolute_discrepancy':max(errors),'contexts':4,'horizon_checks':checks,'crossvalidation_patient_split_checks_passed':True,
        'scope':'Recomputes KM censoring product and all Brier/net-benefit values in Python from retained R out-of-fold predictions. Does not treat repeated CV ranges as confidence intervals.'}
assert report['independent_IPCW_checks_passed']
(OUT/'Independent_clinical_metric_verification_v04.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:report[k] for k in ['independent_IPCW_checks_passed','max_absolute_discrepancy','contexts','crossvalidation_patient_split_checks_passed']},indent=2))
