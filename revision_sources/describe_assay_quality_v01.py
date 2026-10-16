"""Post-analysis descriptive quality audit; no outcome-dependent filtering/refit."""
from pathlib import Path
from datetime import datetime,timezone
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
OUT=HERE/'results_v01'
r=pd.read_csv(OUT/'all_candidate_glioma_response_rows.csv')
m=pd.read_csv(HERE/'frozen_candidate_compound_screen_mapping.csv')
rows=[]
for _,d in m.iterrows():
 x=r[r.broad_id.eq(d.broad_id)&r.screen_id.eq(d.screen_id)]
 rows.append({'original_candidate':d.original_candidate,'broad_id':d.broad_id,'screen_id':d.screen_id,
  'preferred_screen':bool(d.primary_selected),'n_model_rows':len(x),'median_fitted_AUC':x.auc.median(),
  'AUC_ge1_count':int(x.auc.ge(1).sum()),'AUC_lt0_count':int(x.auc.lt(0).sum()),
  'median_source_curve_r2':x.r2.median(),'r2_le0_count':int(x.r2.le(0).sum()),
  'missing_r2_count':int(x.r2.isna().sum()),
  'interpretation':'Descriptive fit/dynamic-range audit added after viewing the new result; no observations removed, winsorized or recoded.'})
q=pd.DataFrame(rows)
q.to_csv(OUT/'assay_quality_descriptive_postanalysis.csv',index=False)
(OUT/'Assay_quality_note_v01.json').write_text(json.dumps({
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'timing':'Post-analysis descriptive audit, transparently identified; primary/sensitivity estimands and all observations unchanged.',
 'finding':'Several preferred MEK/BRAF compound curves have median fitted viability AUC near1 and low source fit r2, limiting dynamic range and interpretation of null score associations.',
 'handling':'Fitted AUC above1 is retained without clipping. No new r2 threshold or favorable-screen substitution was introduced. Lack of support is not proof of no pharmacologic effect.',
 'dataset_window':'Official PRISM secondary screen uses8 four-fold dilutions starting from10 micromolar; compound-specific coverage and source protocols remain part of interpretation.'},indent=2),encoding='utf-8')
print(q[q.preferred_screen][['original_candidate','n_model_rows','median_fitted_AUC','median_source_curve_r2','AUC_ge1_count']].to_string(index=False))
