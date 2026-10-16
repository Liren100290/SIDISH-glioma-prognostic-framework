"""Recompute every drug test from source-derived expression/response subsets.

Does not consume retained scores, correlations, P values or confidence intervals.
The separate full-source preparation script records the original extraction.
"""
from pathlib import Path
from importlib.util import spec_from_file_location,module_from_spec
from datetime import datetime,timezone
import json,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
sys.path.insert(0,str(HERE.parent/'_pydeps'))
import numpy as np
import pandas as pd
s=spec_from_file_location('original_prism_methods',HERE/'prepare_and_analyze_prism_v01.py')
method=module_from_spec(s);s.loader.exec_module(method)
OUT=HERE/'results_v01';OUT.mkdir(exist_ok=True)
IN=HERE/'replay_inputs_v05'
def read(p):return pd.read_csv(p,float_precision='round_trip')

weights=pd.read_csv(ROOT/'05_Analyses/external_validation/locked_nine_gene_coefficients.csv')
genes=weights.gene.tolist();beta=weights.coefficient.to_numpy(float)
assert len(genes)==9 and weights.gene.is_unique
parent=read(IN/'source_glioma_expression_subset.csv')
assert len(parent)==40 and parent.depmap_id.is_unique and parent.candidate_before_expression.all()
assert not any(c in parent for c in ['score','score_per_parent_SD','expression_eligible','all9_finite'])
parent['all9_finite']=np.isfinite(parent[genes].to_numpy(float)).all(axis=1)
parent['expression_eligible']=parent.expression_merge.eq('both')&parent.all9_finite
parent.to_csv(OUT/'glioma_expression_eligibility.csv',index=False)
d=parent[parent.expression_eligible].copy().sort_values('depmap_id').reset_index(drop=True)
assert len(d)==39
x=d[genes].to_numpy(float);mean=x.mean(axis=0);sd=x.std(axis=0,ddof=1)
assert np.isfinite(sd).all() and (sd>0).all()
d['score']=((x-mean)/sd)@beta
d['score_per_parent_SD']=(d.score-d.score.mean())/d.score.std(ddof=1)
names=d.stripped_cell_line_name.str.replace(r'[^A-Za-z0-9]','',regex=True).str.upper()
d['documented_BRAF_model_excluded_in_sensitivity']=names.isin(['AM38','DBTRG05MG'])
d['documented_GBM']=d.lineage_sub_subtype.eq('glioblastoma')
d['documented_adult']=pd.to_numeric(d.age,errors='coerce').ge(18)
d.to_csv(OUT/'glioma_parent_expression_and_scores.csv',index=False)
source_columns=json.loads((IN/'source_gene_columns.json').read_text(encoding='utf-8'))
pd.DataFrame({'gene':genes,'source_column':[source_columns[g] for g in genes],
              'parent_mean_log2_TPM1':mean,'parent_sample_SD':sd,'locked_coefficient':beta}).to_csv(
                  OUT/'frozen_expression_reference_and_weights.csv',index=False)

mapping=pd.read_csv(HERE/'frozen_candidate_compound_screen_mapping.csv')
raw=read(IN/'source_observed_response_subset.csv')
assert len(raw)==687 and 'score' not in raw and 'eligible_response' not in raw
raw=raw.merge(mapping,on=['broad_id','screen_id'],how='left',validate='many_to_one')
model_names=dict(zip(d.depmap_id,d.ccle_name))
raw['metadata_name_consistent']=raw.ccle_name.eq(raw.depmap_id.map(model_names))
raw['response_identity_QC_passed']=raw.passed_str_profiling.astype(str).str.lower().eq('true')
raw['auc_finite']=np.isfinite(pd.to_numeric(raw.auc,errors='coerce'))
raw['eligible_response']=raw.metadata_name_consistent&raw.response_identity_QC_passed&raw.auc_finite
raw.to_csv(OUT/'all_candidate_glioma_response_rows.csv',index=False)
eligible=raw[raw.eligible_response].copy()
keycols=['original_candidate','broad_id','screen_id','depmap_id']
response=eligible.groupby(keycols,as_index=False).agg(auc=('auc','median'),source_rows=('auc','size'),
    source_auc_min=('auc','min'),source_auc_max=('auc','max'),source_r2_median=('r2','median'))
response=response.merge(d[['depmap_id','score','score_per_parent_SD','documented_GBM','documented_adult',
    'documented_BRAF_model_excluded_in_sensitivity']],on='depmap_id',how='left',validate='many_to_one')
assert len(response)==687
response.to_csv(OUT/'model_level_observed_responses_and_scores.csv',index=False)
contexts={
 'all_glioma':set(d.depmap_id),
 'documented_GBM':set(d.loc[d.documented_GBM,'depmap_id']),
 'documented_adult':set(d.loc[d.documented_adult,'depmap_id']),
 'excluding_two_documented_BRAF_models':set(d.loc[~d.documented_BRAF_model_excluded_in_sensitivity,'depmap_id'])}
tests=[]
for context,ids in contexts.items():
    for _,row in mapping.iterrows():
        sub=response[response.broad_id.eq(row.broad_id)&response.screen_id.eq(row.screen_id)&response.depmap_id.isin(ids)].sort_values('depmap_id')
        n=len(sub);ok=n>=10 and sub.score.nunique()>1 and sub.auc.nunique()>1
        rho=lo=hi=p=np.nan;boot=0
        if ok:rho,lo,hi,p,boot=method.association(sub.score.to_numpy(float),sub.auc.to_numpy(float),20261015+len(tests))
        tests.append({'context':context,'original_candidate':row.original_candidate,'broad_id':row.broad_id,
          'screen_id':row.screen_id,'n_models':n,'parent_reference_models':len(d),
          'primary_selected_screen':bool(row.primary_selected),'inferentially_eligible':ok,
          'Spearman_rho':rho,'bootstrap_lower95':lo,'bootstrap_upper95':hi,'permutation_p_two_sided':p,
          'bootstrap_valid':boot,'median_AUC':sub.auc.median() if n else np.nan,
          'min_AUC':sub.auc.min() if n else np.nan,'max_AUC':sub.auc.max() if n else np.nan,
          'source_r2_median':sub.source_r2_median.median() if n else np.nan,
          'hypothesized_direction':'negative rho: higher original score, lower measured viability AUC',
          'analysis_status':'Exploratory public pharmacologic screen; selected compounds and overlapping cell models prevent independent clinical-validation claims'})
tests=pd.DataFrame(tests);tests['BH_q_sensitivity_family']=method.bh(tests.permutation_p_two_sided)
primary=tests[tests.context.eq('all_glioma')&tests.primary_selected_screen].copy()
assert len(primary)==16
primary['BH_q_primary_16']=method.bh(primary.permutation_p_two_sided)
original=json.loads((HERE/'protocol_amendment_candidate_table_v01.json').read_text(encoding='utf-8'))['correct_original_candidates']
primary['original_order']=primary.original_candidate.map({name:i for i,name in enumerate(original)})
primary=primary.sort_values('original_order').drop(columns='original_order')
tests.to_csv(OUT/'all_context_screen_tests.csv',index=False)
primary.to_csv(OUT/'primary_16_candidate_results.csv',index=False)
pd.DataFrame([{'context':c,'models_before_drug_specific_missingness':len(ids)} for c,ids in contexts.items()]).to_csv(OUT/'predeclared_context_counts.csv',index=False)
receipt={'completed_utc':datetime.now(timezone.utc).isoformat(),'models':39,'primary_compounds':16,'tests':len(tests),
 'primary_significant':int(primary.BH_q_primary_16.lt(.05).sum()),
 'inputs':'Source-derived 40-model nine-gene expression and 687 original measured response rows; no retained scores or association results consumed.',
 'scope':'Numerical same-data replay; no new wet-laboratory experiment and no independent clinical-validation claim.'}
(OUT/'Subset_replay_receipt_v05.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt,indent=2),flush=True)
