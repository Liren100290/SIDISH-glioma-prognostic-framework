"""Additional fixed-marker assay matching and sampling-stage description."""
from pathlib import Path
import json,sys
import pandas as pd
import numpy as np
BASE=Path(__file__).resolve().parent.parent;OUT=BASE/'05_Analyses/patient_robustness'
sys.path.insert(0,str(BASE/'05_Analyses/_pydeps'))
from scipy.stats import wilcoxon
NINE=['GFPT2','KAZN','NTN1','PTPRO','SLC22A23','STK32A','STK32B','TCF7L2','TMEM132B']

def statistics(delta):
    x=delta.to_numpy().copy();x[np.abs(x)<=1e-12]=0.;n=len(x)
    rng=np.random.default_rng(20261015);boots=[]
    for _ in range(10):boots.append(np.median(x[rng.integers(0,n,size=(100,n)),:],axis=1))
    lo,hi=np.quantile(np.vstack(boots),[.025,.975],axis=0)
    p=np.asarray([1. if np.all(v==0) else wilcoxon(v,zero_method='wilcox',method='approx').pvalue for v in x.T])
    order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return pd.DataFrame({'gene':delta.columns,'n_donors':n,'median_delta':np.median(x,axis=0),
      'positive_fraction':(x>0).mean(axis=0),'bootstrap_ci_low':lo,'bootstrap_ci_high':hi,
      'wilcoxon_p_exploratory':p,'bh_q_exploratory':q}).set_index('gene')

df=pd.read_csv(OUT/'cleaned_assay_group_expression_means.csv.gz')
genes=[c for c in df.columns if c not in ['patient_id','celltype','assay','state','n_cells','author']]
whole=pd.read_csv(OUT/'cleaned_full_whole_donor_delta.csv',index_col=0)
assert set(genes)==set(whole.columns) and len(genes)==245
keys=['patient_id','celltype','assay']
nc=df.pivot_table(index=keys,columns='state',values='n_cells',fill_value=0)
eligible=nc.index[(nc.h>=10)&(nc.b>=10)&nc.index.get_level_values('patient_id').isin(whole.index)]
expr=df.set_index(keys+['state'])[genes]
delta=expr.xs('h',level='state').reindex(eligible)-expr.xs('b',level='state').reindex(eligible)
delta.to_csv(OUT/'assay_matched_stratum_delta.csv')
paired=delta.groupby('patient_id').mean();paired.to_csv(OUT/'assay_matched_donor_delta.csv')
matched=statistics(paired);matched.to_csv(OUT/'assay_matched_statistics.csv')
same=statistics(whole.reindex(paired.index));same.to_csv(OUT/'assay_same_donors_unadjusted_statistics.csv')
meta=pd.read_csv(OUT/'donor_source_assay_stage.csv').rename(columns={'analysis_patient_id':'patient_id'})
loo=[]
for assay in sorted(meta.assay.unique()):
    excluded=meta.loc[meta.assay==assay,'patient_id'].unique()
    remain=whole.loc[~whole.index.isin(excluded)]
    if len(remain)==len(whole):continue  # No eligible donor carries this assay; not an exclusion test.
    if len(remain)<10:continue
    tab=pd.DataFrame({'gene':whole.columns,'excluded_assay':assay,
       'excluded_eligible_donors':len(whole)-len(remain),'remaining_donors':len(remain),
       'median_delta':remain.median().to_numpy(),'positive_fraction':(remain>0).mean().to_numpy()})
    loo.append(tab)
loo=pd.concat(loo,ignore_index=True);loo.to_csv(OUT/'leave_assay_out.csv',index=False)
summary=loo.groupby('gene').agg(comparisons=('excluded_assay','size'),median_min=('median_delta','min'),
   median_max=('median_delta','max'),positive_median_fraction=('median_delta',lambda x:float((x>1e-12).mean())))
summary.to_csv(OUT/'leave_assay_out_summary.csv')
na=meta.groupby('patient_id').assay.nunique();single=meta.loc[meta.patient_id.isin(na.index[na==1])].set_index('patient_id').assay
stratified=[]
for assay in sorted(single.unique()):
    ids=whole.index.intersection(single.index[single==assay])
    if len(ids)<10:continue
    tab=statistics(whole.reindex(ids));tab.insert(0,'assay',assay);stratified.append(tab.reset_index())
stratified=pd.concat(stratified,ignore_index=True);stratified.to_csv(OUT/'single_assay_donor_statistics.csv',index=False)
stage=meta[['patient_id','stage']].drop_duplicates().set_index('patient_id').stage
stages=[]
stage_counts=stage.reindex(whole.index).value_counts()
for label,n in stage_counts.items():
    if n<10:continue
    tab=statistics(whole.loc[stage.reindex(whole.index)==label]);tab.insert(0,'sampling_stage',label);stages.append(tab.reset_index())
pd.concat(stages,ignore_index=True).to_csv(OUT/'sampling_stage_donor_statistics.csv',index=False)
result={'parent_donors':len(whole),'source_donors_with_multiple_assays':int((na>1).sum()),
 'assays_with_eligible_donors':int(loo.excluded_assay.nunique()),
 'assay_matched_donors':len(paired),'assay_matched_strata':len(eligible),
 'nine_positive_medians':int((matched.loc[NINE].median_delta>1e-12).sum()),
 'nine_ci_above_zero':int((matched.loc[NINE].bootstrap_ci_low>1e-12).sum()),
 'nine_all_assay_leave_out_positive':int((summary.loc[NINE].positive_median_fraction==1).sum()),
 'sampling_stage_eligible_donor_counts':{str(k):int(v) for k,v in stage_counts.items()},
 'single_assay_eligible_strata':stratified.groupby('assay').n_donors.first().to_dict(),
 'nine_gene_assay_matched_statistics':matched.loc[NINE].reset_index().to_dict(orient='records'),
 'interpretation':'Exploratory post-selection sensitivity with fixed genes/labels. Pair within donor, level-3 type and assay; average eligible stratum deltas equally per donor. Leave-assay-out removes every donor carrying the excluded assay, including mixed-assay donors. These checks do not retrain SIDISH or prove the risk labels free of batch effects. Sampling-stage analysis concerns discovery cells and is not treatment adjustment of the downstream bulk survival model.'}
(OUT/'Assay_sensitivity_results_v01.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='nine_gene_assay_matched_statistics'},ensure_ascii=False,indent=2))
print(matched.loc[NINE].to_string())
