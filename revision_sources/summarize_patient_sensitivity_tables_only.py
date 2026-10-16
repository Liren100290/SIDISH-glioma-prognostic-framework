"""Donor-level, cell-type-matched and leave-study-out exploratory sensitivity.

The fixed 245 markers and cell risk labels are never selected again. All new
inference is post-selection sensitivity, not independent biological validation.
"""
from pathlib import Path
import json, sys
import numpy as np
import pandas as pd
import openpyxl
BASE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(BASE/'05_Analyses/_pydeps'))
from scipy.stats import wilcoxon
OUT=BASE/'05_Analyses/patient_robustness'
NINE=['GFPT2','KAZN','NTN1','PTPRO','SLC22A23','STK32A','STK32B','TCF7L2','TMEM132B']
wb=openpyxl.load_workbook(BASE/'00_Submitted/Supplementary_Data.xlsx',read_only=True,data_only=True)
it=iter(wb['ST2_245Markers'].values);cols=next(it)
ref=pd.DataFrame(list(it),columns=cols).dropna(subset=['gene']).set_index('gene');wb.close()
genes=ref.index.tolist()
for k in ['patient_median_delta','patient_support_fraction']:ref[k]=pd.to_numeric(ref[k],errors='coerce')

def bh(p):
    p=np.asarray(p);n=len(p);order=np.argsort(p);q=np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1]
    out=np.empty(n);out[order]=np.minimum(q,1);return out

def stats(delta,seed=20261015,zero_tolerance=0.):
    x=delta.to_numpy().copy()
    if zero_tolerance:x[np.abs(x)<=zero_tolerance]=0.
    n=len(x);rng=np.random.default_rng(seed);boots=[]
    for a in range(0,1000,100):
        rows=rng.integers(0,n,size=(min(100,1000-a),n))
        boots.append(np.median(x[rows,:],axis=1))
    lo,hi=np.quantile(np.vstack(boots),[.025,.975],axis=0)
    p=[]
    for j in range(x.shape[1]):
        v=x[:,j]
        p.append(1. if np.all(v==0) else float(wilcoxon(v,zero_method='wilcox',alternative='two-sided',method='approx').pvalue))
    return pd.DataFrame({'gene':delta.columns,'n_donors':n,'median_delta':np.median(x,axis=0),
      'positive_fraction':(x>0).mean(axis=0),'bootstrap_ci_low':lo,'bootstrap_ci_high':hi,
      'wilcoxon_p_exploratory':p,'bh_q_exploratory':bh(p)}).set_index('gene')

reports={};tables={}
for mode in ['cleaned_full','source_raw','trained']:
    df=pd.read_csv(OUT/(mode+'_group_expression_means.csv.gz'))
    weighted=df[genes].multiply(df.n_cells,axis=0)
    keys=['patient_id','state']
    nc=df.groupby(keys).n_cells.sum();nwide=nc.unstack(fill_value=0)
    eligible=nwide.index[(nwide.h>=10)&(nwide.b>=10)]
    avg=pd.concat([df[keys],weighted],axis=1).groupby(keys).sum().div(nc,axis=0)
    whole=avg.xs('h',level='state').reindex(eligible)-avg.xs('b',level='state').reindex(eligible)
    whole.to_csv(OUT/(mode+'_whole_donor_delta.csv'))
    tab=stats(whole).join(ref[['robustness_tier','patient_median_delta','patient_support_fraction']])
    tab['strict']=(tab.median_delta>.03)&(tab.positive_fraction>.65)
    tab['supportive_ge055']=(~tab.strict)&(tab.median_delta>.01)&(tab.positive_fraction>=.55)
    tab['supportive_gt055']=(~tab.strict)&(tab.median_delta>.01)&(tab.positive_fraction>.55)
    tab.to_csv(OUT/(mode+'_whole_donor_statistics.csv'));tables[mode+'_whole']=tab.reset_index()
    keys2=['patient_id','celltype'];sizes=df.pivot_table(index=keys2,columns='state',values='n_cells',fill_value=0)
    # Keep donor eligibility fixed at the parent whole-donor gate.
    pairs=sizes.index[(sizes.h>=10)&(sizes.b>=10)&sizes.index.get_level_values('patient_id').isin(eligible)]
    by=df.set_index(['patient_id','celltype','state'])[genes]
    contrasts=by.xs('h',level='state').reindex(pairs)-by.xs('b',level='state').reindex(pairs)
    contrasts.to_csv(OUT/(mode+'_patient_celltype_delta.csv'))
    matched=contrasts.groupby('patient_id').mean()  # Equal eligible-type weight within donor.
    matched.to_csv(OUT/(mode+'_equal_celltype_donor_delta.csv'))
    mt=stats(matched,zero_tolerance=1e-12);mt.to_csv(OUT/(mode+'_equal_celltype_statistics.csv'));tables[mode+'_matched']=mt.reset_index()
    same=stats(whole.reindex(matched.index));same.to_csv(OUT/(mode+'_same_donors_unadjusted_statistics.csv'))
    tables[mode+'_same_donors']=same.reset_index()
    counts=contrasts.groupby('celltype').size();type_tables=[]
    for ct,n in counts.items():
        if n<10:continue
        t=stats(contrasts.xs(ct,level='celltype'),zero_tolerance=1e-12);t.insert(0,'celltype',ct);type_tables.append(t.reset_index())
    tstats=pd.concat(type_tables,ignore_index=True);tstats.to_csv(OUT/(mode+'_celltype_statistics.csv'),index=False)
    tables[mode+'_types']=tstats
    author=df[['patient_id','author']].drop_duplicates().set_index('patient_id').author
    loo=[]
    for study,n in author.reindex(whole.index).value_counts().items():
        remaining=whole.loc[author.reindex(whole.index)!=study]
        if len(remaining)<10:continue
        t=pd.DataFrame({'gene':genes,'excluded_study':study,'excluded_donors':int(n),
             'remaining_donors':len(remaining),'median_delta':remaining.median().reindex(genes).to_numpy(),
             'positive_fraction':(remaining>0).mean().reindex(genes).to_numpy()})
        loo.append(t)
    loo=pd.concat(loo,ignore_index=True);loo.to_csv(OUT/(mode+'_leave_study_out.csv'),index=False)
    lsum=loo.groupby('gene').agg(leave_out_comparisons=('excluded_study','size'),
        median_delta_min=('median_delta','min'),median_delta_max=('median_delta','max'),
        positive_median_fraction=('median_delta',lambda x:float((x>0).mean())),
        positive_patient_fraction_min=('positive_fraction','min'),positive_patient_fraction_max=('positive_fraction','max'))
    lsum.to_csv(OUT/(mode+'_leave_study_out_summary.csv'));tables[mode+'_study']=lsum.reset_index()
    nine=tab.loc[NINE,['n_donors','median_delta','positive_fraction','bootstrap_ci_low','bootstrap_ci_high','bh_q_exploratory']].copy()
    nine=nine.join(mt.add_prefix('matched_')).join(lsum.add_prefix('study_'))
    nine.to_csv(OUT/(mode+'_nine_gene_summary.csv'));tables[mode+'_nine']=nine.reset_index()
    reports[mode]={'whole_donors':len(whole),'matched_donors':len(matched),'matched_donor_type_pairs':len(pairs),
      'studies_in_whole_donors':int(author.reindex(whole.index).nunique()),
      'type_pair_counts':{str(k):int(v) for k,v in counts.items()},
      'strict':int(tab.strict.sum()),'supportive_ge055':int(tab.supportive_ge055.sum()),'supportive_gt055':int(tab.supportive_gt055.sum()),
      'submission_max_median_difference':float((tab.median_delta-tab.patient_median_delta).abs().max()),
      'submission_max_fraction_difference':float((tab.positive_fraction-tab.patient_support_fraction).abs().max()),
      'all_marker_positive_median_whole':int((tab.median_delta>0).sum()),
      'all_marker_positive_median_matched':int((mt.median_delta>0).sum()),
      'nine_positive_median_whole':int((tab.loc[NINE].median_delta>0).sum()),
      'nine_positive_median_matched':int((mt.loc[NINE].median_delta>0).sum()),
      'nine_matched_ci_above_zero':int((mt.loc[NINE].bootstrap_ci_low>0).sum()),
      'nine_study_out_positive_all':int((lsum.loc[NINE].positive_median_fraction==1).sum())}
    print(mode,json.dumps(reports[mode],ensure_ascii=False),flush=True)

reports['interpretation']={
 'original_gate':'Reconstruction uses full post-QC/decontX transcriptome before selecting the fixed 245 markers; HVG and source-raw alternatives are sensitivity analyses, not alternative versions of the submitted analysis.',
 'composition':'Pair HR and background cells within the same true donor and level-3 annotation; average eligible type-specific mean-log-expression differences equally within donor. The 81-donor matched analysis has a different eligible population than the original 100-donor analysis; the same-donor unadjusted summary is supplied.',
 'study':'Author/study is a source-atlas annotation. Leave-study-out summaries measure dominance by one source study and do not prove removal of batch effects.',
 'selection':'Markers and SIDISH labels originated from these data. Bootstrap intervals, paired tests and BH-adjusted values are exploratory post-selection sensitivity, not independent confirmatory evidence.',
 'inference':'Patient is the unit. Bootstrap resamples whole patients 1000 times, seed 20261015; intervals are percentile intervals for median paired differences; Wilcoxon is two-sided with zero_method=wilcox and normal approximation; BH is across 245 markers within each reported analysis or cell type.',
 'numerical_ties':'The submitted whole-donor gate retains its exact original >0 definition for reproducibility. In new matched and type-specific analyses, absolute paired differences <=1e-12 log-expression units are treated as zero to prevent floating-point/decontX residuals around 1e-16 from being counted as direction evidence. This numerical safeguard was documented after inspecting residual-scale values; it does not change any locked clinical score coefficient, cutoff or original gene tier.'}
(OUT/'Patient_sensitivity_results_v01.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2),encoding='utf-8')
notes=pd.DataFrame([{'item':k,'description':v} for k,v in reports['interpretation'].items()])
print('Completed all fixed-marker CSV/JSON sensitivity summaries',flush=True)
