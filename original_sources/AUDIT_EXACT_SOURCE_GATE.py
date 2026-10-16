"""Check the exact recovered Step1 inequalities and zero-fraction criterion.

Original saved sources and prior revision results remain untouched.
"""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
BASE=Path(__file__).resolve().parent.parent
ROOT=BASE/'05_Analyses/patient_robustness'
OUT=BASE/'05_Analyses/recovered_original_audit_v01'
if not ROOT.is_dir():
    ROOT=BASE/'retained_revision_inputs/source_gate'
    OUT=BASE/'verification_rerun/exact_source_gate'
OUT.mkdir(parents=True,exist_ok=True)
reports={};all_rows=[]
for mode in ['cleaned_full','source_raw']:
    p=ROOT/(mode+'_group_expression_means.csv.gz')
    data=pd.read_csv(p)
    genes=[c for c in data if c not in ['patient_id','celltype','state','n_cells','author']]
    assert len(genes)==245,genes[:10]
    n=data.groupby(['patient_id','state']).n_cells.sum()
    sums=data[genes].multiply(data.n_cells,axis=0)
    means=pd.concat([data[['patient_id','state']],sums],axis=1).groupby(['patient_id','state']).sum().div(n,axis=0)
    count=n.unstack(fill_value=0)
    eligible=count.index[(count.h>=10)&(count.b>=10)]
    h=means.xs('h',level='state').loc[eligible]
    b=means.xs('b',level='state').loc[eligible]
    zero=((h<=0)|(b<=0)).mean()
    old=pd.read_csv(ROOT/(mode+'_whole_donor_statistics.csv'),index_col=0)
    exact_strict=(old.median_delta>=.03)&(old.positive_fraction>=.65)&(zero.reindex(old.index)<.60)&(old.n_donors>=5)
    exact_support=(~exact_strict)&(old.median_delta>.01)&(old.positive_fraction>=.55)&(zero.reindex(old.index)<.60)&(old.n_donors>=5)
    zero_support=(old.median_delta>.01)&(old.positive_fraction>=.55)&(zero.reindex(old.index)>=.60)&(old.n_donors>=5)
    table=pd.DataFrame(dict(gene=old.index,mode=mode,n_donors=old.n_donors,
      median_delta=old.median_delta,positive_fraction=old.positive_fraction,
      zero_fraction=zero.reindex(old.index),prior_revision_strict=old.strict,
      exact_source_strict=exact_strict,prior_revision_supportive=old.supportive_ge055,
      exact_source_supportive=exact_support,source_zero_inflated_supportive=zero_support))
    table['source_exact_tier']=np.select([exact_strict,exact_support,zero_support],['strict','supportive','zero_inflated_supportive'],default='hold')
    table['prior_submitted_tier']=old.robustness_tier.to_numpy()
    table.to_csv(OUT/(mode+'_exact_source_gate_audit.csv'),index=False)
    report=dict(donors=len(eligible),genes=len(genes),
      exact_source_strict=int(exact_strict.sum()),exact_source_supportive=int(exact_support.sum()),
      zero_inflated_supportive=int(zero_support.sum()),
      strict_differences_from_previous_revision=int((old.strict!=exact_strict).sum()),
      supportive_differences_from_previous_revision=int((old.supportive_ge055!=exact_support).sum()),
      strict_boundary_genes=old.index[(old.median_delta==.03)|(old.positive_fraction==.65)].tolist(),
      gene_zero_fraction_range=[float(zero.min()),float(zero.max())],
      input_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    if mode=='cleaned_full':
        source_strict=old.robustness_tier.str.contains('2A_strict')
        source_support=old.robustness_tier.str.contains('2B_direction')
        report['exact_gate_matches_all245_submitted_tiers']=bool((source_strict==exact_strict).all() and (source_support==exact_support).all())
    reports[mode]=report;all_rows.append(table)
pd.concat(all_rows).to_csv(OUT/'exact_original_patient_gate_audit_all_modes.csv',index=False)
result=dict(source_rule=dict(min_paired_donors=5,min_cells_per_state=10,
  strict='median_delta >=0.03 and positive_fraction >=0.65',
  supportive='not strict; median_delta >0.01 and positive_fraction >=0.55',
  primary_zero_filter='Fraction of donors with HR mean <=0 or BG mean <=0 must be <0.60; the source also permits missing zero_fraction.',
  zero_inflated_supportive='Supportive expression thresholds with zero_fraction >=0.60 are a separate sensitivity tier, not primary metaprogram input.'),
  results=reports,numerical_or_label_outputs_changed=any(
    r['strict_differences_from_previous_revision'] or r['supportive_differences_from_previous_revision'] for r in reports.values()),
  prior_revision_implementation='Used strict > boundaries and omitted the zero filter. The actual fixed inputs are checked explicitly here; equivalence is not asserted for arbitrary new inputs.')
(OUT/'Exact_original_patient_gate_audit_v02.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
