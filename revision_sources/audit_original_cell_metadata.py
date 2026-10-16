"""Read-only source audit. All derived outputs are separate from author originals."""
from pathlib import Path
from collections import Counter
import sys, json, re, hashlib
import numpy as np
import pandas as pd

BASE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(BASE/'05_Analyses'/'_pydeps'))
import h5py
HANDOFF=BASE.parent/'SIDISH'/'SIDISH_全量交接包_20261015'
PROJECT=HANDOFF/'01_项目原件'/'01_SIDISH'/'GBmap_Extend_scRNA_TCGAbulk'
OUT=BASE/'05_Analyses'/'patient_robustness'
OUT.mkdir(exist_ok=True)

def column(group,key):
    obj=group[key]
    if isinstance(obj,h5py.Group):
        vals=obj['categories'].asstr()[:] if obj['categories'].dtype.kind in ['O','S','U'] else obj['categories'][:]
        codes=obj['codes'][:]
        return pd.Categorical.from_codes(codes,vals)
    return obj.asstr()[:] if obj.dtype.kind in ['O','S','U'] else obj[:]

source=HANDOFF/'02_原始数据及预处理'/'scRNA_SIDISH'/'GBMap_extend.h5ad'
trained=PROJECT/'runs'/'os'/'trained_adata.h5ad'
with h5py.File(trained,'r') as f:
    obs=f['obs'];ids=column(obs,obs.attrs['_index'])
    old_pat=np.asarray(column(obs,'patient_id_sc')).astype(str)
    sample=np.asarray(column(obs,'sample_id_sc')).astype(str)
    labels=np.asarray(column(obs,'SIDISH')).astype(str)
    old_types=np.asarray(column(obs,'celltype_major')).astype(str)
    assert len(ids)==len(set(ids))
    prefixes=sample+'_'  # exact prefix defined by the saved exporter metadata
    prefix_ok=np.fromiter((x.startswith(p) for x,p in zip(ids,prefixes)),dtype=bool,count=len(ids))
    assert prefix_ok.all(), 'Do not strip arbitrary strings when matching cell IDs'
    source_ids=np.asarray([x[len(p):] for x,p in zip(ids,prefixes)],dtype=object)
    assert len(source_ids)==len(set(source_ids))

with h5py.File(source,'r') as f:
    obs=f['obs'];source_index=pd.Index(column(obs,obs.attrs['_index']))
    assert source_index.is_unique
    # Reproduce the exact forward identifier transformation in the supplied
    # preprocessing R script (lines 597-607). Signs and parentheses can change.
    # Requiring a unique forward map avoids guessing reverse string edits.
    source_donors=np.asarray(column(obs,'donor_id'))
    def safe(s,pattern):return re.sub(pattern,'-',re.sub(r'\s+','-',str(s)))
    forward=pd.Index(['Disease.'+safe(d,'[^A-Za-z0-9_.-]+')+'_'+safe(s,'[^A-Za-z0-9_.:-]+')
                      for d,s in zip(source_donors,source_index)])
    assert forward.is_unique,'Sanitization collisions require explicit make.unique reconstruction'
    literal_coverage=float(np.mean(source_index.get_indexer(source_ids)>=0))
    row_index=forward.get_indexer(ids)
    assert np.all(row_index>=0),'Some trained cells are absent from the source atlas'
    source_ids=source_index.to_numpy()[row_index]
    md={'trained_cell_id':ids,'source_cell_id':source_ids,'source_row':row_index,
        'old_patient_id':old_pat,'old_sample_id':sample,'SIDISH':labels,'old_celltype':old_types}
    for key in ['donor_id','author','assay','gbmap','annotation_level_1','annotation_level_2','annotation_level_3','stage','sector','location','disease']:
        md[key]=np.asarray(column(obs,key))[row_index]
    d=pd.DataFrame(md)
    raw=f['raw/X/data'][:100000];normalized=f['X/data'][:100000]
    source_expression={'raw_shape':list(map(int,f['raw/X'].attrs['shape'])),
        'raw_sample_integer_fraction':float(np.mean(np.isclose(raw,np.round(raw)))),
        'X_sample_integer_fraction':float(np.mean(np.isclose(normalized,np.round(normalized)))),
        'raw_sample_range':[float(raw.min()),float(raw.max())],
        'X_sample_range':[float(normalized.min()),float(normalized.max())],
        'raw_and_X_structures_have_equal_nnz':len(f['X/data'])==len(f['raw/X/data']),
        'source_citation':f['uns/citation'].asstr()[()]}

assert not d[['donor_id','author','SIDISH']].isna().any().any()
old_to_donor=d.groupby('old_patient_id',observed=True)['donor_id'].nunique()
donor_to_old=d.groupby('donor_id',observed=True)['old_patient_id'].nunique()
donor_to_author=d.groupby('donor_id',observed=True)['author'].nunique()
cross_author=donor_to_author[donor_to_author>1]
d['analysis_patient_id']=d['author'].astype(str)+'::'+d['donor_id'].astype(str)
if len(cross_author)==0:d['analysis_patient_id']=d['donor_id'].astype(str)
counts=d.groupby(['analysis_patient_id','SIDISH'],observed=True).size().unstack(fill_value=0)
counts=counts.reindex(columns=['h','b'],fill_value=0)
counts.columns=['n_HR','n_BG']
counts['eligible_ge10']=(counts.n_HR>=10)&(counts.n_BG>=10)
pt=d[['analysis_patient_id','donor_id','author','assay','old_patient_id']].drop_duplicates()
counts=counts.join(pt.groupby('analysis_patient_id',observed=True).agg(donor_id=('donor_id','first'),author=('author','first'),n_assays=('assay','nunique'),n_old_groups=('old_patient_id','nunique')))
counts.to_csv(OUT/'patient_eligibility.csv')
ct=d.groupby(['analysis_patient_id','annotation_level_3','SIDISH'],observed=True).size().unstack(fill_value=0).reindex(columns=['h','b'],fill_value=0)
ct.columns=['n_HR','n_BG'];ct['eligible_ge10']=(ct.n_HR>=10)&(ct.n_BG>=10)
ct.to_csv(OUT/'patient_celltype_eligibility.csv')
summaries={}
for col in ['author','assay','annotation_level_1','annotation_level_2','annotation_level_3','stage']:
    t=pd.crosstab(d[col],d['SIDISH']);t.to_csv(OUT/(col+'_risk_counts.csv'))
    summaries[col]={'levels':len(t),'counts':t.to_dict(orient='index')}
d.to_csv(OUT/'audited_cell_metadata.csv.gz',index=False,compression='gzip')
np.save(OUT/'source_row_index_in_trained_order.npy',row_index)
report={'source_h5ad':str(source),'trained_h5ad':str(trained),
    'cells':len(d),'risk_counts':dict(Counter(labels)),
    'exact_source_cell_coverage':float(np.mean(row_index>=0)),
    'cell_id_rule':'Reproduce the supplied R script forward rule: Disease.<sanitized donor_id>_<sanitized original source cell ID>; require a unique forward map and 100% trained-cell coverage',
    'cell_id_transform_source':'GBMap author-annotation preprocessing R script, lines 597-607',
    'literal_prefix_removal_coverage':literal_coverage,
    'old_patient_groups':d.old_patient_id.nunique(),'source_donors':d.donor_id.nunique(),
    'analysis_patient_groups':d.analysis_patient_id.nunique(),
    'old_groups_mapping_to_multiple_donors':old_to_donor[old_to_donor>1].to_dict(),
    'donors_mapping_to_multiple_old_groups':donor_to_old[donor_to_old>1].to_dict(),
    'donors_found_in_multiple_studies':cross_author.to_dict(),
    'eligible_ge10_each_state':int(counts.eligible_ge10.sum()),
    'patients_with_no_HR':int((counts.n_HR==0).sum()),
    'patients_with_1_to_9_HR':int(((counts.n_HR>0)&(counts.n_HR<10)).sum()),
    'patients_with_BG_lt10':int((counts.n_BG<10).sum()),
    'eligible_celltype_pairs':int(ct.eligible_ge10.sum()),
    'source_expression':source_expression,'summaries':summaries}
(OUT/'cell_metadata_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,default=lambda x:int(x) if isinstance(x,np.integer) else str(x)),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k not in ['summaries','source_expression']},ensure_ascii=False,indent=2))
