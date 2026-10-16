"""Revision-stage patient robustness reconstruction from read-only author data.

Trained-matrix reconstruction and source-raw sensitivity remain separate.
Original score coefficients, genes and cutoff are never reselected here.
"""
from pathlib import Path
import argparse, sys, json, time
import numpy as np
import pandas as pd
import openpyxl
BASE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(BASE/'05_Analyses'/'_pydeps'))
import h5py
from scipy.sparse import csr_matrix
ROOT=BASE.parent/'SIDISH'/'SIDISH_全量交接包_20261015'
PROJECT=ROOT/'01_项目原件'/'01_SIDISH'/'GBmap_Extend_scRNA_TCGAbulk'
OUT=BASE/'05_Analyses'/'patient_robustness'

def readcol(g,k):
    x=g[k]
    if isinstance(x,h5py.Group):
        v=x['categories'].asstr()[:] if x['categories'].dtype.kind in ['O','S','U'] else x['categories'][:]
        return np.asarray(pd.Categorical.from_codes(x['codes'][:],v))
    return x.asstr()[:] if x.dtype.kind in ['O','S','U'] else x[:]

def make_unique(names):
    reserved=set(names);used=set();number={};out=[]
    for name in names:
        candidate=name
        if candidate in used:
            j=number.get(name,0)+1
            while name+'.'+str(j) in used or name+'.'+str(j) in reserved:j+=1
            candidate=name+'.'+str(j);number[name]=j
        used.add(candidate);out.append(candidate)
    return out

book=openpyxl.load_workbook(BASE/'00_Submitted'/'Supplementary_Data.xlsx',read_only=True,data_only=True)
it=iter(book['ST2_245Markers'].values);hdr=next(it)
reference=pd.DataFrame(list(it),columns=hdr).dropna(subset=['gene']).set_index('gene')
for col in ['patient_median_delta','patient_support_fraction']:
    reference[col]=pd.to_numeric(reference[col],errors='coerce')
book.close();genes=reference.index.astype(str).to_list()
assert len(genes)==245 and len(set(genes))==245
md=pd.read_csv(OUT/'audited_cell_metadata.csv.gz',usecols=['source_row','analysis_patient_id','annotation_level_3','SIDISH','author'])
mi=pd.MultiIndex.from_frame(md[['analysis_patient_id','annotation_level_3','SIDISH']])
group_codes,unique=pd.factorize(mi,sort=True)
groups=unique.to_frame(index=False);groups.columns=['patient_id','celltype','state']
counts=np.bincount(group_codes,minlength=len(groups));groups['n_cells']=counts
author_map=md[['analysis_patient_id','author']].drop_duplicates().set_index('analysis_patient_id')['author']
assert author_map.index.is_unique
groups['author']=groups.patient_id.map(author_map)
groups.to_csv(OUT/'expression_group_metadata.csv',index=False)

def compute(mode):
    start=time.perf_counter()
    path=PROJECT/'runs/os/trained_adata.h5ad' if mode=='trained' else ROOT/'02_原始数据及预处理/scRNA_SIDISH/GBMap_extend.h5ad'
    key='X' if mode=='trained' else 'raw/X'
    with h5py.File(path,'r') as f:
        x=f[key];nrows,ncols=map(int,x.attrs['shape']);ptr=x['indptr'][:]
        if mode=='trained':
            var=readcol(f['var'],f['var'].attrs['_index']).astype(str)
            rowmap=np.arange(nrows,dtype=np.int64)
        else:
            vg=f['raw/var'];ens=readcol(vg,vg.attrs['_index']).astype(str)
            names=readcol(vg,'feature_name').astype(str)
            names=np.where((names=='')|(names=='nan'),ens,names)
            var=np.asarray(make_unique(names.tolist()))
            rowmap=np.full(nrows,-1,dtype=np.int64)
            rowmap[md.source_row.to_numpy()]=np.arange(len(md))
        colidx=pd.Index(var).get_indexer(genes)
        assert np.all(colidx>=0),[genes[i] for i in np.flatnonzero(colidx<0)]
        assert len(np.unique(colidx))==245
        sums=np.zeros((len(groups),245),dtype=np.float64)
        linear=np.zeros_like(sums) if mode=='trained' else None
        block=20000 if mode=='source_raw' else 50000
        for a in range(0,nrows,block):
            b=min(a+block,nrows);selected=rowmap[a:b]>=0
            if not selected.any():continue
            l,u=int(ptr[a]),int(ptr[b])
            data=x['data'][l:u];indices=x['indices'][l:u]
            mat=csr_matrix((data,indices,ptr[a:b+1]-l),shape=(b-a,ncols))
            if not selected.all():mat=mat[selected]
            dest=rowmap[a:b][selected];gc=group_codes[dest]
            total=np.asarray(mat.sum(axis=1)).ravel().astype(np.float64)
            assert np.all(total>0) and np.all(mat.data>=0)
            small=mat[:,colidx].astype(np.float64)
            design=csr_matrix((np.ones(len(dest)),(gc,np.arange(len(dest)))),shape=(len(groups),len(dest)))
            if linear is not None:linear+=(design@small).toarray()
            small.data*=np.repeat(10000./total,np.diff(small.indptr))
            small.data=np.log1p(small.data)
            sums+=(design@small).toarray()
            if a==0 or b==nrows or (a//block)%10==0:
                print(f'{mode}: {b}/{nrows} source rows processed',flush=True)
    payload={'normlog_sums':sums,'genes':np.asarray(genes)}
    if linear is not None:payload['linear_sums']=linear
    np.savez_compressed(OUT/(mode+'_group_expression_sums.npz'),**payload)
    means=sums/counts[:,None]
    pd.concat([groups,pd.DataFrame(means,columns=genes)],axis=1).to_csv(OUT/(mode+'_group_expression_means.csv.gz'),index=False,compression='gzip')
    result={'source':str(path),'matrix_key':key,'mode':mode,'cells_included':int(counts.sum()),
       'genes':245,'patient_celltype_state_groups':len(groups),'elapsed_seconds':time.perf_counter()-start,
       'normalization':'Within-cell sum over all genes available in the specified matrix, target 10000, natural log1p',
       'scope':'Reconstruction from the trained post-decontX/HVG matrix' if mode=='trained' else 'Independent source-raw count sensitivity using the same trained cells and labels'}
    (OUT/(mode+'_expression_computation.json')).write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False),flush=True)
    return payload

def summarize(mode):
    z=np.load(OUT/(mode+'_group_expression_sums.npz'));reports={}
    for value in ['normlog_sums']+(['linear_sums'] if 'linear_sums' in z else []):
        summed=pd.DataFrame(z[value],columns=genes)
        key=groups[['patient_id','state']]
        agg=pd.concat([key,summed],axis=1).groupby(['patient_id','state'],sort=True).sum()
        nc=groups.groupby(['patient_id','state']).n_cells.sum()
        avg=agg.div(nc,axis=0)
        eligible=nc.unstack(fill_value=0).reindex(columns=['h','b'],fill_value=0)
        patients=eligible.index[(eligible.h>=10)&(eligible.b>=10)]
        delta=avg.xs('h',level='state').reindex(patients)-avg.xs('b',level='state').reindex(patients)
        delta.to_csv(OUT/(mode+'_'+value+'_paired_patient_delta.csv'))
        median=delta.median();fraction=(delta>0).mean()
        strict=(median>.03)&(fraction>.65)
        supportive_ex=(~strict)&(median>.01)&(fraction>.55)
        supportive_in=(~strict)&(median>.01)&(fraction>=.55)
        tab=pd.DataFrame({'median_delta':median,'support_fraction':fraction,'n_paired_patients':len(patients),
             'strict':strict,'supportive_exclusive':supportive_ex,'supportive_inclusive':supportive_in})
        tab=tab.join(reference[['patient_median_delta','patient_support_fraction','robustness_tier']])
        tab.to_csv(OUT/(mode+'_'+value+'_gate_comparison.csv'))
        reports[value]={'paired_patients':len(patients),'strict':int(strict.sum()),
            'supportive_exclusive_gt055':int(supportive_ex.sum()),'supportive_inclusive_ge055':int(supportive_in.sum()),
            'strict_genes':strict.index[strict].to_list(),
            'median_absolute_difference_from_submission':float((tab.median_delta-tab.patient_median_delta).abs().median()),
            'max_absolute_difference_from_submission':float((tab.median_delta-tab.patient_median_delta).abs().max()),
            'max_support_fraction_difference':float((tab.support_fraction-tab.patient_support_fraction).abs().max())}
    (OUT/(mode+'_gate_reconstruction_audit.json')).write_text(json.dumps(reports,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(reports,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['trained','source_raw'],required=True);p.add_argument('--summarize-only',action='store_true');a=p.parse_args()
    if not a.summarize_only:compute(a.mode)
    summarize(a.mode)
