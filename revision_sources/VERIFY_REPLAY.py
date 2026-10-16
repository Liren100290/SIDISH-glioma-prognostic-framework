"""Compare all retained numeric result tables, allowing only floating point noise."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
checks=[]
for ref in sorted((ROOT/'REFERENCE_TABLES').rglob('*.csv')):
    rel=ref.relative_to(ROOT/'REFERENCE_TABLES');p=ROOT/rel
    a=pd.read_csv(ref);b=pd.read_csv(p)
    assert list(a.columns)==list(b.columns) and a.shape==b.shape,str(rel)
    gap=0.
    for col in a.columns:
        av=pd.to_numeric(a[col],errors='coerce').astype(float);bv=pd.to_numeric(b[col],errors='coerce').astype(float)
        numeric=(pd.api.types.is_numeric_dtype(a[col]) or av.notna().sum()==a[col].notna().sum())
        if numeric:
            assert np.array_equal(av.isna(),bv.isna()),(str(rel),col,'missingness')
            valid=av.notna();assert np.allclose(av[valid],bv[valid],rtol=1e-10,atol=1e-12),(str(rel),col)
            if valid.any():gap=max(gap,float(np.max(abs(av[valid]-bv[valid]))))
        else:
            assert a[col].fillna('<MISSING>').astype(str).equals(b[col].fillna('<MISSING>').astype(str)),(str(rel),col)
    checks.append({'path':str(rel).replace('\\','/'),'rows':len(a),'maximum_numeric_difference':gap})
out={'all_reference_tables_match':True,'tables':len(checks),'checks':checks,'scope':'Same-data numerical replay, not independent validation.'}
(ROOT/'RUN_LOGS').mkdir(exist_ok=True)
(ROOT/'RUN_LOGS/Reference_table_comparison.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({'all_reference_tables_match':True,'tables':len(checks),'maximum_numeric_difference':max(c['maximum_numeric_difference'] for c in checks)},indent=2))
