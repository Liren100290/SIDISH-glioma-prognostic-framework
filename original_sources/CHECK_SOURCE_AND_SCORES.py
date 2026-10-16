"""Validate source files and retained-score arithmetic; does not retrain models."""
from pathlib import Path
import csv,hashlib,math,statistics
base=Path(__file__).resolve().parent.parent
with (base/'FILE_MANIFEST.csv').open(encoding='utf-8-sig',newline='') as handle:manifest=list(csv.DictReader(handle))
for row in manifest:
    path=base/row['path']
    assert path.is_file() and path.stat().st_size==int(row['bytes']),row['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==row['sha256'],row['path']
root=base/'original_project'
stage2=next(p for p in root.iterdir() if p.is_dir() and p.name.startswith('SIDISH_Step2_'))
stage3=next(p for p in root.iterdir() if p.is_dir() and p.name.startswith('SIDISH_Step3_'))
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as handle:return list(csv.DictReader(handle,delimiter=chr(9)))
coef={r['gene']:float(r['locked_coefficient']) for r in read(stage2/'99_locked_model/SIDISH_Step2_locked_model_coefficients.tsv')}
summary={r['metric']:r['value'] for r in read(stage2/'99_locked_model/SIDISH_Step2_locked_model_summary.tsv')}
cutoff=float(summary['locked_development_cutoff']);assert len(coef)==9
for label,matrix,scores in [('development',stage2/'01_development_data/Step2_development_matrix_zscored.tsv',stage2/'99_locked_model/SIDISH_Step2_locked_development_risk_scores.tsv'),('historical external',stage3/'01_external_validation_data/Step3_external_locked_gene_matrix_zscored.tsv',stage3/'99_final_external_validation/SIDISH_Step3_GSE16011_no_refit_risk_scores.tsv')]:
    data=read(matrix);old=read(scores)
    key=lambda r:((r.get('cohort') or r.get('external_cohort') or 'GSE16011'),r['sample_id'])
    lookup={key(r):r for r in old};assert len(lookup)==len(data)
    delta=0;groups=0;computed=[]
    for row in data:
        score=math.fsum(float(row[g])*c for g,c in coef.items());record=lookup[key(row)]
        delta=max(delta,abs(score-float(record['locked_risk_score'])))
        groups+=('High' if score>=cutoff else 'Low')!=record['locked_risk_group'];computed.append(score)
    assert delta<1e-11 and groups==0
    if label=='development':assert abs(statistics.median(computed)-cutoff)<1e-11
    print(label,'n=',len(data),'max score difference=',delta,'group differences=',groups)
print('Verified',len(manifest),'registered files. Original preprocessing, full model screening, bootstrap and drug training were not run by this check.')
