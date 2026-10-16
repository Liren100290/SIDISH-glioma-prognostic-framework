"""Recompute revision evidence from the included aggregate/clinical inputs.

Does not run the historical TCGA/CGGA model development or original drug engine.
Run with Python 3.12 and the recorded compatible R/package versions.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, os, subprocess, sys
sys.dont_write_bytecode=True

ROOT = Path(__file__).resolve().parent

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''):h.update(block)
    return h.hexdigest()

def manifest_check():
    from ASSEMBLE_INPUT import ensure_input
    ensure_input()
    m=json.loads((ROOT/'FILE_MANIFEST.json').read_text(encoding='utf-8'))
    errors=[]
    for rec in m['files']:
        p=ROOT/rec['path']
        if not p.is_file() or sha(p)!=rec['sha256']:errors.append(rec['path'])
    if errors:raise RuntimeError('Package files differ from manifest: '+repr(errors))
    print('Manifest verified:',len(m['files']),'files',flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--rscript',default='Rscript')
    ap.add_argument('--graphics-rscript')
    ap.add_argument('--check-only',action='store_true')
    ap.add_argument('--skip-hallmark',action='store_true')
    ap.add_argument('--skip-manifest-check',action='store_true',help='Only for an already rerun working copy whose expected outputs have changed.')
    args=ap.parse_args()
    if not args.skip_manifest_check:manifest_check()
    if args.check_only:return
    steps=[
      ('patient_sensitivities',[sys.executable,'05_Analyses/summarize_patient_sensitivity_tables_only.py']),
      ('assay_sensitivities',[sys.executable,'05_Analyses/summarize_assay_sensitivity.py']),
      ('independent_assay_check',[args.rscript,'--vanilla','05_Analyses/verify_assay_statistics.R']),
      ('fixed_score_evaluation',[args.rscript,'--vanilla','05_Analyses/external_validation/evaluate_source_cohorts_v04.R']),
      ('clinical_comparison',[args.rscript,'--vanilla','05_Analyses/external_validation/compare_clinical_models_v04.R']),
      ('clinical_coefficients',[args.rscript,'--vanilla','05_Analyses/external_validation/export_clinical_coefficients_v04.R']),
      ('descriptive_time_auc',[args.rscript,'--vanilla','05_Analyses/external_validation/evaluate_time_auc_v04.R']),
      ('independent_score_check',[sys.executable,'05_Analyses/external_validation/verify_source_results_v04.py']),
      ('independent_clinical_metric_check',[sys.executable,'05_Analyses/external_validation/verify_clinical_metrics_v04.py']),
      ('measured_glioma_replay',[sys.executable,'05_Analyses/glioma_drug_audit/replay_prism_from_subsets_v05.py']),
      ('measured_assay_quality',[sys.executable,'05_Analyses/glioma_drug_audit/describe_assay_quality_v01.py']),
    ]
    if not args.skip_hallmark:
        steps.append(('hallmark_aggregate_replay',[args.rscript,'--vanilla','05_Analyses/functional_enrichment/rerun_hallmark_from_aggregates_v04.R']))
    if args.graphics_rscript:
        steps.append(('scientific_figures',[args.graphics_rscript,'--vanilla','06_Figures/draw_revision_figures_v04.R']))
        steps.append(('measured_drug_figures',[args.graphics_rscript,'--vanilla','06_Figures/draw_prism_figure_v05.R']))
        steps.append(('figure_metadata_v04',[args.graphics_rscript,'--vanilla','06_Figures/audit_figure_exports_v04.R']))
        steps.append(('figure_metadata_v05',[args.graphics_rscript,'--vanilla','06_Figures/audit_figure_exports_v05.R']))
    steps.append(('compare_retained_reference_tables',[sys.executable,'VERIFY_REPLAY.py']))
    folder=ROOT/'RUN_LOGS';folder.mkdir(exist_ok=True)
    receipts=[]
    env=os.environ.copy();env['PYTHONIOENCODING']='utf-8';env['PYTHONDONTWRITEBYTECODE']='1'
    for label,command in steps:
        print('RUN',label,flush=True)
        start=datetime.now(timezone.utc).isoformat()
        with (folder/(label+'.log')).open('w',encoding='utf-8') as log:
            p=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env=env)
        receipts.append({'step':label,'command':command,'started_utc':start,'finished_utc':datetime.now(timezone.utc).isoformat(),'exit_code':p.returncode})
        (folder/'run_receipt.json').write_text(json.dumps(receipts,indent=2),encoding='utf-8')
        if p.returncode:raise RuntimeError(label+' failed; inspect its RUN_LOGS log')
        print('PASS',label,flush=True)
    print('Revision replay complete. Numerical reproduction is not independent clinical or biological validation.',flush=True)

if __name__=='__main__':main()
