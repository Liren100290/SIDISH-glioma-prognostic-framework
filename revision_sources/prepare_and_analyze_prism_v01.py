"""Exploratory measured glioma-screen evaluation of all frozen drug candidates.

Scores and subgroup references retain the submitted nine weights. PRISM screening
is separate from the reported GDSC2 prediction resource; cell models can overlap.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE.parent / '_pydeps'))
import numpy as np
import pandas as pd
from scipy.stats import rankdata

SRC = HERE / 'primary_sources'
OUT = HERE / 'results_v01'
OUT.mkdir(exist_ok=True)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def bh(values):
    p = np.asarray(values, dtype=float)
    q = np.full(len(p), np.nan)
    good = np.flatnonzero(np.isfinite(p))
    if not len(good):
        return q
    order = good[np.argsort(p[good])]
    adjusted = np.minimum.accumulate((p[order] * len(order) / np.arange(1, len(order) + 1))[::-1])[::-1]
    q[order] = np.minimum(adjusted, 1)
    return q


def association(x, y, seed):
    n = len(x)
    assert n >= 10 and np.isfinite(x).all() and np.isfinite(y).all()
    rx = rankdata(x)
    ry = rankdata(y)
    rx -= rx.mean()
    ry -= ry.mean()
    den = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    observed = float(rx @ ry / den)
    rng = np.random.default_rng(seed)
    # Random permutations are two-sided; the score direction is never flipped.
    perms = np.asarray([rng.permutation(ry) for _ in range(10000)])
    simulated = (perms @ rx) / den
    p = float((1 + np.count_nonzero(np.abs(simulated) >= abs(observed) - 1e-14)) / 10001)
    ix = rng.integers(0, n, size=(1000, n))
    bx = rankdata(x[ix], axis=1)
    by = rankdata(y[ix], axis=1)
    bx -= bx.mean(axis=1, keepdims=True)
    by -= by.mean(axis=1, keepdims=True)
    bd = np.sqrt((bx * bx).sum(axis=1) * (by * by).sum(axis=1))
    boots = np.divide((bx * by).sum(axis=1), bd, out=np.full(1000, np.nan), where=bd > 0)
    valid = boots[np.isfinite(boots)]
    lo, hi = np.quantile(valid, [.025, .975])
    return observed, float(lo), float(hi), p, len(valid)


def main():
    for name in ['analysis_protocol_v01.json', 'protocol_amendment_candidate_table_v01.json',
                 'schema_and_identity_amendment_v01.json']:
        assert (HERE / name).exists(), name
    for stage in ['metadata_files', 'large_data']:
        ledger = json.loads((SRC / ('retrieval_ledger_' + stage + '.json')).read_text(encoding='utf-8'))
        assert all(r['publisher_size_and_hash_verified'] for r in ledger['records'])
    mapping = pd.read_csv(HERE / 'frozen_candidate_compound_screen_mapping.csv')
    original_candidates = json.loads((HERE / 'protocol_amendment_candidate_table_v01.json').read_text(encoding='utf-8'))['correct_original_candidates']
    assert len(original_candidates) == 16
    weights = pd.read_csv(ROOT / '05_Analyses/external_validation/locked_nine_gene_coefficients.csv')
    genes = weights.gene.tolist()
    beta = weights.coefficient.to_numpy(float)
    assert len(genes) == 9
    meta = pd.read_csv(HERE / 'model_metadata_eligibility_before_expression.csv')
    candidates = meta[meta.candidate_before_expression].copy()
    assert len(candidates) == 40 and candidates.depmap_id.is_unique

    header = pd.read_csv(SRC / 'CCLE_19Q4_expression.csv', nrows=0).columns.tolist()
    gene_columns = {}
    for gene in genes:
        found = [s for s in header if re.sub(r'\s*\([^)]*\)\s*$', '', s) == gene]
        assert len(found) == 1, (gene, found)
        gene_columns[gene] = found[0]
    expr = pd.read_csv(SRC / 'CCLE_19Q4_expression.csv', usecols=[header[0]] + list(gene_columns.values()))
    expr = expr.rename(columns={header[0]: 'depmap_id', **{v:k for k,v in gene_columns.items()}})
    assert expr.depmap_id.is_unique
    parent = candidates.merge(expr, on='depmap_id', how='left', validate='one_to_one', indicator='expression_merge')
    parent['all9_finite'] = np.isfinite(parent[genes].to_numpy(float)).all(axis=1)
    parent['expression_eligible'] = parent.expression_merge.eq('both') & parent.all9_finite
    parent.to_csv(OUT / 'glioma_expression_eligibility.csv', index=False)
    d = parent[parent.expression_eligible].copy().sort_values('depmap_id').reset_index(drop=True)
    assert len(d) >= 10 and d.depmap_id.is_unique
    x = d[genes].to_numpy(float)
    mean = x.mean(axis=0)
    sd = x.std(axis=0, ddof=1)
    assert np.isfinite(sd).all() and (sd > 0).all(), dict(zip(genes, sd))
    # The source is already log2(TPM+1); no second logarithm is applied.
    z = (x - mean) / sd
    d['score'] = z @ beta
    d['score_per_parent_SD'] = (d.score - d.score.mean()) / d.score.std(ddof=1)
    stripped = d.stripped_cell_line_name.str.replace(r'[^A-Za-z0-9]', '', regex=True).str.upper()
    d['documented_BRAF_model_excluded_in_sensitivity'] = stripped.isin(['AM38', 'DBTRG05MG'])
    d['documented_GBM'] = d.lineage_sub_subtype.eq('glioblastoma')
    d['documented_adult'] = pd.to_numeric(d.age, errors='coerce').ge(18)
    d.to_csv(OUT / 'glioma_parent_expression_and_scores.csv', index=False)
    pd.DataFrame({'gene': genes, 'source_column': [gene_columns[g] for g in genes],
                  'parent_mean_log2_TPM1': mean, 'parent_sample_SD': sd,
                  'locked_coefficient': beta}).to_csv(OUT / 'frozen_expression_reference_and_weights.csv', index=False)

    selected = []
    total = 0
    source_cols = pd.read_csv(SRC / 'PRISM_19Q4_secondary_dose_response.csv', nrows=0).columns.tolist()
    assert 'convergence' not in source_cols, 'Inspect the documented convergence flag before any analysis.'
    columns = ['broad_id','depmap_id','ccle_name','screen_id','auc','r2','name','passed_str_profiling','row_name']
    keys = set(zip(mapping.broad_id, mapping.screen_id))
    all_cancer_candidate_rows = 0
    for chunk in pd.read_csv(SRC / 'PRISM_19Q4_secondary_dose_response.csv', usecols=columns, chunksize=50000):
        total += len(chunk)
        wanted = pd.Series(list(zip(chunk.broad_id, chunk.screen_id)), index=chunk.index).isin(keys)
        all_cancer_candidate_rows += int(wanted.sum())
        chosen = chunk[wanted & chunk.depmap_id.isin(d.depmap_id)].copy()
        selected.append(chosen)
    raw = pd.concat(selected, ignore_index=True)
    raw = raw.merge(mapping, on=['broad_id','screen_id'], how='left', validate='many_to_one')
    model_names = dict(zip(d.depmap_id, d.ccle_name))
    raw['metadata_name_consistent'] = raw.ccle_name.eq(raw.depmap_id.map(model_names))
    raw['response_identity_QC_passed'] = raw.passed_str_profiling.astype(str).str.lower().eq('true')
    raw['auc_finite'] = np.isfinite(pd.to_numeric(raw.auc, errors='coerce'))
    raw['eligible_response'] = raw.metadata_name_consistent & raw.response_identity_QC_passed & raw.auc_finite
    raw.to_csv(OUT / 'all_candidate_glioma_response_rows.csv', index=False)
    eligible = raw[raw.eligible_response].copy()
    keycols = ['original_candidate','broad_id','screen_id','depmap_id']
    response = eligible.groupby(keycols, as_index=False).agg(
        auc=('auc','median'), source_rows=('auc','size'),
        source_auc_min=('auc','min'), source_auc_max=('auc','max'),
        source_r2_median=('r2','median'))
    response = response.merge(d[['depmap_id','score','score_per_parent_SD','documented_GBM',
                     'documented_adult','documented_BRAF_model_excluded_in_sensitivity']],
                     on='depmap_id', how='left', validate='many_to_one')
    assert response.groupby(['broad_id','screen_id','depmap_id']).size().eq(1).all()
    response.to_csv(OUT / 'model_level_observed_responses_and_scores.csv', index=False)

    contexts = {
        'all_glioma': set(d.depmap_id),
        'documented_GBM': set(d.loc[d.documented_GBM,'depmap_id']),
        'documented_adult': set(d.loc[d.documented_adult,'depmap_id']),
        'excluding_two_documented_BRAF_models': set(d.loc[~d.documented_BRAF_model_excluded_in_sensitivity,'depmap_id']),
    }
    tests = []
    for context, ids in contexts.items():
        for i, row in mapping.iterrows():
            sub = response[response.broad_id.eq(row.broad_id) & response.screen_id.eq(row.screen_id)
                           & response.depmap_id.isin(ids)].sort_values('depmap_id')
            n = len(sub)
            eligible_test = n >= 10 and sub.score.nunique() > 1 and sub.auc.nunique() > 1
            rho = lo = hi = p = np.nan
            boot = 0
            if eligible_test:
                seed = 20261015 + len(tests)
                rho,lo,hi,p,boot = association(sub.score.to_numpy(float), sub.auc.to_numpy(float),seed)
            tests.append({'context':context,'original_candidate':row.original_candidate,
                'broad_id':row.broad_id,'screen_id':row.screen_id,'n_models':n,
                'parent_reference_models':len(d),'primary_selected_screen':bool(row.primary_selected),
                'inferentially_eligible':eligible_test,'Spearman_rho':rho,'bootstrap_lower95':lo,
                'bootstrap_upper95':hi,'permutation_p_two_sided':p,'bootstrap_valid':boot,
                'median_AUC':sub.auc.median() if n else np.nan,
                'min_AUC':sub.auc.min() if n else np.nan,'max_AUC':sub.auc.max() if n else np.nan,
                'source_r2_median':sub.source_r2_median.median() if n else np.nan,
                'hypothesized_direction':'negative rho: higher original score, lower measured viability AUC',
                'analysis_status':'Exploratory public pharmacologic screen; selected compounds and overlapping cell models prevent independent clinical-validation claims'})
    tests = pd.DataFrame(tests)
    tests['BH_q_sensitivity_family'] = bh(tests.permutation_p_two_sided)
    primary = tests[tests.context.eq('all_glioma') & tests.primary_selected_screen].copy()
    assert len(primary) == 16 and primary.original_candidate.nunique() == 16
    primary['BH_q_primary_16'] = bh(primary.permutation_p_two_sided)
    primary['original_order'] = primary.original_candidate.map({n:i for i,n in enumerate(original_candidates)})
    primary = primary.sort_values('original_order').drop(columns='original_order')
    tests.to_csv(OUT / 'all_context_screen_tests.csv', index=False)
    primary.to_csv(OUT / 'primary_16_candidate_results.csv', index=False)
    pd.DataFrame([{'context':c,'models_before_drug_specific_missingness':len(ids)} for c,ids in contexts.items()]).to_csv(
        OUT / 'predeclared_context_counts.csv', index=False)
    audit = {
        'completed_utc':datetime.now(timezone.utc).isoformat(),
        'source_response_rows':total,'all_cancer_frozen_candidate_rows':all_cancer_candidate_rows,
        'metadata_glioma_candidates':len(candidates),'all9_expression_eligible_parent':len(d),
        'missing_expression_models':parent.loc[~parent.expression_eligible,'depmap_id'].tolist(),
        'raw_glioma_candidate_rows':len(raw),'identity_name_disagreements':int((~raw.metadata_name_consistent).sum()),
        'failed_response_STR_rows':int((~raw.response_identity_QC_passed).sum()),
        'nonfinite_AUC_rows':int((~raw.auc_finite).sum()),'model_level_responses':len(response),
        'technical_duplicate_rows_collapsed':len(eligible)-len(response),
        'primary_candidates':16,'all_prespecified_context_screen_tests':len(tests),
        'primary_models_per_drug_min':int(primary.n_models.min()),'primary_models_per_drug_max':int(primary.n_models.max()),
        'primary_BH_significant_count':int(primary.BH_q_primary_16.lt(.05).sum()),
        'primary_negative_direction_count':int(primary.Spearman_rho.lt(0).sum()),
        'nine_weights_refit':False,'clinical_cutoff_optimized_or_reused_as_drug_threshold':False,
        'subgroups_restandardized':False,'source_AUC_clipped':False,'outcome_r2_filter_applied':False,
        'expression_representation':'Official CCLE19Q4 log2(TPM+1), already transformed; one exact source column per original gene.',
        'parent_reference_bootstrap_interpretation':'Bootstrap intervals conditional on observed full-parent expression reference; resampling unique model observations, not drug concentrations or wells.',
        'BRAF_sensitivity_interpretation':'Excludes two model identities documented as BRAF V600E in DOI10.1186/s40478-019-0774-7; other models are not inferred to be BRAF-wildtype.',
        'independence_limit':'PRISM measures response in a separate screen from reported GDSC2 predictions. The exact original GDSC2 training membership is not source-matched, and model identities can overlap. Candidates were selected before this post-selection evaluation. No causal, BBB, prospective or clinical validation claim.',
        'input_sha256':{n:sha(SRC/n) for n in ['PRISM_19Q4_secondary_dose_response.csv','CCLE_19Q4_expression.csv']},
        'script_sha256':sha(Path(__file__)),
    }
    (OUT / 'PRISM_glioma_audit_v01.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(audit,ensure_ascii=False,indent=2),flush=True)
    print(primary[['original_candidate','screen_id','n_models','Spearman_rho','bootstrap_lower95','bootstrap_upper95','permutation_p_two_sided','BH_q_primary_16']].to_string(index=False),flush=True)


if __name__ == '__main__':
    main()
