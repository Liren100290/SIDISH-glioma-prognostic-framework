from pathlib import Path
from datetime import datetime,timezone
import json
p=Path(__file__).resolve().parent/'clinical_comparison_protocol_v04.json'
if p.exists():raise SystemExit('Preserve existing protocol')
d={
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'stage':'After locked-score associations were calculated; before any new paired clinical-model performance, resampling, calibration, Brier or decision-curve calculation',
 'population':'GSE16011 revision-source reconstruction,254 patients; contexts age/sex/grade, plus KPS, molecular variables, and recorded-treatment context, each with explicit complete cases',
 'models':'Clinical Cox model versus the same terms plus the fixed nine-gene score as one predictor; new within-cohort fitting, not external validation of a development-trained combined model',
 'apparent_comparison':'Same complete cases for every pair; 1000 paired nonparametric patient bootstrap replicates with both models refitted, percentile CI for apparent C-index difference; no claim of optimism-corrected validation',
 'internal_validation':'10 repeats of 5-fold cross-validation, seed20261015+repeat, folds balanced within event/grade strata. All clinical coefficients, the score multiplier and baseline hazards estimated in training folds. Nine internal gene weights and full parent expression reference remain fixed.',
 'crossfit_concordance':'Compare out-of-fold predicted log cumulative hazard at36 months; report mean and range over the10 repeats as resampling stability, not a confidence interval. Average out-of-fold event probabilities across repeats for calibration and prediction error.',
 'horizons_months':[12,36,60],
 'censoring':'Kaplan-Meier estimate of censoring survival. IPCW Brier and net benefit retain the full evaluation denominator; deaths by horizon weighted by1/G(T-), survivors beyond horizon by1/G(t), earlier censored observations contribute zero weight. Skip horizons with censoring survival<0.1.',
 'calibration':'Five ordered predicted-risk groups; observed risk from Kaplan-Meier and Greenwood95%intervals, with group size and at-risk counts',
 'decision_analysis':'Hypothetical thresholds0.05-0.80 by0.05; no clinical action or optimal threshold is established, no treatment-effect claim',
 'limits':['New within-cohort models and cross-validation do not supply missing development-trained external clinical models','Complete-case selection, historical diagnosis, treatment timing and small molecular subsets limit interpretation','Report every context/time point and unfavorable differences; no analysis/threshold selection by performance']
}
p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8');print(p)
