from pathlib import Path
from datetime import datetime,timezone
import json
p=Path(__file__).resolve().parent/'source_scoring_plan_v04.json'
if p.exists():raise SystemExit('Preserve existing scoring plan')
d={
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'stage':'After input/eligibility audits, before any new GSE16011 reconstruction or REMBRANDT score-survival calculation',
 'parents':['GSE16011 revision-source reconstruction (254 adults,230 events)','GSE108474 timing-uncertain adult grade-II-IV evaluation (277 adults,235 events)','GSE43378 all-50 timing-uncertain case-mix context plus source-supported primary-GBM subset'],
 'scoring':'Original nine weights and cutoff; sample gene SD/mean from the eligible full parent cohort, frozen in every subgroup and complete-case analysis. Continuous-score HR unit is one full-parent score SD.',
 'metrics':'Harrell concordance, 1000 patient bootstrap percentile intervals seed20261015; fixed-cutoff HR/logrank; continuous score associations with CI. Higher original score always indicates predicted higher risk.',
 'subgroups':'Every grade II/III/IV and measured IDH1 group with >=20 patients and >=10 events; smaller groups retained in attrition counts. No subgroup rescaling.',
 'covariate_contexts':['Unadjusted','Age per10years (or five-year-band midpoint),sex,grade when variable','Add KPS per10points where complete','GSE16011 molecular context additionally includes source IDH1(R132),1p/19q codeletion','GSE16011 treatment-record context additionally includes KPS,radiotherapy,any chemotherapy and complete-vs-other resection'],
 'fitting':'Complete cases separately identified per context. Remove only invariant covariate terms; record those omissions. No inferential claim from sparse fits (require >=30 observations and >=20 events). Cox proportional hazards checked for score, exploratory.',
 'clinical_comparison':'Within-cohort clinical models are newly fitted and explicitly labeled. If conducted, apparent paired-bootstrap increments and cross-validation are distinct from external validation of a development-trained combined model.',
 'retention':'Keep all directions/intervals/pvalues. Do not optimize cutoff, alter nine internal weights, select a favorable validation cohort, or flip outcome/score orientation.',
 'clinical_limits':'Treatment timing is incomplete and historical grade/IDH1 data cannot establish modern integrated diagnoses or treatment independence.'
}
p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8');print(p)
