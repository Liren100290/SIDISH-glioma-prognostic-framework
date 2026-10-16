from pathlib import Path
from datetime import datetime,timezone
import json
p=Path(__file__).resolve().parent/'additional_candidate_schema_amendment_v04.json'
if p.exists():raise SystemExit('Preserve existing amendment')
d={
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'stage':'After downloading and inspecting primary-source column definitions and metadata, before calculating any REMBRANDT risk score or score-survival result',
 'age':'Only five-year age ranges are supplied. Include age bands whose lower limit is >=18 (in practice >=20). Exclude 15-19 as adult status unresolved, plus younger/missing ages; use band midpoint for exploratory covariate adjustment and label it approximate.',
 'event':'The primary clinical field is EVENT_OS, a formal event indicator. The primary dictionary defines its event as death and overall survival as time to death or last follow-up. Interpret 1=death and 0=no death recorded; retain the source definitions and flag any contradictory source evidence rather than recoding by model performance.',
 'histology':'Require WHO II/III/IV and DISEASE_TYPE in GBM, ASTROCYTOMA, OLIGODENDROGLIOMA or MIXED. Historical categories do not establish modern integrated molecular diagnosis.',
 'specimen_timing':'TIMEPOINT is blank for all 671 records. Thus the available adult diffuse-glioma evaluation is sampling-timing uncertain. It must not be described as confirmed primary pretreatment validation. Prior therapy fields contain recorded exposure or blanks, with blanks treated as unknown.',
 'unique_patients':'The official biospecimen map contains 541 distinct subjects and 541 distinct specimens. Therefore no specimen averaging is required in these downloaded inputs.',
 'clinical_models':'Exploratory score association adjusted for age-band midpoint, sex, grade where complete; KPS-complete adjustment separately. Do not fit ambiguous blank treatment records as untreated controls.',
 'subgroups':'Grade II/III/IV checks retain the full eligible reference; report all groups with >=20 patients and >=10 events, otherwise descriptive counts only.',
 'reporting':'All eligible results retained regardless of direction. These protocol clarifications precede scoring; they do not resolve source recruitment overlap or establish clinical utility.'
}
p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
print(p)
