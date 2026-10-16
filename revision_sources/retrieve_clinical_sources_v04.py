from audit_external_sources_v04 import fetch,OUT
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import json,urllib.request

sources={
 'GSE16011_supplement_tables_record.json':'https://api.figshare.com/v2/articles/22383048',
 'GSE16011_supplement_methods_record.json':'https://api.figshare.com/v2/articles/22383051',
 'GSE108474_clinical.data.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108474/suppl/GSE108474_REMBRANDT_clinical.data.txt.gz',
 'GSE108474_biospecimen_mapping_GEO.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108474/suppl/GSE108474_REMBRANDT_biospecimen_mapping_GEO.txt.gz',
 'GSE108474_data.dictionary.xlsx':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108474/suppl/GSE108474_REMBRANDT_data.dictionary.xlsx',
 'GSE108474_GeneExpression.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108474/suppl/GSE108474_REMBRANDT_GeneExpression.txt.gz',
 'GSE108474_pubmed.xml':'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=30106394&retmode=xml',
}
with ThreadPoolExecutor(max_workers=3) as pool: records=list(pool.map(fetch,sources.items()))
for n in ['GSE16011_supplement_tables_record.json','GSE16011_supplement_methods_record.json']:
 p=OUT/n
 if not p.exists():continue
 d=json.loads(p.read_text(encoding='utf-8'))
 if d.get('resource_doi','').lower()!='10.1158/0008-5472.can-09-2307':
  raise ValueError('Figshare primary article DOI mismatch: '+str(d.get('resource_doi')))
 for f in d.get('files',[]):
  if f['name'].endswith('.pdf'):records.append(fetch(('GSE16011_'+f['name'],f['download_url'])))
(OUT/'retrieval_ledger_clinical.json').write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),'records':records},ensure_ascii=False,indent=2),encoding='utf-8')
for r in records: print(json.dumps({k:r[k] for k in ['name','usable','bytes','error'] if k in r},ensure_ascii=False))
