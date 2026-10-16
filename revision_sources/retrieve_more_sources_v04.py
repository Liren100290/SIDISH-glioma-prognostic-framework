"""Additional public source retrieval; records access failures without changing outcomes."""
from audit_external_sources_v04 import fetch, OUT
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import json, urllib.request, hashlib

sources={
 'GSE16011_publisher.html':'https://aacrjournals.org/cancerres/article/69/23/9065/553005/Intrinsic-Gene-Expression-Profiles-of-Gliomas-Are',
 'GSE7696_gds_search.json':'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSE7696%5BACCN%5D&retmode=json',
 'GSE108474_series.soft.txt':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108474&targ=self&form=text&view=full',
 'GSE7696_paper_unige_search.json':'https://archive-ouverte.unige.ch/rest/search?query=10.1200%2FJCO.2007.15.7164',
 'GPL8542_platform.soft.txt':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8542&targ=self&form=text&view=full',
}
with ThreadPoolExecutor(max_workers=3) as pool:
 records=list(pool.map(fetch,sources.items()))
url='https://api.figshare.com/v2/articles/search'
payload=json.dumps({'resource_doi':'10.1158/0008-5472.CAN-09-2307','page_size':100}).encode()
rec={'name':'GSE16011_figshare_search.json','url':url,'method':'POST','request':payload.decode()}
try:
 req=urllib.request.Request(url,data=payload,headers={'Content-Type':'application/json','User-Agent':'AcademicRevisionAudit/1.0'})
 with urllib.request.urlopen(req,timeout=45) as r:blob=r.read()
 parsed=json.loads(blob)
 (OUT/rec['name']).write_bytes(blob)
 rec.update(usable=True,bytes=len(blob),sha256=hashlib.sha256(blob).hexdigest(),results=len(parsed))
except Exception as ex:rec.update(usable=False,error=type(ex).__name__+': '+str(ex))
records.append(rec)
(OUT/'retrieval_ledger_additional.json').write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),'records':records},ensure_ascii=False,indent=2),encoding='utf-8')
for r in records:print(json.dumps({k:r[k] for k in ['name','usable','bytes','error','results'] if k in r},ensure_ascii=False))
