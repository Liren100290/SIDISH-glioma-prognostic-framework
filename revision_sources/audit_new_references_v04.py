from audit_external_sources_v04 import fetch,OUT
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from pathlib import Path
import json,urllib.parse

dois=['10.1111/cas.12214','10.1038/sdata.2018.158','10.1093/nar/gks461','10.1016/j.cels.2015.12.004','10.1002/bimj.200610301','10.1186/1472-6947-8-53']
sources={f'new_reference_{i+40}_crossref.json':'https://api.crossref.org/works/'+urllib.parse.quote(doi,safe='') for i,doi in enumerate(dois)}
sources['GSE108474_article_ncbi.xml']='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=6091243'
sources={n:u for n,u in sources.items() if not (OUT/n).exists()}
with ThreadPoolExecutor(max_workers=3) as pool:records=list(pool.map(fetch,sources.items()))
refs=[]
for i,doi in enumerate(dois):
 p=OUT/f'new_reference_{i+40}_crossref.json'
 if not p.exists():continue
 d=json.loads(p.read_text(encoding='utf-8'))['message']
 assert d['DOI'].lower()==doi.lower()
 authors=[{'family':a.get('family',''),'given':a.get('given','')} for a in d.get('author',[])]
 refs.append({'number':i+40,'doi':doi,'title':d['title'][0],'authors':authors,'journal':d.get('container-title',[''])[0],
              'year':d.get('published',{}).get('date-parts',[[None]])[0][0],'volume':d.get('volume',''),'issue':d.get('issue',''),
              'pages':d.get('page',''),'article_number':d.get('article-number',''),'source':str(p)})
here=Path(__file__).resolve().parent
(here/'New_references_verified_v04.json').write_text(json.dumps(refs,ensure_ascii=False,indent=2),encoding='utf-8')
lines=[]
for r in refs:
 lines+=['TY  - JOUR','TI  - '+r['title']]
 lines+=['AU  - '+a['family']+', '+a['given'] for a in r['authors']]
 lines+=['JO  - '+r['journal'],'PY  - '+str(r['year']),'VL  - '+r['volume'],'IS  - '+r['issue'],'SP  - '+(r['pages'] or r['article_number']),'DO  - '+r['doi'],'UR  - https://doi.org/'+r['doi'],'ER  - ','']
(here/'New_references_verified_v04.ris').write_text('\n'.join(lines),encoding='utf-8')
ledger_path=OUT/'retrieval_ledger_new_references.json'
old=json.loads(ledger_path.read_text(encoding='utf-8'))['records'] if ledger_path.exists() else []
byname={r['name']:r for r in old}
for r in records:byname[r['name']]=r
ledger_path.write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),'records':list(byname.values())},ensure_ascii=False,indent=2),encoding='utf-8')
for r in refs:print(r['number'],r['title'],r['year'],r['doi'])
