"""Retrieve primary records for revision-stage validation; never alter inputs/results.

Network failures are evidence of retrieval failure, not proof of source facts.
The fetched originals and a hash/URL ledger are retained for inspection.
"""
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import hashlib, json, urllib.request, urllib.error, re, sys
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
OUT = HERE / 'primary_sources_v04'
OUT.mkdir(exist_ok=True)
SOURCES = {
 'GSE43378_article.xml': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7657127/fullTextXML',
 'GSE43378_series.soft.txt': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43378&targ=self&form=text&view=full',
 'GSE7696_series.soft.txt': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7696&targ=self&form=text&view=full',
 'GSE16011_series.soft.txt': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16011&targ=self&form=text&view=full',
 'GSE7696_crossref.json': 'https://api.crossref.org/works/10.1200/JCO.2007.15.7164',
 'GSE7696_biostudies_info.json': 'https://www.ebi.ac.uk/biostudies/api/v1/studies/E-GEOD-7696/info',
 'GSE7696_biostudies_study.json': 'https://www.ebi.ac.uk/biostudies/api/v1/studies/E-GEOD-7696',
}
FOLLOWUP = {
 'GSE43378_article_ncbi.xml': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=7657127',
 'GSE43378_article_bioc.json': 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/PMC7657127/unicode',
 'GSE7696_publisher.pdf': 'https://ascopubs.org/doi/pdfdirect/10.1200/JCO.2007.15.7164',
 'E-GEOD-7696.idf.txt': 'https://ftp.ebi.ac.uk/biostudies/fire/E-GEOD-/696/E-GEOD-7696/Files/E-GEOD-7696.idf.txt',
 'E-GEOD-7696.sdrf.txt': 'https://ftp.ebi.ac.uk/biostudies/fire/E-GEOD-/696/E-GEOD-7696/Files/E-GEOD-7696.sdrf.txt',
 'GSE55918_series.soft.txt': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55918&targ=self&form=text&view=full',
 'GSE16011_series_matrix.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16011/matrix/GSE16011_series_matrix.txt.gz',
 'GSE16011_pubmed.xml': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=19920198&retmode=xml',
}

def fetch(item):
    name, url = item
    rec = {'name': name, 'url': url, 'retrieved_utc': datetime.now(timezone.utc).isoformat()}
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'AcademicRevisionAudit/1.0', 'Accept': '*/*'})
        with urllib.request.urlopen(req, timeout=45) as r:
            blob = r.read()
            rec.update(http_status=r.status, final_url=r.url, content_type=r.headers.get('Content-Type', ''))
        # A challenge page is not a usable primary record.
        challenge = any(x in blob[:6000].lower() for x in [b'recaptcha', b'checking your browser', b'captcha-delivery'])
        rec.update(bytes=len(blob), sha256=hashlib.sha256(blob).hexdigest(), challenge_page=challenge)
        dest = OUT / name
        if dest.exists() and dest.read_bytes() != blob:
            dest = OUT / (Path(name).stem + '_retrieved_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + Path(name).suffix)
        if not dest.exists():
            dest.write_bytes(blob)
        rec['saved_as'] = str(dest)
        rec['usable'] = not challenge
        if blob.lstrip().startswith(b'[Error]'):
            rec.update(usable=False,source_error_response=True)
        if name.endswith('.pdf') and not blob.startswith(b'%PDF'):
            rec.update(usable=False,unexpected_pdf_content=True)
        if name.endswith('.json') and rec['usable']:
            try: json.loads(blob)
            except Exception as ex: rec.update(usable=False,invalid_json=str(ex))
        if name.endswith('.xml') and not challenge:
            root = ET.fromstring(blob)
            rec['fulltext_body_present'] = root.find('.//body') is not None
            txt = '\n'.join(' '.join(e.itertext()) for e in root.iter() if e.tag in ('p', 'title', 'table-wrap', 'supplementary-material'))
            (OUT / (Path(name).stem + '_text.txt')).write_text(txt, encoding='utf-8')
            rec['evidence_snippets'] = [s[:1600] for s in txt.splitlines() if re.search(r'Niigata|untreated|32 |surviv|supplement|Table S|primary|patient', s, re.I)][:28]
        elif name.endswith('.soft.txt') and not challenge:
            txt = blob.decode('utf-8', errors='replace')
            rec['evidence_snippets'] = [s for s in txt.splitlines() if re.search(r'summary|overall_design|pubmed|supplementary_file|relation|sample_id', s, re.I)][:20]
        elif name == 'GSE7696_crossref.json':
            msg = json.loads(blob)['message']
            rec['publisher_links'] = msg.get('link', [])
    except Exception as ex:
        rec.update(usable=False, error=type(ex).__name__ + ': ' + str(ex))
    return rec

if __name__ == '__main__':
    followup = '--followup' in sys.argv
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(fetch, (FOLLOWUP if followup else SOURCES).items()))
    ledger = {'created_utc': datetime.now(timezone.utc).isoformat(), 'purpose': 'Primary-source audit; no event recoding or score modification performed', 'records': results}
    (OUT / ('retrieval_ledger_followup.json' if followup else 'retrieval_ledger.json')).write_text(json.dumps(ledger, ensure_ascii=False, indent=2), encoding='utf-8')
    for r in results:
        print(json.dumps({k: r[k] for k in ['name', 'usable', 'bytes', 'error', 'publisher_links'] if k in r}, ensure_ascii=False))
