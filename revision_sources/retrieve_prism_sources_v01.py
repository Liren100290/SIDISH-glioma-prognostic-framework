"""Download immutable public PRISM/CCLE source files and verify publisher hashes.

No score-response associations are calculated here. Originals are retained.
"""
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import argparse
import hashlib
import json
import urllib.request

HERE = Path(__file__).resolve().parent
OUT = HERE / 'primary_sources'
OUT.mkdir(exist_ok=True)


def checksum(p, method):
    h = hashlib.new(method)
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def download(job):
    target, record, source_doi = job
    p = OUT / target
    rec = {'file': target, 'original_file_name': record['name'], 'file_id': record['id'],
           'url': record['download_url'], 'source_doi': source_doi,
           'expected_bytes': record.get('size'), 'expected_md5': record.get('computed_md5'),
           'retrieval_started_utc': datetime.now(timezone.utc).isoformat()}
    try:
        if p.exists():
            assert p.stat().st_size == record['size']
            if record.get('computed_md5'):
                assert checksum(p, 'md5') == record['computed_md5']
            rec['cache_verified'] = True
        else:
            partial = p.with_name(p.name + '.download_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
            req = urllib.request.Request(record['download_url'], headers={'User-Agent': 'AcademicRevisionAudit/1.0'})
            with urllib.request.urlopen(req, timeout=45) as r, partial.open('xb') as f:
                rec.update(http_status=r.status, final_url=r.url, content_type=r.headers.get('Content-Type'))
                while True:
                    block = r.read(1024 * 1024)
                    if not block:
                        break
                    f.write(block)
            assert partial.stat().st_size == record['size'], (target, partial.stat().st_size, record['size'])
            if record.get('computed_md5'):
                assert checksum(partial, 'md5') == record['computed_md5'], target
            assert not p.exists()
            partial.rename(p)
        rec.update(bytes=p.stat().st_size, sha256=checksum(p, 'sha256'), md5=checksum(p, 'md5'),
                   publisher_size_and_hash_verified=True, completed_utc=datetime.now(timezone.utc).isoformat())
    except Exception as e:
        rec.update(publisher_size_and_hash_verified=False, error=type(e).__name__ + ': ' + str(e))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['metadata_files', 'large_data'])
    a = ap.parse_args()
    assert (HERE / 'analysis_protocol_v01.json').exists()
    assert (HERE / 'protocol_amendment_candidate_table_v01.json').exists()
    pr = json.loads((OUT / 'PRISM_19Q4_figshare_v4.json').read_text(encoding='utf-8'))
    ex = json.loads((OUT / 'DepMap_19Q4_figshare.json').read_text(encoding='utf-8'))
    assert pr['doi'] == '10.6084/m9.figshare.9393293.v4'
    assert ex['doi'] == '10.6084/m9.figshare.11384241.v3'
    specs = [
        ('PRISM_secondary_readme.txt', pr, 'secondary-screen-readme.txt'),
        ('PRISM_secondary_cell_line_info.csv', pr, 'secondary-screen-cell-line-info.csv'),
        ('PRISM_secondary_treatment_info.csv', pr, 'secondary-screen-replicate-collapsed-treatment-info.csv'),
        ('CCLE_19Q4_sample_info.csv', ex, 'sample_info.csv'),
        ('DepMap_19Q4_README.txt', ex, 'README'),
    ] if a.stage == 'metadata_files' else [
        ('PRISM_19Q4_secondary_dose_response.csv', pr, 'secondary-screen-dose-response-curve-parameters.csv'),
        ('CCLE_19Q4_expression.csv', ex, 'CCLE_expression.csv'),
    ]
    if a.stage == 'large_data':
        assert (HERE / 'schema_and_identity_amendment_v01.json').exists()
    jobs = []
    for target, source, name in specs:
        found = [f for f in source['files'] if f['name'] == name]
        assert len(found) == 1
        jobs.append((target, found[0], source['doi']))
    with ThreadPoolExecutor(max_workers=2) as pool:
        records = list(pool.map(download, jobs))
    (OUT / ('retrieval_ledger_' + a.stage + '.json')).write_text(json.dumps({
        'stage': a.stage, 'created_utc': datetime.now(timezone.utc).isoformat(),
        'records': records}, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps([{k: r[k] for k in ['file', 'bytes', 'publisher_size_and_hash_verified', 'error'] if k in r}
                      for r in records], ensure_ascii=False, indent=2), flush=True)
    if not all(r['publisher_size_and_hash_verified'] for r in records):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
