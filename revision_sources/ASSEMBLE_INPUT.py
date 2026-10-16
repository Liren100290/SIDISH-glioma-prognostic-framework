"""Reassemble one fixed research input; verify every byte by SHA-256."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def ensure_input():
    spec=json.loads((ROOT/'PACKAGED_INPUT_PARTS.json').read_text(encoding='utf-8'))
    target=ROOT/spec['target']
    assert target.resolve().is_relative_to(ROOT.resolve())
    if target.exists():
        if target.stat().st_size!=spec['target_bytes'] or digest(target)!=spec['target_sha256']:
            raise RuntimeError('Existing assembled input differs from the documented SHA-256; no file was overwritten.')
        return target
    for rec in spec['parts']:
        p=ROOT/rec['path']
        assert p.resolve().is_relative_to(ROOT.resolve())
        if not p.is_file():
            raise RuntimeError('Missing input part: '+rec['path']+'. Extract all four Supplementary Code 1 data ZIPs into the same parent as the core ZIP first.')
        if p.stat().st_size!=rec['bytes'] or digest(p)!=rec['sha256']:
            raise RuntimeError('Input part checksum mismatch: '+rec['path'])
    target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('xb') as out:
        for rec in spec['parts']:
            with (ROOT/rec['path']).open('rb') as inp:
                for block in iter(lambda:inp.read(1024*1024),b''):out.write(block)
    if target.stat().st_size!=spec['target_bytes'] or digest(target)!=spec['target_sha256']:
        raise RuntimeError('Assembled input checksum mismatch.')
    print('Assembled whole-transcriptome input; SHA-256 verified.',flush=True)
    return target
if __name__=='__main__':ensure_input()
