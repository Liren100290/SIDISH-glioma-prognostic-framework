"""Verify the exact original/revision source mirrors and syntax of Python sources."""
from pathlib import Path
import ast
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parent

def main():
    with (ROOT / 'SOURCE_MAP.csv').open(encoding='utf-8', newline='') as stream:
        records = list(csv.DictReader(stream))
    assert len(records) == 49, 'Expected 49 preserved source files'
    for record in records:
        path = ROOT / record['repository_path']
        assert path.resolve().is_relative_to(ROOT.resolve())
        payload = path.read_bytes()
        assert hashlib.sha256(payload).hexdigest() == record['sha256'], path.name
        if path.suffix.lower() == '.py':
            ast.parse(payload.decode('utf-8-sig'), filename=str(path))
    model = json.loads((ROOT / 'locked_model.json').read_text(encoding='utf-8'))
    assert len(model['genes']) == 9
    assert model['cutoff'] == -0.0269027400037899
    print('Verified 49 exact source mirrors, Python syntax, nine-gene model and fixed cutoff.')

if __name__ == '__main__':
    main()
