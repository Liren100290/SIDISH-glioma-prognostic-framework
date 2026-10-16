"""Download/verify the six fixed release assets; optionally extract without overwriting."""
from pathlib import Path
from urllib.request import Request, urlopen
from zipfile import ZipFile
import argparse
import hashlib

ROOT = Path(__file__).resolve().parent
REPOSITORY = 'Liren100290/SIDISH-glioma-prognostic-framework'
TAG = 'v1.0.0'

def sha(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()

def expected_assets():
    entries = {}
    for line in (ROOT / 'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
        checksum, name = line.split(None, 1)
        name = name.strip()
        if len(checksum) != 64 or Path(name).name != name:
            raise ValueError('Invalid checksum entry')
        entries[name] = checksum
    return entries

def extract_safely(archive, destination):
    with ZipFile(archive) as z:
        if z.testzip() is not None:
            raise RuntimeError('ZIP CRC mismatch: ' + archive.name)
        for info in z.infolist():
            if info.is_dir():
                continue
            target = destination / info.filename
            if not target.resolve().is_relative_to(destination.resolve()):
                raise RuntimeError('Unsafe archive member')
            data = z.read(info)
            if target.exists():
                if hashlib.sha256(data).hexdigest() != sha(target):
                    raise RuntimeError('Existing file differs; no overwrite: ' + str(target))
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open('xb') as stream:
                stream.write(data)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('download'))
    parser.add_argument('--verify-only', action='store_true', help='Do not use the network; require all six ZIPs locally.')
    parser.add_argument('--extract', action='store_true', help='Extract verified archives into the output directory; never replace differing files.')
    options = parser.parse_args()
    destination = options.output.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    for name, expected in expected_assets().items():
        path = destination / name
        if not path.exists():
            if options.verify_only:
                raise FileNotFoundError(path)
            url = f'https://github.com/{REPOSITORY}/releases/download/{TAG}/{name}'
            request = Request(url, headers={'User-Agent': 'SIDISH-glioma-reproducibility/1.0.0'})
            with urlopen(request, timeout=60) as response, path.open('xb') as stream:
                for block in iter(lambda: response.read(1024 * 1024), b''):
                    stream.write(block)
        if sha(path) != expected:
            raise RuntimeError('SHA-256 mismatch: ' + name + '; file retained for inspection, nothing deleted.')
        print('SHA-256 verified:', name)
        if options.extract:
            extract_safely(path, destination)
    print('All release assets verified. No analytical code was automatically executed.')

if __name__ == '__main__':
    main()
