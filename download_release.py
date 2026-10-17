"""Download/verify seven ZIPs and the README from v1.0.1; safely extract ZIPs only."""
from pathlib import Path
from urllib.request import Request, urlopen
from zipfile import ZipFile
import argparse
import hashlib

ROOT = Path(__file__).resolve().parent
REPOSITORY = 'Liren100290/SIDISH-glioma-prognostic-framework'
TAG = 'v1.0.1'
# SHA-256 of the release manifest using LF line endings.
MANIFEST_SHA256 = '1fc1d4024059e2e02bb86100a18498a65123ec8ae9b54cb282d227fbf85fd181'

def sha(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()

def expected_assets():
    manifest = (ROOT / 'SHA256SUMS.txt').read_bytes().replace(b'\r\n', b'\n')
    if hashlib.sha256(manifest).hexdigest() != MANIFEST_SHA256:
        raise RuntimeError('Checksum manifest does not match release ' + TAG + '; use the matching downloader and manifest.')
    entries = {}
    for line in manifest.decode('utf-8').splitlines():
        checksum, name = line.split(None, 1)
        name = name.strip()
        if (len(checksum) != 64 or any(c not in '0123456789abcdef' for c in checksum)
                or Path(name).name != name or name in entries):
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
    parser.add_argument('--verify-only', action='store_true', help='Do not use the network; require all seven ZIPs and the release README locally.')
    parser.add_argument('--extract', action='store_true', help='Extract verified archives into the output directory; never replace differing files.')
    options = parser.parse_args()
    destination = options.output.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    verified_archives = []
    for name, expected in expected_assets().items():
        path = destination / name
        if not path.exists():
            if options.verify_only:
                raise FileNotFoundError(path)
            url = f'https://github.com/{REPOSITORY}/releases/download/{TAG}/{name}'
            request = Request(url, headers={'User-Agent': 'SIDISH-glioma-reproducibility/1.0.1'})
            with urlopen(request, timeout=60) as response, path.open('xb') as stream:
                for block in iter(lambda: response.read(1024 * 1024), b''):
                    stream.write(block)
        if sha(path) != expected:
            raise RuntimeError('SHA-256 mismatch: ' + name + '; file retained for inspection, nothing deleted.')
        print('SHA-256 verified:', name)
        if path.suffix.lower() == '.zip':
            verified_archives.append(path)
    if options.extract:
        for archive in verified_archives:
            extract_safely(archive, destination)
    print('All release assets verified. No analytical code was automatically executed.')

if __name__ == '__main__':
    main()
