#!/usr/bin/env python3
"""Download public experiment files with resumable, verified HTTP ranges."""
import concurrent.futures as cf
import hashlib
import json
import os
import pathlib
import time
import threading
import urllib.error
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[2] / 'input/antkowiak2020'
FILES = [
    (22305708, 'I16_S2_R1_001.fastq.gz', 994447128, '6d48dee09e7cbe3f11bc749c9a7839c6'),
    (22305990, 'I18_S3_R1_001.fastq.gz', 3150366942, '3e76a04053fde5c0ee43e81d2c50c807'),
    (49503285, 'File3.fq.gz', 5723104027, 'b08fb297825ee97adee126c900bd50c6'),
]


def download(spec):
    fid, name, size, expected = spec
    raw = ROOT / 'raw'
    raw.mkdir(parents=True, exist_ok=True)
    chunk_mb = int(os.environ.get('ANTK_CHUNK_MB', '8'))
    parts = raw / (name + ('.parts' if chunk_mb == 8 else f'.parts{chunk_mb}M'))
    parts.mkdir(exist_ok=True)
    url = (f'https://ndownloader.figshare.com/files/{fid}' if fid == 49503285
           else f'https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/{fid}/{name}')
    chunk = chunk_mb * 1024 * 1024
    if chunk_mb == 1:
        for old in (raw / (name + '.parts')).glob('[0-9]*'):
            data = old.read_bytes()
            for offset in range(0, len(data), chunk):
                part_id = int(old.name)*8 + offset//chunk
                target = parts / f'{part_id:05d}'
                if not target.exists():
                    target.write_bytes(data[offset:offset+chunk])
    start = time.monotonic()
    # The public downloader issues a short-lived S3 link. Share each link among
    # range requests instead of requesting a new API redirect for every MB.
    class NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, req, fp, code, msg, headers, newurl):
            return None
    opener = urllib.request.build_opener(NoRedirect())
    link_lock = threading.Lock()
    signed = [None, 0.]

    def request_url():
        if fid != 49503285:
            return url
        with link_lock:
            if signed[0] is not None and time.monotonic()-signed[1] < 3:
                return signed[0]
            req = urllib.request.Request(url+'?fresh='+str(time.time_ns()), method='HEAD')
            try:
                opener.open(req, timeout=30)
                raise RuntimeError('expected public download redirect')
            except urllib.error.HTTPError as err:
                if err.code not in (301, 302, 303, 307, 308):
                    raise
                signed[:] = [err.headers['Location'], time.monotonic()]
            return signed[0]

    def part(i):
        lo, hi = i * chunk, min(size, (i + 1) * chunk) - 1
        p = parts / f'{i:05d}'
        if p.exists() and p.stat().st_size == hi - lo + 1:
            return
        for attempt in range(12):
            try:
                req = urllib.request.Request(request_url(), headers={'Range': f'bytes={lo}-{hi}'})
                with urllib.request.urlopen(req, timeout=120) as r:
                    assert r.status == 206, (r.status, r.headers)
                    assert r.headers['Content-Range'] == f'bytes {lo}-{hi}/{size}'
                    b = r.read()
                assert len(b) == hi - lo + 1
                p.write_bytes(b)
                return
            except Exception:
                if attempt == 11:
                    raise
                time.sleep(2 + attempt)

    n = (size + chunk - 1) // chunk
    with cf.ThreadPoolExecutor(int(os.environ.get('ANTK_DOWNLOAD_WORKERS', '64'))) as pool:
        futures = [pool.submit(part, i) for i in range(n)]
        for done, future in enumerate(cf.as_completed(futures), 1):
            future.result()
            if done % 16 == 0 or done == n:
                print(name, f'{done}/{n} chunks', round(time.monotonic() - start, 1), flush=True)
    md5, sha = hashlib.md5(), hashlib.sha256()
    target = raw / name
    with target.open('wb') as out:
        for i in range(n):
            b = (parts / f'{i:05d}').read_bytes()
            out.write(b)
            md5.update(b)
            sha.update(b)
    assert target.stat().st_size == size
    assert expected is None or md5.hexdigest() == expected
    record = dict(file_id=fid, name=name, bytes=size, md5=md5.hexdigest(),
                  expected_md5=expected, sha256=sha.hexdigest(), url=url,
                  seconds=time.monotonic() - start)
    (raw / (name + '.verified.json')).write_text(json.dumps(record, indent=2))
    print('VERIFIED', json.dumps(record), flush=True)
    for p in parts.iterdir():
        p.unlink()
    parts.rmdir()


if __name__ == '__main__':
    import sys
    for spec in FILES:
        if not sys.argv[1:] or str(spec[0]) in sys.argv[1:]:
            download(spec)
