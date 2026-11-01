#!/usr/bin/env python3
"""Build the exact MUSCLE 3.8.31 source on the local architecture."""
import concurrent.futures
import hashlib
import pathlib
import re
import subprocess
import tarfile
import urllib.request

ROOT=pathlib.Path(__file__).resolve().parent
ARCHIVE=ROOT/'sources/muscle3.8.31_src.tar.gz'
EXPECTED='43c5966a82133bd7da5921e8142f2f592c2b5f53d802f0527a2801783af809ad'

def main():
    if not ARCHIVE.exists():
        ARCHIVE.write_bytes(urllib.request.urlopen('https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_src.tar.gz',timeout=120).read())
    assert hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()==EXPECTED
    destination=ROOT/'sources/muscle3';destination.mkdir(exist_ok=True)
    with tarfile.open(ARCHIVE) as archive:
        for member in archive.getmembers():
            if not member.isfile():continue
            p=(destination/member.name).resolve()
            assert destination.resolve() in p.parents
            p.parent.mkdir(parents=True,exist_ok=True)
            p.write_bytes(archive.extractfile(member).read())
    source=destination/'muscle3.8.31/src'
    names=re.search(r"CPPNames='([^']+)'",(source/'mk').read_text()).group(1).split()
    def compile_one(name):
        subprocess.run(['c++','-std=c++11','-O3','-DNDEBUG=1','-D_FILE_OFFSET_BITS=64','-c',name+'.cpp','-o',name+'.o'],cwd=source,check=True)
    with concurrent.futures.ThreadPoolExecutor(8) as pool:list(pool.map(compile_one,names))
    subprocess.run(['c++','-O3','-o',str((ROOT/'bin/muscle3').resolve())]+[str(source/(n+'.o')) for n in names],check=True)
    subprocess.run([str(ROOT/'bin/muscle3'),'-version'],check=True)

if __name__=='__main__':main()
