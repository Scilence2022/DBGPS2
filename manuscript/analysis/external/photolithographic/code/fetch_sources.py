#!/usr/bin/env python3
"""Fetch the pinned public code and design files used by this benchmark."""
import concurrent.futures
import hashlib
import json
import pathlib
import subprocess
import tarfile
import urllib.request

ROOT=pathlib.Path(__file__).resolve().parents[2]/'input/antkowiak2020'
SOURCES=ROOT/'sources'
PINS={
    'noisy_dna_data_storage':('MLI-lab/noisy_dna_data_storage','41422be448b3c1da9041c3a2967d2754d2ad8b1a'),
    'dna_rs_coding':('reinhardh/dna_rs_coding','455e1a5182bf64281f54d434a4ed5ca456c2373e'),
}
DT='861d3b95b52d1b121e9d85040f20c8ecf0cda2d4'

def main():
    for d in ['sources','raw','processed','results','bin']:(ROOT/d).mkdir(parents=True,exist_ok=True)
    for name,(repo,rev) in PINS.items():
        target=SOURCES/name
        if not (target/'.git').exists():
            subprocess.run(['git','clone',f'https://github.com/{repo}.git',str(target)],check=True)
        actual=subprocess.check_output(['git','-C',str(target),'rev-parse','HEAD'],text=True).strip()
        if actual!=rev:subprocess.run(['git','-C',str(target),'checkout','--detach',rev],check=True)
    def fetch(i):
        u=f'https://raw.githubusercontent.com/fml-ethz/dt4dds-challenges_notebooks/{DT}/data_experimental/Photolithographic_Antkowiak/File{i}/design_files.fasta'
        p=SOURCES/f'File{i}_design.fasta';b=urllib.request.urlopen(u,timeout=120).read();p.write_bytes(b)
        return dict(file=p.name,url=u,sha256=hashlib.sha256(b).hexdigest(),bytes=len(b))
    with concurrent.futures.ThreadPoolExecutor(3) as pool:manifest=list(pool.map(fetch,[1,2,3]))
    for article in [12128994,12129090,27134316]:
        b=urllib.request.urlopen(f'https://api.figshare.com/v2/articles/{article}',timeout=60).read()
        (SOURCES/f'figshare_{article}.json').write_bytes(b)
    (ROOT/'results/source_pins.json').write_text(json.dumps(dict(repositories=PINS,dt4dds_commit=DT,files=manifest),indent=2))
    print('Sources fetched. MUSCLE 3.8.31 source is available at https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_src.tar.gz')

if __name__=='__main__':main()
