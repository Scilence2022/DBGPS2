#!/usr/bin/env python3
"""Give an isolated experimental binary an explicit read-counting length cap."""
import hashlib
import json
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'input/antkowiak2020';BENCH=ROOT/'benchmarks/antkowiak2020'
BUILD=DATA/'build_block_rs';BUILD.mkdir(exist_ok=True)
base=DATA/'build_inner_rs/graph_inner_rs.generated.c';source=base.read_text()
old='int n=atoi(argv[2]),k=atoi(argv[3]);'
assert source.count(old)==1
source=source.replace(old,old+'\n    const char *limit_env=getenv("ANTK_READ_LIMIT");\n    int read_limit=limit_env?atoi(limit_env):60;\n    if(read_limit<60||read_limit>300)return 8;')
for old,new in [('10000000,8,60,filter_shift','10000000,8,read_limit,filter_shift'),('10000000,8,60);','10000000,8,read_limit);')]:
    assert source.count(old)==1,old;source=source.replace(old,new)
generated=BUILD/'graph_block_rs.generated.c';generated.write_text(source);objects=[]
for path in [generated]+[ROOT/p for p in ['dbgps_codec.c','dbgps_redundancy.c','dbgps_strand.c','kthread.c']]:
    obj=BUILD/(path.stem+'.o');objects.append(str(obj))
    subprocess.run(['cc','-O3','-Wall','-I'+str(ROOT),'-I'+str(BENCH),'-c',str(path),'-o',str(obj)],check=True)
subprocess.run(['c++','-O3','-std=c++11','-I'+str(BENCH),str(BENCH/'graph_inner_rs.cpp')]+objects+['-lz','-lpthread','-o',str(DATA/'bin/graph_block_rs')],check=True)
(BUILD/'source_manifest.json').write_text(json.dumps(dict(base_generated_sha256=hashlib.sha256(base.read_bytes()).hexdigest(),
    generated_sha256=hashlib.sha256(generated.read_bytes()).hexdigest(),
    change='Only the read-counting length cap is configurable through ANTK_READ_LIMIT; default 60',
    binary_sha256=hashlib.sha256((DATA/'bin/graph_block_rs').read_bytes()).hexdigest()),indent=2))
