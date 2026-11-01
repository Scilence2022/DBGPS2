#!/usr/bin/env python3
"""Add candidate callbacks to an ignored copy; leave the strict adapter intact."""
import hashlib
import json
import pathlib
import subprocess

ROOT=pathlib.Path(__file__).resolve().parents[2]
BENCH=ROOT/'benchmarks/antkowiak2020'
DATA=ROOT/'input/antkowiak2020'
BUILD=DATA/'build_inner_rs';BUILD.mkdir(exist_ok=True)
source=(BENCH/'graph_adapter.c').read_text()
base_sha=hashlib.sha256(source.encode()).hexdigest()
def replace(old,new):
    global source
    assert source.count(old)==1,old
    source=source.replace(old,new)
replace('#include "../../DBGPS2.c"','#include '+json.dumps(str(ROOT/'DBGPS2.c')))
replace('#include <sys/time.h>','#include <sys/time.h>\n#include "graph_inner_rs.h"\nstatic const char *evaluation_reference;')
replace('char fname[1024];','char fname[1024],inner_prefix[1024];\n    snprintf(inner_prefix,sizeof(inner_prefix),"%s.k%d.c%d.r%.2f",prefix,k,cut,ratio);\n    inner_begin(n,inner_prefix);')
replace('if(!overflow)for(int j=0;j<count;++j)if(rs_valid(a[j].dna)) {',
        'if(!overflow)for(int j=0;j<count;++j)if(inner_add(a[j].dna,id)==0) {')
replace('fclose(out);fclose(status);free(a);free(b);',
        'fclose(out);fclose(status);free(a);free(b);\n    inner_finish(evaluation_reference);')
replace('int n=atoi(argv[2]),k=atoi(argv[3]);','int n=atoi(argv[2]),k=atoi(argv[3]);\n    evaluation_reference=argv[5];')
replace('if(argc>6) {','''const char *grid=getenv("ANTK_GRID");
    if(grid) {
        char *pairs=strdup(grid),*save=NULL;
        for(char *p=strtok_r(pairs,",",&save);p;p=strtok_r(NULL,",",&save)) {
            char *colon=strchr(p,':');if(!colon)return 3;*colon=0;
            run(h,seeds,ns,n,k,atoi(colon+1),atof(p),argv[4]);
        }
        free(pairs);
    } else if(argc>6) {''')
generated=BUILD/'graph_inner_rs.generated.c'
generated.write_text('/* Generated from graph_adapter.c SHA256 '+base_sha+' */\n'+source)
objects=[]
for path in [generated]+[ROOT/p for p in ['dbgps_codec.c','dbgps_redundancy.c','dbgps_strand.c','kthread.c']]:
    obj=BUILD/(path.stem+'.o');objects.append(str(obj))
    subprocess.run(['cc','-O3','-Wall','-I'+str(ROOT),'-I'+str(BENCH),'-c',str(path),'-o',str(obj)],check=True)
subprocess.run(['c++','-O3','-std=c++11','-I'+str(BENCH),str(BENCH/'graph_inner_rs.cpp')]+objects+['-lz','-lpthread','-o',str(DATA/'bin/graph_inner_rs')],check=True)
(BUILD/'source_manifest.json').write_text(json.dumps(dict(base_adapter_sha256=base_sha,generated_sha256=hashlib.sha256(generated.read_bytes()).hexdigest()),indent=2))
