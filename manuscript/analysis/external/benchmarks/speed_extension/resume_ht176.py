from pathlib import Path
import subprocess,sys,time,json
r=Path(__file__).resolve().parent
while not (r/'data/song_ht176.preparation.json').exists():time.sleep(2)
for script,args in [('smoke_test.py',[]),('fragment_smoke_test.py',[]),('provenance.py',[]),('run_benchmark.py',['--pilot','--only','song_ht176'])]:
 print('PREPARE',script,flush=True);subprocess.run([sys.executable,str(r/script),*args],check=True)
print('START REMAINING REPEATS WITH HT176',flush=True)
subprocess.run([sys.executable,str(r/'run_benchmark.py')],check=True)
