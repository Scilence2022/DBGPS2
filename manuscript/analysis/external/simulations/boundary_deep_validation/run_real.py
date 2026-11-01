"""External behavior check on the already-inspected full raw HT176 library."""
from pathlib import Path
import csv
import json
import subprocess
from run_deep import digest

R=Path(__file__).resolve().parent
raw=['/Users/song/Documents/123/HT176_FDDP220002706-1a_1.fq.gz',
     '/Users/song/Documents/123/HT176_FDDP220002706-1a_2.fq.gz']
ref='/Users/song/dbgps2_ht/ref.fa';prefix=R/'real_HT176_raw'
cmd=[str(R/'bin/blind_profile'),'-x',ref,'-M',str(R/'src/legacy-2base.map'),
     '-a','101010102','-b','101295684','-k','31','-t','3','-N','-v','4093',
     '-c','1,3,5','-r','1.25,1.5,2,2.55,3,4,8,0','-o',str(prefix),*raw]
lock=dict(command=cmd,script_sha256=digest(Path(__file__)),binary_sha256=digest(R/'bin/blind_profile'),
          policies_sha256=digest(R/'blind_policy.py'),interpretation='Previously inspected real library; not an unseen test set. All queried IDs enter blind objectives, including absent designs.')
lp=R/'real_design_and_source_lock.json'
if lp.exists():assert json.loads(lp.read_text())==lock
else:lp.write_text(json.dumps(lock,indent=2))
if not Path(str(prefix)+'.json').exists():
    with Path(str(prefix)+'.stdout').open('w') as out,Path(str(prefix)+'.stderr').open('w') as err:
        p=subprocess.run(cmd,stdout=out,stderr=err)
    assert p.returncode==0,p.returncode
    grid=list(csv.DictReader(Path(str(prefix)+'.grid.tsv').open(),delimiter='\t'));assert len(grid)==24
    prior=list(csv.DictReader((R.parent/'reconstruction_boundary/real_HT176_raw.grid.tsv').open(),delimiter='\t'))
    checks=0
    for row in grid:
        old=next((a for a in prior if a['c']==row['c'] and a['r']==row['r']),None)
        if old:
            for field in row:
                if field not in ['seconds','native_checks']:assert row[field]==old[field],(field,row,old)
            checks+=1
    assert checks==16
    Path(str(prefix)+'.json').write_text(json.dumps(dict(command=cmd,k=31,
        settings=dict(cohort='real_HT176_raw',n=210000,seed=0,geometry='real',length=0,depth=0,recovery=0,error=0,channel='real',sampling='real',abundance_sigma=0),
        grid=grid,prior_points_identical=checks),indent=2))
print('REAL_BLIND_ALL_COMPLETE',flush=True)
