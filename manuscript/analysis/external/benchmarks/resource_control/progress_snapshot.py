from pathlib import Path
import json
from datetime import datetime,timezone
R=Path(__file__).resolve().parent
try:rows=json.loads((R/'resource_results.json').read_text())
except FileNotFoundError:rows=[]
lines=(R/'resource_matrix.log').read_text().splitlines();start=next((x for x in reversed(lines) if x.startswith('START ')),None)
print('Completed',len(rows),'/162')
if start:
 print(start);_,mode,name,stamp=start.split(' ',3);elapsed=(datetime.now(timezone.utc)-datetime.fromisoformat(stamp)).total_seconds();print('Current measurement elapsed',round(elapsed),'seconds');p=R/mode/'runs'/(name+'.stdout')
 if p.exists():
  output=p.read_text().splitlines();progress=[x for x in output if x.startswith(('CLUSTERED','MSA '))]
  for line in progress[-2:]:print(line)
for line in lines[-3:]:
 if line.startswith('DONE'):print(line)

if any('Traceback' in line for line in lines[-20:]):
 print('MEASUREMENT_ERROR')
 print('\n'.join(lines[-8:]))
