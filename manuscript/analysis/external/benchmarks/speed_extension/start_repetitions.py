from pathlib import Path
import json,time,subprocess,sys
r=Path(__file__).resolve().parent
while True:
 try:
  settings=json.loads((r/'selected_settings.json').read_text());tests=json.loads((r/'smoke_validation.json').read_text());inputs=json.loads((r/'input_provenance.json').read_text())
  assert len(settings)==3 and len(tests)==3 and len(inputs)>=10
  break
 except (FileNotFoundError,json.JSONDecodeError,AssertionError):time.sleep(2)
print('All preparation and smoke-test gates complete',flush=True)
subprocess.run([sys.executable,str(r/'run_benchmark.py')],check=True)
