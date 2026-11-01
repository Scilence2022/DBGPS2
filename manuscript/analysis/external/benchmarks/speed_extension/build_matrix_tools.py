"""Rebuild native graph callers; NGmerge/MUSCLE sources are pinned separately."""
from pathlib import Path
import json,platform,subprocess
R=Path(__file__).resolve().parent;S=R/'src';(R/'bin').mkdir(exist_ok=True)
commands=[]
for stem in ['format_graph','format_graph_bloom']:
 cmd=['cc','-O3','-o',str(R/'bin'/stem),str(S/(stem+'.c')),*[str(S/n) for n in ['dbgps_codec.c','dbgps_redundancy.c','dbgps_strand.c','kthread.c']],'-lz','-lpthread']
 if platform.system()=='Linux':cmd.append('-lm')
 subprocess.run(cmd,check=True);commands.append(cmd)
(R/'sources/rebuild_commands.json').write_text(json.dumps(commands,indent=2))
