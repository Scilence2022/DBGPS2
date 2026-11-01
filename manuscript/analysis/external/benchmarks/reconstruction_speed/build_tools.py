"""Rebuild the graph adapter and inner-code selector from the saved snapshot."""
from pathlib import Path
import subprocess

root=Path(__file__).resolve().parent
source=root/'src'
(root/'bin').mkdir(exist_ok=True)
commands=[
    ['cc','-O3','-Wall','-o',str(root/'bin/graph_adapter'),
     str(source/'benchmarks/antkowiak2020/graph_adapter.c'),
     *[str(source/n) for n in ['dbgps_codec.c','dbgps_redundancy.c','dbgps_strand.c','kthread.c']],
     '-lz','-lpthread'],
    ['c++','-O3','-std=c++11','-o',str(root/'bin/inner_select'),str(source/'inner_select.cpp')],
    ['cc','-O3','-o',str(root/'bin/prepare_reads'),
     str(source/'benchmarks/antkowiak2020/prepare_reads.c'),'-lz']
]
for command in commands:
    subprocess.run(command,check=True)
subprocess.run([str(root/'bin/muscle3'),'-version'],check=True)
