#!/usr/bin/env python3
"""Build benchmark adapters while preserving the downloaded upstream sources."""
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[2]
BENCH = ROOT / 'benchmarks/antkowiak2020'
DATA = ROOT / 'input/antkowiak2020'
REPO = DATA / 'sources/dna_rs_coding'
PORTABLE = DATA / 'sources/portable_rs'
PORTABLE.mkdir(exist_ok=True)

for source in (REPO / 'include').glob('*.hpp'):
    text = subprocess.check_output(['git', '-C', str(REPO), 'show',
                                    '8e7e3b5:include/' + source.name], text=True)
    text = '\n'.join(line for line in text.splitlines()
                     if 'boost/' not in line and 'using namespace boost' not in line)
    if source.name == 'encodedecode.hpp':
        text = text.replace('uniform_int_distribution<uint>', 'LegacyUniform<uint>')
        # The deposited physical File 1 uses the spectral, non-systematic outer
        # representation. Select the upstream spectral decoder, not its default.
        text = text.replace('outercode.RS_shortened_decode(infrec[nb],cw[nb])',
                            'outercode.RS_decode_spec(infrec[nb],cw[nb])')
        # Fail clearly when decoding has failed, instead of allocating gigabytes
        # from an erroneous length header in the original demonstration program.
        text = text.replace('str.resize(filesize);',
                            'if(filesize > str.size()) throw std::runtime_error("invalid decoded file length");\n\tstr.resize(filesize);')
    (PORTABLE / source.name).write_text(text)

def run(args):
    subprocess.run(args, cwd=ROOT, check=True)

common = ['dbgps_codec.c', 'dbgps_redundancy.c', 'dbgps_strand.c', 'kthread.c']
run(['cc', '-O3', '-Wall', '-I.', str(BENCH/'graph_adapter.c')] + common +
    ['-lz', '-lpthread', '-o', str(DATA/'bin/graph_adapter')])
run(['cc', '-O3', str(BENCH/'prepare_reads.c'), '-lz', '-o', str(DATA/'bin/prepare_reads')])
run(['c++', '-O3', '-std=c++11', '-I'+str(BENCH), '-I'+str(PORTABLE),
     str(BENCH/'rs_driver.cpp'), '-o', str(DATA/'bin/rs_driver')])
run(['c++', '-O3', '-std=c++11', '-DFAST_GF14', '-I'+str(BENCH), '-I'+str(PORTABLE),
     str(BENCH/'rs_driver.cpp'), '-o', str(DATA/'bin/rs_driver_fast')])
run(['c++', '-O3', '-std=c++11', str(BENCH/'inner_select.cpp'),
     '-o', str(DATA/'bin/inner_select')])
