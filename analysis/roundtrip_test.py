#!/usr/bin/env python3
"""Round-trip through the REAL decode_fountain.py path: encode a synthetic
file with our PRNG port (random 32-bit seeds, screened like the original
encoder), then decode via decode_fountain.py and check md5."""
import hashlib
import random
import subprocess
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from fountain_common import PRNG, ReedSolomon2, bytes_to_dna, screen_dna

K, CS = 2000, 32
random.seed(11)
chunks = [random.randrange(256) for _ in range(K * CS)]
file_bytes = bytes(chunks)

prng = PRNG(K=K, delta=0.001, c=0.025)
rs = ReedSolomon2()

droplets = []
while len(droplets) < int(K * 1.15):
    seed = random.randrange(1 << 32)
    _, d, nums = prng.get_src_blocks(seed)
    data = bytearray(CS)
    for n in nums:
        base = n * CS
        for i in range(CS):
            data[i] ^= file_bytes[base + i]
    word = rs.encode(list(seed.to_bytes(4, 'big')) + list(data))
    dna = bytes_to_dna(bytes(word))
    if screen_dna(dna[:136]):
        droplets.append(dna)
print("encoded droplets:", len(droplets), flush=True)
random.shuffle(droplets)

with open('/tmp/rt_droplets.txt', 'w') as f:
    f.write('\n'.join(droplets) + '\n')

out = subprocess.run(
    [sys.executable, __file__.rsplit('/', 1)[0] + '/decode_fountain.py',
     '/tmp/rt_droplets.txt', '--out', '/tmp/rt_decoded.bin',
     '-n', str(K)],
    capture_output=True, text=True)
print(out.stderr)
print(out.stdout)
got = open('/tmp/rt_decoded.bin', 'rb').read()
want = file_bytes
print("size:", len(got), "expected:", len(want))
print("md5 match:", hashlib.md5(got).hexdigest() == hashlib.md5(want).hexdigest())
