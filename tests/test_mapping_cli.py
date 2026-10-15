"""Generate strands independently and test configurable codecs through the graph."""

import binascii
from pathlib import Path
import random
import subprocess
import sys
import tempfile


PRIMERS = ("CCTGCAGAGTAGCATGTC", "CTGACACTGATGCATCCG")


def word(value, bases):
    return "".join("ACGT"[(value >> shift) & 3] for shift in range(2 * bases - 2, -1, -2))


def encode(raw, codewords, bits, padding_value=0):
    binary = "".join(f"{byte:08b}" for byte in raw)
    padding = -len(binary) % bits
    binary += f"{padding_value:0{padding}b}" if padding else ""
    return "".join(codewords[int(binary[pos:pos + bits], 2)] for pos in range(0, len(binary), bits))


def main():
    binary = str(Path(sys.argv[1]).resolve())
    repo = Path(__file__).resolve().parent.parent
    legacy = {}
    for line in (repo / "examples/legacy-2base.map").read_text().splitlines():
        fields = line.split("#", 1)[0].split()
        if len(fields) == 2 and set(fields[0]) <= set("ACGT"):
            legacy[int(fields[1], 2)] = fields[0]
    assert len(legacy) == 16
    rng = random.Random(182)
    # Includes partial packed DNA bytes, CRC padding >= 8 bits, reads > 200 bp,
    # nonseparable permutations, and index values beyond the former uint32 limit.
    cases = [(1, 2, 16, 140), (2, 4, 16, 140), (4, 8, 20, 140),
             (3, 6, 24, 120), (3, 5, 24, 120), (2, 3, 16, 80),
             (1, 1, 16, 80), (5, 8, 15, 250), (7, 8, 14, 98),
             (8, 16, 16, 144), (6, 12, 24, 96), (8, 15, 0, 64),
             (2, 4, 16, 0), (2, 4, 16, 140)]
    with tempfile.TemporaryDirectory(prefix="dbgps-mapping-") as directory:
        directory = Path(directory)
        for case, (bases, bits, index_bases, payload_bases) in enumerate(cases):
            codewords = [word(x, bases) for x in rng.sample(range(4 ** bases), 2 ** bits)]
            if case == len(cases) - 1:
                codewords = [legacy[value] for value in range(16)]
            mapping_file = directory / f"case{case}.map"
            mapping_file.write_text(f"bases {bases}\nbits {bits}\n" + "".join(
                f"{dna} {value:0{bits}b}\n" for value, dna in enumerate(codewords)))
            index_bits = index_bases // bases * bits
            payload_bytes = payload_bases // bases * bits // 8
            first = (1 << (index_bits - 2)) + 1000 if index_bits else 0
            indices = (first, first + 1) if index_bits else (0,)
            expected = {}
            raw_by_index = {}
            for index in indices:
                raw = index.to_bytes(index_bits // 8, "big") + bytes(rng.randrange(256) for _ in range(payload_bytes))
                raw += binascii.crc_hqx(raw, 0).to_bytes(2, "big")
                raw_by_index[index] = raw
                expected[index] = encode(raw, codewords, bits)
            reads = directory / f"case{case}.fa"

            def write_reads(strands):
                reads.write_text("".join(f">{index}_{copy}\n{PRIMERS[0]}{dna}{PRIMERS[1]}\n"
                                         for index, dna in strands.items() for copy in range(8)))

            write_reads(expected)
            k = 31 if index_bases else 17
            common = [binary, "-M", str(mapping_file), "-k", str(k), "-i", str(index_bases),
                      "-l", str(payload_bases), "-a", str(indices[0]), "-b", str(indices[-1])]

            def run(extra=(), repeats=1, wanted=None, use_mapping=True):
                prefix = directory / f"output{case}"
                options = common if use_mapping else common[:1] + common[3:]
                args = options + list(extra) + ["-o", str(prefix), str(reads)]
                result = subprocess.run(args, capture_output=True, text=True, timeout=30)
                assert result.returncode == 0, (args, result.stderr)
                rows = [line.split("\t") for line in prefix.with_suffix(".strands").read_text().splitlines()]
                target = expected if wanted is None else wanted
                assert len(rows) == repeats * len(target), (args, rows, result.stdout)
                for index, dna in rows:
                    assert dna == target[int(index)], (args, index, dna)
                for index in target:
                    assert sum(int(row[0]) == index for row in rows) == repeats

            run(["-t", "1"])
            run(["-B", "12", "-F", "23", "-t", "3"])
            if case in (1, 4, 7):
                run(["-t", "3", "-C", "2", "-R", "2.65"], repeats=6)
            if bits == 5:
                # A legal DNA word can carry nonzero padding while the CRC itself
                # remains correct. It must still be rejected.
                write_reads({i: encode(raw, codewords, bits, 1) for i, raw in raw_by_index.items()})
                run(wanted={})
                write_reads({i: encode(raw[:-1] + bytes([raw[-1] ^ 1]), codewords, bits)
                             for i, raw in raw_by_index.items()})
                run(wanted={})
            if case == 5:
                forbidden = next(word(x, bases) for x in range(4 ** bases) if word(x, bases) not in codewords)
                write_reads({i: dna[:index_bases] + forbidden + dna[index_bases + bases:]
                             for i, dna in expected.items()})
                run(wanted={})
            if case == len(cases) - 1:
                # The legacy permutation is opt-in; default canonical seeds
                # must not silently reconstruct strands using that codebook.
                run(wanted={}, use_mapping=False)
            print(f"mapping CLI: {bases} bases/{bits} bits, index={index_bases}, payload={payload_bases} passed")

        invalid_options = [("-B", "-1"), ("-B", "65"), ("-B", "12x"),
                           ("-B", "999999999999999999999999"), ("-F", "22"),
                           ("-F", "41"), ("-F", "27x"), ("-i", "15"), ("-l", "139"), ("-k", "32"), ("-i", "32"),
                           ("-a", "-1"), ("-b", "18446744073709551616"),
                           ("-a", "1001", "-b", "1000"), ("-b", str(1 << 32)),
                           ("-M", str(directory / "missing.map"))]
        for options in invalid_options:
            result = subprocess.run([binary, "-a", "1000", "-b", "1000"] + list(options) + [str(reads)],
                                    capture_output=True, text=True, timeout=10)
            assert result.returncode == 1 and "ERROR:" in result.stderr, (options, result)
            assert "Counting k-mers" not in result.stderr
        print("mapping CLI: incompatible layouts and invalid index ranges rejected before counting")


if __name__ == "__main__":
    main()
