"""Independent byte encoders exercise both graph engines and their wire format."""

import binascii
from pathlib import Path
import random
import subprocess
import sys
import tempfile

from test_mapping_cli import encode, word, PRIMERS


def gf_mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 256:
            a ^= 0x11d
    return result


def canonical_spec(spec):
    if not spec.startswith("rs"):
        return spec
    amount = spec.split(":")[1] if ":" in spec else "16"
    return f"rs:{int(amount[:-1] if amount.endswith('B') else amount)}B"


def protect(raw, spec):
    if spec == "none":
        return raw
    if spec == "crc16":
        return raw + binascii.crc_hqx(raw, 0).to_bytes(2, "big")
    if spec == "crc32":
        return raw + binascii.crc32(raw).to_bytes(4, "big")
    if spec == "crc32c":
        crc = 0xffffffff
        for byte in raw:
            crc ^= byte
            for _ in range(8):
                crc = (crc >> 1) ^ (0x82f63b78 if crc & 1 else 0)
        return raw + (crc ^ 0xffffffff).to_bytes(4, "big")
    if spec.startswith("rs"):
        parity = int(canonical_spec(spec).split(":")[1][:-1])
        generator = [1]
        root = 1
        for _ in range(parity):
            product = [0] * (len(generator) + 1)
            for i, coefficient in enumerate(generator):
                product[i] ^= coefficient
                product[i + 1] ^= gf_mul(coefficient, root)
            generator = product
            root = gf_mul(root, 2)
        remainder = list(raw) + [0] * parity
        for i in range(len(raw)):
            factor = remainder[i]
            for j, coefficient in enumerate(generator):
                remainder[i + j] ^= gf_mul(factor, coefficient)
        return raw + bytes(remainder[len(raw):])
    if spec == "hamming":
        trailer = ""
        for byte in raw:
            positions = dict(zip((3, 5, 6, 7, 9, 10, 11, 12), map(int, f"{byte:08b}")))
            for position in (1, 2, 4, 8):
                positions[position] = sum(value for at, value in positions.items()
                                          if at & position) % 2
            overall = sum(positions.values()) % 2
            trailer += "".join(str(positions[pos]) for pos in (8, 4, 2, 1)) + str(overall)
        trailer += "0" * (-len(trailer) % 8)
        return raw + bytes(int(trailer[i:i + 8], 2) for i in range(0, len(trailer), 8))
    raise AssertionError(spec)


def main():
    binaries = [str(Path(arg).resolve()) for arg in sys.argv[1:]]
    rng = random.Random(20261015)
    cases = [(1, 2, 16, 80), (2, 4, 16, 80), (3, 5, 24, 48),
             (3, 6, 24, 48), (8, 15, 0, 64), (1, 2, 16, 0),
             (1, 2, 16, 224)]  # long canonical reads must not be truncated at 200
    with tempfile.TemporaryDirectory(prefix="dbgps-ecc-") as temporary:
        directory = Path(temporary)
        for case, (bases, bits, index_bases, payload_bases) in enumerate(cases):
            mapping_args = []
            words = [word(i, bases) for i in range(1 << bits)]
            if bases != 1:
                words = [word(i, bases) for i in rng.sample(range(4 ** bases), 1 << bits)]
                mapping = directory / "mapping.map"
                mapping.write_text(f"bases {bases}\nbits {bits}\n" + "".join(
                    f"{dna} {i:0{bits}b}\n" for i, dna in enumerate(words)))
                mapping_args = ["-M", str(mapping)]
            index_bytes = index_bases // bases * bits // 8
            index = 1000 if index_bytes else 0
            raw = index.to_bytes(index_bytes, "big") + bytes(rng.randrange(256) for _ in
                                                          range(payload_bases // bases * bits // 8))
            specs = ["none", "crc16", "crc32", "crc32c", "rs:4B", "rs:8B", "rs:4", "rs", "hamming"]
            if case == len(cases) - 1:
                specs = ["rs:4B"]  # fits DBGPS2's total <300 bp limit
            for spec in specs:
                protected = protect(raw, spec)
                expected = encode(protected, words, bits)
                reads = directory / "reads.fa"

                def write_reads(strands):
                    reads.write_text("".join(f">{i}_{copy}\n{PRIMERS[0]}{dna}{PRIMERS[1]}\n"
                                             for i, dna in enumerate(strands) for copy in range(8)))

                def run(wanted, extra=()):
                    for binary in binaries:
                        prefix = directory / Path(binary).name
                        args = [binary, *mapping_args, "-E", spec, "-k", "31" if index_bytes else "17",
                                "-i", str(index_bases), "-l", str(payload_bases), "-a", str(index),
                                "-b", str(index), "-t", "2", "-o", str(prefix), *extra, str(reads)]
                        result = subprocess.run(args, capture_output=True, text=True, timeout=30)
                        assert result.returncode == 0, (args, result.stderr)
                        label = f"Redundancy (-E): {canonical_spec(spec)};"
                        assert label in result.stdout and label in prefix.with_suffix(".log").read_text(), result.stdout
                        rows = prefix.with_suffix(".strands").read_text().splitlines()
                        assert rows == ([f"{index}\t{wanted}"] if wanted is not None else []), (
                            binary, case, spec, rows, wanted, result.stdout, result.stderr)

                write_reads([expected])
                run(expected, ("-B", "12", "-F", "23"))
                if payload_bases and case != len(cases) - 1:
                    damaged = bytearray(protected)
                    damaged[index_bytes] ^= 0x80
                    damaged_dna = encode(damaged, words, bits)
                    write_reads([damaged_dna])
                    correcting = spec.startswith("rs") or spec == "hamming"
                    run(expected if correcting else damaged_dna if spec == "none" else None)
                    if correcting:
                        write_reads([expected, damaged_dna])
                        run(expected)  # equivalent corrected candidates are one result
                        damaged[-1] ^= 0x80  # payload and parity errors in distinct codewords
                        write_reads([encode(damaged, words, bits)])
                        run(expected)
                if -len(protected) * 8 % bits:
                    write_reads([encode(protected, words, bits, 1)])
                    run(None)  # mapping padding is never repaired or ignored
                print(f"redundancy CLI: {spec}, {bases} bases/{bits} bits, payload={payload_bases} passed")

        for binary in binaries:
            for spec in ("unknown", "rs:", "rs:0", "rs:255", "rs:-1", "rs:4x", "crc16:4", "rs:254",
                         "rs:8b", "rs:B", "rs:0B", "rs:255B", "rs:8BB", "rs:8Bjunk"):
                result = subprocess.run([binary, "-E", spec, "-a", "0", "-b", "0", str(reads)],
                                        capture_output=True, text=True, timeout=10)
                assert result.returncode == 1 and "ERROR:" in result.stderr, (spec, result)
                assert "Counting k-mers" not in result.stderr


if __name__ == "__main__":
    main()
