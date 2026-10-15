"""Check strand reconstruction and reproducibility through the real counters."""

import binascii
import gzip
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def main():
    repo = Path(__file__).resolve().parent.parent
    def encode(data):
        return "".join("ACGT"[(byte >> shift) & 3] for byte in data for shift in (6, 4, 2, 0))

    rng = random.Random(17)
    primers = ("CCTGCAGAGTAGCATGTC", "CTGACACTGATGCATCCG")
    raw_strands = {}
    for index in (1000, 1001):
        data = index.to_bytes(4, "big") + bytes(rng.randrange(256) for _ in range(35))
        crc = binascii.crc_hqx(data, 0).to_bytes(2, "big")
        raw_strands[index] = data + crc

    with tempfile.TemporaryDirectory(prefix="dbgps-regression-") as directory:
        directory = Path(directory)
        for binary in sys.argv[1:]:
            has_crc = Path(binary).name == "DBGPS2"
            expected = {index: encode(raw) for index, raw in raw_strands.items()}
            inputs = []
            for mate in (1, 2):
                path = directory / f"{Path(binary).name}-reads{mate}.fq.gz"
                with gzip.open(path, "wt") as out:
                    for index, strand in expected.items():
                        sequence = primers[0] + strand + primers[1]
                        assert len(sequence) == 200
                        for copy in range(10):
                            out.write(f"@{index}_{copy}/{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n")
                inputs.append(str(path))
            expected_rows = {index: sequence if has_crc else sequence[:-8]
                             for index, sequence in expected.items()}
            for threads in (1, 3):
                for repeat in (0, 1):
                    for sweep in (False, True):
                        prefix = directory / f"{Path(binary).name}-{threads}-{repeat}-{sweep}"
                        args = [binary, "-a", "1000", "-b", "1001", "-t", str(threads), "-o", str(prefix)]
                        copies = 1
                        if sweep:
                            args += ["-C", "2", "-R", "2.65"]
                            copies = 6  # two cutoffs at each of three coverage ratios
                        result = subprocess.run(args + inputs, capture_output=True, text=True, timeout=30)
                        assert result.returncode == 0, result.stderr
                        assert "Mapping: none" in result.stdout and "direct 2-bit packing" in result.stdout
                        rows = [line.split("\t") for line in prefix.with_suffix(".strands").read_text().splitlines()]
                        assert len(rows) == 2 * copies, (args, rows, result.stdout, result.stderr)
                        for index, sequence in rows:
                            assert sequence == expected_rows[int(index)], (args, index, sequence)
                        for index in expected_rows:
                            assert sum(int(row[0]) == index for row in rows) == copies
            prefix = directory / "explicit-canonical"
            args = [binary, "-M", str(repo / "examples/acgt-1base.map"),
                    "-a", "1000", "-b", "1001", "-o", str(prefix)]
            result = subprocess.run(args + inputs, capture_output=True, text=True, timeout=30)
            assert result.returncode == 0, result.stderr
            rows = [line.split("\t") for line in prefix.with_suffix(".strands").read_text().splitlines()]
            assert len(rows) == len(expected_rows)
            assert {int(index): dna for index, dna in rows} == expected_rows
            print(f"{Path(binary).name}: exact CLI strands, repeatability, thread counts and ratio sweep passed")


if __name__ == "__main__":
    main()
