"""Check that DBGPS2-RPC reproduces DBGPS2 batch output and serves the same strands."""

import binascii
import gzip
import json
from pathlib import Path
import random
import re
import socket
import subprocess
import sys
import tempfile
import time

PRIMERS = ("CCTGCAGAGTAGCATGTC", "CTGACACTGATGCATCCG")
FIRST_INDEX, STRANDS = 1000, 40
RANGE = ["-a", str(FIRST_INDEX), "-b", str(FIRST_INDEX + STRANDS - 1)]
TIMING = re.compile(r"^(K-mer counting|Strand decoding) time: .*\n", re.M)
REPO = Path(__file__).resolve().parent.parent
LEGACY_MAP = str(REPO / "examples" / "legacy-2base.map")
CANONICAL_MAP = str(REPO / "examples" / "acgt-1base.map")

CASES = (
    [],
    ["-c", "2"],
    ["-C", "3"],
    ["-r", "1.5", "-R", "1.7"],
    ["-s", "8"],
    ["-k", "25"],
    ["-B", "12"],
    ["-F", "23"],
    ["-F", "23", "-B", "12"],
    ["-E", "none", "-l", "156"],  # the following primer base is the lookahead
    ["-E", "none", "-l", "140"],
    ["-E", "crc32"],
    ["-M", CANONICAL_MAP],
    ["-M", LEGACY_MAP],
)
SERVED_CASES = ([], ["-E", "none", "-l", "156"], ["-F", "23", "-B", "12"])


def encode(data):
    return "".join("ACGT"[(byte >> shift) & 3] for byte in data for shift in (6, 4, 2, 0))


def write_reads(directory):
    """Write noisy reads whose lookahead base after half of the strands branches."""
    rng = random.Random(23)
    paths = [directory / "reads1.fq.gz", directory / "reads2.fq.gz"]
    outputs = [gzip.open(path, "wt") for path in paths]
    for offset in range(STRANDS):
        data = (FIRST_INDEX + offset).to_bytes(4, "big") + bytes(rng.randrange(256) for _ in range(35))
        strand = encode(data + binascii.crc_hqx(data, 0).to_bytes(2, "big"))
        for copy in range(24):
            tail = PRIMERS[1]
            if offset % 2 and copy % 2:
                tail = "G" + tail[1:]
            read = list(PRIMERS[0] + strand + tail)
            for position in range(len(read)):
                if rng.random() < 0.004:
                    read[position] = rng.choice("ACGT")
            read = "".join(read)
            outputs[copy % 2].write(f"@{offset}_{copy}\n{read}\n+\n{'I' * len(read)}\n")
    for output in outputs:
        output.close()
    return [str(path) for path in paths]


def batch(binary, prefix, options, reads, to_files=True):
    args = [binary, *RANGE, *options]
    if to_files:
        args += ["-o", str(prefix)]
    result = subprocess.run(args + reads, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, (args, result.stderr)
    stdout = TIMING.sub("", result.stdout)
    if not to_files:
        return stdout.replace(str(prefix), "PREFIX"), "", ""
    log = TIMING.sub("", prefix.with_suffix(".log").read_text()).replace(str(prefix), "PREFIX")
    return stdout.replace(str(prefix), "PREFIX"), log, prefix.with_suffix(".strands").read_text()


def served_strands(binary, options, reads, log_path):
    with socket.socket() as reservation:
        reservation.bind(("127.0.0.1", 0))
        port = reservation.getsockname()[1]
    with log_path.open("w") as log:
        process = subprocess.Popen([binary, "-S", "-P", str(port), *options, *reads], stdout=log, stderr=log)
        try:
            deadline = time.monotonic() + 30
            while "Listening on port" not in log_path.read_text():
                assert process.poll() is None, log_path.read_text()
                assert time.monotonic() < deadline, log_path.read_text()
                time.sleep(0.02)
            rows = []
            with socket.create_connection(("127.0.0.1", port), timeout=10) as connection:
                with connection.makefile("rwb") as stream:
                    for index in range(FIRST_INDEX, FIRST_INDEX + STRANDS):
                        request = {"jsonrpc": "2.0", "id": index, "method": "path/find",
                                   "params": {"index": index}}
                        stream.write(json.dumps(request).encode() + b"\n")
                        stream.flush()
                        result = json.loads(stream.readline())["result"]
                        if result["found"]:
                            rows.append(f"{index}\t{result['sequence']}\n")
            return "".join(rows)
        finally:
            process.terminate()
            process.wait(timeout=10)


def main():
    assembler, rpc = sys.argv[1:3]
    with tempfile.TemporaryDirectory(prefix="dbgps-rpc-parity-") as directory:
        directory = Path(directory)
        reads = write_reads(directory)
        decoded = {}
        for number, options in enumerate(CASES):
            prefix = directory / f"case{number}"
            expected = batch(assembler, prefix, options, reads)
            observed = batch(rpc, prefix, options, reads)
            assert observed[0] == expected[0], (options, observed[0], expected[0])
            assert observed[1] == expected[1], (options, observed[1], expected[1])
            assert observed[2] == expected[2], (options, observed[2], expected[2])
            decoded[number] = expected[2].count("\n")
            # CRC32 and the legacy mapping do not match these reads: both must find nothing.
            assert decoded[number] > 0 or "crc32" in options or LEGACY_MAP in options, options

        options = CASES[0]
        prefix = directory / "stdout"
        assert batch(rpc, prefix, options, reads, False) == batch(assembler, prefix, options, reads, False)

        for options in SERVED_CASES:
            prefix = directory / "served"
            _, _, expected = batch(assembler, prefix, options, reads)
            assert expected
            assert served_strands(rpc, options, reads, directory / "server.log") == expected, options

        for payload, status in ((272, 0), (276, 1), (284, 1)):  # total strand must be < 300 bp
            for binary in (assembler, rpc):
                result = subprocess.run([binary, "-l", str(payload), *RANGE, reads[0]],
                                        capture_output=True, text=True, timeout=60)
                assert result.returncode == status, (binary, payload, result.stderr)

    print(f"rpc parity: {len(CASES)} batch and {len(SERVED_CASES)} served configurations match "
          f"({max(decoded.values())} strands at most)")


if __name__ == "__main__":
    main()
