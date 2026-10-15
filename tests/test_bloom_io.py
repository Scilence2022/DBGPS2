"""Both Bloom modes reject malformed/unreplayable input without hanging."""
import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile


with tempfile.TemporaryDirectory(prefix="dbgps-bloom-io-") as folder:
    folder = Path(folder)
    good = folder / "good.fa"
    good.write_text(">good\n" + "A" * 31 + "\n")
    malformed = folder / "bad.fq"
    malformed.write_text("@bad\n" + "C" * 40 + "\n+\nIII\n")
    truncated = folder / "truncated.fa.gz"
    truncated.write_bytes(gzip.compress(good.read_bytes())[:-8])
    fifo = folder / "reads.fifo"
    os.mkfifo(fifo)
    for binary in sys.argv[1:]:
        common = [binary, "-a", "1000", "-b", "1000"]
        for options in ([], ["-B", "12"], ["-F", "23"], ["-B", "12", "-F", "23"]):
            for bad in (malformed, truncated, folder / "missing.fa"):
                result = subprocess.run(common + options + [str(good), str(bad)],
                                        capture_output=True, text=True, timeout=10)
                assert result.returncode == 1 and "ERROR:" in result.stderr, (options, bad, result)
        result = subprocess.run(common + ["-F", "23", str(fifo)],
                                capture_output=True, text=True, timeout=3)
        assert result.returncode == 1 and "ERROR:" in result.stderr, result
        for options in ([], ["-B", "0", "-F", "0"]):
            result = subprocess.run(common + options + [str(good)], capture_output=True, text=True, timeout=10)
            assert result.returncode == 0 and "0 bitmap bytes" in result.stdout, result
        print(f"{Path(binary).name}: malformed FASTQ, truncated gzip, missing files, FIFO and explicit disable passed")
