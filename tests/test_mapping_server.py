"""Exercise the running JSON-RPC server and client with independent codebooks."""

import json
from pathlib import Path
import random
import socket
import subprocess
import sys
import tempfile
import time

from test_mapping_cli import encode, word, PRIMERS


def main():
    binary, client = (str(Path(arg).resolve()) for arg in sys.argv[1:])
    repo = Path(__file__).resolve().parent.parent
    rng = random.Random(193)
    cases = [(None, 1, 2, 16, 80, 1000),
             ("legacy", 2, 4, 16, 80, 1000),
             ("custom", 3, 5, 24, 120, (1 << 35) + 123),
             ("custom", 5, 8, 15, 100, (1 << 22) + 31),
             (None, 1, 2, 28, 80, (1 << 54) + 123)]
    with tempfile.TemporaryDirectory(prefix="dbgps-rpc-") as directory:
        directory = Path(directory)
        for case, (kind, bases, bits, index_bases, payload_bases, index) in enumerate(cases):
            words = [word(i, bases) for i in range(2 ** bits)]
            mapping_args = []
            if kind == "legacy":
                mapping = repo / "examples/legacy-2base.map"
                records = [line.split() for line in mapping.read_text().splitlines()
                           if line and line[0] in "ACGT"]
                words = [dna for dna, value in sorted(records, key=lambda row: int(row[1], 2))]
                mapping_args = ["-M", str(mapping)]
            elif kind:
                words = [word(i, bases) for i in rng.sample(range(4 ** bases), 2 ** bits)]
                mapping = directory / f"case{case}.map"
                mapping.write_text(f"bases {bases}\nbits {bits}\n" + "".join(
                    f"{dna} {value:0{bits}b}\n" for value, dna in enumerate(words)))
                mapping_args = ["-M", str(mapping)]
            index_bits = index_bases // bases * bits
            raw = index.to_bytes(index_bits // 8, "big") + bytes(
                rng.randrange(256) for _ in range(payload_bases // bases * bits // 8))
            expected = encode(raw, words, bits)
            reads = directory / f"case{case}.fa"
            reads.write_text("".join(f">{copy}\n{PRIMERS[0]}{expected}{PRIMERS[1]}\n"
                                     for copy in range(10)))
            with socket.socket() as reservation:
                reservation.bind(("127.0.0.1", 0))
                port = reservation.getsockname()[1]
            log_path = directory / f"server{case}.log"
            with log_path.open("w") as log:
                server_args = (["-S", "-P", str(port)] if case % 3 == 0 else
                               ["-S", str(port)] if case % 3 == 1 else [f"-S{port}"])
                proc = subprocess.Popen([binary, *server_args, "-i", str(index_bases),
                                         "-l", str(payload_bases), *mapping_args, str(reads)],
                                        stdout=log, stderr=log)
                try:
                    deadline = time.monotonic() + 10
                    while "Listening on port" not in log_path.read_text():
                        assert proc.poll() is None, log_path.read_text()
                        assert time.monotonic() < deadline, log_path.read_text()
                        time.sleep(0.02)
                    with socket.create_connection(("127.0.0.1", port), timeout=5) as connection:
                        with connection.makefile("rwb") as stream:
                            request_id = 0

                            def rpc(method, params=None):
                                nonlocal request_id
                                request_id += 1
                                request = {"jsonrpc": "2.0", "id": request_id,
                                           "method": method, "params": params or {}}
                                stream.write(json.dumps(request).encode() + b"\n")
                                stream.flush()
                                response = json.loads(stream.readline())
                                assert response["id"] == request_id, response
                                return response

                            info = rpc("graph/info")["result"]
                            assert info["index_byte_num"] == index_bits // 8
                            assert info["mapping"]["custom"] == bool(kind)
                            assert (info["mapping"]["bases"], info["mapping"]["bits"]) == (bases, bits)
                            for value in (str(index), index) if index < 2 ** 53 else (str(index),):
                                result = rpc("path/find", {"index": value})["result"]
                                assert result == {"index": value, "found": True, "sequence": expected}, result
                            root = (PRIMERS[0] + expected[:index_bases])[-31:]
                            result = rpc("path/findKmer", {"kmer": root})["result"]
                            assert result["found"] and result["sequence"] == expected[index_bases:], result
                            assert rpc("kmer/cov", {"kmer": root})["result"]["coverage"] > 0
                            for invalid in ({"indexLen": 32}, {"indexLen": -1},
                                            {"pathLen": payload_bases - 1}, {"pathLen": "24"},
                                            {"pathLen": 2 ** 33}, {"index": str(1 << index_bits)},
                                            {"index": "18446744073709551616"}, {"index": "1x"},
                                            {"index": -1}, {"index": 1.5}, {"index": True},
                                            {"index": 2 ** 53}, {"maxCovRatio": 0},
                                            {"skipRatios": -1}):
                                response = rpc("path/find", {"index": str(index), **invalid})
                                assert response["error"]["code"] == -32602, (invalid, response)
                            response = rpc("path/findKmer", {"kmer": root, "pathLen": payload_bases - 1})
                            assert response["error"]["code"] == -32602
                            result = subprocess.run([client, "-p", str(port), "-c", f"find {index}"],
                                                    capture_output=True, text=True, timeout=10)
                            assert result.returncode == 0, result.stderr
                            reply = json.loads(result.stdout)["result"]
                            assert reply["index"] == str(index) and reply["sequence"] == expected, reply
                            assert rpc("server/shutdown")["result"]["status"] == "shutting down"
                    # Wake the existing blocking accept loop after the shutdown reply.
                    try:
                        with socket.create_connection(("127.0.0.1", port), timeout=1):
                            pass
                    except ConnectionRefusedError:
                        pass
                    assert proc.wait(timeout=5) == 0, log_path.read_text()
                finally:
                    if proc.poll() is None:
                        proc.kill()
                        proc.wait(timeout=5)
            print(f"mapping server: {kind or 'canonical'}, {bases} bases/{bits} bits, exact index={index} passed")


if __name__ == "__main__":
    main()
