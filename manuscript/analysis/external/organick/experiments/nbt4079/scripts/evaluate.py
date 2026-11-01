#!/usr/bin/env python3
"""Evaluate complete 110 nt outputs against all 607150 author references."""
from collections import Counter
import csv
import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def evaluate():
    references = {}
    with gzip.open(ROOT / "data/id20.refs.txt.gz", "rt") as handle:
        for line in handle:
            assert line.startswith("5'-") and line.rstrip().endswith("-3'")
            seq = line.strip()[3:-3]
            assert len(seq) == 150
            assert seq[:20] == "AGTGCAACAAGTCAATCCGT" and seq[-20:] == "AATTGAATGCTTGCTTGCCG"
            body = seq[20:130]
            assert body[:14] not in references
            references[body[:14]] = body
    n = len(references)
    results = []
    for path in sorted((ROOT / "results").glob("*.tsv")):
        if ".summary." in path.name or path.name == "summary.tsv":
            continue
        counts = Counter()
        found, wrong = set(), set()
        nohp_found, nohp_wrong = set(), set()
        nohp_unknown = 0
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["seed", "support", "payload"]:
                continue
            for row in reader:
                prefix, candidate = row["seed"][-14:], row["payload"]
                counts["returned_outputs"] += 1
                assert len(candidate) == 110 and candidate[:14] == prefix
                homopolymer = any(a == b for a, b in zip(candidate, candidate[1:]))
                if prefix not in references:
                    counts["outputs_unknown_address"] += 1
                    nohp_unknown += not homopolymer
                elif candidate == references[prefix]:
                    found.add(prefix)
                    if not homopolymer:
                        nohp_found.add(prefix)
                else:
                    wrong.add(prefix)
                    if not homopolymer:
                        nohp_wrong.add(prefix)
                counts["outputs_with_homopolymers"] += homopolymer
        assert not (found & wrong)
        missing = n - len(found) - len(wrong)
        result = {"run": path.stem, "designed_strands": n, **counts,
                  "exact": len(found), "wrong_known_address": len(wrong),
                  "missing": missing, "exact_rate": len(found)/n,
                  "error_erasure_fraction": (missing+2*len(wrong))/n,
                  "output_precision": len(found)/counts["returned_outputs"] if counts["returned_outputs"] else 0}
        results.append(result)
        print(json.dumps(result), flush=True)
        filtered = dict(result)
        filtered.update(run=path.stem+".nohp", exact=len(nohp_found),
                        wrong_known_address=len(nohp_wrong), outputs_unknown_address=nohp_unknown,
                        returned_outputs=len(nohp_found)+len(nohp_wrong)+nohp_unknown,
                        missing=n-len(nohp_found)-len(nohp_wrong), outputs_with_homopolymers=0)
        filtered["exact_rate"] = filtered["exact"]/n
        filtered["output_precision"] = filtered["exact"]/filtered["returned_outputs"] if filtered["returned_outputs"] else 0
        filtered["error_erasure_fraction"] = (filtered["missing"]+2*filtered["wrong_known_address"])/n
        results.append(filtered)
    fields = ["run", "designed_strands", "returned_outputs", "exact", "wrong_known_address", "missing",
              "outputs_unknown_address", "outputs_with_homopolymers", "exact_rate", "output_precision",
              "error_erasure_fraction"]
    with (ROOT / "results/summary.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    (ROOT / "results/summary.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    evaluate()
