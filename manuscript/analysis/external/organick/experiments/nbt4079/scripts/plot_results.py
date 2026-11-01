#!/usr/bin/env python3
"""Standalone scientific comparison figure; run with matplotlib installed."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
results = {row["run"]: row for row in json.loads((ROOT / "results/summary.json").read_text())}
pre = json.loads((ROOT / "results/preprocessing.json").read_text())
samples = ["5x_raw", "10x_raw", "full"]
configs = ["c0_r2.55", "c1_r2.55", "modal"]
labels = ["DBGPS2: count >=1, ratio 2.55", "DBGPS2: count >=2, ratio 2.55",
          "Modal observed read (baseline)"]
colors = ["#2878B5", "#179C7D", "#7A7E87"]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.7), layout="constrained")
x = np.arange(3)
for j, (cfg, label, color) in enumerate(zip(configs, labels, colors)):
    rows = [results[sample+"."+cfg] for sample in samples]
    axes[0].bar(x + (j-1)*0.25, [row["exact_rate"]*100 for row in rows], width=0.24,
                color=color, label=label)
    axes[1].bar(x + (j-1)*0.25,
                [row["wrong_known_address"] + row.get("outputs_unknown_address", 0) for row in rows],
                width=0.24, color=color)
for ax in axes:
    ax.set_xticks(x, [f'{pre["samples"][sample]["mean_raw_coverage"]:.2f}x' for sample in samples])
    ax.set_xlabel("Mean raw FASTQ records / designed strand")
    ax.spines[["right", "top"]].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#E7E7E7")
axes[0].set_ylabel("Exact 110 nt reconstruction (% of 607,150)")
axes[0].set_ylim(0, 102)
axes[1].set_ylabel("Incorrect output count (log scale)")
axes[1].set_yscale("log")
axes[1].set_ylim(50, 2e6)
axes[0].set_title("Exact strand recovery")
axes[1].set_title("Errors accepted without a strand checksum")
fig.suptitle("Organick et al. 2018 - File 16, author reads and references", fontsize=14)
handles, legend_labels = axes[0].get_legend_handles_labels()
fig.legend(handles, legend_labels, loc="outside lower center", ncol=2, frameon=False, fontsize=9)
fig.savefig(ROOT / "results/recovery_comparison.png", dpi=200)
fig.savefig(ROOT / "results/recovery_comparison.svg")
