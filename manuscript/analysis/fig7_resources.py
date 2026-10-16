#!/usr/bin/env python3
"""Figure 7: run time and memory of the DBGPS2 executable, with and without Bloom filters."""
import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import figstyle as fs

fs.setup()
B = os.path.join(fs.RESULTS, "bloom")
SAMP = ["HT_con", "HT3", "HT5"]
MODES = [("exact", "Exact table"), ("B12", "Query Bloom (-B 12)"), ("F34", "Singleton Bloom (-F 34)"),
         ("F34B12", "Both (-F 34 -B 12)")]
MCOL = {"exact": "#2a78d6", "B12": "#86b6ef", "F34": "#1baf7a", "F34B12": "#6250d6", "dbgps": fs.DBGPS_C}


def parse(s, mode):
    t = open(f"{B}/{s}_{mode}.time").read()
    out = open(f"{B}/{s}_{mode}.stdout").read()
    real = float(re.search(r"([\d.]+) real", t).group(1))
    rss = int(re.search(r"(\d+)\s+maximum resident set size", t).group(1)) / 1024 ** 3
    dec = float(re.findall(r"Strand decoding time: ([\d.]+)", out)[-1])
    n = int(re.findall(r"Strands decoded: (\d+)", out)[-1])
    strands = set(open(f"{B}/{s}_{mode}.strands").read().split("\n")) - {""}
    return dict(sample=s, mode=mode, wall=real, rss_gb=rss, decode_s=dec, count_wall=real - dec, decoded=n,
                strands=strands)


rows = [parse(s, m) for s in SAMP for m, _ in MODES + [("dbgps", "")]
        if os.path.exists(f"{B}/{s}_{m}.time")]
df = pd.DataFrame(rows)
for s in SAMP:
    ref = df[(df["sample"] == s) & (df["mode"] == "exact")].strands
    if len(ref):
        ref = ref.iloc[0]
        df.loc[df["sample"] == s, "identical_to_exact"] = df[df["sample"] == s].strands.apply(lambda x: x == ref)
df.drop(columns="strands").to_csv(os.path.join(fs.DATADIR, "resources_dbgps2_executable.csv"), index=False)
print(df.drop(columns="strands").round(2).to_string())

fig = plt.figure(figsize=(fs.DOUBLE, 72 * fs.MM))
gs = fig.add_gridspec(1, 3, wspace=0.42, left=0.065, right=0.985, top=0.76, bottom=0.12)
x = np.arange(len(SAMP)); w = 0.2

ax = fig.add_subplot(gs[0, 0])
for i, (m, lab) in enumerate(MODES):
    v = [df[(df["sample"] == s) & (df["mode"] == m)].rss_gb.max() for s in SAMP]
    ax.bar(x + (i - 1.5) * w, v, w, color=MCOL[m], label=lab, edgecolor="white", lw=0.3)
ax.set_xticks(x, [fs.LABEL[s] for s in SAMP]); ax.set_ylabel("Peak memory (GB)")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=5.1, columnspacing=0.8)
fs.panel(ax, "a", y=1.22)

ax = fig.add_subplot(gs[0, 1])
for i, (m, lab) in enumerate(MODES):
    d = [df[(df["sample"] == s) & (df["mode"] == m)] for s in SAMP]
    cw = [dd.count_wall.max() for dd in d]; dw = [dd.decode_s.max() for dd in d]
    ax.bar(x + (i - 1.5) * w, cw, w, color=MCOL[m], edgecolor="white", lw=0.3)
    ax.bar(x + (i - 1.5) * w, dw, w, bottom=cw, color=MCOL[m], alpha=0.45, edgecolor="white", lw=0.3, hatch="////")
ax.set_xticks(x, [fs.LABEL[s] for s in SAMP]); ax.set_ylabel("Wall-clock time (s)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor="#9a9a9a", label="k-mer counting (16 threads)"),
                   Patch(facecolor="#9a9a9a", alpha=0.45, hatch="////", label="Path search (1 thread)")],
          loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=1, fontsize=5.1)
fs.panel(ax, "b", y=1.22)

ax = fig.add_subplot(gs[0, 2])
w2 = 0.36
v1 = [df[(df["sample"] == s) & (df["mode"] == "dbgps")].decode_s.max() for s in SAMP]
v2 = [df[(df["sample"] == s) & (df["mode"] == "exact")].decode_s.max() for s in SAMP]
ax.bar(x - w2 / 2, v1, w2, color=fs.DBGPS_C, label="DBGPS (r → ∞, c = 10)")
ax.bar(x + w2 / 2, v2, w2, color=fs.DBGPS2_C, label="DBGPS2 (r = 2.55, c = 1)")
for i in range(len(SAMP)):
    if v1[i] and v2[i]:
        ax.text(x[i], max(v1[i], v2[i]) * 1.04, f"{v1[i] / v2[i]:.1f}×", ha="center", fontsize=5.2, color=fs.INK2)
ax.set_xticks(x, [fs.LABEL[s] for s in SAMP]); ax.set_ylabel("Path search time, 285,583 indices (s)")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=1, fontsize=5.1)
fs.panel(ax, "c", y=1.22)

fs.save(fig, "Fig7_resources")
