#!/usr/bin/env python3
"""Supplementary Fig. 2: parameter sweep of the two DBGPS2 Bloom filters (DBGPS has none).

-B (query prefilter): bits per key. -F (two-pass singleton filter): log2 of the total bitmap bits.
DBGPS2, r = 2.55, c = 1, 3 threads, full depth, one run at a time.
"""
import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import figstyle as fs

fs.setup()
B = os.path.join(fs.RESULTS, "bloom")
SW = os.path.join(fs.RESULTS, "bloom_sweep")
SAMP = ["HT_con", "HT5"]


def parse(path_stub):
    t = open(path_stub + ".time").read(); out = open(path_stub + ".stdout").read()
    real = float(re.search(r"([\d.]+) real", t).group(1))
    rss = int(re.search(r"(\d+)\s+maximum resident set size", t).group(1)) / 1024 ** 3
    dec = float(re.findall(r"Strand decoding time: ([\d.]+)", out)[-1])
    strands = set(open(path_stub + ".strands").read().split("\n")) - {""}
    m = re.search(r"Bloom query filter \(-B\): \d+ bits/key; (\d+) bitmap bytes", out)
    return dict(wall=real, rss_gb=rss, decode_s=dec, count_wall=real - dec, bitmap_gib=int(m.group(1)) / 1024 ** 3 if m else 0.0,
                strands=strands)


rows = []
for s in SAMP:
    ref = parse(f"{B}/{s}_exact")
    rows.append(dict(sample=s, option="exact", value=0, **{k: v for k, v in ref.items() if k != "strands"},
                     identical=True))
    for mode, path in [("B12", f"{B}/{s}_B12"), ("F34", f"{B}/{s}_F34")] + \
                      [(m, f"{SW}/{s}_{m}") for m in ("B8", "B16", "B24", "F32", "F36")]:
        if os.path.exists(path + ".time"):
            r = parse(path)
            rows.append(dict(sample=s, option=mode[0], value=int(mode[1:]),
                             **{k: v for k, v in r.items() if k != "strands"}, identical=r["strands"] == ref["strands"]))
df = pd.DataFrame(rows)
df.to_csv(os.path.join(fs.DATADIR, "bloom_parameter_sweep.csv"), index=False)
print(df.round(2).to_string())

fig, axs = plt.subplots(1, 4, figsize=(fs.DOUBLE, 52 * fs.MM))
fig.subplots_adjust(left=0.065, right=0.99, top=0.84, bottom=0.2, wspace=0.55)
col = {"HT_con": fs.COLOR["HT_con"], "HT5": fs.COLOR["HT5"]}


def series(s, opt, col_name, rel=False):
    d = df[(df["sample"] == s) & (df.option == opt)].sort_values("value")
    ex = df[(df["sample"] == s) & (df.option == "exact")].iloc[0]
    y = d[col_name].values - (ex[col_name] if rel == "diff" else 0)
    if rel == "ratio":
        y = d[col_name].values / ex[col_name]
    return d.value.values, y


for s in SAMP:
    x, y = series(s, "B", "decode_s", "ratio")
    axs[0].plot(x, y, "-o", color=col[s], ms=3, lw=1, label=fs.LABEL[s])
    x, y = series(s, "B", "bitmap_gib")
    axs[1].plot(x, y, "-o", color=col[s], ms=3, lw=1)
    x, y = series(s, "F", "rss_gb", "ratio")
    axs[2].plot(x, y, "-o", color=col[s], ms=3, lw=1)
    x, y = series(s, "F", "count_wall", "ratio")
    axs[3].plot(x, y, "-o", color=col[s], ms=3, lw=1)
axs[0].axhline(1, color=fs.INK2, lw=0.6, ls="--")
axs[0].set_xlabel("Query filter, bits per key (-B)"); axs[0].set_ylabel("Path search time\nrelative to exact table")
axs[0].set_ylim(0.5, 1.1)
axs[1].set_xlabel("Query filter, bits per key (-B)"); axs[1].set_ylabel("Bloom bitmap size (GiB)")
axs[2].axhline(1, color=fs.INK2, lw=0.6, ls="--")
axs[2].set_xlabel("Singleton filter, log$_2$ bits (-F)"); axs[2].set_ylabel("Peak memory\nrelative to exact table")
axs[3].axhline(1, color=fs.INK2, lw=0.6, ls="--")
axs[3].set_xlabel("Singleton filter, log$_2$ bits (-F)"); axs[3].set_ylabel("Counting time\nrelative to exact table")
axs[2].set_xticks([32, 34, 36]); axs[3].set_xticks([32, 34, 36])
axs[0].legend(loc="lower left")
for a, l in zip(axs, "abcd"):
    fs.panel(a, l, x=-0.3, y=1.08)
fs.save(fig, "FigS5_bloom")
