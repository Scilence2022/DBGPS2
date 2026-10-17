#!/usr/bin/env python3
"""Figure 7: run time and memory of DBGPS (no Bloom filter) and DBGPS2 with its two Bloom filters.

DBGPS2 has two independent Bloom options: -B (query prefilter, for lookup speed) and -F (two-pass
singleton filter, for memory). DBGPS has neither; it is emulated with a very large ratio bound (-r).
All timings were measured one at a time with 3 threads (the default) on an idle machine.
"""
import os, re, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import figstyle as fs

fs.setup()
B = os.path.join(fs.RESULTS, "bloom")      # exact / -B / -F / -F -B at c = 1 (DBGPS at c = 10)
SP1 = os.path.join(fs.RESULTS, "speed3")   # DBGPS vs DBGPS2 at c = 1
SP2 = os.path.join(fs.RESULTS, "speed2")   # c = 10, 11, 12 replicates: -r and -B contributions
SAMP = ["HT_con", "HT3", "HT5"]
CONF = [("dbgps", "DBGPS (no Bloom filter)", fs.DBGPS_C),
        ("exact", "DBGPS2", "#2a78d6"),
        ("B12", "DBGPS2 -B 12 (query)", "#86b6ef"),
        ("F34", "DBGPS2 -F 34 (memory)", "#1baf7a"),
        ("F34B12", "DBGPS2 -F 34 -B 12", "#6250d6")]


def dec_times(path):
    return [float(x) for x in re.findall(r"Strand decoding time: ([\d.]+)", open(path).read())]


def parse(s, mode):
    t = open(f"{B}/{s}_{mode}.time").read()
    out = open(f"{B}/{s}_{mode}.stdout").read()
    real = float(re.search(r"([\d.]+) real", t).group(1))
    rss = int(re.search(r"(\d+)\s+maximum resident set size", t).group(1)) / 1024 ** 3
    dec = dec_times(f"{B}/{s}_{mode}.stdout")[-1]
    strands = set(open(f"{B}/{s}_{mode}.strands").read().split("\n")) - {""}
    return dict(sample=s, mode=mode, wall=real, rss_gb=rss, decode_s=dec, count_wall=real - dec,
                decoded=int(re.findall(r"Strands decoded: (\d+)", out)[-1]), strands=strands)


df = pd.DataFrame([parse(s, m) for s in SAMP for m, _, _ in CONF if os.path.exists(f"{B}/{s}_{m}.time")])
for s in SAMP:
    ref = df[(df["sample"] == s) & (df["mode"] == "exact")].strands.iloc[0]
    df.loc[df["sample"] == s, "identical_to_dbgps2_exact"] = df[df["sample"] == s].strands.apply(
        lambda x, r=ref: x == r)
df.drop(columns="strands").to_csv(os.path.join(fs.DATADIR, "resources_dbgps_dbgps2_bloom.csv"), index=False)

# DBGPS (r = 10^6) versus DBGPS2 (r = 2.55): path-search CPU time at c = 1 (one run) and c = 10 (mean of c = 10-12)
sp = []
for s in SAMP:
    for alg in ("dbgps", "dbgps2"):
        f = f"{SP1}/{s}_{alg}.stdout"
        if os.path.exists(f):
            out = open(f).read()
            reps = [dec_times(f)[-1]]
            if alg == "dbgps2":  # extra repeats of the control run (a first run was interrupted by a system sleep)
                reps += [dec_times(g)[-1] for g in sorted(glob.glob(f"{SP1}_ctrl_repeats/{s}_dbgps2_rep*.stdout"))]
            sp.append(dict(sample=s, c=1, mode=alg, decode_s=float(np.median(reps)), n_runs=len(reps),
                           decoded=int(re.findall(r"Strands decoded: (\d+)", out)[-1])))
        f = f"{SP2}/{s}_{alg}.stdout"
        if os.path.exists(f):
            out = open(f).read()
            sp.append(dict(sample=s, c=10, mode=alg, decode_s=float(np.mean(dec_times(f))),
                           decoded=int(re.findall(r"Strands decoded: (\d+)", out)[0])))
sp = pd.DataFrame(sp)
sp.to_csv(os.path.join(fs.DATADIR, "speed_dbgps_vs_dbgps2.csv"), index=False)

# contributions of -r and -B to the path-search time, relative to DBGPS (c = 10, 11, 12)
rep = []
for s in SAMP:
    t = {m: dec_times(f"{SP2}/{s}_{m}.stdout") for m in ("dbgps", "dbgps_B12", "dbgps2", "dbgps2_B12")
         if os.path.exists(f"{SP2}/{s}_{m}.stdout")}
    if len(t) == 4:
        for i, c in enumerate((10, 11, 12)):
            for m, lab in (("dbgps2", "r"), ("dbgps_B12", "B"), ("dbgps2_B12", "rB")):
                rep.append(dict(sample=s, c=c, option=lab, ratio=t[m][i] / t["dbgps"][i], time_s=t[m][i],
                                dbgps_s=t["dbgps"][i]))
rep = pd.DataFrame(rep)
rep.to_csv(os.path.join(fs.DATADIR, "speed_r_and_B_contributions.csv"), index=False)
print(df.drop(columns="strands").round(2).to_string()); print(sp.round(1).to_string()); print(rep.round(2).to_string())

fig = plt.figure(figsize=(fs.DOUBLE, 108 * fs.MM))
gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.34, left=0.075, right=0.985, top=0.88, bottom=0.085)
x = np.arange(len(SAMP)); w = 0.16
lab = [fs.LABEL[s] for s in SAMP]
off = lambda i: (i - 2) * w

# a: memory
ax = fig.add_subplot(gs[0, 0])
for i, (m, l, col) in enumerate(CONF):
    v = [df[(df["sample"] == s) & (df["mode"] == m)].rss_gb.max() for s in SAMP]
    ax.bar(x + off(i), v, w, color=col, label=l, edgecolor="white", lw=0.3)
ax.set_xticks(x, lab); ax.set_ylabel("Peak memory (GB)")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=5.0, columnspacing=0.8)
fs.panel(ax, "a", y=1.2)

# b: wall-clock time
ax = fig.add_subplot(gs[0, 1])
for i, (m, l, col) in enumerate(CONF):
    d = [df[(df["sample"] == s) & (df["mode"] == m)] for s in SAMP]
    cw = [dd.count_wall.max() for dd in d]; dw = [dd.decode_s.max() for dd in d]
    ax.bar(x + off(i), cw, w, color=col, edgecolor="white", lw=0.3)
    ax.bar(x + off(i), dw, w, bottom=cw, color=col, alpha=0.45, edgecolor="white", lw=0.3, hatch="////")
ax.set_xticks(x, lab); ax.set_ylabel("Wall-clock time (s)")
ax.legend(handles=[Patch(facecolor="#9a9a9a", label="k-mer counting (3 threads)"),
                   Patch(facecolor="#9a9a9a", alpha=0.45, hatch="////", label="Path search (1 thread)")],
          loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=1, fontsize=5.1)
fs.panel(ax, "b", y=1.2)

# c: path search, DBGPS versus DBGPS2
ax = fig.add_subplot(gs[1, 0]); w2 = 0.19
spec = [(1, "dbgps", fs.DBGPS_C, 0.45), (1, "dbgps2", fs.DBGPS2_C, 0.45),
        (10, "dbgps", fs.DBGPS_C, 1.0), (10, "dbgps2", fs.DBGPS2_C, 1.0)]
for i, (c, m, col, al) in enumerate(spec):
    v = [sp[(sp["sample"] == s) & (sp.c == c) & (sp["mode"] == m)].decode_s.max() for s in SAMP]
    ax.bar(x + (i - 1.5) * w2, v, w2, color=col, alpha=al, edgecolor="white", lw=0.3)
for j, s in enumerate(SAMP):
    for c, k0 in ((1, 0), (10, 2)):
        a = sp[(sp["sample"] == s) & (sp.c == c) & (sp["mode"] == "dbgps")].decode_s.max()
        b = sp[(sp["sample"] == s) & (sp.c == c) & (sp["mode"] == "dbgps2")].decode_s.max()
        ax.text(x[j] + (k0 + 0.5 - 1.5) * w2, max(a, b) * 1.15, f"{a / b:.1f}×", ha="center", fontsize=4.9,
                color=fs.INK2)
ax.set_yscale("log"); ax.set_ylim(5, 3000)
ax.set_xticks(x, lab); ax.set_ylabel("Path search time (s)")
ax.legend(handles=[Patch(facecolor=fs.DBGPS_C, label="DBGPS (r → ∞)"), Patch(facecolor=fs.DBGPS2_C, label="DBGPS2 (r = 2.55)"),
                   Patch(facecolor="#9a9a9a", alpha=0.45, label="c = 1"), Patch(facecolor="#9a9a9a", label="c = 10")],
          loc="upper right", ncol=2, fontsize=4.9, columnspacing=0.8)
fs.panel(ax, "c", y=1.06)

# d: how much of the speed-up comes from the ratio bound (-r) and how much from the query Bloom filter (-B)
ax = fig.add_subplot(gs[1, 1]); w3 = 0.24
opts = [("r", "-r 2.55", "#2a78d6"), ("B", "-B 12 (r → ∞)", "#86b6ef"), ("rB", "-r 2.55 and -B 12", "#184f95")]
for i, (o, l, col) in enumerate(opts):
    for j, s in enumerate(SAMP):
        d = rep[(rep["sample"] == s) & (rep.option == o)]
        if len(d):
            ax.bar(j + (i - 1) * w3, d.ratio.mean(), w3, color=col, edgecolor="white", lw=0.3,
                   label=l if j == 0 else None)
            ax.scatter(np.full(len(d), j + (i - 1) * w3) + np.linspace(-0.05, 0.05, len(d)), d.ratio, s=5,
                       color=fs.INK, zorder=3, lw=0)
ax.axhline(1, color=fs.INK2, lw=0.6, ls="--")
ax.set_xticks(x, lab); ax.set_ylim(0, 1.25)
ax.set_ylabel("Path search time relative to\nDBGPS (r → ∞, no Bloom filter)")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=3, fontsize=5.0, columnspacing=0.8)
fs.panel(ax, "d", y=1.2)

fs.save(fig, "Fig7_resources")
