#!/usr/bin/env python3
"""Figure 6: k-mer size and sequencing depth define the operating envelope."""
import os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as fs

fs.setup()
FULL = os.path.join(fs.RESULTS, "full")
KS = os.path.join(fs.RESULTS, "ksweep")
DEP = os.path.join(fs.RESULTS, "depth")
N = 210000
R_DEF = 2.55


def grid_for(s, k):
    f = f"{FULL}/{s}.grid.tsv" if k == 31 else f"{KS}/{s}_k{k}.grid.tsv"
    if not os.path.exists(f):
        return None
    g = pd.read_csv(f, sep="\t")
    return g if len(g) >= 90 else None


def stats_for(s, k):
    f = f"{FULL}/{s}.strandstats.tsv" if k == 31 else f"{KS}/{s}_k{k}.strandstats.tsv"
    return pd.read_csv(f, sep="\t") if os.path.exists(f) else None


def info_for(s, k):
    f = f"{FULL}/{s}.info.tsv" if k == 31 else f"{KS}/{s}_k{k}.info.tsv"
    if not os.path.exists(f):
        return None
    return dict(l.rstrip("\n").split("\t") for l in open(f) if "\t" in l)


KLIST = [17, 19, 21, 23, 25, 27, 31]
rows = []
for s in fs.SAMPLES:
    for k in KLIST:
        g = grid_for(s, k)
        if g is None:
            continue
        st = stats_for(s, k); inf = info_for(s, k)
        d1 = g[g.r == 0]; d2 = g[g.r == R_DEF]; d3 = g[g.r > 0]
        rows.append(dict(sample=s, k=k, dbgps=d1.correct.max(), dbgps_c=int(d1.loc[d1.correct.idxmax()].c),
                         dbgps2=d2.correct.max(), dbgps2_c=int(d2.loc[d2.correct.idxmax()].c),
                         dbgps2_any=d3.correct.max(),
                         wrong=int(g.wrong.sum()), wrong_dbgps2=int(d2.wrong.sum()),
                         wrong_r_inf=int(d1.wrong.sum()), wrong_r3=int(g[g.r == 3].wrong.sum()),
                         wrong_r4=int(g[g.r == 4].wrong.sum()),
                         settings_with_wrong_r_inf=int((d1.wrong > 0).sum()),
                         median_min_cov=float(st.min_cov.median()) if st is not None else np.nan,
                         intact=float((st.min_cov >= 2).mean()) if st is not None else np.nan,
                         distinct=float(inf["distinct_kmers"]) if inf else np.nan,
                         rss=float(inf.get("peak_rss_mb_final", np.nan)) if inf else np.nan))
kt = pd.DataFrame(rows)
kt.to_csv(os.path.join(fs.DATADIR, "ksize_sweep.csv"), index=False)
print(kt.to_string())

fig = plt.figure(figsize=(fs.DOUBLE, 112 * fs.MM))
gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.42, left=0.07, right=0.985, top=0.95, bottom=0.085)

# a: DBGPS2 recovery vs k
ax = fig.add_subplot(gs[0, 0])
for s in fs.SAMPLES:
    d = kt[kt["sample"] == s].sort_values("k")
    ax.plot(d.k, d.dbgps2 / N * 100, color=fs.COLOR[s], lw=1, marker="o", ms=2.4, mec="white", mew=0.3,
            label=fs.LABEL[s])
ax.set_xlabel("k-mer size k"); ax.set_ylabel("Strands recovered by DBGPS2 (%)")
ax.set_xticks(KLIST)
ax.legend(loc="lower center", ncol=3, columnspacing=0.8, fontsize=5.2)
ax.set_title("DBGPS2 (r = 2.55, best c)", fontsize=6.3, pad=2)
fs.panel(ax, "a")

# b: DBGPS recovery vs k
ax = fig.add_subplot(gs[0, 1])
for s in fs.SAMPLES:
    d = kt[kt["sample"] == s].sort_values("k")
    ax.plot(d.k, d.dbgps / N * 100, color=fs.COLOR[s], lw=1, marker="s", ms=2.4, mec="white", mew=0.3)
ax.set_xlabel("k-mer size k"); ax.set_ylabel("Strands recovered by DBGPS (%)")
ax.set_xticks(KLIST); ax.set_ylim(0, 102)
ax.set_title("DBGPS (r → ∞, best c)", fontsize=6.3, pad=2)
fs.panel(ax, "b")

# c: minimum path coverage vs k
ax = fig.add_subplot(gs[0, 2])
for s in fs.SAMPLES:
    d = kt[kt["sample"] == s].sort_values("k")
    ax.plot(d.k, d.median_min_cov, color=fs.COLOR[s], lw=1, marker="o", ms=2.4, mec="white", mew=0.3)
ax.set_yscale("log")
ax.set_xlabel("k-mer size k"); ax.set_ylabel("Median minimum true-path coverage")
ax.set_xticks(KLIST)
fs.panel(ax, "c")

# d-f: sequencing depth
dep = []
for f in glob.glob(f"{DEP}/*_p*.grid.tsv"):
    tag = os.path.basename(f)[:-9]
    s, kk, pp = tag.rsplit("_", 2)
    if os.path.getsize(f) == 0:
        continue
    g = pd.read_csv(f, sep="\t")
    if len(g) < 39:
        continue
    for r in (0, R_DEF, 3):
        d = g[g.r == r]
        dep.append(dict(sample=s, k=int(kk[1:]), pct=int(pp[1:]), r=r, best=d.correct.max()))
for s in fs.SAMPLES:  # 100% from the full-depth runs
    for k in (31, 21):
        g = grid_for(s, k)
        if g is None:
            continue
        for r in (0, R_DEF, 3):
            dep.append(dict(sample=s, k=k, pct=100, r=r, best=g[g.r == r].correct.max()))
dep = pd.DataFrame(dep)
dep.to_csv(os.path.join(fs.DATADIR, "depth_series.csv"), index=False)

for j, s in enumerate(("HT_con", "HT5")):
    ax = fig.add_subplot(gs[1, j])
    if len(dep):
        for (r, k), col, ls, lab in (((0, 31), fs.DBGPS_C, "-", "DBGPS, k = 31"),
                                     ((R_DEF, 31), fs.DBGPS2_C, "-", "DBGPS2, k = 31"),
                                     ((R_DEF, 21), fs.DBGPS2_C, (0, (2, 1)), "DBGPS2, k = 21")):
            d = dep[(dep["sample"] == s) & (dep.r == r) & (dep.k == k)].sort_values("pct")
            if len(d):
                ax.plot(d.pct * 2 * {"HT_con": 61.66, "HT5": 58.64}[s] / 100 * 1e6 / N, d.best / N * 100,
                        color=col, ls=ls, lw=1.1, marker="o", ms=2.2, mec="white", mew=0.3, label=lab)
    ax.set_xscale("log")
    ax.set_xlabel("Reads per designed strand")
    ax.set_ylabel("Strands recovered (%)"); ax.set_ylim(0, 102)
    ax.set_title(fs.LABEL[s], fontsize=6.3, pad=2)
    if j == 0:
        ax.legend(loc="lower right")
    fs.panel(ax, "de"[j])

# f: reads per strand needed for 95% recovery
ax = fig.add_subplot(gs[1, 2])
need = []
reads_m = {"HT_con": 61.66, "HT1": 57.51, "HT2": 58.97, "HT3": 66.22, "HT4": 78.24, "HT5": 58.64}
for s in fs.SAMPLES:
    for (r, k), lab in (((0, 31), "DBGPS k31"), ((R_DEF, 31), "DBGPS2 k31"), ((R_DEF, 21), "DBGPS2 k21")):
        d = dep[(dep["sample"] == s) & (dep.r == r) & (dep.k == k)].sort_values("pct") if len(dep) else []
        val = np.nan
        if len(d):
            x = d.pct.values * 2 * reads_m[s] * 1e6 / 100 / N; y = d.best.values / N
            ok = np.where(y >= 0.95)[0]
            if len(ok):
                i = ok[0]
                if i == 0:
                    val = x[0]
                else:  # log-linear interpolation between the bracketing depths
                    t = (0.95 - y[i - 1]) / (y[i] - y[i - 1])
                    val = float(np.exp(np.log(x[i - 1]) + t * (np.log(x[i]) - np.log(x[i - 1]))))
        need.append(dict(sample=fs.LABEL[s], method=lab, reads_per_strand_for_95pct=val))
need = pd.DataFrame(need)
need.to_csv(os.path.join(fs.DATADIR, "depth_needed_95.csv"), index=False)
x = np.arange(len(fs.SAMPLES)); w = 0.27
for i, (lab, col, hatch) in enumerate((("DBGPS k31", fs.DBGPS_C, None), ("DBGPS2 k31", fs.DBGPS2_C, None),
                                       ("DBGPS2 k21", "#86b6ef", None))):
    v = need[need.method == lab].reads_per_strand_for_95pct.values if len(need) else np.full(6, np.nan)
    ax.bar(x + (i - 1) * w, np.nan_to_num(v, nan=0), w, color=col, label=lab, edgecolor="white", lw=0.3)
    for xi, vi in zip(x + (i - 1) * w, v):
        if np.isnan(vi):
            ax.text(xi, 1.2, "n.r.", ha="center", va="bottom", fontsize=4.5, rotation=90, color=fs.INK2)
ax.set_yscale("log"); ax.set_ylim(1, 20000)
ax.set_xticks(x, [fs.LABEL[s] for s in fs.SAMPLES])
ax.set_ylabel("Reads per strand for 95% recovery")
ax.legend(loc="upper center", fontsize=5.1, ncol=3, columnspacing=0.6, handlelength=1)
ax.spines["left"].set_bounds(1, 1000); ax.yaxis.set_minor_locator(plt.NullLocator())
ax.set_yticks([1, 10, 100, 1000], ["1", "10", "100", "1,000"])
fs.panel(ax, "f")

fs.save(fig, "Fig6_k_depth")
print(need.pivot(index="sample", columns="method", values="reads_per_strand_for_95pct").round(1))
