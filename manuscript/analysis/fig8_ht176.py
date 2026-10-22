#!/usr/bin/env python3
"""Figure 8: DBGPS2 on a sample aged for 176 days at 70 degrees C (HT176)."""
import os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import figstyle as fs

fs.setup()
R = fs.RESULTS
FULL, KS, DEP, QC = (os.path.join(R, d) for d in ("full", "ksweep", "depth", "qc"))
N = 210000
S = "HT176"
READS = {"HT_con": 61656962, "HT5": 58635539, "HT176": 60642860}
R_DEF = 2.55


def grid(k):
    f = f"{FULL}/{S}.grid.tsv" if k == 31 else f"{KS}/{S}_k{k}.grid.tsv"
    return pd.read_csv(f, sep="\t") if os.path.exists(f) and os.path.getsize(f) else None


ks_list = [k for k in (17, 19, 21, 23, 25, 27, 31) if grid(k) is not None and len(grid(k)) >= 90]
fig = plt.figure(figsize=(fs.DOUBLE, 118 * fs.MM))
gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.42, left=0.07, right=0.985, top=0.94, bottom=0.085)
CMP = [("HT_con", fs.COLOR["HT_con"], "Control"), ("HT5", fs.COLOR["HT5"], "HT5"), (S, fs.COLOR[S], "HT176")]

# a: fragment length
ax = fig.add_subplot(gs[0, 0])
for s, col, lab in CMP:
    h = np.loadtxt(f"{QC}/{s}.frag.tsv", skiprows=1)[:, 1]
    h = np.r_[h, np.zeros(11)]
    ax.step(np.arange(len(h)), np.cumsum(h) / h.sum() * 100, where="post", color=col, lw=1.1, label=lab)
ax.set_xlim(0, 210); ax.set_ylim(0, 102); ax.set_xticks([0, 50, 100, 150, 200])
ax.set_xlabel("Fragment length (nt)"); ax.set_ylabel("Aligned pairs ≤ length (%)")
ax.legend(loc="upper left", handlelength=1.2)
fs.panel(ax, "a")

# b: minimum true-path coverage
ax = fig.add_subplot(gs[0, 1])
for s, col, lab in CMP:
    st = pd.read_csv(f"{FULL}/{s}.strandstats.tsv", sep="\t")
    v = np.sort(np.maximum(st.min_cov.values, 0.5))
    ax.plot(v, np.arange(1, len(v) + 1) / len(v) * 100, color=col, lw=1.1)
ax.set_xscale("log"); ax.set_xlim(0.5, 2000)
ax.set_xlabel("Minimum k-mer coverage along strand path"); ax.set_ylabel("Strands (cumulative %)")
fs.panel(ax, "b")

# c: recovery versus cutoff
ax = fig.add_subplot(gs[0, 2])
g31 = grid(31)
for (g, r, col, ls, lab) in [(g31, 0, fs.DBGPS_C, "-", "DBGPS, k = 31"), (g31, R_DEF, fs.DBGPS2_C, "-", "DBGPS2, k = 31")] + \
        ([(grid(21), R_DEF, fs.DBGPS2_C, (0, (2, 1)), "DBGPS2, k = 21")] if 21 in ks_list else []):
    d = g[g.r == r].sort_values("c")
    ax.plot(d.c, d.correct / N * 100, color=col, ls=ls, lw=1.1, marker="o", ms=2, mec="white", mew=0.3, label=lab)
ax.set_xscale("log"); ax.set_ylim(0, 102)
ax.set_xlabel("Coverage cutoff c"); ax.set_ylabel("Strands recovered (%)")
ax.legend(loc="lower right", fontsize=5.2)
fs.panel(ax, "c")

# d: strands that are not recovered, by cause
ax = fig.add_subplot(gs[1, 0])
cats = [("no_path", "Path broken", "#b9b8b3"), ("overflow", "Path overflow", "#eb6834"),
        ("no_valid", "No CRC-valid path", "#eda100"), ("ambiguous", "Ambiguous", "#6250d6"), ("wrong", "Wrong strand", "#e34948")]
bars = [("DBGPS\nk = 31", g31[g31.r == 0]), ("DBGPS2\nk = 31", g31[g31.r == R_DEF])] + \
       ([("DBGPS2\nk = 21", grid(21)[grid(21).r == R_DEF])] if 21 in ks_list else [])
for i, (lab, d) in enumerate(bars):
    row = d.loc[d.correct.idxmax()]
    bottom = 0
    for k, l, col in cats:
        v = row[k] / N * 100
        ax.bar(i, v, 0.7, bottom=bottom, color=col, edgecolor="white", lw=0.3, label=l if i == 0 else None)
        bottom += v
    ax.text(i, bottom + 0.05, f"{bottom:.2f}%", ha="center", fontsize=5.2, color=fs.INK2)
ax.set_xticks(range(len(bars)), [b[0] for b in bars]); ax.set_ylim(0, 3.2)
ax.set_ylabel("Strands not recovered (%)")
ax.legend(loc="upper right", fontsize=4.9, handlelength=0.9, labelspacing=0.3)
fs.panel(ax, "d")

# e: sequencing depth
ax = fig.add_subplot(gs[1, 1])
rows = []
for f in glob.glob(f"{DEP}/{S}_k*_p*.grid.tsv"):
    if os.path.getsize(f) == 0:
        continue
    tag = os.path.basename(f)[:-9]; _, kk, pp = tag.rsplit("_", 2)
    g = pd.read_csv(f, sep="\t")
    if len(g) < 39:
        continue
    for r in (0, R_DEF):
        rows.append(dict(k=int(kk[1:]), pct=int(pp[1:]), r=r, best=g[g.r == r].correct.max()))
for k in (31, 21):
    for r in (0, R_DEF):
        g = grid(k)
        if g is not None and len(g) >= 90:
            rows.append(dict(k=k, pct=100, r=r, best=g[g.r == r].correct.max()))
dep = pd.DataFrame(rows)
dep.to_csv(os.path.join(fs.DATADIR, "ht176_depth_series.csv"), index=False)
for (r, k), col, ls, lab in [((0, 31), fs.DBGPS_C, "-", "DBGPS, k = 31"), ((R_DEF, 31), fs.DBGPS2_C, "-", "DBGPS2, k = 31"),
                              ((R_DEF, 21), fs.DBGPS2_C, (0, (2, 1)), "DBGPS2, k = 21")]:
    d = dep[(dep.r == r) & (dep.k == k)].sort_values("pct") if len(dep) else []
    if len(d):
        ax.plot(d.pct * 2 * READS[S] / 100 / N, d.best / N * 100, color=col, ls=ls, lw=1.1, marker="o", ms=2.2,
                mec="white", mew=0.3, label=lab)
ax.axhline(95, color=fs.GRID, lw=0.7, zorder=0)
ax.set_xscale("log"); ax.set_ylim(0, 102)
ax.set_xlabel("Reads per designed strand"); ax.set_ylabel("Strands recovered (%)")
fs.panel(ax, "e")

# f: k-mer size
ax = fig.add_subplot(gs[1, 2])
kt = []
for k in ks_list:
    g = grid(k)
    kt.append(dict(k=k, dbgps=g[g.r == 0].correct.max(), dbgps2=g[g.r == R_DEF].correct.max(),
                   dbgps2_any=g[g.r > 0].correct.max(), wrong_dbgps2=int(g[g.r == R_DEF].wrong.sum()),
                   wrong_dbgps=int(g[g.r == 0].wrong.sum())))
kt = pd.DataFrame(kt)
kt.to_csv(os.path.join(fs.DATADIR, "ht176_ksize.csv"), index=False)
if len(kt):
    ax.plot(kt.k, kt.dbgps / N * 100, color=fs.DBGPS_C, marker="s", ms=2.6, mec="white", mew=0.3, lw=1.1, label="DBGPS")
    ax.plot(kt.k, kt.dbgps2 / N * 100, color=fs.DBGPS2_C, marker="o", ms=2.6, mec="white", mew=0.3, lw=1.1, label="DBGPS2")
    ax.set_xticks(kt.k)
ax.set_ylim(0, 102)
ax.set_xlabel("k-mer size k"); ax.set_ylabel("Best strand recovery (%)")
ax.legend(loc="center right", fontsize=5.2)
fs.panel(ax, "f")

fs.save(fig, "Fig8_HT176")
print(kt.round(0).to_string() if len(kt) else "no k-sweep yet")
