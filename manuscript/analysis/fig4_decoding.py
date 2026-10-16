#!/usr/bin/env python3
"""Figure 4: strand recovery of DBGPS (ratio constraint disabled) vs DBGPS2 across aging."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import figstyle as fs

fs.setup()
FULL = os.path.join(fs.RESULTS, "full")
N = 210000
R_DEF = 2.55
grids = {s: pd.read_csv(f"{FULL}/{s}.grid.tsv", sep="\t") for s in fs.SAMPLES}
allg = pd.concat([g.assign(sample=fs.LABEL[s]) for s, g in grids.items()])
allg.to_csv(os.path.join(fs.DATADIR, "decoding_grid_full_depth.csv"), index=False)

best = []
for s, g in grids.items():
    d1 = g[g.r == 0]; d2 = g[g.r == R_DEF]
    b1 = d1.loc[d1.correct.idxmax()]; b2 = d2.loc[d2.correct.idxmax()]
    ball = g[g.r > 0].loc[g[g.r > 0].correct.idxmax()]
    best.append(dict(sample=fs.LABEL[s], dbgps_best_c=int(b1.c), dbgps_best=int(b1.correct),
                     dbgps2_best_c=int(b2.c), dbgps2_best=int(b2.correct),
                     dbgps2_c1=int(d2[d2.c == 1].correct.iloc[0]),
                     dbgps2_any_best=int(ball.correct), dbgps2_any_r=ball.r, dbgps2_any_c=int(ball.c)))
best = pd.DataFrame(best)
best.to_csv(os.path.join(fs.DATADIR, "best_recovery.csv"), index=False)
print(best.to_string())

fig = plt.figure(figsize=(fs.DOUBLE, 128 * fs.MM))
outer = fig.add_gridspec(2, 1, height_ratios=[0.8, 1.15], hspace=0.5, left=0.065, right=0.985,
                         top=0.955, bottom=0.075)
top = outer[0].subgridspec(1, 6, wspace=0.12)
bot = outer[1].subgridspec(1, 3, wspace=0.45, width_ratios=[1, 1.05, 1.15])

# a: recovery vs c for each sample
for i, s in enumerate(fs.SAMPLES):
    ax = fig.add_subplot(top[0, i])
    g = grids[s]
    for r, col, lab in ((0, fs.DBGPS_C, "DBGPS (r → ∞)"), (R_DEF, fs.DBGPS2_C, "DBGPS2 (r = 2.55)")):
        d = g[g.r == r].sort_values("c")
        ax.plot(d.c, d.correct / N * 100, color=col, lw=1.1, marker="o", ms=2, mec="white", mew=0.3, label=lab)
    ax.set_xscale("log")
    ax.set_ylim(0, 102)
    ax.set_xticks([1, 10, 100], ["1", "10", "100"])
    ax.set_title(fs.LABEL[s], fontsize=6.3, pad=2)
    ax.axhline(100, color=fs.GRID, lw=0.6, zorder=0)
    if i == 0:
        ax.set_ylabel("Strands recovered (%)")
        fs.panel(ax, "a", x=-0.42, y=1.04)
    else:
        ax.set_yticklabels([])
    if i == 2:
        ax.set_xlabel("Coverage cutoff c", x=1.05)
    if i == 0:
        ax.legend(loc="lower right", fontsize=5.2, handlelength=1.2)

# b: best recovery per sample
ax = fig.add_subplot(bot[0, 0])
x = np.arange(len(fs.SAMPLES)); w = 0.38
ax.bar(x - w / 2, best.dbgps_best / N * 100, w, color=fs.DBGPS_C, label="DBGPS, best c")
ax.bar(x + w / 2, best.dbgps2_best / N * 100, w, color=fs.DBGPS2_C, label="DBGPS2, best c")
for i in range(len(x)):
    ax.text(x[i] - w / 2, best.dbgps_best[i] / N * 100 + 1.2, f"{best.dbgps_best[i] / N * 100:.1f}",
            ha="center", fontsize=4.6, rotation=90, va="bottom", color=fs.INK2)
    ax.text(x[i] + w / 2, best.dbgps2_best[i] / N * 100 + 1.2, f"{best.dbgps2_best[i] / N * 100:.1f}",
            ha="center", fontsize=4.6, rotation=90, va="bottom", color=fs.INK2)
ax.set_xticks(x, best["sample"])
ax.set_ylabel("Strands recovered (%)")
ax.set_ylim(0, 125); ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.legend(loc="upper right", ncol=1)
fs.panel(ax, "b")

# c: r x c landscape for the most aged sample
ax = fig.add_subplot(bot[0, 1])
s = "HT5"
g = grids[s]
rs = sorted(g.r.unique(), key=lambda v: (v == 0, v))
cs = sorted(g.c.unique())
mat = np.array([[g[(g.r == r) & (g.c == c)].correct.iloc[0] / N * 100 for c in cs] for r in rs])
cmap = LinearSegmentedColormap.from_list("blu", ["#f5f9fe", "#9ec5f4", "#3987e5", "#184f95", "#0d366b"])
im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=100, origin="lower")
ax.set_xticks(range(len(cs)), [str(c) for c in cs], rotation=90, fontsize=5)
ax.set_yticks(range(len(rs)), [("∞ (DBGPS)" if r == 0 else f"{r:g}") for r in rs], fontsize=5.5)
ax.set_xlabel("Coverage cutoff c"); ax.set_ylabel("Ratio bound r")
bi = np.unravel_index(np.argmax(mat), mat.shape)
ax.scatter(bi[1], bi[0], marker="*", s=22, color="white", edgecolor=fs.INK, linewidth=0.4)
ax.spines[:].set_visible(False); ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
cb.set_label("Strands recovered (%)", fontsize=5.8); cb.ax.tick_params(labelsize=5.3, length=1.5)
cb.outline.set_linewidth(0.4)
ax.set_title(f"{fs.LABEL[s]}: recovery landscape", fontsize=6.3, pad=3)
fs.panel(ax, "c", x=-0.3)

# d: outcome composition at each method's best setting
ax = fig.add_subplot(bot[0, 2])
cats = [("correct", "Recovered", "#2a78d6"), ("no_path", "Path broken", "#b9b8b3"),
        ("overflow", "Path overflow", "#eb6834"), ("no_valid", "No CRC-valid path", "#eda100"),
        ("ambiguous", "Ambiguous", "#6250d6"), ("wrong", "Wrong strand", "#e34948")]
xs, labs = [], []
for i, s in enumerate(fs.SAMPLES):
    g = grids[s]
    for j, (r, tag) in enumerate(((0, "D1"), (R_DEF, "D2"))):
        d = g[g.r == r]
        rowb = d.loc[d.correct.idxmax()]
        xpos = i * 2.6 + j * 1.05
        bottom = 0
        for k, lab, col in cats:
            v = rowb[k] / N * 100
            ax.bar(xpos, v, 0.95, bottom=bottom, color=col, edgecolor="white", lw=0.3,
                   label=lab if (i == 0 and j == 0) else None)
            bottom += v
        xs.append(xpos); labs.append("D" if r == 0 else "D2")
ax.set_xticks(xs, ["1" if l == "D" else "2" for l in labs], fontsize=5)
for i, s in enumerate(fs.SAMPLES):
    ax.text(i * 2.6 + 0.52, -14, fs.LABEL[s], ha="center", fontsize=5.5, color=fs.INK)
ax.text(-1.2, -6.2, "DBGPS", ha="right", fontsize=5, color=fs.INK2)
ax.set_ylim(0, 100)
ax.set_ylabel("Indexed strands (%)")
ax.legend(loc="lower center", bbox_to_anchor=(0.47, 1.0), ncol=3, fontsize=4.9, handlelength=0.9, columnspacing=0.6, handletextpad=0.4)
fs.panel(ax, "d", x=-0.2, y=1.13)

fs.save(fig, "Fig4_decoding")
