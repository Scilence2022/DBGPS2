#!/usr/bin/env python3
"""Figure 5: why the ratio constraint works - search cost, rescued strands, recovery ceiling."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as fs

fs.setup()
FULL = os.path.join(fs.RESULTS, "full")
N = 210000
R_DEF = 2.55
grids, stats, outc = {}, {}, {}
for s in fs.SAMPLES:
    grids[s] = pd.read_csv(f"{FULL}/{s}.grid.tsv", sep="\t")
    stats[s] = pd.read_csv(f"{FULL}/{s}.strandstats.tsv", sep="\t")
    o = pd.read_csv(f"{FULL}/{s}.outcomes.tsv", sep="\t", dtype={"id": np.int64}, header=None, skiprows=1,
                    names=["id", "o"])
    hdr = open(f"{FULL}/{s}.outcomes.tsv").readline().rstrip("\n").split("\t")[1:]
    outc[s] = (o, hdr)


def outcome_col(s, c, r):
    o, hdr = outc[s]
    key = f"c{int(c)}_r{r:g}"
    j = hdr.index(key)
    return o.o.str[j].values


fig = plt.figure(figsize=(fs.DOUBLE, 115 * fs.MM))
gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.42, left=0.07, right=0.985, top=0.95, bottom=0.085)

# a: search cost
ax = fig.add_subplot(gs[0, 0])
for s, mk in (("HT_con", "o"), ("HT5", "s")):
    g = grids[s]
    for r, col in ((0, fs.DBGPS_C), (R_DEF, fs.DBGPS2_C)):
        d = g[g.r == r].sort_values("c")
        ax.plot(d.c, d.mean_nodes_ref, color=col, lw=1, marker=mk, ms=2.4, mec="white", mew=0.3)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Coverage cutoff c"); ax.set_ylabel("Graph nodes visited per index")
ax.axhline(150, color=fs.GRID, lw=0.6, zorder=0)
ax.text(1.0, 163, "single path (150 nodes)", ha="left", va="bottom", fontsize=4.8, color=fs.INK2)
ax.legend(handles=[Line2D([], [], color=fs.DBGPS_C, lw=1, label="DBGPS"),
                   Line2D([], [], color=fs.DBGPS2_C, lw=1, label="DBGPS2"),
                   Line2D([], [], color=fs.INK2, marker="o", lw=0, ms=2.4, label="Control"),
                   Line2D([], [], color=fs.INK2, marker="s", lw=0, ms=2.4, label="HT5")],
          loc="upper right", ncol=2, columnspacing=0.8)
fs.panel(ax, "a")

# b: path overflow under DBGPS
ax = fig.add_subplot(gs[0, 1])
for s in fs.SAMPLES:
    d = grids[s][grids[s].r == 0].sort_values("c")
    ax.plot(d.c, d.overflow / N * 100, color=fs.COLOR[s], lw=1, label=fs.LABEL[s])
ax.set_xscale("log")
ax.set_xlabel("Coverage cutoff c"); ax.set_ylabel("Indices with path overflow (%)")
ax.set_title("DBGPS (r → ∞)", fontsize=6.3, pad=2)
ax.legend(loc="upper right")
fs.panel(ax, "b")

# c: rescued strands in HT5 by minimum true-path coverage
ax = fig.add_subplot(gs[0, 2])
s = "HT5"
g = grids[s]
d1 = g[g.r == 0].loc[g[g.r == 0].correct.idxmax()]
o1 = outcome_col(s, d1.c, 0); o2 = outcome_col(s, 1, R_DEF)
st = stats[s].set_index("id").loc[outc[s][0].id]
mc = st.min_cov.values
groups = [("Both", (o1 == "S") & (o2 == "S"), "#9ec5f4"),
          ("DBGPS2 only", (o1 != "S") & (o2 == "S"), fs.DBGPS2_C),
          ("Neither", (o1 != "S") & (o2 != "S"), "#b9b8b3")]
for lab, m, col in groups:
    v = np.sort(np.maximum(mc[m], 0.5))
    ax.plot(v, np.arange(1, len(v) + 1) / len(v) * 100, color=col, lw=1.2, label=lab)
ax.axvline(d1.c + 0.5, color=fs.DBGPS_C, lw=0.7, ls="--", label=f"DBGPS c = {int(d1.c)}")
ax.set_xscale("log"); ax.set_xlim(0.5, 300); ax.set_ylim(0, 101)
ax.set_xlabel("Minimum true-path k-mer coverage"); ax.set_ylabel("Strands in group (cumulative %)")
ax.set_title("HT5: strands recovered by", fontsize=6.3, pad=2)
ax.legend(loc="lower right", fontsize=5.0, handlelength=1.2, labelspacing=0.3, borderaxespad=0.1)
fs.panel(ax, "c")
only1 = ((o1 == "S") & (o2 != "S")).sum()

# d: recovery relative to the coverage-limited ceiling
ax = fig.add_subplot(gs[1, 0])
rows = []
for i, s in enumerate(fs.SAMPLES):
    st = stats[s]
    intact = (st.min_cov >= 2).mean() * 100
    intact_r = ((st.min_cov >= 2) & (st.r_needed <= R_DEF)).mean() * 100
    g = grids[s]
    b1 = g[g.r == 0].correct.max() / N * 100
    b2 = g[g.r == R_DEF].correct.max() / N * 100
    rows.append(dict(sample=fs.LABEL[s], intact_path=intact, intact_within_r=intact_r, dbgps=b1, dbgps2=b2))
    ax.plot([i, i], [b1, intact], color=fs.GRID, lw=2.2, zorder=0, solid_capstyle="round")
    ax.scatter(i, intact, marker="_", s=60, color=fs.INK, lw=1.1, zorder=3)
    ax.scatter(i, b2, s=16, color=fs.DBGPS2_C, edgecolor="white", lw=0.4, zorder=4)
    ax.scatter(i, b1, s=16, color=fs.DBGPS_C, edgecolor="white", lw=0.4, zorder=4)
ceil = pd.DataFrame(rows)
ceil.to_csv(os.path.join(fs.DATADIR, "recovery_ceiling.csv"), index=False)
ax.set_xticks(range(len(fs.SAMPLES)), [fs.LABEL[s] for s in fs.SAMPLES])
ax.set_ylabel("Strands (%)"); ax.set_ylim(40, 101)
ax.legend(handles=[Line2D([], [], color=fs.INK, marker="_", lw=0, ms=7, label="Intact true path (all k-mers ≥2)"),
                   Line2D([], [], color=fs.DBGPS2_C, marker="o", lw=0, ms=3.5, label="DBGPS2 recovered"),
                   Line2D([], [], color=fs.DBGPS_C, marker="o", lw=0, ms=3.5, label="DBGPS recovered (best c)")],
          loc="lower left")
fs.panel(ax, "d")

# e: why DBGPS2 still loses strands
ax = fig.add_subplot(gs[1, 1])
cats = [("Coverage gap (min cov ≤1)", "#b9b8b3"), ("Ratio > r on true path", "#eda100")]
for i, s in enumerate(fs.SAMPLES):
    o, hdr = outc[s]
    o2 = outcome_col(s, 1, R_DEF)
    st = stats[s].set_index("id").loc[o.id]
    lost = o2 != "S"
    gap = lost & (st.min_cov.values <= 1)
    rat = lost & ~gap & (st.r_needed.values > R_DEF)
    oth = lost & ~gap & ~rat
    bottom = 0
    assert oth.sum() == 0  # every lost strand has a coverage gap or a steeper-than-r true transition
    for (lab, col), m in zip(cats, (gap, rat)):
        v = m.sum() / N * 100
        ax.bar(i, v, 0.7, bottom=bottom, color=col, edgecolor="white", lw=0.3, label=lab if i == 0 else None)
        bottom += v
ax.set_xticks(range(len(fs.SAMPLES)), [fs.LABEL[s] for s in fs.SAMPLES])
ax.set_ylabel("Strands not recovered by DBGPS2 (%)")
ax.set_ylim(0, 6.5)
ax.legend(loc="upper left", fontsize=5.2)
fs.panel(ax, "e")

# f: ratio bound needed to admit the true path
ax = fig.add_subplot(gs[1, 2])
for s in fs.SAMPLES:
    st = stats[s]
    v = np.sort(st.r_needed[st.min_cov >= 2].values)
    ax.plot(v, np.arange(1, len(v) + 1) / N * 100, color=fs.COLOR[s], lw=1, label=fs.LABEL[s])
for rv in (1.5, 2.55):
    ax.axvline(rv, color=fs.DBGPS2_C, lw=0.6, ls="--" if rv == 1.5 else "-")
ax.set_ylim(0, 108)
ax.text(2.65, 102, "r = 2.55", fontsize=5, color=fs.INK2, va="center")
ax.text(1.45, 102, "1.5", fontsize=5, color=fs.INK2, va="center", ha="right")
ax.set_xscale("log"); ax.set_xlim(1, 20)
ax.set_xticks([1, 2, 5, 10, 20], ["1", "2", "5", "10", "20"])
ax.set_xlabel("Smallest r admitting the true path, r$_{needed}$")
ax.set_ylabel("Strands (cumulative %)")
ax.legend(loc="lower right")
fs.panel(ax, "f")

fs.save(fig, "Fig5_mechanism")
print(ceil.round(2).to_string())
print("HT5 DBGPS-only strands:", only1)
