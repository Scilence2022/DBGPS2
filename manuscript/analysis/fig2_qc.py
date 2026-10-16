#!/usr/bin/env python3
"""Figure 2: sequencing-level signatures of thermal aging (alignment QC)."""
import json, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import figstyle as fs

fs.setup()
QC = os.path.join(fs.RESULTS, "qc")
TYPES = ["C>T", "G>A", "A>G", "T>C", "C>A", "G>T", "A>C", "T>G", "C>G", "G>C", "A>T", "T>A"]

rows, subs, pos, ins_len, per_strand, edit = [], {}, {}, {}, {}, {}
for s in fs.SAMPLES:
    js = [json.load(open(f"{QC}/{s}_{r}.json")) for r in (1, 2)]
    agg = {k: sum(j[k] for j in js) for k in ("reads", "mapped", "aligned_bases", "mismatches",
                                               "ins_bases", "del_bases")}
    sub = {t: sum(j["subs"].get(t, 0) for j in js) for t in TYPES}
    subs[s] = sub
    n_sub = sum(sub.values())
    adapter = np.mean([j["adapter_frac"] for j in js])
    polyg = np.mean([j["polyg_frac"] for j in js])
    # insert length: adapter start position in the read (reads without adapter: >= 150)
    h = np.zeros(151)
    for j in js:
        for k, v in j["adapter_pos"].items():
            h[min(int(k), 150)] += v
        h[150] += j["reads"] - sum(j["adapter_pos"].values())
    ins_len[s] = h
    e = np.zeros(21)
    for j in js:
        for k, v in j["edit_hist"].items():
            e[int(k)] += v
    edit[s] = e
    p = sum(np.loadtxt(f"{QC}/{s}_{r}.pos.tsv", skiprows=1) for r in (1, 2))
    p[:, 0] /= 2  # position column was summed over the two read files
    pos[s] = p
    ps = pd.concat([pd.read_csv(f"{QC}/{s}_{r}.per_strand.tsv", sep="\t", header=None, names=["id", "n"])
                    for r in (1, 2)]).groupby("id").n.sum()
    per_strand[s] = ps
    ct = sub["C>T"] + sub["G>A"]
    rows.append(dict(sample=fs.LABEL[s], reads_qc=agg["reads"], mapped_frac=agg["mapped"] / agg["reads"],
                     adapter_frac=adapter, polyG_frac=polyg,
                     mean_q=np.mean([j["mean_q"] for j in js]), q30_frac=np.mean([j["q30_frac"] for j in js]),
                     sub_rate=agg["mismatches"] / agg["aligned_bases"],
                     ins_rate=agg["ins_bases"] / agg["aligned_bases"],
                     del_rate=agg["del_bases"] / agg["aligned_bases"],
                     CT_GA_frac_of_subs=ct / n_sub,
                     CT_GA_rate=ct / agg["aligned_bases"],
                     perfect_read_frac=float(e[0] / e.sum()),
                     median_insert=int(np.searchsorted(np.cumsum(h) / h.sum(), 0.5)),
                     strands_hit=len(ps), reads_per_strand_cv=float(ps.std() / ps.mean())))
qc = pd.DataFrame(rows)
qc.to_csv(os.path.join(fs.DATADIR, "qc_summary.csv"), index=False)
print(qc.round(5).to_string())

fig = plt.figure(figsize=(fs.DOUBLE, 112 * fs.MM))
gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.42, left=0.07, right=0.985, top=0.95, bottom=0.09)
x = np.arange(len(fs.SAMPLES))
cols = [fs.COLOR[s] for s in fs.SAMPLES]
labels = [fs.LABEL[s] for s in fs.SAMPLES]

# a: insert length cumulative distribution
ax = fig.add_subplot(gs[0, 0])
for s in fs.SAMPLES:
    h = ins_len[s]
    c = np.cumsum(h) / h.sum()
    ax.step(np.arange(151), c * 100, where="post", color=fs.COLOR[s], lw=1, label=fs.LABEL[s])
ax.set_xlabel("Insert length (nt, from adapter position)")
ax.set_ylabel("Reads with insert ≤ length (%)")
ax.set_xlim(0, 150); ax.set_ylim(0, 100)
ax.text(148, 4, "≥150 nt\n(full-length\nreads)", ha="right", va="bottom", fontsize=5.3, color=fs.INK2)
ax.legend(loc="upper left", ncol=1, handlelength=1.0)
fs.panel(ax, "a")

# b: fraction short inserts and mapping
ax = fig.add_subplot(gs[0, 1])
w = 0.38
ax.bar(x - w / 2, qc.adapter_frac * 100, w, color=cols, edgecolor="white", lw=0.4)
ax.bar(x + w / 2, (1 - qc.mapped_frac) * 100, w, color=cols, edgecolor="white", lw=0.4, hatch="////", alpha=0.55)
for i, v in enumerate(qc.adapter_frac * 100):
    ax.text(i - w / 2, v + 1.5, f"{v:.0f}", ha="center", fontsize=5.2, color=fs.INK2)
ax.set_xticks(x, labels, rotation=0)
ax.set_ylabel("Reads (%)")
ax.set_ylim(0, 112)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor="#9a9a9a", label="Adapter read-through (insert <150 nt)"),
                   Patch(facecolor="#9a9a9a", alpha=0.55, hatch="////", label="Unmapped")],
          loc="upper left", bbox_to_anchor=(0, 1.13))
fs.panel(ax, "b")

# c: per-base error rates
ax = fig.add_subplot(gs[0, 2])
mk = {"sub_rate": ("o", "Substitution"), "del_rate": ("s", "Deletion"), "ins_rate": ("^", "Insertion")}
for k, (m, lab) in mk.items():
    ax.plot(x, qc[k] * 1e3, color=fs.INK2, lw=0.6, zorder=1)
    ax.scatter(x, qc[k] * 1e3, c=cols, marker=m, s=14, edgecolor="white", linewidth=0.4, zorder=2, label=lab)
ax.set_xticks(x, labels)
ax.set_ylabel("Errors per 1,000 aligned bases")
leg = ax.legend(loc="upper left", markerscale=1)
for hnd in leg.legend_handles:
    hnd.set_facecolor("#9a9a9a")
fs.panel(ax, "c")

# d: substitution spectrum heat map
ax = fig.add_subplot(gs[1, 0])
mat = np.array([[subs[s][t] / sum(subs[s].values()) * 100 for t in TYPES] for s in fs.SAMPLES])
cmap = LinearSegmentedColormap.from_list("seq", ["#fbf3ee", "#f2a07b", "#d1502f", "#771b10"])
im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=30)
ax.set_xticks(range(len(TYPES)), TYPES, rotation=90)
ax.set_yticks(range(len(fs.SAMPLES)), labels)
for i in range(mat.shape[0]):
    for j in range(2):
        ax.text(j, i, f"{mat[i, j]:.0f}", ha="center", va="center", fontsize=5,
                color="white" if mat[i, j] > 18 else fs.INK)
ax.spines[:].set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
cb.set_label("Share of substitutions (%)", fontsize=5.8)
cb.ax.tick_params(labelsize=5.5, length=1.5)
cb.outline.set_linewidth(0.4)
ax.set_xlabel("Substitution (designed strand → read)")
fs.panel(ax, "d", x=-0.22)

# e: positional depth along the strand
ax = fig.add_subplot(gs[1, 1])
for s in fs.SAMPLES:
    p = pos[s]
    depth = p[:, 1]
    ax.plot(p[:, 0], depth / depth[18:182].max(), color=fs.COLOR[s], lw=1, label=fs.LABEL[s])
ax.axvspan(0, 18, color=fs.GRID, lw=0); ax.axvspan(182, 200, color=fs.GRID, lw=0)
ax.text(9, 0.04, "P1", ha="center", fontsize=5.5, color=fs.INK2)
ax.text(191, 0.04, "P2", ha="center", fontsize=5.5, color=fs.INK2)
ax.set_xlim(0, 200); ax.set_ylim(0, 1.05)
ax.set_xlabel("Position along 200-nt strand")
ax.set_ylabel("Relative read depth")
fs.panel(ax, "e")

# f: substitution rate along the strand
ax = fig.add_subplot(gs[1, 2])
for s in fs.SAMPLES:
    p = pos[s]
    rate = p[:, 2] / np.maximum(p[:, 1], 1) * 1e3
    k = 5
    sm = np.convolve(rate, np.ones(k) / k, mode="same")
    ax.plot(p[2:-2, 0], sm[2:-2], color=fs.COLOR[s], lw=1)
ax.axvspan(0, 18, color=fs.GRID, lw=0); ax.axvspan(182, 200, color=fs.GRID, lw=0)
ax.set_xlim(0, 200); ax.set_ylim(bottom=0)
ax.set_xlabel("Position along 200-nt strand")
ax.set_ylabel("Substitutions per 1,000 bases")
fs.panel(ax, "f")

fs.save(fig, "Fig2_sequencing_damage")
