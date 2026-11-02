#!/usr/bin/env python3
"""Figure 1: DBGPS reconstruction and the local coverage ratio added by DBGPS2.

All counts and coverage profiles are schematic. The worked window in panel f uses
the integer bounds of the search: max(floor(p/r), c+1) .. min(floor(p*r), 16383).
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import figstyle as fs

SHARED, SHARED_EC = "#f3f2ef", "#b9b8b3"
ADDED, ADDED_EC = "#e3eefb", fs.DBGPS2_C
REJECT = "#a9a8a3"

fs.setup()
fig = plt.figure(figsize=(fs.DOUBLE, 124 * fs.MM))
gs = fig.add_gridspec(3, 4, height_ratios=[0.58, 1.0, 1.08], hspace=0.3, wspace=0.42,
                      left=0.055, right=0.985, top=0.975, bottom=0.07)


def box(ax, x, y, w, h, text, fc=SHARED, ec=SHARED_EC, fontsize=5.6, bold=False, color=fs.INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.012",
                                fc=fc, ec=ec, lw=0.6))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize,
            fontweight="bold" if bold else "normal", color=color, linespacing=1.25)


def arrow(ax, x0, y0, x1, y1, color=fs.INK2, ls="-", lw=0.7, zorder=1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=6, lw=lw,
                                 color=color, ls=ls, shrinkA=0, shrinkB=0, zorder=zorder))


# a: the shared framework, with the DBGPS2 additions attached to the stage they act on
ax = fig.add_subplot(gs[0, :])
ax.set_xlim(-0.004, 1.004); ax.set_ylim(0, 1); ax.axis("off")
steps = ["Sequencing\nreads", "Count $k$-mers\ninto a graph with\n$k$-mer coverage", "Build a seed\nfor each index",
         "Extend candidate paths\nthrough admissible\nchild $k$-mers", "Validate\ncomplete paths\nwith the code",
         "Accept the strand\nif one valid\nsequence remains"]
widths = [0.085, 0.145, 0.105, 0.31, 0.11, 0.135]
gap = (1 - sum(widths)) / (len(steps) - 1)
xs = []
x = 0.0
for i, (t, w) in enumerate(zip(steps, widths)):
    box(ax, x, 0.50, w, 0.44, t)
    xs.append(x)
    if i < len(steps) - 1:
        arrow(ax, x + w + 0.006, 0.72, x + w + gap - 0.006, 0.72)
    x += w + gap
half = widths[3] / 2 - 0.012
tags = [(xs[1] + 0.002, widths[1] - 0.004, "-F two-pass singleton filter\n(optional, less memory)", False),
        (xs[3] + 0.002, half, "-B query prefilter\n(optional, fewer lookups)", False),
        (xs[3] + widths[3] - 0.002 - half, half, "-r local coverage ratio\n(core algorithmic change)", True)]
for tx, tw, t, bold in tags:
    box(ax, tx, 0.03, tw, 0.28, t, fc=ADDED, ec=ADDED_EC, fontsize=5.1, bold=bold, color=fs.INK)
    arrow(ax, tx + tw / 2, 0.33, tx + tw / 2, 0.475, color=fs.DBGPS2_C, lw=0.6)
lx = xs[4] + 0.01
ax.add_patch(Rectangle((lx, 0.21), 0.016, 0.09, fc=SHARED, ec=SHARED_EC, lw=0.5))
ax.text(lx + 0.022, 0.255, "Shared with DBGPS", fontsize=5.3, va="center")
ax.add_patch(Rectangle((lx, 0.05), 0.016, 0.09, fc=ADDED, ec=ADDED_EC, lw=0.5))
ax.text(lx + 0.022, 0.095, "Added in DBGPS2", fontsize=5.3, va="center", color=fs.DBGPS2_C)
fs.panel(ax, "a", x=-0.05, y=0.82)

# b: DBGPS path search on a schematic graph
ax = fig.add_subplot(gs[1, :])
ax.set_xlim(0, 100); ax.set_ylim(0, 1); ax.axis("off")
segs = [(4, 7, "Primer", "#d9d8d3"), (11, 7, "Index", "#9ec5f4"), (18, 47, "Payload", "#e3eefb"),
        (65, 6, "Code", "#86b6ef"), (71, 7, "Primer", "#d9d8d3")]
for x0, L, t, c in segs:
    ax.add_patch(Rectangle((x0, 0.87), L, 0.11, fc=c, ec="white", lw=0.8))
    ax.text(x0 + L / 2, 0.925, t, ha="center", va="center", fontsize=5.3)
ax.annotate("", xy=(7.5, 0.825), xytext=(18, 0.825),
            arrowprops=dict(arrowstyle="|-|,widthA=0.18,widthB=0.18", lw=0.6, color=fs.INK2, shrinkA=0, shrinkB=0))
ax.text(12.75, 0.75, "seed $k$-mer", ha="center", va="center", fontsize=5.2, color=fs.INK2)
arrow(ax, 19, 0.825, 71, 0.825, color=fs.INK2, lw=0.6)
ax.text(45, 0.755, "one base per step; four possible children, kept if cov(child) > $c$ (here $c$ = 2)",
        ha="center", va="center", fontsize=5.3, color=fs.INK2)


def node(ax, x, y, n, color=fs.INK, fc="white", tc=None):
    ax.scatter([x], [y], s=62, fc=fc, ec=color, lw=0.8, zorder=3)
    ax.text(x, y - 0.004, str(n), ha="center", va="center", fontsize=5.1, zorder=4, color=tc or fs.INK)


true_x = [12.75, 22, 31, 40, 49, 58, 67]
true_n = [42, 40, 43, 38, 36, 39, 41]
yt, ye, yd = 0.53, 0.25, 0.25
for a, b in zip(true_x[:-1], true_x[1:]):
    arrow(ax, a + 1.6, yt, b - 1.6, yt, color=fs.INK, lw=0.9)
for i, (x, n) in enumerate(zip(true_x, true_n)):
    node(ax, x, yt, n, fc="#ececea" if i == 0 else "white")
# a low-count error child is dropped
arrow(ax, 23.3, yt - 0.05, 29.6, yd + 0.05, color=REJECT, ls="--", lw=0.7)
node(ax, 31, yd, 1, color=REJECT, tc=fs.INK2)
ax.text(31, yd - 0.135, "cov ≤ $c$: dropped", ha="center", va="center", fontsize=5.2, color=fs.INK2)
# a stronger error child passes the cutoff and becomes a second candidate
ex = [49, 58, 67]
arrow(ax, 41.3, yt - 0.05, 47.6, yd + 0.05, color=fs.DBGPS_C, lw=0.8)
for a, b in zip(ex[:-1], ex[1:]):
    arrow(ax, a + 1.6, ye, b - 1.6, ye, color=fs.DBGPS_C, lw=0.8)
for x, n in zip(ex, [6, 5, 6]):
    node(ax, x, ye, n, color=fs.DBGPS_C)
ax.text(58, ye - 0.135, "error branch with cov > $c$: kept as a second candidate", ha="center", va="center",
        fontsize=5.2, color=fs.DBGPS_C)
arrow(ax, 68.8, yt, 74.8, yt, color=fs.INK, lw=0.7)
arrow(ax, 68.8, ye, 74.8, ye, color=fs.DBGPS_C, lw=0.7)
box(ax, 75, yt - 0.075, 8.2, 0.15, "Code valid", fc="#e7f3ec", ec="#3f9a6b", fontsize=5.3)
box(ax, 75, ye - 0.075, 8.2, 0.15, "Code invalid", fc="#fbe3d6", ec=fs.DBGPS_C, fontsize=5.3)
arrow(ax, 83.6, yt, 87.6, yt, color=fs.INK, lw=0.7)
box(ax, 88, yt - 0.105, 11.6, 0.21, "One distinct valid\nsequence: strand\naccepted", fontsize=5.1)
ax.text(93.8, ye, "Too many candidates:\nbudget exhausted,\nindex lost", ha="center", va="center",
        fontsize=5.0, color=fs.INK2, linespacing=1.25)
fs.panel(ax, "b", x=-0.05, y=0.84)

# c-e: one coverage profile under the two pruning rules
rng = np.random.default_rng(3)
steps_ = np.arange(150)
true = 40 * (0.12 + 0.88 * ((steps_ - 75) / 75) ** 2) * np.exp(rng.normal(0, 0.06, 150))
true[40:44] *= 0.85
noise_x = np.sort(rng.choice(150, 38, replace=False))
noise_y = np.exp(rng.normal(np.log(3.2), 0.55, noise_x.size))
noise_y[noise_x < 25] *= 3.2  # strongly covered regions spawn stronger error branches
noise_y[noise_x > 125] *= 3.2


def cartoon(ax, title, sub, ylabel=True):
    ax.plot(steps_, true, color=fs.INK, lw=1.0, label="True path")
    ax.scatter(noise_x, noise_y, s=4, color=fs.DBGPS_C, lw=0, zorder=3, label="Error branches")
    ax.set_yscale("log"); ax.set_ylim(1, 80); ax.set_xlim(0, 149)
    ax.set_xticks([0, 50, 100, 149]); ax.set_xticklabels(["0", "50", "100", "149"])
    ax.set_xlabel("Extension step")
    if ylabel:
        ax.set_ylabel("$k$-mer coverage")
    ax.set_title(title, fontsize=6.2, pad=11, fontweight="bold")
    ax.text(0.5, 1.02, sub, transform=ax.transAxes, ha="center", va="bottom", fontsize=5.1, color=fs.INK2)


ax = fig.add_subplot(gs[2, 0])
cartoon(ax, "DBGPS, low cutoff", "error branches kept, search overflows")
ax.axhline(2, color=fs.DBGPS_C, lw=0.8, ls="--")
ax.text(147, 2.15, "$c$", ha="right", va="bottom", fontsize=6)
ax.legend(loc="upper center", fontsize=5.0, handletextpad=0.4)
fs.panel(ax, "c", x=-0.3, y=1.1)

ax = fig.add_subplot(gs[2, 1])
cartoon(ax, "DBGPS, high cutoff", "true path broken, strand lost", ylabel=False)
ax.axhline(9, color=fs.DBGPS_C, lw=0.8, ls="--")
ax.text(147, 9.6, "$c$", ha="right", va="bottom", fontsize=6)
lowm = true < 9
ax.scatter(steps_[lowm], true[lowm], s=3, color=fs.INK, marker="x", lw=0.5, zorder=4)
fs.panel(ax, "d", x=-0.2, y=1.1)

ax = fig.add_subplot(gs[2, 2])
cartoon(ax, "DBGPS2, local window", "child must stay within $r$-fold of its parent", ylabel=False)
r = 2.55
band = ax.fill_between(steps_[1:], np.maximum(true[:-1] / r, 2), np.minimum(true[:-1] * r, 80),
                       color=fs.DBGPS2_C, alpha=0.18, lw=0)
ax.axhline(2, color=fs.DBGPS_C, lw=0.6, ls="--")
ax.text(147, 2.15, "$c$", ha="right", va="bottom", fontsize=6)
ax.legend([band], ["Window"], loc="upper center", fontsize=5.0, handletextpad=0.4)
fs.panel(ax, "e", x=-0.2, y=1.1)

# f: worked integer window at r = 2.55, c = 1
ax = fig.add_subplot(gs[2, 3])
r, c, child = 2.55, 1, 3
rows = [(1.0, 50, "Parent coverage 50", "rejected", fs.DBGPS_C), (0.0, 5, "Parent coverage 5", "kept", fs.DBGPS2_C)]
for y, p, label, verdict, col in rows:
    lo, hi = max(math.floor(p / r), c + 1), min(math.floor(p * r), 16383)
    ax.add_patch(Rectangle((lo, y - 0.13), hi - lo, 0.26, fc=fs.DBGPS2_C, alpha=0.18, lw=0))
    ax.plot([p, p], [y - 0.13, y + 0.13], color=fs.INK, lw=1.1, solid_capstyle="butt")
    ax.scatter([child], [y], s=14, fc="white" if verdict == "rejected" else col, ec=col, lw=0.8, zorder=3)
    ax.text(1.1, y + 0.3, label, fontsize=5.3, va="center", bbox=dict(fc="white", ec="none", pad=0.6))
    ax.text(hi * 1.1, y, f"{lo}–{hi}", fontsize=5.2, va="center", color=fs.DBGPS2_C)
    ax.text(child, y - 0.3, f"child 3: {verdict}", fontsize=5.2, va="center", ha="center", color=col,
            bbox=dict(fc="white", ec="none", pad=0.6))
ax.axvline(child, color=fs.INK2, lw=0.4, ls=":", zorder=0)
ax.set_xscale("log"); ax.set_xlim(1, 400); ax.set_ylim(-0.55, 1.55)
ax.set_xticks([1, 10, 100]); ax.set_xticklabels(["1", "10", "100"])
ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlabel("$k$-mer coverage")
ax.set_title("Same count, opposite decision", fontsize=6.2, pad=11, fontweight="bold")
ax.text(0.5, 1.02, "integer window at $r$ = 2.55, $c$ = 1", transform=ax.transAxes, ha="center", va="bottom",
        fontsize=5.1, color=fs.INK2)
fs.panel(ax, "f", x=-0.09, y=1.1)

fs.save(fig, "Fig1_concept")
