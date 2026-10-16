#!/usr/bin/env python3
"""Figure 1: study design and the principle of ratio-constrained path search."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import figstyle as fs

fs.setup()
fig = plt.figure(figsize=(fs.DOUBLE, 98 * fs.MM))
gs = fig.add_gridspec(3, 3, height_ratios=[0.42, 0.26, 1.25], hspace=0.32, wspace=0.34,
                      left=0.06, right=0.985, top=0.975, bottom=0.085)


def box(ax, x, y, w, h, text, fc="#f3f2ef", ec="#b9b8b3", fontsize=5.8, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.02",
                                fc=fc, ec=ec, lw=0.6))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize,
            fontweight="bold" if bold else "normal", color=fs.INK, linespacing=1.25)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=6,
                                 lw=0.7, color=fs.INK2))


# a: workflow
ax = fig.add_subplot(gs[0, :])
ax.set_xlim(-0.01, 1.01); ax.set_ylim(0, 1); ax.axis("off")
steps = [
    ("6.5 MB digital file\n(DBGPS encoding)", "#f3f2ef"),
    ("210,000 oligos\n200 nt each", "#f3f2ef"),
    ("Thermal aging\ncontrol, HT1–HT5", "#fbe3d6"),
    ("Illumina PE150\n58–78 M pairs/sample", "#f3f2ef"),
    ("k-mer counting\n(optional Bloom filters)", "#e3eefb"),
    ("Path search from\neach index seed", "#e3eefb"),
    ("CRC-validated\nstrands", "#e3eefb"),
]
n = len(steps); w = 0.118; gap = (1 - n * w) / (n - 1)
for i, (t, c) in enumerate(steps):
    x = i * (w + gap)
    box(ax, x, 0.36, w, 0.5, t, fc=c)
    if i < n - 1:
        arrow(ax, x + w + 0.004, 0.61, x + w + gap - 0.004, 0.61)
# heat-aging thermometer of samples
for j, s in enumerate(fs.SAMPLES):
    x = 2 * (w + gap) + 0.006 + j * (w - 0.012) / 6
    ax.add_patch(Rectangle((x, 0.2), (w - 0.012) / 6 - 0.003, 0.1, color=fs.COLOR[s], lw=0))
ax.text(2 * (w + gap) + w / 2, 0.02, "increasing heat exposure →", ha="center", fontsize=5, color=fs.INK2)
ax.text(4 * (w + gap), 0.92, "DBGPS2 (this work)", fontsize=6, color=fs.DBGPS2_C, fontweight="bold")
fs.panel(ax, "a", x=-0.045, y=0.85)

# b: strand layout
ax = fig.add_subplot(gs[1, :])
ax.set_xlim(-2, 202); ax.set_ylim(0, 1); ax.axis("off")
segs = [(0, 18, "P1 (18)", "#d9d8d3"), (18, 16, "Index (16)", "#9ec5f4"),
        (34, 140, "Payload (140 nt = 35 bytes)", "#e3eefb"), (174, 8, "CRC", "#86b6ef"),
        (182, 18, "P2 (18)", "#d9d8d3")]
for x0, L, t, c in segs:
    ax.add_patch(Rectangle((x0, 0.35), L, 0.4, fc=c, ec="white", lw=0.8))
    ax.text(x0 + L / 2, 0.55, t, ha="center", va="center", fontsize=5.5, color=fs.INK)
# seed k-mer and extension
ax.annotate("", xy=(3, 0.25), xytext=(34, 0.25), arrowprops=dict(arrowstyle="|-|", lw=0.6, color=fs.INK2,
                                                                  mutation_scale=2))
ax.text(18.5, 0.04, "seed 31-mer", ha="center", fontsize=5.2, color=fs.INK2)
ax.annotate("", xy=(183, 0.25), xytext=(34, 0.25), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=fs.DBGPS2_C,
                                                                    mutation_scale=5))
ax.text(108, 0.04, "149 single-base extensions, each choosing among 4 children", ha="center", fontsize=5.2,
        color=fs.INK2)
fs.panel(ax, "b", x=-0.045, y=0.75)

# c-e: cartoon coverage along a degraded strand
rng = np.random.default_rng(3)
steps_ = np.arange(150)
true = 40 * (0.12 + 0.88 * ((steps_ - 75) / 75) ** 2) * np.exp(rng.normal(0, 0.06, 150))
true[40:44] *= 0.85
noise_x = np.sort(rng.choice(150, 38, replace=False))
noise_y = np.exp(rng.normal(np.log(3.2), 0.55, noise_x.size))
noise_y[noise_x < 25] *= 3.2  # high-coverage regions spawn stronger noise branches
noise_y[noise_x > 125] *= 3.2


def cartoon(ax, title, sub):
    ax.plot(steps_, true, color=fs.INK, lw=1.1, label="True path")
    ax.scatter(noise_x, noise_y, s=6, color=fs.DBGPS_C, lw=0, zorder=3, label="Noise branches")
    ax.set_yscale("log"); ax.set_ylim(1, 80); ax.set_xlim(0, 149)
    ax.set_xlabel("Extension step")
    ax.set_ylabel("k-mer coverage")
    ax.set_title(title, fontsize=6.3, pad=12)
    ax.text(0.5, 1.015, sub, transform=ax.transAxes, ha="center", va="bottom", fontsize=5.4, color=fs.INK2)


ax = fig.add_subplot(gs[2, 0])
cartoon(ax, "Low global cutoff (DBGPS, small c)", "noise branches retained → path explosion")
ax.axhline(2, color=fs.DBGPS_C, lw=0.8, ls="--")
ax.text(148, 2.15, "c", ha="right", va="bottom", fontsize=6, color=fs.INK)
ax.legend(loc="lower left", fontsize=5.2)
fs.panel(ax, "c", x=-0.2, y=1.1)

ax = fig.add_subplot(gs[2, 1])
cartoon(ax, "High global cutoff (DBGPS, large c)", "true path broken where coverage is low → strand lost")
ax.axhline(9, color=fs.DBGPS_C, lw=0.8, ls="--")
ax.text(148, 9.6, "c", ha="right", va="bottom", fontsize=6, color=fs.INK)
lowm = true < 9
ax.scatter(steps_[lowm], true[lowm], s=4, color=fs.INK, marker="x", lw=0.6, zorder=4)
fs.panel(ax, "d", x=-0.2, y=1.1)

ax = fig.add_subplot(gs[2, 2])
cartoon(ax, "Local ratio window (DBGPS2)", "children accepted only within r-fold of parent coverage")
r = 2.55
ax.fill_between(steps_[1:], np.maximum(true[:-1] / r, 2), np.minimum(true[:-1] * r, 80),
                color=fs.DBGPS2_C, alpha=0.18, lw=0, label="cov(parent)/r … cov(parent)·r")
ax.axhline(2, color=fs.DBGPS_C, lw=0.6, ls="--")
ax.legend(loc="lower left", fontsize=5.2, handles=[ax.collections[-1]])
fs.panel(ax, "e", x=-0.2, y=1.1)

fs.save(fig, "Fig1_concept")
