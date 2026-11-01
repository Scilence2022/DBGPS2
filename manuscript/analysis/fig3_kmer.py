#!/usr/bin/env python3
"""Figure 3: thermal aging erodes the coverage separation between true and noise k-mers."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import figstyle as fs

fs.setup()
FULL = os.path.join(fs.RESULTS, "full")
DQ = os.path.join(fs.RESULTS, "dstoreqc")
samples = [s for s in fs.SAMPLES if os.path.exists(f"{FULL}/{s}.strandstats.tsv")]


def spectrum(s):
    d = pd.read_csv(f"{FULL}/{s}.spectrum.tsv", sep="\t", comment="#")
    head = open(f"{FULL}/{s}.spectrum.tsv").readline()
    meta = dict(x.split("=") for x in head[2:].split())
    return d, int(meta["ref_kmers_total"]), int(meta["ref_kmers_absent"])


rows = []
fig = plt.figure(figsize=(fs.DOUBLE, 118 * fs.MM))
gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.42, left=0.07, right=0.985, top=0.95, bottom=0.08)

# a: coverage spectra, true vs noise (control and most aged)
ax = fig.add_subplot(gs[0, 0])
show = [s for s in ("HT_con", "HT5") if s in samples]
for s in show:
    d, tot, absent = spectrum(s)
    ax.plot(d["cov"], d.ref, color=fs.COLOR[s], lw=1.1, label=f"{fs.LABEL[s]} true")
    ax.plot(d["cov"], d.noise, color=fs.COLOR[s], lw=0.8, ls=(0, (2, 1.2)), label=f"{fs.LABEL[s]} noise")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1, 3000); ax.set_ylim(1, 2e9)
ax.set_xlabel("k-mer coverage")
ax.set_ylabel("Distinct k-mers")
ax.legend(loc="upper right", ncol=1)
fs.panel(ax, "a")

# b: per-sample Kd and Kn as a function of cutoff c (from the spectra)
ax = fig.add_subplot(gs[0, 1])
cs = np.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80, 100])
for s in samples:
    d, tot, absent = spectrum(s)
    cov = d["cov"].values; ref = d.ref.values; noise = d.noise.values
    kd, kn = [], []
    for c in cs:
        ref_above = ref[cov > c].sum()
        kd.append(1 - ref_above / tot)
        kn.append(noise[cov > c].sum() / max(ref_above, 1))
        rows.append(dict(sample=fs.LABEL[s], c=c, Kd=kd[-1], Kn=kn[-1]))
    ax.plot(np.array(kd) * 100, kn, color=fs.COLOR[s], lw=1, marker="o", ms=2.2, mec="white", mew=0.3,
            label=fs.LABEL[s])
    for c, xk, yk in zip(cs, kd, kn):
        if s in ("HT_con", "HT5") and c in (1, 10, 50):
            ax.annotate(f"c={c}", (xk * 100, yk), fontsize=5, color=fs.INK2, xytext=(3, 2),
                        textcoords="offset points")
ax.set_xscale("symlog", linthresh=0.1); ax.set_yscale("log")
ax.set_xlabel("True k-mers lost, Kd (%)")
ax.set_ylabel("Noise / true k-mers retained, Kn")
ax.legend(loc="upper right", ncol=2, columnspacing=0.8)
fs.panel(ax, "b")
pd.DataFrame(rows).to_csv(os.path.join(fs.DATADIR, "kmer_KdKn_by_cutoff.csv"), index=False)

# c: distribution of the minimum true-path coverage per strand
ax = fig.add_subplot(gs[0, 2])
stats = {s: pd.read_csv(f"{FULL}/{s}.strandstats.tsv", sep="\t") for s in samples}
for s in samples:
    v = np.sort(stats[s].min_cov.values)
    ax.plot(np.maximum(v, 0.5), np.arange(1, len(v) + 1) / len(v) * 100, color=fs.COLOR[s], lw=1)
ax.set_xscale("log")
ax.set_xlim(0.5, 2000)
ax.set_xlabel("Minimum k-mer coverage along strand path")
ax.set_ylabel("Strands (cumulative %)")
fs.panel(ax, "c")

# d: coverage ratio distributions, true transitions vs noise branches
ax = fig.add_subplot(gs[1, 0])
for s in show:
    h = pd.read_csv(f"{FULL}/{s}.ratiohist.tsv", sep="\t")
    x = h.log2_ratio_lo + 0.025
    ax.plot(x, h["true"] / h["true"].sum() * 100, color=fs.COLOR[s], lw=1.1, label=f"{fs.LABEL[s]} true → true")
    ax.plot(x, h["alt_ge2"] / max(h["alt_ge2"].sum(), 1) * 100, color=fs.COLOR[s], lw=0.8, ls=(0, (2, 1.2)),
            label=f"{fs.LABEL[s]} true → noise (cov ≥2)")
for v in (-np.log2(2.55), np.log2(2.55)):
    ax.axvline(v, color=fs.DBGPS2_C, lw=0.7)
ax.set_yscale("log"); ax.set_ylim(1e-4, 3e5)
ax.set_yticks([1e-4, 1e-2, 1, 100], ["10$^{-4}$", "10$^{-2}$", "1", "100"])
ax.spines["left"].set_bounds(1e-4, 100); ax.yaxis.set_minor_locator(plt.NullLocator())
ax.set_xlim(-10, 3)
ax.text(1.5, 300, "r = 2.55\nwindow", ha="left", va="center", fontsize=5, color=fs.INK2)
ax.set_xlabel("log$_2$(child coverage / parent coverage)")
ax.set_ylabel("Transitions (%)")
ax.legend(loc="upper left", fontsize=5.0, bbox_to_anchor=(0, 1.02))
fs.panel(ax, "d")

# e: noise branches surviving a cutoff vs a ratio window
ax = fig.add_subplot(gs[1, 1])
x = np.arange(len(samples))
w = 0.38
a = [stats[s].n_alt_ge2.mean() for s in samples]
b = [stats[s].n_alt_ge2_r255.mean() for s in samples]
ax.bar(x - w / 2, a, w, color=fs.DBGPS_C, label="cov ≥ 2 (c = 1)")
ax.bar(x + w / 2, b, w, color=fs.DBGPS2_C, label="cov ≥ 2 and within r = 2.55")
ax.set_xticks(x, [fs.LABEL[s] for s in samples])
ax.set_yscale("log"); ax.set_ylim(0.005, 5000)
ax.set_yticks([0.01, 0.1, 1, 10, 100], ["0.01", "0.1", "1", "10", "100"])
ax.spines["left"].set_bounds(0.005, 100); ax.yaxis.set_minor_locator(plt.NullLocator())
for i in range(len(samples)):
    ax.text(x[i] - w / 2, a[i] * 1.2, f"{a[i]:.0f}", ha="center", fontsize=5, color=fs.INK2)
    ax.text(x[i] + w / 2, b[i] * 1.2, f"{b[i]:.2f}", ha="center", fontsize=5, color=fs.INK2)
ax.set_ylabel("Noise branches per strand path")
ax.legend(loc="upper left", ncol=1)
fs.panel(ax, "e")

# f: an example coverage trace from the most aged sample
ax = fig.add_subplot(gs[1, 2])
s = show[-1]
tr = pd.read_csv(f"{FULL}/{s}.trace.tsv", sep="\t")
st = stats[s].set_index("id")
cand = [i for i in tr.id.unique() if st.loc[i, "min_cov"] >= 3 and st.loc[i, "n_alt_ge2"] >= 3]
ex = cand[0] if cand else tr.id.unique()[0]
t = tr[tr.id == ex]
ax.plot(t.step, t.true_cov, color=fs.INK, lw=1, label="True path")
lo = t.parent_cov / 2.55; hi = t.parent_cov * 2.55
ax.fill_between(t.step, lo, hi, color=fs.DBGPS2_C, alpha=0.15, lw=0, label="Ratio window (r = 2.55)")
alts = t[["alt1", "alt2", "alt3"]].values
for j in range(3):
    m = alts[:, j] > 0
    ax.scatter(t.step[m], alts[m, j], s=5, color=fs.DBGPS_C, lw=0, label="Noise branch" if j == 0 else None)
ax.set_yscale("log"); ax.set_ylim(0.7, 3000)
ax.set_xlabel("Extension step")
ax.set_ylabel("k-mer coverage")
ax.set_title(f"{fs.LABEL[s]}, strand {ex}", fontsize=6)
ax.legend(loc="upper right", fontsize=5.2)
fs.panel(ax, "f")

fs.save(fig, "Fig3_kmer_graph")
