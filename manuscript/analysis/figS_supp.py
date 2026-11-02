#!/usr/bin/env python3
"""Supplementary figures and tables."""
import os, glob, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import figstyle as fs

fs.setup()
FULL = os.path.join(fs.RESULTS, "full")
N = 210000

# ---- Supplementary Fig. 1: recovery landscapes (r x c), all samples, k = 31 ----
fig, axs = plt.subplots(2, 4, figsize=(fs.DOUBLE, 92 * fs.MM), sharex=True, sharey=True)
fig.subplots_adjust(left=0.075, right=0.915, top=0.94, bottom=0.14, wspace=0.08, hspace=0.3)
axs.flat[-1].axis("off")
for ax, s in zip(axs.flat, fs.ALL):
    g = pd.read_csv(f"{FULL}/{s}.grid.tsv", sep="\t")
    rs = sorted(g.r.unique(), key=lambda v: (v == 0, v)); cs = sorted(g.c.unique())
    mat = np.array([[g[(g.r == r) & (g.c == c)].correct.iloc[0] / N * 100 for c in cs] for r in rs])
    im = ax.imshow(mat, aspect="auto", cmap="Blues", vmin=0, vmax=100, origin="lower")
    ax.set_title(fs.LABEL[s], fontsize=6.3, pad=2)
    ax.set_xticks(range(len(cs)), [str(c) for c in cs], rotation=90, fontsize=5)
    ax.set_yticks(range(len(rs)), [("∞" if r == 0 else f"{r:g}") for r in rs], fontsize=5.5)
    ax.spines[:].set_visible(False); ax.tick_params(length=0)
for ax in axs[1]:
    ax.set_xlabel("Coverage cutoff $c$")
for ax in axs[:, 0]:
    ax.set_ylabel("Ratio bound $r$")
cax = fig.add_axes([0.93, 0.25, 0.01, 0.5])
cb = fig.colorbar(im, cax=cax); cb.set_label("Strands recovered (%)", fontsize=6)
cb.ax.tick_params(labelsize=5.5, length=1.5); cb.outline.set_linewidth(0.4)
fs.save(fig, "FigS1_landscapes")

# ---- Supplementary Table 1: per-sample summary ----
qc = pd.read_csv(os.path.join(fs.DATADIR, "qc_summary.csv"))
reads = {"HT176": 60642860, "HT_con": 61656962, "HT1": 57505049, "HT2": 58970885, "HT3": 66224448, "HT4": 78239619, "HT5": 58635539}
rows = []
for s in fs.ALL:
    inf = dict(l.rstrip("\n").split("\t") for l in open(f"{FULL}/{s}.info.tsv") if "\t" in l)
    head = open(f"{FULL}/{s}.spectrum.tsv").readline()
    meta = dict(x.split("=") for x in head[2:].split())
    st = pd.read_csv(f"{FULL}/{s}.strandstats.tsv", sep="\t")
    q = qc[qc["sample"] == fs.LABEL[s]].iloc[0]
    rows.append({"Sample": fs.LABEL[s], "Read pairs (M)": reads[s] / 1e6,
                 "Read pairs per strand": reads[s] / N,
                 "Adapter read-through (%)": q.adapter_frac * 100, "Median fragment (nt)": q.median_fragment, "Pairs with both mates aligned (%)": q.pairs_aligned_frac * 100,
                 "Mapped (%)": q.mapped_frac * 100, "Substitutions per kb": q.sub_rate * 1e3,
                 "C>T+G>A per kb": q.CT_GA_rate * 1e3, "Q30 (%)": q.q30_frac * 100,
                 "Distinct 31-mers (M)": int(inf["distinct_kmers"]) / 1e6,
                 "True 31-mers absent (%)": int(meta["ref_kmers_absent"]) / int(meta["ref_kmers_total"]) * 100,
                 "Median strand coverage": st.median_cov.median(),
                 "Median minimum path coverage": st.min_cov.median(),
                 "Peak memory (GB)": float(inf["peak_rss_mb_final"]) / 1024})
t1 = pd.DataFrame(rows)
t1.to_csv(os.path.join(fs.DATADIR, "SupplementaryTable1_samples.csv"), index=False)
print(t1.round(2).T.to_string())

# ---- Supplementary Table 2: decoding fidelity over the full grid ----
rows = []
for s in fs.ALL:
    g = pd.read_csv(f"{FULL}/{s}.grid.tsv", sep="\t")
    for lab, d in (("DBGPS (r → ∞)", g[g.r == 0]), ("DBGPS2 (r = 1.5–6)", g[g.r > 0]),
                   ("DBGPS2 (r = 2.55)", g[g.r == 2.55])):
        rows.append({"Sample": fs.LABEL[s], "Algorithm": lab, "Settings": len(d),
                     "Strands decoded (sum)": int(d.decoded_ref.sum()),
                     "Wrong strands": int(d.wrong.sum()), "Ambiguous indices": int(d.ambiguous.sum()),
                     "Spurious decodes of unused indices": int(d.spurious_nonref.sum()),
                     "Path overflows": int(d.overflow.sum())})
t2 = pd.DataFrame(rows)
t2.to_csv(os.path.join(fs.DATADIR, "SupplementaryTable2_fidelity.csv"), index=False)
print(t2.to_string())

# ---- Supplementary Table 3: k-mer size sweep ----
kt = pd.read_csv(os.path.join(fs.DATADIR, "ksize_sweep.csv"))
t3 = pd.DataFrame({"Sample": kt["sample"].map(fs.LABEL), "k": kt.k,
                   "Distinct k-mers (M)": (kt.distinct / 1e6).round(1),
                   "Median minimum path coverage": kt.median_min_cov.round().astype(int),
                   "Intact true paths (%)": (kt.intact * 100).round(2),
                   "DBGPS best (%)": (kt.dbgps / N * 100).round(2), "DBGPS best c": kt.dbgps_c,
                   "DBGPS2 r=2.55 best (%)": (kt.dbgps2 / N * 100).round(2), "DBGPS2 best c": kt.dbgps2_c,
                   "Wrong strands, r → ∞ (sum over c)": kt.wrong_r_inf,
                   "Wrong strands, r = 2.55": kt.wrong_dbgps2, "Wrong strands, r = 3": kt.wrong_r3,
                   "Wrong strands, r = 4": kt.wrong_r4})
t3.to_csv(os.path.join(fs.DATADIR, "SupplementaryTable4_ksize.csv"), index=False)
print(t3.to_string())
