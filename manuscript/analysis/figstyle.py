"""Shared matplotlib style for the manuscript figures (Nature-family layout)."""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt

MM = 1 / 25.4
SINGLE = 89 * MM
DOUBLE = 183 * MM

SAMPLES = ["HT_con", "HT1", "HT2", "HT3", "HT4", "HT5"]
LABEL = {"HT_con": "Control", "HT1": "HT1", "HT2": "HT2", "HT3": "HT3", "HT4": "HT4", "HT5": "HT5"}
# control: neutral ink; heat series: validated single-hue ordinal ramp (light -> dark)
COLOR = {"HT_con": "#3b3b3b", "HT1": "#f2a07b", "HT2": "#e57850", "HT3": "#d1502f",
         "HT4": "#a8321c", "HT5": "#771b10"}
# methods: categorical slots 1-2 of the validated palette
DBGPS2_C = "#2a78d6"
DBGPS_C = "#eb6834"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"

RESULTS = os.environ.get("DBGPS_RESULTS", "/Users/song/dbgps2_ht/results")
HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "..", "figures")
DATADIR = os.path.join(HERE, "..", "data")


def setup():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6.5,
        "axes.titlesize": 6.5,
        "axes.labelsize": 6.5,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 5.8,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.size": 2.2,
        "ytick.major.size": 2.2,
        "xtick.minor.size": 1.2,
        "ytick.minor.size": 1.2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "lines.linewidth": 1.0,
        "lines.markersize": 3,
        "legend.frameon": False,
        "legend.handlelength": 1.4,
        "legend.borderaxespad": 0.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": 600,
        "figure.dpi": 150,
    })
    os.makedirs(FIGDIR, exist_ok=True)
    os.makedirs(DATADIR, exist_ok=True)


def panel(ax, letter, x=-0.18, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="left", color=INK)


def save(fig, name):
    fig.savefig(os.path.join(FIGDIR, name + ".pdf"))
    fig.savefig(os.path.join(FIGDIR, name + ".png"), dpi=300)
    plt.close(fig)
