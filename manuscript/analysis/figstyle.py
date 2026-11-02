"""Shared matplotlib style for the manuscript figures (Nature-family layout)."""
import os
import matplotlib as mpl
import matplotlib._mathtext as mathtext
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

MM = 1 / 25.4
SINGLE = 89 * MM
DOUBLE = 183 * MM

SAMPLES = ["HT_con", "HT1", "HT2", "HT3", "HT4", "HT5"]
ALL = SAMPLES + ["HT176"]   # HT176: 176 days at 70 C (not part of the HT1-HT5 series)
LABEL = {"HT_con": "Control", "HT1": "HT1", "HT2": "HT2", "HT3": "HT3", "HT4": "HT4", "HT5": "HT5", "HT176": "HT176"}
# control: neutral ink; heat series: validated single-hue ordinal ramp (light -> dark)
COLOR = {"HT176": "#CC79A7", "HT_con": "#3b3b3b", "HT1": "#f2a07b", "HT2": "#e57850", "HT3": "#d1502f",
         "HT4": "#a8321c", "HT5": "#771b10"}
# colour-blind safe palette (Wong, Nat. Methods 8, 441; 2011), shared by every figure
BLUE, ORANGE, GREEN, PURPLE = "#0072B2", "#D55E00", "#009E73", "#CC79A7"
YELLOW, SKY = "#E69F00", "#56B4E9"
GRAY, LIGHT, RED = "#7B8389", "#D9DEE2", "#9E2A2B"
DBGPS2_C, DBGPS_C = BLUE, ORANGE
# outcome classes keep one colour across figures
OUTCOME = {"correct": BLUE, "no_path": LIGHT, "overflow": ORANGE, "no_valid": YELLOW, "ambiguous": SKY, "wrong": RED}
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"

RESULTS = os.environ.get("DBGPS_RESULTS", "/Users/song/dbgps2_ht/results")
HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "..", "figures")
DATADIR = os.path.join(HERE, "data")


def setup():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica"],
        "font.size": 6.5,
        "axes.titlesize": 6.5,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6,
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
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "text.color": INK,
        "lines.linewidth": 1.0,
        "lines.markersize": 3,
        "legend.frameon": False,
        "legend.handlelength": 1.4,
        "legend.borderaxespad": 0.2,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": 600,
        "savefig.facecolor": "white",
        "figure.dpi": 150,
    })
    mathtext.SHRINK_FACTOR = 0.8  # keeps axis exponents at 5 pt or more
    os.makedirs(FIGDIR, exist_ok=True)
    os.makedirs(DATADIR, exist_ok=True)


def panel(ax, letter, x=-0.18, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="left", color=INK, gid="panel")


def plabel(ax, letter, dx=-28, dy=5):
    """Panel letter at a fixed offset from the top left of the axes, in points."""
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(dx, dy), textcoords="offset points",
                fontsize=8, fontweight="bold", va="bottom", ha="left", color=INK, gid="panel")


def finalize(fig):
    """Apply the text rules: 5-7 pt, no coloured text except white on dark fills."""
    for t in fig.findobj(mpl.text.Text):
        if t.get_gid() == "panel":
            continue
        t.set_fontsize(min(7, max(5, t.get_fontsize())))
        if min(to_rgb(t.get_color())) < 0.95:
            t.set_color(INK)


def save(fig, name):
    finalize(fig)
    fig.savefig(os.path.join(FIGDIR, name + ".pdf"))
    fig.savefig(os.path.join(FIGDIR, name + ".png"), dpi=450)
    plt.close(fig)
