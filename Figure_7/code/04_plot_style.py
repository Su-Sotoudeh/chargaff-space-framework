"""
04_plot_style.py
------------------
Shared matplotlib/seaborn style settings for publication-grade figures
(Nature Communications / NAR conventions):
  - sans-serif font (Arial/Helvetica fallback to DejaVu Sans)
  - colorblind-safe palette (Okabe-Ito)
  - 300+ dpi raster, vector PDF also saved
  - consistent panel labeling helper
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

# ---- Okabe-Ito colorblind-safe palette ----
OKABE_ITO = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "black": "#000000",
    "gray": "#999999",
    "cyan": '#66CCEE',
    "red": "red",
    'Teal':'#009E73',
    "Reddish":'#CC79A7',
    '1':"#D55E00",
    '2':"#0072B2"
}

BASE_COLORS = {
    "A": OKABE_ITO["vermillion"],
    "T": OKABE_ITO["blue"],
    "G": OKABE_ITO["bluish_green"],
    "C": OKABE_ITO["orange"],
}

METHOD_TYPE_COLORS = {
    "formula": OKABE_ITO["blue"],
    "ml": OKABE_ITO["vermillion"],
}


def set_publication_style():
    sns.set_theme(style="ticks", context="paper")
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.titlesize": 12,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "savefig.dpi": 300,
        "figure.dpi": 150,
        "pdf.fonttype": 42,   # editable text in PDF (Illustrator-friendly)
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.bbox": "tight",
    })


def add_panel_label(ax, label, x=-0.15, y=1.08, fontsize=13):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=fontsize,
             fontweight="bold", va="top", ha="left")


def save_fig(fig, path_no_ext):
    fig.savefig(f"{path_no_ext}-300.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{path_no_ext}-600.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{path_no_ext}.pdf", bbox_inches="tight")
    fig.savefig(f"{path_no_ext}.svg", bbox_inches="tight")
    print(f"Saved: {path_no_ext}.png / .pdf")
