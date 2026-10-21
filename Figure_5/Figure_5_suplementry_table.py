"""
Supplementary figure — per-class summary of the Panel B Spearman results
=========================================================================
One heatmap: rows = pooled data + seven feature classes, columns = the key
pairwise Spearman correlations from each class's 20x20 Panel B matrix.

Values are read from the CSVs written by Figure_5_ABC.py
(Figure5_panelB_spearman_<tag>.csv). If a CSV is not found, the built-in
value for that row (taken from the printed matrices) is used instead, and a
notice is printed, so the script always runs.
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =====================================================================
# CONFIG
# =====================================================================
CSV_DIR = str(config.output_dir(5))   # where Figure5_panelB_spearman_<tag>.csv were written
CSV_PATTERN = "Figure5_panelB_spearman_{tag}.csv"

# (row label on the figure, tag used in the CSV file name)
CLASSES = [
    ("Pooled (all classes)", "all"),      # change "all" to match your pooled CSV name
    ("CDS", "CDS"),
    ("UTR", "UTR"),
    ("Intron", "INT"),
    ("SINE", "SINE"),
    ("LINE", "LINE"),
    ("LTR", "LTR"),
    ("DNA transposon", "DNA"),
]

# Columns: (row condition, column condition) looked up in each matrix
COLUMNS = [
    ("+M", "-RC"),   # opposite-sign reverse complement
    ("+M", "+RC"),   # same-sign reverse complement
    ("+M", "+C"),    # plain complement
    ("+M", "+R"),    # plain reverse
    ("UM", "0M"),    # unannotated vs combined
    ("UM", "+M"),    # unannotated vs +annotation
    ("UM", "WM"),    # unannotated vs whole sequence
]
GROUP_SIZES = [1, 3, 3]
GROUP_TITLES = ["Opposite-sign\nreverse complement",
                "Same-sign reverse complement\nand plain transforms",
                "Unannotated sequence\nvs. other conditions"]

USE_K_LABELS = True     # show native transformation M as K, as in Figure 5

SAVE_NAME = str(config.output_dir(5) / "Supp_Fig5_class_summary_heatmap")
SAVE_DPI = 600
CMAP = "RdBu_r"         # standard diverging map so moderate values stay visible

# Built-in fallback values, same order as COLUMNS
BUILTIN = {
    "Pooled (all classes)": [0.9999, 0.8018, -0.9826, -0.8171, 0.9858, 0.9322, 0.9943],
    "CDS":                  [0.9999, 0.6251, -0.9011, -0.7179, 0.9383, 0.8420, 0.9998],
    "UTR":                  [0.9998, 0.8440, -0.9844, -0.8568, 0.9916, 0.9494, 1.0000],
    "Intron":               [0.9999, 0.7750, -0.9902, -0.7824, 0.9939, 0.9314, 0.9992],
    "SINE":                 [1.0000, 0.5311, -0.9127, -0.6111, 0.8710, 0.7795, 0.9957],
    "LINE":                 [1.0000, 0.3460, -0.9594, -0.3755, 0.9881, 0.7828, 0.9980],
    "LTR":                  [0.9999, 0.7355, -0.9295, -0.7966, 0.9765, 0.9065, 0.9997],
    "DNA transposon":       [0.9999, 0.9552, -0.9964, -0.9584, 0.9939, 0.9824, 1.0000],
}


def disp(label):
    return label.replace("M", "K") if USE_K_LABELS else label


def get_row(label, tag):
    path = os.path.join(CSV_DIR, CSV_PATTERN.format(tag=tag))
    if os.path.exists(path):
        df = pd.read_csv(path, index_col=0)
        df.index = [str(i).replace("K", "M") for i in df.index]      # accept K or M labels
        df.columns = [str(c).replace("K", "M") for c in df.columns]
        print(f"[csv]      {label:<22} <- {path}")
        return [float(df.loc[a, b]) for a, b in COLUMNS]
    print(f"[built-in] {label:<22} (no file {path})")
    return BUILTIN[label]


def build_table():
    rows = {label: get_row(label, tag) for label, tag in CLASSES}
    col_names = [f"{disp(a)} vs {disp(b)}" for a, b in COLUMNS]
    return pd.DataFrame(rows, index=col_names).T


def draw(table):
    data = table.values
    n_rows, n_cols = data.shape

    fig, ax = plt.subplots(figsize=(11, 6.4))
    im = ax.imshow(data, cmap=CMAP, vmin=-1, vmax=1, aspect="auto")

    # cell values
    for i in range(n_rows):
        for j in range(n_cols):
            v = data[i, j]
            ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=10,
                    color="white" if abs(v) > 0.7 else "black")

    # thin white cell borders
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", length=0)

    # group separators (columns) and pooled/classes separator (rows)
    edge = -0.5
    for size in GROUP_SIZES[:-1]:
        edge += size
        ax.axvline(edge, color="black", linewidth=2.0)
    ax.axhline(0.5, color="black", linewidth=2.0)

    # labels
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([c.replace(" vs ", "\nvs ") for c in table.columns], fontsize=10)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(table.index, fontsize=11)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # group titles above the columns
    start = 0
    for size, title in zip(GROUP_SIZES, GROUP_TITLES):
        center = start + (size - 1) / 2
        ax.text(center, -0.75, title, ha="center", va="bottom", fontsize=11,
                fontweight="bold")
        start += size

    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label("Spearman's \u03c1", fontsize=11)
    cb.outline.set_visible(False)

    fig.tight_layout()
    return fig


def main():
    table = build_table()
    print("\n" + table.round(4).to_string() + "\n")
    table.round(4).to_csv(f"{SAVE_NAME}_values.csv")

    fig = draw(table)
    fig.savefig(f"{SAVE_NAME}.png", dpi=SAVE_DPI, facecolor="white")
    fig.savefig(f"{SAVE_NAME}.pdf", facecolor="white")
    fig.savefig(f"{SAVE_NAME}.tiff", dpi=SAVE_DPI, format="tiff", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.show()


if __name__ == "__main__":
    main()