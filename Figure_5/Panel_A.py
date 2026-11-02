
"""
Figure 5 — Panel A ONLY
========================
20 per-condition k-mer heatmaps (5 sign-rows: +, -, 0, U, W;
4 type-columns: M, C, R, RC), sharing one color scale with a single
colorbar.

Color scaling:
- Shared global percentile-based scale across all 20 heatmaps.
- Default upper limit: 98th percentile.
- Values above the limit are saturated at the maximum color.
- Original normalized data remain unchanged.
"""

import os
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# =====================================================================
# CONFIG — DATA
# =====================================================================

FILE_PATH = str(config.figure_input(5, '4k_kmers_normalised_MCRRC_PC0UW.csv'))
# Figure files are written into Figure_5/outputs/ (the inputs above are absolute paths)
os.chdir(config.output_dir(5))


CONDITIONS = [
    "+M", "+C", "+R", "+RC",
    "-M", "-C", "-R", "-RC",
    "0M", "0C", "0R", "0RC",
    "UM", "UC", "UR", "URC",
    "WM", "WC", "WR", "WRC"
]

ROW_PREFIXES = ["+", "-", "0", "U", "W"]
COL_SUFFIXES = ["M", "C", "R", "RC"]

ROW_DISPLAY = ["+", "\u2212", "0", "U", "W"]
COL_DISPLAY = ["M", "C", "R", "RC"]


# =====================================================================
# CONFIG — COLOR SCALING
# =====================================================================

COLOR_VMIN = 0

# Try 95, 98, 99, or 100
COLOR_VMAX_PERCENTILE = 99.93

COLOR_CMAP = "viridis"


# =====================================================================
# CONFIG — LAYOUT (INCHES)
# =====================================================================

GRID_SIDE_IN = 16.0

LEFT_MARGIN_IN = 0.8
GAP_BEFORE_CBAR_IN = 0.35
CBAR_WIDTH_IN = 0.35
RIGHT_MARGIN_IN = 1.3

TOP_MARGIN_IN = 0.9
BOTTOM_MARGIN_IN = 0.7

WSPACE = 0.03
HSPACE = 0.06

COLUMN_LABEL_SIZE = 20
ROW_LABEL_SIZE = 20
SUPTITLE_SIZE = 18
XLABEL_SIZE = 13
CBAR_LABEL_SIZE = 13

PANEL_BORDER_WIDTH = 1.6

SAVE_NAME = 'Fig_5_panel_A_4mer-pattern'
SAVE_DPI = 600


# =====================================================================
# STEP 1 — LOAD & NORMALIZE
# =====================================================================

def load_normalized(file_path):

    df = pd.read_csv(file_path, index_col=0)

    n_cols = df.shape[1]

    for i in range(0, n_cols, 20):

        block_cols = df.columns[i:i + 20]

        for g in range(5):

            cols = block_cols[g * 4:(g + 1) * 4]

            group_sum = df.loc[:, cols].sum(axis=1)

            df.loc[:, cols] = (
                    df.loc[:, cols]
                    .div(group_sum, axis=0)
                    * 100
            )

    # Exclude Drosophila
    df = df[
        ~df.index.astype(str).str.contains(
            "Drosophila",
            case=False
        )
    ]

    return df.reset_index(drop=True)


# =====================================================================
# STEP 2 — SORT / BUILD PIPELINE
# =====================================================================

def build_condition_matrices(df, conditions=CONDITIONS):

    cols = df.columns.tolist()

    parts = str(cols[0]).split('_')

    separator = '__' if len(parts) >= 3 else '_'

    unique_kmers = sorted(
        set(c.split(separator)[0] for c in cols)
    )

    chunks = {
        cond: [] for cond in conditions
    }

    kmer_scores = []

    for kmer in unique_kmers:

        kmer_cols = [
            f"{kmer}{separator}{cond}"
            for cond in conditions
        ]

        if not all(col in df.columns for col in kmer_cols):
            continue

        sort_col = f"{kmer}{separator}+M"

        df_chunk = (
            df[kmer_cols]
            .sort_values(
                by=sort_col,
                ascending=True
            )
            .reset_index(drop=True)
        )

        kmer_scores.append(
            df_chunk[sort_col].mean()
        )

        for cond in conditions:

            chunks[cond].append(
                df_chunk[[f"{kmer}{separator}{cond}"]]
            )

    if not kmer_scores:
        raise ValueError(
            "No complete k-mer blocks were found. "
            "Check the column names and separator."
        )

    order = np.argsort(
        kmer_scores,
        kind="stable"
    )

    matrices = {
        cond: pd.concat(
            [chunks[cond][i] for i in order],
            axis=1
        )
        for cond in conditions
    }

    return matrices


# =====================================================================
# STEP 3 — LAYOUT
# =====================================================================

def compute_layout():

    fig_w = (
            LEFT_MARGIN_IN
            + GRID_SIDE_IN
            + GAP_BEFORE_CBAR_IN
            + CBAR_WIDTH_IN
            + RIGHT_MARGIN_IN
    )

    fig_h = (
            TOP_MARGIN_IN
            + GRID_SIDE_IN
            + BOTTOM_MARGIN_IN
    )

    grid_left = LEFT_MARGIN_IN / fig_w

    grid_right = (
                         LEFT_MARGIN_IN + GRID_SIDE_IN
                 ) / fig_w

    grid_top = 1 - TOP_MARGIN_IN / fig_h

    grid_bottom = (
            1 - (TOP_MARGIN_IN + GRID_SIDE_IN) / fig_h
    )

    cbar_left = (
            grid_right + GAP_BEFORE_CBAR_IN / fig_w
    )

    cbar_right = (
            cbar_left + CBAR_WIDTH_IN / fig_w
    )

    return dict(
        fig_w=fig_w,
        fig_h=fig_h,
        grid=(
            grid_left,
            grid_right,
            grid_bottom,
            grid_top
        ),
        cbar=(
            cbar_left,
            cbar_right,
            grid_bottom,
            grid_top
        )
    )


# =====================================================================
# STEP 4 — GLOBAL COLOR SCALE
# =====================================================================

def calculate_color_scale(matrices):

    all_values = np.concatenate([
        m.to_numpy(dtype=float).ravel()
        for m in matrices.values()
    ])

    # Ignore NaN and infinite values
    all_values = all_values[
        np.isfinite(all_values)
    ]

    if len(all_values) == 0:
        raise ValueError("No finite heatmap values found.")

    vmin = COLOR_VMIN

    vmax = np.percentile(
        all_values,
        COLOR_VMAX_PERCENTILE
    )

    if vmax <= vmin:
        raise ValueError(
            "Invalid color scale. Check COLOR_VMIN "
            "and the selected percentile."
        )

    clipped = np.sum(all_values > vmax)

    print("\n" + "=" * 50)
    print("GLOBAL HEATMAP COLOR SCALE")
    print("=" * 50)

    print(f"Total values: {len(all_values):,}")
    print(f"Minimum: {all_values.min():.5f}")
    print(f"Maximum: {all_values.max():.5f}")
    print(f"Selected percentile: {COLOR_VMAX_PERCENTILE}")
    print(f"Color minimum: {vmin:.5f}")
    print(f"Color maximum: {vmax:.5f}")

    print(
        f"Values saturated: {clipped:,} "
        f"({100 * clipped / len(all_values):.3f}%)"
    )

    print("=" * 50 + "\n")

    return vmin, vmax


# =====================================================================
# STEP 5 — DRAW FIGURE 5A
# =====================================================================

def draw_panel_A(matrices):

    layout = compute_layout()

    fig = plt.figure(
        figsize=(
            layout["fig_w"],
            layout["fig_h"]
        )
    )

    g_left, g_right, g_bottom, g_top = layout["grid"]

    gs = gridspec.GridSpec(
        len(ROW_PREFIXES),
        len(COL_SUFFIXES),
        left=g_left,
        right=g_right,
        bottom=g_bottom,
        top=g_top,
        wspace=WSPACE,
        hspace=HSPACE,
        figure=fig
    )

    # One shared global scale across all 20 heatmaps
    vmin, vmax = calculate_color_scale(matrices)

    im = None

    for r, rp in enumerate(ROW_PREFIXES):

        for c, cs in enumerate(COL_SUFFIXES):

            ax = fig.add_subplot(gs[r, c])

            data = matrices[f"{rp}{cs}"].to_numpy(
                dtype=float
            )

            im = ax.imshow(
                data,
                aspect="auto",
                cmap=COLOR_CMAP,
                vmin=vmin,
                vmax=vmax,
                interpolation="nearest"
            )

            ax.set_xticks([])
            ax.set_yticks([])

            # Panel borders
            for spine in ax.spines.values():

                spine.set_visible(True)
                spine.set_color("black")
                spine.set_linewidth(PANEL_BORDER_WIDTH)

            # Column labels
            if r == 0:

                ax.set_title(
                    COL_DISPLAY[c],
                    fontsize=COLUMN_LABEL_SIZE,
                    fontweight="bold",
                    pad=8
                )

            # Row labels
            if c == 0:

                ax.set_ylabel(
                    ROW_DISPLAY[r],
                    fontsize=ROW_LABEL_SIZE,
                    fontweight="bold",
                    rotation=0,
                    labelpad=15,
                    va="center"
                )

            else:

                ax.set_ylabel("")

    # -------------------------------------------------------------
    # SHARED COLORBAR
    # -------------------------------------------------------------

    c_left, c_right, c_bottom, c_top = layout["cbar"]

    gs_cbar = gridspec.GridSpec(
        1,
        1,
        left=c_left,
        right=c_right,
        bottom=c_bottom,
        top=c_top,
        figure=fig
    )

    cax = fig.add_subplot(gs_cbar[0, 0])

    cbar = fig.colorbar(
        im,
        cax=cax
    )

    cbar.set_label(
        "Share of 4-mer type within condition (%)",
        fontsize=CBAR_LABEL_SIZE
    )

    # -------------------------------------------------------------
    # X-AXIS DESCRIPTION
    # -------------------------------------------------------------

    fig.text(
        (g_left + g_right) / 2,
        g_bottom - 0.03,
        "256 possible 4-mers, ordered by mean +M composition "
        "(ascending, left \u2192 right)",
        ha="center",
        va="top",
        fontsize=XLABEL_SIZE
    )

    return fig


# =====================================================================
# STEP 6 — MAIN
# =====================================================================

def main():

    print("Loading data...")

    df = load_normalized(FILE_PATH)

    print("Data shape:", df.shape)

    print("Building condition matrices...")

    matrices = build_condition_matrices(df)

    print("Drawing Figure 5A...")

    fig = draw_panel_A(matrices)

    # TIFF
    fig.savefig(
        f'{SAVE_NAME}-P{COLOR_VMAX_PERCENTILE}-{SAVE_DPI}.tiff',
        dpi=SAVE_DPI,
        format='tiff',
        facecolor='white',
        edgecolor='none',
        pil_kwargs={"compression": "tiff_lzw"}
    )

    # PNG
    fig.savefig(
        f'{SAVE_NAME}-P{COLOR_VMAX_PERCENTILE}-{SAVE_DPI}.png',
        dpi=SAVE_DPI,
        format='png',
        facecolor='white',
        edgecolor='none'
    )

    print("\nFigure saved successfully.")

    plt.show()


if __name__ == "__main__":
    main()
