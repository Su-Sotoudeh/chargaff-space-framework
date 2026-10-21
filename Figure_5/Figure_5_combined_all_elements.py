"""
Figure 5 — Combined Panel A/B/C/D generator
============================================
Panel A (top-left)   : 20 per-condition k-mer heatmaps (5 sign-rows x 4
                        transformation-columns), grid forced to an exact square.
Panel B (top-right)  : 20x20 Spearman correlation matrix, forced to the SAME
                        exact square size as Panel A.
Panel C (middle)     : Mean +K/-K/+RC/-RC/+C/-C/+R/-R percentage per k-mer
                        (256 k-mers), row-gapped strip heatmap.
Panel D (bottom)     : 7 repeat-element rows (CDS, UTR, INT, SINE, LINE, LTR,
                        DNA) showing mean +K percentage per k-mer (one
                        independent colorbar per row). X-order matches Panel C.

Display convention: native/measured transformation "M" is labeled "K",
and the zero-sign "0" is labeled "N" on the figure.
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import spearmanr

# =====================================================================
# CONFIG — data
# =====================================================================

_OUT5 = config.output_dir(5)          # Figure_5/outputs (written by count_from_masked_file.py / Figure_5_<class>.py)
FILE_PATH_AB = str(_OUT5 / '4k_kmers_normalised_MCRRC_PC0UW.csv')
FILE_PATH_C  = FILE_PATH_AB   # same file, loaded with different normalisation

# ---- Panel D inputs -------------------------------------------------
MASTER_FILE_D = str(_OUT5 / '4k_kmers_by_bumber_MCRRC_PC0UW.csv')
ELEMENT_FILES_D = {
    "+M-CDS":  '4k_kmers_normalised_MCRRC_PC0UW_CDS.csv',
    "+M-UTR":  '4k_kmers_normalised_MCRRC_PC0UW_UTR.csv',
    "+M-INT":  '4k_kmers_normalised_MCRRC_PC0UW_INT.csv',
    "+M-SINE": '4k_kmers_normalised_MCRRC_PC0UW_SINE.csv',
    "+M-LINE": '4k_kmers_normalised_MCRRC_PC0UW_LINE.csv',
    "+M-LTR":  '4k_kmers_normalised_MCRRC_PC0UW_LTR.csv',
    "+M-DNA":  '4k_kmers_normalised_MCRRC_PC0UW_DNA.csv',
}
ELEMENT_FILES_D = {k: str(_OUT5 / v) for k, v in ELEMENT_FILES_D.items()}
PANEL_D_ROW_LABELS = ["+K-CDS", "+K-UTR", "+K-INT",
                      "+K-SINE", "+K-LINE", "+K-LTR", "+K-DNA"]
PANEL_D_ROWS_USED  = None   # None -> use all rows of every D file

CONDITIONS = ["+M", "+C", "+R", "+RC", "-M", "-C", "-R", "-RC",
              "0M", "0C", "0R", "0RC", "UM", "UC", "UR", "URC",
              "WM", "WC", "WR", "WRC"]

ORDERED_CONDITIONS_B = [
    "+M", "-M", "0M", "UM", "WM",
    "+C", "-C", "0C", "UC", "WC",
    "+R", "-R", "0R", "UR", "WR",
    "+RC", "-RC", "0RC", "URC", "WRC",
]

ROW_PREFIXES_A = ["+", "-", "0", "U", "W"]
COL_SUFFIXES_A = ["M", "C", "R", "RC"]
ROW_DISPLAY_A  = ["+", "\u2212", "N", "U", "W"]   # '0' shown as 'N'
COL_DISPLAY_A  = ["K", "C", "R", "RC"]            # 'M' shown as 'K'


def disp(label):
    """Display label: M -> K, and zero-sign '0' -> 'N' (e.g. '0M' -> 'NK')."""
    return label.replace("M", "K").replace("0", "N")


PANEL_C_ROWS_USED   = range(0, 2443)   # use None to keep all rows after filter
PANEL_C_ROW_ORDER   = ["+M", "-RC", "+RC", "-M", "+C", "-R", "+R", "-C"]
PANEL_C_ROW_COLORS  = ["red", "blue", "red", "blue", "red", "blue", "red", "blue"]
PANEL_C_GAP_HEIGHTS = [0.04, 0.35, 0.04, 0.35, 0.04, 0.35, 0.04]

SAVE_B_CSV       = str(_OUT5 / "Figure5_panelB_spearman.csv")
CSV_USE_K_LABELS = False

# =====================================================================
# CONFIG — layout (all sizes in inches)
# =====================================================================
SQUARE_SIDE_IN    = 9.0
MARGIN_LEFT_IN    = 1.3
GAP_A_TO_CBAR_IN  = 0.35
CBAR_WIDTH_IN     = 0.30
GAP_CBAR_TO_B_IN  = 2.5
GAP_B_TO_CBAR_IN  = 0.40
MARGIN_RIGHT_IN   = 1.8
MARGIN_TOP_IN     = 1.6
GAP_AB_TO_C_IN    = 1.7
PANEL_C_HEIGHT_IN = 3.4
GAP_C_TO_D_IN     = 1.1      # slightly more room for Panel C x-label
PANEL_D_HEIGHT_IN = 3.4      # taller so the 7 colorbars have breathing room
MARGIN_BOTTOM_IN  = 0.9

# ---- Species to exclude from every input ----------------------------
EXCEPTION_SPECIES = ['C_elegans', 'Drosophila', 'Gallus_Gallus',
                     'Zebrafish', 'Turtle']


def compute_layout():
    s = SQUARE_SIDE_IN
    fig_w = (MARGIN_LEFT_IN + s + GAP_A_TO_CBAR_IN + CBAR_WIDTH_IN
             + GAP_CBAR_TO_B_IN + s + GAP_B_TO_CBAR_IN + CBAR_WIDTH_IN
             + MARGIN_RIGHT_IN)
    fig_h = (MARGIN_TOP_IN + s + GAP_AB_TO_C_IN + PANEL_C_HEIGHT_IN
             + GAP_C_TO_D_IN + PANEL_D_HEIGHT_IN
             + MARGIN_BOTTOM_IN)

    def rect(x_in, y_from_top_in, w_in, h_in):
        return (x_in / fig_w, 1 - (y_from_top_in + h_in) / fig_h,
                w_in / fig_w, h_in / fig_h)

    x_A     = MARGIN_LEFT_IN
    x_cbarA = x_A + s + GAP_A_TO_CBAR_IN
    x_B     = x_cbarA + CBAR_WIDTH_IN + GAP_CBAR_TO_B_IN
    x_cbarB = x_B + s + GAP_B_TO_CBAR_IN

    y_top = MARGIN_TOP_IN
    y_C   = MARGIN_TOP_IN + s + GAP_AB_TO_C_IN
    y_D   = y_C + PANEL_C_HEIGHT_IN + GAP_C_TO_D_IN

    panel_width = x_B + s - x_A       # same horizontal extent as C

    return dict(
        fig_w=fig_w, fig_h=fig_h,
        A=rect(x_A, y_top, s, s),
        cbarA=rect(x_cbarA, y_top, CBAR_WIDTH_IN, s),
        B=rect(x_B, y_top, s, s),
        cbarB=rect(x_cbarB, y_top, CBAR_WIDTH_IN, s),
        C=rect(x_A, y_C, panel_width, PANEL_C_HEIGHT_IN),
        cbarC=rect(x_cbarB, y_C, CBAR_WIDTH_IN, PANEL_C_HEIGHT_IN),
        D=rect(x_A, y_D, panel_width, PANEL_D_HEIGHT_IN),
        x_cbarD_in=x_cbarB,
        letters=dict(
            A=(0.1 / fig_w,            1 - (y_top - 0.75) / fig_h),
            B=((x_B - 1.5) / fig_w,    1 - (y_top - 0.75) / fig_h),
            C=(0.1 / fig_w,            1 - (y_C   - 0.75) / fig_h),
            D=(0.1 / fig_w,            1 - (y_D   - 0.75) / fig_h),
        ),
    )


# =====================================================================
# STEP 1 — SHARED LOADING
# =====================================================================

def load_normalized(file_path):
    """Panels A & B: within-group percentage per chromosome."""
    df = pd.read_csv(file_path, index_col=0)
    n_cols = df.shape[1]
    for i in range(0, n_cols, 20):
        block_cols = df.columns[i:i + 20]
        for g in range(5):
            cols = block_cols[g * 4:(g + 1) * 4]
            group_sum = df.loc[:, cols].sum(axis=1)
            df.loc[:, cols] = df.loc[:, cols].div(group_sum, axis=0) * 100
    return _drop_exceptions_and_reset(df)


def load_raw_counts(file_path):
    """Panel C / Panel D: raw counts, no per-row normalisation."""
    df = pd.read_csv(file_path, index_col=0)
    return _drop_exceptions_and_reset(df)


def _drop_exceptions_and_reset(df):
    """Drop any chromosome whose index starts with an exception species."""
    mask = ~df.index.astype(str).str.startswith(tuple(EXCEPTION_SPECIES))
    return df[mask].reset_index(drop=True)


# =====================================================================
# STEP 2 — SHARED K-MER SORT/BUILD PIPELINE (used by A, B, C and D)
# =====================================================================

def build_condition_matrices(df, conditions=CONDITIONS):
    cols = df.columns.tolist()
    parts = str(cols[0]).split('_')
    separator = '__' if len(parts) >= 3 else '_'

    unique_kmers = sorted(set(c.split(separator)[0] for c in cols))

    chunks = {cond: [] for cond in conditions}
    kmer_scores, kmers_used = [], []

    for kmer in unique_kmers:
        kmer_cols = [f"{kmer}{separator}{cond}" for cond in conditions]
        if not all(col in df.columns for col in kmer_cols):
            continue

        sort_col = f"{kmer}{separator}+M"
        df_chunk = (df[kmer_cols]
                    .sort_values(by=sort_col, ascending=True)
                    .reset_index(drop=True))
        kmer_scores.append(df_chunk[sort_col].mean())
        kmers_used.append(kmer)

        for cond in conditions:
            chunks[cond].append(df_chunk[[f"{kmer}{separator}{cond}"]])

    order = np.argsort(kmer_scores, kind="stable")
    kmer_order_sorted = [kmers_used[i] for i in order]

    matrices = {}
    for cond in conditions:
        mat = pd.concat([chunks[cond][i] for i in order], axis=1)
        mat.columns = kmer_order_sorted
        matrices[cond] = mat
    return matrices


# =====================================================================
# STEP 3 — PANEL B
# =====================================================================

def compute_spearman_matrix(matrices, ordered_conditions=ORDERED_CONDITIONS_B):
    def flat_corr(a, b):
        r, _ = spearmanr(a.values.flatten(), b.values.flatten())
        return r

    n = len(ordered_conditions)
    sim = np.zeros((n, n))
    for i, ci in enumerate(ordered_conditions):
        for j, cj in enumerate(ordered_conditions):
            sim[i, j] = flat_corr(matrices[ci], matrices[cj])
    return pd.DataFrame(sim, index=ordered_conditions, columns=ordered_conditions)


def report_and_save_spearman(sim_df, csv_path=SAVE_B_CSV):
    out = sim_df.copy()
    if CSV_USE_K_LABELS:
        # CSV: M -> K only, keep 0 as 0 for downstream compatibility
        out.index   = [c.replace("M", "K") for c in out.index]
        out.columns = [c.replace("M", "K") for c in out.columns]

    with pd.option_context("display.width", 300,
                           "display.max_columns", None,
                           "display.max_rows", None,
                           "display.float_format", "{:.4f}".format):
        print("\n" + "=" * 70)
        print("PANEL B — Spearman correlation matrix (20 x 20)")
        print("=" * 70)
        print(out.to_string())
        print("=" * 70)

    out.to_csv(csv_path, float_format="%.6f")
    print(f"Panel B matrix saved to: {csv_path}\n")


def build_custom_diverging_cmap():
    colors = [
        (0.00, "#0000A8"), (0.06, "#4169E1"), (0.20, "#D6E4FF"),
        (0.20, "#FFFFFF"), (0.80, "#FFFFFF"),
        (0.80, "#FFE0E0"), (0.94, "#E94B4B"), (1.00, "#B00000"),
    ]
    return LinearSegmentedColormap.from_list("Chargaff_Spearman_Custom", colors)


# =====================================================================
# STEP 4 — PANEL C
# =====================================================================

def compute_panel_C_heatmap(raw_matrices, rows=PANEL_C_ROWS_USED,
                            row_order=PANEL_C_ROW_ORDER,
                            return_order=False):
    sl = slice(None) if rows is None else slice(rows.start, rows.stop)

    mpM, mmM   = raw_matrices["+M"].iloc[sl].mean(),  raw_matrices["-M"].iloc[sl].mean()
    mpRC, mmRC = raw_matrices["+RC"].iloc[sl].mean(), raw_matrices["-RC"].iloc[sl].mean()
    mpC, mmC   = raw_matrices["+C"].iloc[sl].mean(),  raw_matrices["-C"].iloc[sl].mean()
    mpR, mmR   = raw_matrices["+R"].iloc[sl].mean(),  raw_matrices["-R"].iloc[sl].mean()

    total_MRC = mpM + mmM + mpRC + mmRC
    total_CR  = mpC + mmC + mpR  + mmR

    perc = {
        "+M":  mpM  / total_MRC * 100, "-M":  mmM  / total_MRC * 100,
        "+RC": mpRC / total_MRC * 100, "-RC": mmRC / total_MRC * 100,
        "+C":  mpC  / total_CR  * 100, "-C":  mmC  / total_CR  * 100,
        "+R":  mpR  / total_CR  * 100, "-R":  mmR  / total_CR  * 100,
    }

    sorted_idx = perc["+M"].sort_values(ascending=False).index
    for key in perc:
        perc[key] = perc[key].loc[sorted_idx]

    heatmap = np.vstack([perc[cond].values for cond in row_order])
    if return_order:
        return heatmap, sorted_idx
    return heatmap


# =====================================================================
# STEP 4b — PANEL D
# =====================================================================

def compute_panel_D_heatmap(element_matrices_list,
                            kmer_order=None,
                            master_matrices=None,
                            rows=PANEL_D_ROWS_USED):
    """
    Compute the 7 x N heatmap for Panel D.
    Each row = one repeat element, value = mean +M percentage relative to
    (+M + -M + +RC + -RC), per k-mer.

    kmer_order: shared x-ordering (usually Panel C's) to force column alignment.
    If None, ordering comes from master_matrices (its +M percentage, descending).
    """
    sl = slice(None) if rows is None else slice(rows.start, rows.stop)

    def plusM_perc(matrices):
        pM  = matrices["+M"].iloc[sl].mean()
        mM  = matrices["-M"].iloc[sl].mean()
        pRC = matrices["+RC"].iloc[sl].mean()
        mRC = matrices["-RC"].iloc[sl].mean()
        total = pM + mM + pRC + mRC
        return (pM / total) * 100

    if kmer_order is None:
        if master_matrices is None:
            raise ValueError("Need kmer_order or master_matrices.")
        kmer_order = plusM_perc(master_matrices).sort_values(ascending=False).index

    rows_out = []
    for em in element_matrices_list:
        perc = plusM_perc(em).reindex(kmer_order)
        rows_out.append(perc.values)

    return np.vstack(rows_out), kmer_order


# =====================================================================
# STEP 5 — PANEL DRAWING (A, B, C)
# =====================================================================

def draw_panel_A(fig, rect, cbar_rect, matrices,
                 row_prefixes=ROW_PREFIXES_A, col_suffixes=COL_SUFFIXES_A,
                 row_display=ROW_DISPLAY_A, col_display=COL_DISPLAY_A,
                 cmap="viridis", wspace=0.04, hspace=0.06):
    left, bottom, width, height = rect
    gs = gridspec.GridSpec(len(row_prefixes), len(col_suffixes),
                           left=left, right=left + width,
                           bottom=bottom, top=bottom + height,
                           wspace=wspace, hspace=hspace, figure=fig)

    vmin = min(m.values.min() for m in matrices.values())
    vmax = max(m.values.max() for m in matrices.values())

    im = None
    for r, rp in enumerate(row_prefixes):
        for c, cs in enumerate(col_suffixes):
            ax = fig.add_subplot(gs[r, c])
            data = matrices[f"{rp}{cs}"].values
            im = ax.imshow(data, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(True); spine.set_color("black")
                spine.set_linewidth(1.2)
            if r == 0:
                ax.set_title(col_display[c], fontsize=15,
                             fontweight="bold", pad=6)
            if c == 0:
                ax.set_ylabel(row_display[r], fontsize=15, fontweight="bold",
                              rotation=0, labelpad=14, va="center")
            if r == len(row_prefixes) - 1:
                ax.set_xlabel("k-mers", fontsize=9)

    cax = fig.add_axes(cbar_rect)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Share of 4-mer type within condition (%)", fontsize=12)


def draw_panel_B(fig, rect, cbar_rect, sim_df,
                 ordered_conditions=ORDERED_CONDITIONS_B):
    sim_plot = sim_df.loc[ordered_conditions[::-1], ordered_conditions]
    ax  = fig.add_axes(rect)
    cax = fig.add_axes(cbar_rect)

    cmap = build_custom_diverging_cmap()
    norm = Normalize(vmin=-1, vmax=1)

    mask_annot = np.zeros(sim_plot.shape, dtype=bool)
    mask_annot[15:20, 0:5] = True
    annot_values = np.where(mask_annot,
                            np.vectorize(lambda v: f"{v:.2f}")(sim_plot.values),
                            "")

    sns.heatmap(sim_plot, cmap=cmap, norm=norm, linecolor="black",
                ax=ax, cbar_ax=cax, annot=annot_values, fmt="",
                annot_kws={"size": 9}, square=False,
                cbar_kws={"label": "Spearman's \u03c1"})
    cax.yaxis.label.set_size(11)

    n = len(ordered_conditions)
    for i in range(n + 1):
        ax.axhline(i, color="black", linewidth=0.4, alpha=0.2)
        ax.axvline(i, color="black", linewidth=0.4, alpha=0.2)
    for i in range(0, n + 1, 5):
        ax.axhline(i, color="black", linewidth=2.2)
        ax.axvline(i, color="black", linewidth=2.2)

    ax.plot([2, 2], [15, 20], color="black", linewidth=1.4)
    ax.plot([0, 5], [18, 18], color="black", linewidth=1.4)

    ax.set_xticklabels([disp(c) for c in sim_plot.columns],
                       rotation=0, fontsize=12)
    ax.set_yticklabels([disp(c) for c in sim_plot.index],
                       rotation=0, fontsize=12)
    ax.set_xlabel("Condition (column)", fontsize=14)
    ax.set_ylabel("Condition (row)", fontsize=14)
    return ax


def draw_panel_C(fig, rect, cbar_rect, heatmap_data,
                 gap_heights=PANEL_C_GAP_HEIGHTS,
                 row_labels=PANEL_C_ROW_ORDER,
                 row_colors=PANEL_C_ROW_COLORS,
                 cmap_name="viridis"):
    n_rows, n_cols = heatmap_data.shape
    assert len(gap_heights) == n_rows - 1

    ax  = fig.add_axes(rect)
    cax = fig.add_axes(cbar_rect)

    row_heights, row_values = [], []
    for i in range(n_rows):
        row_heights.append(1.0)
        row_values.append(heatmap_data[i, :])
        if i != n_rows - 1:
            row_heights.append(gap_heights[i])
            row_values.append(np.full(n_cols, np.nan))

    y_edges = np.concatenate([[0], np.cumsum(row_heights)])
    x_edges = np.arange(n_cols + 1)

    masked = np.ma.masked_invalid(np.array(row_values))
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(color="white")

    mesh = ax.pcolormesh(x_edges, y_edges, masked, cmap=cmap,
                         shading="flat", edgecolors="none")
    ax.invert_yaxis()

    real_row_indices = list(range(0, len(row_heights), 2))
    real_row_centers = [(y_edges[i] + y_edges[i + 1]) / 2
                        for i in real_row_indices]
    ax.set_yticks(real_row_centers)
    ax.set_yticklabels([disp(r) for r in row_labels], fontsize=13)
    for label, color in zip(ax.get_yticklabels(), row_colors):
        label.set_color(color); label.set_fontweight("bold")

    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=0)

    cb = fig.colorbar(mesh, cax=cax)
    cb.outline.set_visible(False)
    cb.set_label("Mean percentage within\nK/RC or C/R group (%)", fontsize=12)
    return ax


# =====================================================================
# STEP 5b — PANEL D DRAWING (own per-row colorbars on the right)
# =====================================================================

def draw_panel_D(fig, rect, heatmap_data, row_labels,
                 fig_w, fig_h, x_cbar_in,
                 cbar_pad_in=0.10, cbar_width_in=0.20,
                 cbar_label_pad_in=0.45,
                 row_cbar_shrink=0.60,
                 gap_heights=None, cmap_name="viridis"):
    n_rows, n_cols = heatmap_data.shape
    if gap_heights is None:
        gap_heights = [0.02] * (n_rows - 1)
    assert len(gap_heights) == n_rows - 1

    ax = fig.add_axes(rect)

    # --- Build layout (rows + gaps) --------------------------------
    row_heights, row_values, row_start = [], [], []
    entry_idx = 0
    for i in range(n_rows):
        row_start.append(entry_idx)
        row_heights.append(1.0)
        row_values.append(heatmap_data[i, :])
        entry_idx += 1
        if i != n_rows - 1:
            row_heights.append(gap_heights[i])
            row_values.append(np.full(n_cols, np.nan))
            entry_idx += 1

    masked = np.ma.masked_invalid(np.array(row_values))
    x_edges = np.arange(n_cols + 1)
    y_edges = np.concatenate([[0], np.cumsum(row_heights)])
    y_total = y_edges[-1]

    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(color="white")

    element_meshes = {}
    for i in range(n_rows):
        y_top = y_edges[row_start[i]]
        y_bot = y_edges[row_start[i] + 1]
        row_masked = masked[row_start[i]:row_start[i] + 1, :]
        mesh_i = ax.pcolormesh(x_edges, [y_top, y_bot], row_masked,
                               cmap=cmap, shading="flat", edgecolors="none")
        element_meshes[i] = mesh_i

    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])
    ax.invert_yaxis()

    real_rows = list(range(0, len(row_heights), 2))
    real_centers = [(y_edges[i] + y_edges[i + 1]) / 2 for i in real_rows]
    ax.set_yticks(real_centers)
    ax.set_yticklabels(row_labels, fontsize=13)
    for lbl in ax.get_yticklabels():
        lbl.set_color("black"); lbl.set_fontweight("bold")

    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=0)

    # --- Per-row colorbars, aligned with cbarB / cbarC --------------
    left, bottom, width, height = rect
    x_cbar_frac = x_cbar_in / fig_w
    cbar_w_frac = cbar_width_in / fig_w

    for i in range(n_rows):
        y_top_d = y_edges[row_start[i]]
        y_bot_d = y_edges[row_start[i] + 1]
        fig_top = bottom + height * (1 - y_top_d / y_total)
        fig_bot = bottom + height * (1 - y_bot_d / y_total)
        slot_h  = fig_top - fig_bot
        c_h     = slot_h * row_cbar_shrink
        margin  = (slot_h - c_h) / 2

        cax_i = fig.add_axes([x_cbar_frac, fig_bot + margin,
                              cbar_w_frac, c_h])
        cbar_i = fig.colorbar(element_meshes[i], cax=cax_i)
        cbar_i.outline.set_visible(False)
        cbar_i.ax.tick_params(labelsize=6, pad=2)

        vmin = element_meshes[i].norm.vmin
        vmax = element_meshes[i].norm.vmax
        cbar_i.set_ticks([vmin, vmax])
        cbar_i.set_ticklabels([f"{vmin:.1f}", f"{vmax:.1f}"])

    # --- Vertical axis label to the right of the colorbars ----------
    label_x = (x_cbar_in + cbar_width_in + cbar_label_pad_in) / fig_w
    label_y = bottom + height / 2
    fig.text(label_x, label_y,
             "Percentage (%) of +K on 4-type transition",
             rotation=90, va="center", ha="center", fontsize=11)

    return ax


# =====================================================================
# STEP 6 — ASSEMBLE THE FULL FIGURE
# =====================================================================

def main():
    # --- Panels A & B -----------------------------------------------
    df_ab       = load_normalized(FILE_PATH_AB)
    print("Rows kept after exception filter:", len(df_ab))
    print("Any exception species still present?",
          bool(df_ab.index.astype(str)
               .str.startswith(tuple(EXCEPTION_SPECIES)).any()))

    matrices_ab = build_condition_matrices(df_ab)
    sim_df      = compute_spearman_matrix(matrices_ab)
    report_and_save_spearman(sim_df)

    # --- Panel C (and its k-mer x-order, shared with D) -------------
    df_c        = load_raw_counts(FILE_PATH_C)
    matrices_c  = build_condition_matrices(df_c)
    heatmap_C, sorted_idx_C = compute_panel_C_heatmap(matrices_c,
                                                      return_order=True)

    # --- Panel D ----------------------------------------------------
    df_master_d       = load_raw_counts(MASTER_FILE_D)
    matrices_master_d = build_condition_matrices(df_master_d)

    element_matrices_list = []
    for label, path in ELEMENT_FILES_D.items():
        df_e = load_raw_counts(path)
        element_matrices_list.append(build_condition_matrices(df_e))

    heatmap_D, _ = compute_panel_D_heatmap(
        element_matrices_list,
        kmer_order=sorted_idx_C,          # lock x-order to Panel C
        master_matrices=matrices_master_d,
        rows=PANEL_D_ROWS_USED,
    )

    # --- Figure -----------------------------------------------------
    L   = compute_layout()
    fig = plt.figure(figsize=(L["fig_w"], L["fig_h"]))

    draw_panel_A(fig, L["A"], L["cbarA"], matrices_ab)
    draw_panel_B(fig, L["B"], L["cbarB"], sim_df)
    axC = draw_panel_C(fig, L["C"], L["cbarC"], heatmap_C)
    axD = draw_panel_D(fig, L["D"], heatmap_D, PANEL_D_ROW_LABELS,
                       fig_w=L["fig_w"], fig_h=L["fig_h"],
                       x_cbar_in=L["x_cbarD_in"])

    # Panel C x-label — handles PANEL_C_ROWS_USED = None gracefully
    if PANEL_C_ROWS_USED is None:
        row_txt = "all chromosomes"
    else:
        row_txt = f"first {PANEL_C_ROWS_USED.stop} chromosomes"

    axC.set_xlabel(
        f"{heatmap_C.shape[1]} 4-mers, ranked by mean +K percentage "
        f"(descending, averaged over {row_txt})",
        fontsize=12,
    )

    # Panel D x-label
    axD.set_xlabel("256 4-mers (same x-order as Panel C)", fontsize=12)

    # Panel letters
    for name, (x, y) in L["letters"].items():
        fig.text(x, y, name, fontsize=24, fontweight="bold")

    # =================================================================
    # SAVE — journal-ready outputs (no bbox_inches="tight")
    # =================================================================
    BASE_NAME = str(_OUT5 / "Figure5_combined_panels_ABCD")

    fig.savefig(f"{BASE_NAME}.png", dpi=300,
                facecolor="white", edgecolor="none")

    fig.savefig(f"{BASE_NAME}.tiff", dpi=300, format="tiff",
                facecolor="white", edgecolor="none",
                pil_kwargs={"compression": "tiff_lzw"})

    fig.savefig(f"{BASE_NAME}_600dpi.tiff", dpi=600, format="tiff",
                facecolor="white", edgecolor="none",
                pil_kwargs={"compression": "tiff_lzw"})

    fig.savefig(f"{BASE_NAME}.pdf", format="pdf",
                facecolor="white", edgecolor="none")

    fig.savefig(f"{BASE_NAME}.svg", format="svg",
                facecolor="white", edgecolor="none")

    print("Saved: "
          f"{BASE_NAME}.png | .tiff | _600dpi.tiff | .pdf | .svg")

    plt.show()


if __name__ == "__main__":
    main()