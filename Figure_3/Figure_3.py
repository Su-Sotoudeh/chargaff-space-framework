import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

# =====================================================================
# 1. CONFIG
# =====================================================================
FOLDERS_PATH = str(config.require(config.MAMMALS_CHR_DIR, "Put the per-chromosome tables in Dataframes/Mammals/<organism>/<chromosome>.csv"))
EXCLUDE_FOLDERS = None  # ["Mus"]
MIN_CHR_MB = 10
OUT_PREFIX = str(config.output_dir(3) / "Figure_3")

PLOT_STYLE = "bar"          # "bar" or "pointrange", for panels A/B
PORTION_STYLE = "bar"       # "bar" or "violin", for panel D (% of chromosome)

YLIM_A2T = (0.5, 1.65)
YLIM_G2C = (0.8, 1.2)
YLIM_PCT_ELEMENTS = None
YLIM_PCT_AGGREGATE = None
HEATMAP_VMAX = 0.1

STRAND_ELEMENTS = ["CDS", "UTR", "Intron", "LINE", "SINE", "LTR", "DNA"]
SUM_ELEMENT = "Annotation"
UNANNOTATED_LABEL = "Unannotated"

ALL_ELEMENTS_AB = STRAND_ELEMENTS + [SUM_ELEMENT, UNANNOTATED_LABEL]
CONDITIONS = ["combined", "plus", "minus"]

# ---------------------------------------------------------------------
# 2. LAYOUT
# ---------------------------------------------------------------------
FIGSIZE = (15, 9.5)

PANEL_RECTS = {
    "A": [0.055, 0.565, 0.42, 0.365],
    "B": [0.055, 0.085, 0.42, 0.365],
    "C": [0.460, 0.500, 0.400, 0.410],
    "D": [0.540, 0.085, 0.350, 0.320],
}

PORTION_SPLIT_RATIO = [1.55, 1.0]
PORTION_WSPACE = 0.018

ELEMENT_GROUP_GAP = 1.6
SECTION_GAP = 0.9

# ---------------------------------------------------------------------
# 3. COLOR PALETTE
# ---------------------------------------------------------------------
COLOR_COMBINED = "#7F8C8D"
COLOR_PLUS = "#E76F51"
COLOR_MINUS = "#2A9D8F"
COLOR_UNANNOTATED = "#7D5BA6"
EDGE_COLOR = "#3D3D3D"

# =====================================================================
# 2. PARSE ALL FILES
# =====================================================================
data = {el: {c: {"A2T": [], "G2C": []} for c in CONDITIONS} for el in ALL_ELEMENTS_AB}
pct_data = {el: {c: [] for c in CONDITIONS} for el in ALL_ELEMENTS_AB}

cross_dev = {(i, j): [] for i in STRAND_ELEMENTS + [None] for j in STRAND_ELEMENTS + [None]}

folders_seen = 0
folders_with_analyzed_chr = 0
chromosomes_analyzed = 0


def safe_div(a, b):
    try:
        return float(a) / float(b)
    except ZeroDivisionError:
        return 1.0


def deviation_score(A, T, C, G):
    """Chargaff-parity deviation: sqrt(AT_skew^2 + GC_skew^2). 0 = perfect parity."""
    at_denom = A + T
    gc_denom = G + C
    at_skew = (A - T) / at_denom if at_denom != 0 else 0.0
    gc_skew = (G - C) / gc_denom if gc_denom != 0 else 0.0
    return (at_skew ** 2 + gc_skew ** 2) ** 0.5


if EXCLUDE_FOLDERS is None:
    excluded = set()
else:
    excluded = set(EXCLUDE_FOLDERS)

folders = [f for f in os.listdir(FOLDERS_PATH) if f not in excluded]

for folder in folders:
    folder_dir = os.path.join(FOLDERS_PATH, folder)
    if not os.path.isdir(folder_dir):
        continue

    folders_seen += 1
    folder_has_analyzed_chr = False

    for chr_file in os.listdir(folder_dir):
        if not chr_file.endswith(".csv"):
            continue

        df = pd.read_csv(os.path.join(folder_dir, chr_file), nrows=18)
        df.columns.values[0] = "Element"
        df[["A", "T", "C", "G"]] = df[["A", "T", "C", "G"]].astype(float)

        nonseq = df.loc[16, "A"]
        df.loc[16, ["A", "T", "C", "G"]] = nonseq / 4

        chr_bp = int(df.iloc[17, 5])
        if chr_bp / 1_000_000 < MIN_CHR_MB:
            continue

        chromosomes_analyzed += 1
        if not folder_has_analyzed_chr:
            folders_with_analyzed_chr += 1
            folder_has_analyzed_chr = True

        plus_counts = {}
        minus_counts = {}

        sum_A_p = sum_T_p = sum_C_p = sum_G_p = sum_bp_p = 0.0
        sum_A_m = sum_T_m = sum_C_m = sum_G_m = sum_bp_m = 0.0

        for i, element in enumerate(STRAND_ELEMENTS):
            A_p, T_p, G_p, C_p, bp_p = df.iloc[i, 1], df.iloc[i, 2], df.iloc[i, 4], df.iloc[i, 3], df.iloc[i, 5]
            A_m, T_m, G_m, C_m, bp_m = df.iloc[7 + i, 1], df.iloc[7 + i, 2], df.iloc[7 + i, 4], df.iloc[7 + i, 3], df.iloc[7 + i, 5]

            plus_counts[element] = (A_p, T_p, C_p, G_p, bp_p)
            minus_counts[element] = (A_m, T_m, C_m, G_m, bp_m)

            # Panels A/B: ratios
            data[element]["combined"]["A2T"].append(safe_div(A_p + A_m, T_p + T_m))
            data[element]["combined"]["G2C"].append(safe_div(G_p + G_m, C_p + C_m))
            data[element]["plus"]["A2T"].append(safe_div(A_p, T_p))
            data[element]["plus"]["G2C"].append(safe_div(G_p, C_p))
            data[element]["minus"]["A2T"].append(safe_div(A_m, T_m))
            data[element]["minus"]["G2C"].append(safe_div(G_m, C_m))

            # Panel D: % of chromosome
            pct_data[element]["plus"].append(100 * bp_p / chr_bp)
            pct_data[element]["minus"].append(100 * bp_m / chr_bp)
            pct_data[element]["combined"].append(100 * (bp_p + bp_m) / chr_bp)

            sum_A_p += A_p
            sum_T_p += T_p
            sum_C_p += C_p
            sum_G_p += G_p
            sum_bp_p += bp_p

            sum_A_m += A_m
            sum_T_m += T_m
            sum_C_m += C_m
            sum_G_m += G_m
            sum_bp_m += bp_m

        # Annotation = sum of the 7 elements
        data[SUM_ELEMENT]["plus"]["A2T"].append(safe_div(sum_A_p, sum_T_p))
        data[SUM_ELEMENT]["plus"]["G2C"].append(safe_div(sum_G_p, sum_C_p))
        data[SUM_ELEMENT]["minus"]["A2T"].append(safe_div(sum_A_m, sum_T_m))
        data[SUM_ELEMENT]["minus"]["G2C"].append(safe_div(sum_G_m, sum_C_m))
        data[SUM_ELEMENT]["combined"]["A2T"].append(safe_div(sum_A_p + sum_A_m, sum_T_p + sum_T_m))
        data[SUM_ELEMENT]["combined"]["G2C"].append(safe_div(sum_G_p + sum_G_m, sum_C_p + sum_C_m))

        # Unannotated (combined only)
        A_u, T_u, C_u, G_u, bp_u = df.iloc[15, 1], df.iloc[15, 2], df.iloc[15, 3], df.iloc[15, 4], df.iloc[15, 5]
        data[UNANNOTATED_LABEL]["combined"]["A2T"].append(safe_div(A_u, T_u))
        data[UNANNOTATED_LABEL]["combined"]["G2C"].append(safe_div(G_u, C_u))
        pct_data[UNANNOTATED_LABEL]["combined"].append(100 * bp_u / chr_bp)

        # Shared (row 14): split 50/50 into annotation plus/minus totals
        bp_shared = df.iloc[14, 5]
        pct_data[SUM_ELEMENT]["plus"].append(100 * (sum_bp_p + bp_shared / 2) / chr_bp)
        pct_data[SUM_ELEMENT]["minus"].append(100 * (sum_bp_m + bp_shared / 2) / chr_bp)
        pct_data[SUM_ELEMENT]["combined"].append(100 * (sum_bp_p + sum_bp_m + bp_shared) / chr_bp)

        # Panel C: pairwise cross-strand Chargaff deviation
        for el_i in STRAND_ELEMENTS:
            A_p, T_p, C_p, G_p, _ = plus_counts[el_i]
            cross_dev[(el_i, None)].append(deviation_score(A_p, T_p, C_p, G_p))
            for el_j in STRAND_ELEMENTS:
                A_m, T_m, C_m, G_m, _ = minus_counts[el_j]
                cross_dev[(el_i, el_j)].append(
                    deviation_score(A_p + A_m, T_p + T_m, C_p + C_m, G_p + G_m)
                )

        for el_j in STRAND_ELEMENTS:
            A_m, T_m, C_m, G_m, _ = minus_counts[el_j]
            cross_dev[(None, el_j)].append(deviation_score(A_m, T_m, C_m, G_m))


# =====================================================================
# 3. SUMMARY STATISTICS
# =====================================================================
def mean_std(d, key1, key2, metric=None):
    vals = d[key1][key2][metric] if metric is not None else d[key1][key2]
    if len(vals) == 0:
        return None, None
    return float(np.mean(vals)), float(np.std(vals))


def clipped_yerr(mean, std):
    """Asymmetric error bar that never lets mean-std dip below 0."""
    lower = min(std, mean)
    upper = std
    return [[lower], [upper]]


def print_ratio_summary(source, metric_name, title):
    print(f"\n{title}")
    print("-" * len(title))
    for el in ALL_ELEMENTS_AB:
        conds = ["combined"] if el == UNANNOTATED_LABEL else CONDITIONS
        for cond in conds:
            mean, std = mean_std(source, el, cond, metric_name)
            n = len(source[el][cond][metric_name])
            if mean is None:
                print(f"{el:12s} {cond:9s} n=0")
            else:
                print(f"{el:12s} {cond:9s} n={n:4d} mean={mean:10.6f} std={std:10.6f}")


def print_pct_summary(source, title):
    print(f"\n{title}")
    print("-" * len(title))
    for el in ALL_ELEMENTS_AB:
        conds = ["combined"] if el == UNANNOTATED_LABEL else CONDITIONS
        for cond in conds:
            mean, std = mean_std(source, el, cond)
            n = len(source[el][cond])
            if mean is None:
                print(f"{el:12s} {cond:9s} n=0")
            else:
                print(f"{el:12s} {cond:9s} n={n:4d} mean={mean:10.6f}% std={std:10.6f}")


def print_cross_dev_summary(cross_dict, title):
    print(f"\n{title}")
    print("-" * len(title))
    row_keys = STRAND_ELEMENTS + [None]
    col_keys = [None] + list(reversed(STRAND_ELEMENTS))

    for rk in row_keys:
        for ck in col_keys:
            if rk is None and ck is None:
                continue
            vals = cross_dict.get((rk, ck), [])
            plus_label = f"{rk}+" if rk is not None else "Alone"
            minus_label = f"{ck}-" if ck is not None else "Alone"
            if len(vals) == 0:
                print(f"{plus_label:12s} vs {minus_label:12s} n=0")
            else:
                print(
                    f"{plus_label:12s} vs {minus_label:12s} "
                    f"n={len(vals):4d} mean={np.mean(vals):10.6f} std={np.std(vals):10.6f}"
                )


# =====================================================================
# 4. STYLE
# =====================================================================
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

BAR_W = 0.26
HATCHED = {SUM_ELEMENT, UNANNOTATED_LABEL}


def build_x_positions(elements):
    centers = {}
    bar_positions = []
    x = 0.0
    for idx, el in enumerate(elements):
        if el in HATCHED and (idx == 0 or elements[idx - 1] not in HATCHED):
            x += SECTION_GAP
        centers[el] = x
        if el == UNANNOTATED_LABEL:
            bar_positions.append((el, "combined", x))
        else:
            bar_positions.append((el, "combined", x - BAR_W))
            bar_positions.append((el, "plus", x))
            bar_positions.append((el, "minus", x + BAR_W))
        x += ELEMENT_GROUP_GAP
    return centers, bar_positions


def draw_panel_bar(ax, source, metric, ylabel, panel_label, ylim=None):
    elements = ALL_ELEMENTS_AB
    centers, bar_positions = build_x_positions(elements)

    for el, cond, xpos in bar_positions:
        mean, std = mean_std(source, el, cond, metric)
        if mean is None:
            continue
        if el == UNANNOTATED_LABEL:
            color = COLOR_UNANNOTATED
        else:
            color = {"combined": COLOR_COMBINED, "plus": COLOR_PLUS, "minus": COLOR_MINUS}[cond]
        width = BAR_W if el != UNANNOTATED_LABEL else BAR_W * 2.6
        hatch = "///" if el in HATCHED else None
        ax.bar(
            xpos, mean, width=width, yerr=clipped_yerr(mean, std), capsize=2,
            color=color, edgecolor=EDGE_COLOR, linewidth=0.5, hatch=hatch,
            error_kw={"elinewidth": 0.6, "capthick": 0.6, "ecolor": "#333333"}
        )

    if metric is not None:
        ax.axhline(1.0, color="black", linestyle="--", linewidth=0.7, zorder=0, alpha=0.6)

    sep = (centers[STRAND_ELEMENTS[-1]] + centers[SUM_ELEMENT]) / 2
    ax.axvline(sep, color="gray", linestyle=":", linewidth=0.6, alpha=0.6)

    ax.set_xticks(list(centers.values()))
    ax.set_xticklabels(list(centers.keys()), rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.margins(x=0.03)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_title(panel_label, loc="left", fontweight="bold", fontsize=12, pad=10)


def draw_panel_pointrange(ax, source, metric, ylabel, panel_label, ylim=None):
    elements = ALL_ELEMENTS_AB
    centers, bar_positions = build_x_positions(elements)

    for el, cond, xpos in bar_positions:
        mean, std = mean_std(source, el, cond, metric)
        if mean is None:
            continue
        if el == UNANNOTATED_LABEL:
            color = COLOR_UNANNOTATED
        else:
            color = {"combined": COLOR_COMBINED, "plus": COLOR_PLUS, "minus": COLOR_MINUS}[cond]
        marker = "D" if el in HATCHED else "o"
        ax.errorbar(
            xpos, mean, yerr=clipped_yerr(mean, std), fmt=marker, color=color,
            markersize=5, markeredgecolor=EDGE_COLOR, markeredgewidth=0.5,
            elinewidth=1.0, capsize=2.5, capthick=1.0, zorder=3
        )

    if metric is not None:
        ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8, zorder=0, alpha=0.7)

    sep = (centers[STRAND_ELEMENTS[-1]] + centers[SUM_ELEMENT]) / 2
    ax.axvline(sep, color="gray", linestyle=":", linewidth=0.6, alpha=0.5)

    ax.set_xticks(list(centers.values()))
    ax.set_xticklabels(list(centers.keys()), rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.margins(x=0.05)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_title(panel_label, loc="left", fontweight="bold", fontsize=12, pad=10)


draw_panel = draw_panel_bar if PLOT_STYLE == "bar" else draw_panel_pointrange


def sorted_strand_elements():
    return sorted(
        STRAND_ELEMENTS,
        key=lambda el: mean_std(pct_data, el, "combined")[0] or 0,
        reverse=True
    )


def _draw_group(ax, elements, source, use_labels, value_suffix=""):
    centers = {}
    bar_positions = []
    x = 0.0
    for el in elements:
        centers[el] = x
        if el == UNANNOTATED_LABEL:
            bar_positions.append((el, "combined", x))
        else:
            bar_positions.append((el, "combined", x - BAR_W))
            bar_positions.append((el, "plus", x))
            bar_positions.append((el, "minus", x + BAR_W))
        x += ELEMENT_GROUP_GAP

    label_queue = []
    for el, cond, xpos in bar_positions:
        mean, std = mean_std(source, el, cond)
        if mean is None:
            continue
        color = COLOR_UNANNOTATED if el == UNANNOTATED_LABEL else \
            {"combined": COLOR_COMBINED, "plus": COLOR_PLUS, "minus": COLOR_MINUS}[cond]
        width = BAR_W if el != UNANNOTATED_LABEL else BAR_W * 2.6
        hatch = "///" if el in HATCHED else None
        ax.bar(
            xpos, mean, width=width, yerr=clipped_yerr(mean, std), capsize=2,
            color=color, edgecolor=EDGE_COLOR, linewidth=0.5, hatch=hatch,
            error_kw={"elinewidth": 0.6, "capthick": 0.6, "ecolor": "#333333"}
        )
        if use_labels and cond == "combined":
            label_queue.append((xpos, mean, std))

    ax.set_xticks(list(centers.values()))
    ax.set_xticklabels(list(centers.keys()), rotation=45, ha="right")
    ax.spines["top"].set_visible(False)
    ax.margins(x=0.08, y=0.15)
    return centers, label_queue


def draw_portion_bar(ax_elem, ax_agg, panel_label):
    elements_ordered = sorted_strand_elements()
    _, labels_e = _draw_group(ax_elem, elements_ordered, pct_data, use_labels=True)
    if YLIM_PCT_ELEMENTS is not None:
        ax_elem.set_ylim(*YLIM_PCT_ELEMENTS)
    else:
        ax_elem.relim()
        ax_elem.autoscale_view()

    for xpos, mean, std in labels_e:
        ax_elem.text(
            xpos, mean + std + ax_elem.get_ylim()[1] * 0.02,
            f"{mean:.1f}%", ha="center", va="bottom", fontsize=6.5
        )

    _, labels_a = _draw_group(ax_agg, [SUM_ELEMENT, UNANNOTATED_LABEL], pct_data, use_labels=True)
    if YLIM_PCT_AGGREGATE is not None:
        ax_agg.set_ylim(*YLIM_PCT_AGGREGATE)
    else:
        ax_agg.relim()
        ax_agg.autoscale_view()

    for xpos, mean, std in labels_a:
        ax_agg.text(
            xpos, mean + std + ax_agg.get_ylim()[1] * 0.02,
            f"{mean:.1f}%", ha="center", va="bottom", fontsize=6.5
        )

    ax_elem.spines["right"].set_visible(False)
    ax_agg.spines["left"].set_visible(False)
    ax_agg.yaxis.tick_right()
    ax_agg.yaxis.set_label_position("right")

    ax_elem.set_ylabel("% of chromosome (elements)")
    ax_agg.set_ylabel("% of chromosome (aggregate)")
    ax_elem.set_title(panel_label, loc="left", fontweight="bold", fontsize=12, pad=10)


def draw_portion_violin(ax_elem, ax_agg, panel_label):
    def _violin_group(ax, elements):
        centers = {}
        x = 0.0
        for el in elements:
            centers[el] = x
            if el == UNANNOTATED_LABEL:
                positions = [(el, "combined", x)]
            else:
                positions = [(el, "combined", x - BAR_W), (el, "plus", x), (el, "minus", x + BAR_W)]

            for e, cond, xpos in positions:
                vals = pct_data[e][cond]
                if len(vals) == 0:
                    continue
                color = COLOR_UNANNOTATED if e == UNANNOTATED_LABEL else \
                    {"combined": COLOR_COMBINED, "plus": COLOR_PLUS, "minus": COLOR_MINUS}[cond]
                width = BAR_W * 1.8 if e != UNANNOTATED_LABEL else BAR_W * 3.5
                parts = ax.violinplot([vals], positions=[xpos], widths=width, showmeans=True, showextrema=False)
                for body in parts["bodies"]:
                    body.set_facecolor(color)
                    body.set_edgecolor(EDGE_COLOR)
                    body.set_linewidth(0.5)
                    body.set_alpha(0.9)
                    if e in HATCHED:
                        body.set_hatch("///")
                parts["cmeans"].set_color(EDGE_COLOR)
                parts["cmeans"].set_linewidth(1.0)
            x += ELEMENT_GROUP_GAP

        ax.set_xticks(list(centers.values()))
        ax.set_xticklabels(list(centers.keys()), rotation=45, ha="right")
        ax.spines["top"].set_visible(False)
        ax.margins(x=0.08)

    _violin_group(ax_elem, sorted_strand_elements())
    if YLIM_PCT_ELEMENTS is not None:
        ax_elem.set_ylim(*YLIM_PCT_ELEMENTS)

    _violin_group(ax_agg, [SUM_ELEMENT, UNANNOTATED_LABEL])
    if YLIM_PCT_AGGREGATE is not None:
        ax_agg.set_ylim(*YLIM_PCT_AGGREGATE)

    ax_elem.spines["right"].set_visible(False)
    ax_agg.spines["left"].set_visible(False)
    ax_agg.yaxis.tick_right()
    ax_agg.yaxis.set_label_position("right")

    ax_elem.set_ylabel("% of chromosome (elements)")
    ax_agg.set_ylabel("% of chromosome (aggregate)")
    ax_elem.set_title(panel_label, loc="left", fontweight="bold", fontsize=12, pad=10)
    ax_agg.set_title("Annotation / Unannotated\n(independent scale)", fontsize=8, pad=10)


draw_portion_panel = draw_portion_bar if PORTION_STYLE == "bar" else draw_portion_violin


def draw_cross_heatmap(ax, cross_dict, panel_label, vmax=None):
    row_keys = STRAND_ELEMENTS + [None]
    col_keys = [None] + list(reversed(STRAND_ELEMENTS))

    mat = np.full((len(row_keys), len(col_keys)), np.nan)
    for r, rk in enumerate(row_keys):
        for c, ck in enumerate(col_keys):
            if rk is None and ck is None:
                continue
            vals = cross_dict.get((rk, ck), [])
            if len(vals) == 0:
                continue
            mat[r, c] = np.mean(vals)

    cmap = mcolors.LinearSegmentedColormap.from_list("parity_target_2c", ["#1a9850", "white"])
    cmap.set_bad("#1a9850")
    cmap.set_over("white")
    vmax_used = vmax if vmax is not None else HEATMAP_VMAX
    norm = mcolors.Normalize(vmin=0, vmax=vmax_used, clip=False)
    im = ax.imshow(mat, cmap=cmap, norm=norm, aspect="equal")

    row_labels = [f"{e}+" for e in STRAND_ELEMENTS] + ["Alone"]
    col_labels = ["Alone"] + [f"{e}-" for e in reversed(STRAND_ELEMENTS)]
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels)

    for r in range(mat.shape[0]):
        for c in range(mat.shape[1]):
            if np.isnan(mat[r, c]):
                continue
            txt_color = "white" if mat[r, c] > vmax_used * 0.55 else "black"
            ax.text(c, r, f"{mat[r, c]:.3f}", ha="center", va="center", fontsize=7, color=txt_color)

    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="max")
    cbar.set_label(r"Deviation from parity   $\sqrt{AT_{skew}^2+GC_{skew}^2}$" + "\n(green = 0 = perfect parity)",
               labelpad=15)  # Increase this value to move label further right


    ax.set_title(panel_label, loc="left", fontweight="bold", fontsize=12, pad=10)


    return mat, row_labels, col_labels


# =====================================================================
# 5. BUILD FIGURE
# =====================================================================
fig = plt.figure(figsize=FIGSIZE)

ax_a = fig.add_axes(PANEL_RECTS["A"])
ax_b = fig.add_axes(PANEL_RECTS["B"])
ax_c = fig.add_axes(PANEL_RECTS["C"])

d_left, d_bottom, d_width, d_height = PANEL_RECTS["D"]
gap = PORTION_WSPACE * d_width
ratio_sum = sum(PORTION_SPLIT_RATIO)
elem_w = (d_width - gap) * PORTION_SPLIT_RATIO[0] / ratio_sum
agg_w = (d_width - gap) * PORTION_SPLIT_RATIO[1] / ratio_sum
ax_d_elem = fig.add_axes([d_left, d_bottom, elem_w, d_height])
ax_d_agg = fig.add_axes([d_left + elem_w + gap, d_bottom, agg_w, d_height])

draw_panel(ax_a, data, "A2T", "A / T", "A", ylim=YLIM_A2T)
draw_panel(ax_b, data, "G2C", "G / C", "B", ylim=YLIM_G2C)
ax_b.set_xlabel("Element")

mat_dev, rlab, clab = draw_cross_heatmap(ax_c, cross_dev, "C")

draw_portion_panel(ax_d_elem, ax_d_agg, "D")
ax_d_elem.set_xlabel("Element", labelpad=24)
ax_d_agg.set_xlabel("")

if PLOT_STYLE == "bar":
    legend_handles = [
        Patch(facecolor=COLOR_COMBINED, edgecolor=EDGE_COLOR, linewidth=0.5, label="Combined (+/-)"),
        Patch(facecolor=COLOR_PLUS, edgecolor=EDGE_COLOR, linewidth=0.5, label="Plus annotation"),
        Patch(facecolor=COLOR_MINUS, edgecolor=EDGE_COLOR, linewidth=0.5, label="Minus annotation"),
        Patch(facecolor=COLOR_COMBINED, edgecolor=EDGE_COLOR, hatch="///", linewidth=0.5, label="Annotation (7 elements)"),
        Patch(facecolor=COLOR_UNANNOTATED, edgecolor=EDGE_COLOR, hatch="///", linewidth=0.5, label="Unannotated"),
    ]
else:
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_COMBINED, markeredgecolor=EDGE_COLOR, markersize=6, label="Combined (+/-)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_PLUS, markeredgecolor=EDGE_COLOR, markersize=6, label="Plus strand"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_MINUS, markeredgecolor=EDGE_COLOR, markersize=6, label="Minus strand"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COLOR_COMBINED, markeredgecolor=EDGE_COLOR, markersize=6, label="Annotation (sum of 7)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COLOR_UNANNOTATED, markeredgecolor=EDGE_COLOR, markersize=6, label="Unannotated"),
    ]

DPI = 300

fig.legend(handles=legend_handles, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 0.995), frameon=False, fontsize=7.5)

#plt.savefig(f'{OUT_PREFIX}_{PLOT_STYLE}-{DPI}.tiff',dpi=DPI,                  # 600 = gold standard for line art + scatterformat='tiff',
#            bbox_inches=None,facecolor='white',edgecolor='none',pil_kwargs={"compression": "tiff_lzw"})  # lossless compression
#fig.savefig(f"{OUT_PREFIX}_{PLOT_STYLE}-{DPI}.pdf", dpi=DPI, bbox_inches="tight")

#fig.savefig(f"{OUT_PREFIX}_{PLOT_STYLE}-{DPI}.png", dpi=DPI, bbox_inches="tight")


# ── 3. SVG — pure vector, no DPI needed ───────────────────────
'''plt.savefig(f'{OUT_PREFIX}_{PLOT_STYLE}.svg',
            format='svg',
            bbox_inches=None,
            facecolor='white',
            edgecolor='none')'''

plt.savefig(f'{OUT_PREFIX}_{PLOT_STYLE}-{DPI}.tiff',dpi=DPI,                  # 600 = gold standard for line art + scatterformat='tiff',
            bbox_inches=None,facecolor='white',edgecolor='none',pil_kwargs={"compression": "tiff_lzw"})  # lossless compression
fig.savefig(f"{OUT_PREFIX}_{PLOT_STYLE}-{DPI}.pdf", dpi=DPI, bbox_inches="tight")

fig.savefig(f"{OUT_PREFIX}_{PLOT_STYLE}-{DPI}.png", dpi=DPI, bbox_inches="tight")


# ── 3. SVG — pure vector, no DPI needed ───────────────────────
plt.savefig(f'{OUT_PREFIX}_{PLOT_STYLE}.svg',
            format='svg',
            bbox_inches=None,
            facecolor='white',
            edgecolor='none')




rows_out = []
row_keys = STRAND_ELEMENTS + [None]
col_keys = [None] + list(reversed(STRAND_ELEMENTS))
for rk in row_keys:
    for ck in col_keys:
        if rk is None and ck is None:
            continue
        vals = cross_dev.get((rk, ck), [])
        rows_out.append({
            "plus_element": (rk + "+") if rk else "Alone",
            "minus_element": (ck + "-") if ck else "Alone",
            "n_chromosomes": len(vals),
            "deviation_mean": np.mean(vals) if vals else np.nan,
            "deviation_std": np.std(vals) if vals else np.nan,
        })
pd.DataFrame(rows_out).to_csv(f"{OUT_PREFIX}_cross_deviation_mean_std.csv", index=False)

plt.show()

print("\nAnalysis counts")
print("---------------")
print(f"Folders found:                     {folders_seen}")
print(f"Folders with analyzed chromosomes: {folders_with_analyzed_chr}")
print(f"Chromosomes analyzed:              {chromosomes_analyzed}")

print_ratio_summary(data, "A2T", "Panel A values: A/T")
print_ratio_summary(data, "G2C", "Panel B values: G/C")
print_cross_dev_summary(cross_dev, "Panel C values: cross-strand Chargaff deviation")
print_pct_summary(pct_data, "Panel D values: % chromosome occupied")

print("\nSaved files")
print("-----------")
print(f"{OUT_PREFIX}_{PLOT_STYLE}.pdf")
print(f"{OUT_PREFIX}_{PLOT_STYLE}.png")
print(f"{OUT_PREFIX}_cross_deviation_mean_std.csv")
