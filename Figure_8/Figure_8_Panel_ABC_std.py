"""
Figure 8 — Three-panel publication figure + full statistics report
(fixed version)

Panel A:
    Equilibrium mixing — deviation-from-parity distribution vs. sampling
    fraction, all seven taxonomic groups pooled (one random draw, Repeat 1).

Panel B:
    Mammals Gene+TE vs. Mammals Gene-only, now truly MATCHED: only
    chromosomes (Chr_ID x Fraction) present, with valid values for all three
    strand conditions, in BOTH annotation sets are plotted and tested.

Panel C:
    Mean deviation per taxonomic group across sampling fraction, split by
    strand condition (Plus / Minus / Mixed). Uses ALL repeats; error bars are
    SD (or SEM) across repeats.

CHANGES IN THIS VERSION
    1. Panel B matching fix (common Chr_IDs + common rows per fraction).
    2. New wide per-chromosome table (Plus/Minus/Mixed in the same row), which
       allows "Mixed lower than both orientations" and threshold percentages.
    3. Extended statistics report printed AND saved as CSV files in
       OUT_DIR/stats (use them as Supplementary tables):
         - dataset size per group
         - pooled equilibrium table (means, medians, fold reduction,
           % chromosomes in equilibrium, ratio to full annotation)
         - per-group equilibrium table
         - Panel C group table (fold reduction, stability, repeat CV)
         - clipped-point counts for Panels A and B (y-axis limits)
         - Panel B counts check + paired statistics (Wilcoxon)
         - sensitivity analysis without the Chordata/Mammals overlap
         - key numbers for the Results text
    4. Panel C warning no longer fires for the 100% fraction (one draw by
       design, so no error bar there).
    5. Panel C x-ticks no longer depend on the last group in the loop.
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.stats import gaussian_kde, wilcoxon


# ============================================================================
# CONFIGURATION
# ============================================================================

GROUP_FILES = [
    "Bacteria",
    "Archaea",
    "Fungi",
    "Arthropoda",
    "Plantae",
    "Chordata",
    "Mammals",
]

DATA_DIR = str(config.FIG8_SAMPLING_DIR)   # Dataframes/figure8_equilibrium_other_groups/<Group>_sampling_results.csv

OUT_DIR = str(config.output_dir(8, "Figs"))

# Statistics tables (CSV) are saved here.
STATS_DIR = os.path.join(OUT_DIR, "stats")

# ----------------------------------------------------------------------------
# Panels A and B
# ----------------------------------------------------------------------------
REPEAT_TO_USE = 1
CHR_ID_COL = "Chr_ID"

# ----------------------------------------------------------------------------
# Panel C
# ----------------------------------------------------------------------------
USE_ALL_REPEATS_PANEL_C = True
PANEL_C_ERROR_TYPE = "std"        # "std" or "sem"
PANEL_C_MIN_REPEATS = 2           # fractions < 1.0 should have at least this many

# ----------------------------------------------------------------------------
# Statistics report settings
# ----------------------------------------------------------------------------
# "In equilibrium" is reported in two ways:
#   (1) Mixed deviation lower than BOTH single orientations (no threshold)
#   (2) Mixed deviation below each of these numeric thresholds
EQUILIBRIUM_THRESHOLDS = (0.01, 0.02, 0.05)

# Fractions printed in the per-group / sensitivity tables
REPORT_FRACTIONS = [0.01, 0.05, 0.10, 0.20, 0.40, 1.00]

# Stability of Plus/Minus is judged from this fraction upward
STABILITY_FROM_FRACTION = 0.05

# ============================================================================
# FRACTIONS TO DISPLAY
# (values that do not exist in the data, e.g. 0.001, are simply skipped)
# ============================================================================

FRACTIONS_TO_SHOW_A = [0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 1.0]

FRACTIONS_TO_SHOW_B = [0.01, 0.05, 0.1, 0.4, 1.0]

FRACTIONS_TO_SHOW_C = [
    0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
]

# ============================================================================
# Y-LIMITS  (points above these limits are counted in the report)
# ============================================================================

Y_LIM_A = (0.0, 0.5)
Y_LIM_B = (0.0, 0.25)
Y_LIM_C = (0.0, 0.13)

# ============================================================================
# PANEL C MODE
# ============================================================================

PANEL_C_MODE = "split_1x3"        # "combined" or "split_1x3"

# ============================================================================
# FIGURE SIZE
# ============================================================================

FIG_WIDTH_MM = 210
FIG_HEIGHT_MM = 260
DPI_RASTER = 300

TICK_LABEL_SIZE_A = 6
TICK_LABEL_SIZE_B = 6
TICK_LABEL_SIZE_C = 6

# ============================================================================
# LAYOUT CONTROLS
# ============================================================================

A_GROUP_SPACING = 1.9
A_STRAND_SPACING = 0.4
A_VIOLIN_WIDTH = 0.20
A_BOX_WIDTH = 0.04
A_JITTER_WIDTH = 0.035
A_JITTER_OFFSET = -0.07
A_BOX_OFFSET = -0.035

B_GROUP_SPACING = 3.60
B_SOURCE_SPACING = 1.1
B_STRAND_SPACING = 0.3
B_VIOLIN_WIDTH = 0.15
B_BOX_WIDTH = 0.040
B_JITTER_WIDTH = 0.025
B_JITTER_OFFSET = -0.075
B_BOX_OFFSET = -0.022

SHOW_CONNECTOR_LINES_A = False
SHOW_CONNECTOR_LINES_B = False

PANEL_C_ERROR_CAPSIZE = 2.5
PANEL_C_ERROR_LINEWIDTH = 0.8
PANEL_C_MARKERSIZE = 3.2
PANEL_C_LINEWIDTH = 1.0

# ============================================================================
# STRAND ORDER AND COLORS
# ============================================================================

STRAND_ORDER = ["Plus", "Minus", "Mixed"]
SKEW_COLS = [f"{s}_skew" for s in STRAND_ORDER]

PALETTE = {
    "Plus": "#D55E00",
    "Minus": "#007C91",
    "Mixed": "#6A5A91",
}

MAMMALS_GENE_TE_SHADES = {
    "Plus": "#66A61E",
    "Minus": "#4D8A17",
    "Mixed": "#2F5D0E",
}

MAMMALS_GENE_SHADES = {
    "Plus": "#E15759",
    "Minus": "#C83E3E",
    "Mixed": "#8F2626",
}

SOURCE_LABELS_SHORT = {
    "mammals_gene_te": "GeneTE",
    "mammals_gene": "Gene",
}

GROUP_ORDER = [
    "Bacteria",
    "Archaea",
    "Fungi",
    "Arthropoda",
    "Plantae",
    "Chordata",
    "Mammals",
]

GROUP_COLORS = {
    "Bacteria": "#0072B2",
    "Archaea": "#E69F00",
    "Fungi": "#009E73",
    "Arthropoda": "#D55E00",
    "Plantae": "#7B61A8",
    "Chordata": "#8C6D4F",
    "Mammals": "#CC79A7",
}

STRAND_LINESTYLE = {"Plus": "-", "Minus": "--", "Mixed": ":"}
STRAND_MARKER = {"Plus": "o", "Minus": "s", "Mixed": "^"}

# ============================================================================
# STYLE
# ============================================================================

mm_to_in = 1 / 25.4

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7.5,
    "axes.linewidth": 0.75,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


# ============================================================================
# SMALL HELPERS
# ============================================================================

def print_section(title):
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def show(df, title=None):
    if title:
        print("\n" + title)
    print(df.round(5).to_string(index=False))


def save_table(df, name):
    os.makedirs(STATS_DIR, exist_ok=True)
    path = os.path.join(STATS_DIR, f"Fig8_{name}.csv")
    df.to_csv(path, index=False)
    print(f"[saved] {path}")


def fraction_in(value, allowed):
    return any(np.isclose(value, a) for a in allowed)


# ============================================================================
# SKEW CALCULATION
# ============================================================================

def calculate_skew_from_arrays(A, T, C, G):
    """
    Normalized Chargaff deviation per row:

        sqrt( ((A-T)/(A+T))^2 + ((G-C)/(G+C))^2 )
    """
    A = np.asarray(A, dtype=float)
    T = np.asarray(T, dtype=float)
    C = np.asarray(C, dtype=float)
    G = np.asarray(G, dtype=float)

    at_den = A + T
    gc_den = G + C

    at_part = np.where(at_den == 0, np.nan, (A - T) / at_den)
    gc_part = np.where(gc_den == 0, np.nan, (G - C) / gc_den)

    return np.sqrt(at_part ** 2 + gc_part ** 2)


def compute_skew_long(df, repeat_value=REPEAT_TO_USE):
    """
    Long table (Fraction, Strand, Skew) for one selected Repeat.
    Used for the Panel A plot (unchanged from the original script).
    """
    if repeat_value is not None:
        df = df[df["Repeat"] == repeat_value]

    records = []
    for fraction_step, mask in df.groupby("Fraction"):
        for strand in STRAND_ORDER:
            skew = calculate_skew_from_arrays(
                mask[f"{strand}_A_raw"].to_numpy(dtype=float),
                mask[f"{strand}_T_raw"].to_numpy(dtype=float),
                mask[f"{strand}_C_raw"].to_numpy(dtype=float),
                mask[f"{strand}_G_raw"].to_numpy(dtype=float),
            )
            valid = skew[~np.isnan(skew)]
            records.extend(
                {"Fraction": fraction_step, "Strand": strand, "Skew": value}
                for value in valid
            )
    return pd.DataFrame.from_records(records)


def compute_skew_wide(df, repeat_value=REPEAT_TO_USE):
    """
    One row per chromosome x fraction, with Plus_skew, Minus_skew and
    Mixed_skew in the SAME row (needed for paired / equilibrium statistics).
    """
    if repeat_value is not None:
        df = df[df["Repeat"] == repeat_value]

    if CHR_ID_COL in df.columns:
        chr_ids = df[CHR_ID_COL].astype(str).to_numpy()
    else:
        chr_ids = np.array(["NA"] * len(df))

    out = pd.DataFrame({
        "Chr_ID": chr_ids,
        "Fraction": df["Fraction"].to_numpy(dtype=float),
    })

    if "Organ_ID" in df.columns:
        out.insert(0, "Organ_ID", df["Organ_ID"].astype(str).to_numpy())

    for strand in STRAND_ORDER:
        out[f"{strand}_skew"] = calculate_skew_from_arrays(
            df[f"{strand}_A_raw"].to_numpy(dtype=float),
            df[f"{strand}_T_raw"].to_numpy(dtype=float),
            df[f"{strand}_C_raw"].to_numpy(dtype=float),
            df[f"{strand}_G_raw"].to_numpy(dtype=float),
        )

    return out.reset_index(drop=True)


def restrict_to_common_rows(wide_a, wide_b):
    """
    Keep only (Chr_ID, Fraction) rows that have valid values for all three
    strand conditions in BOTH tables.
    """
    a = wide_a.dropna(subset=SKEW_COLS)
    b = wide_b.dropna(subset=SKEW_COLS)

    keys = (
        a[["Chr_ID", "Fraction"]].drop_duplicates()
        .merge(
            b[["Chr_ID", "Fraction"]].drop_duplicates(),
            on=["Chr_ID", "Fraction"],
            how="inner",
        )
    )

    a = a.merge(keys, on=["Chr_ID", "Fraction"], how="inner")
    b = b.merge(keys, on=["Chr_ID", "Fraction"], how="inner")

    return a.reset_index(drop=True), b.reset_index(drop=True)


def wide_to_long(wide, source):
    id_vars = [c for c in ("Organ_ID", "Chr_ID", "Fraction") if c in wide.columns]

    long_df = wide.melt(
        id_vars=id_vars,
        value_vars=SKEW_COLS,
        var_name="Strand",
        value_name="Skew",
    )
    long_df["Strand"] = long_df["Strand"].str.replace("_skew", "", regex=False)
    long_df["Source"] = source

    return long_df.dropna(subset=["Skew"]).reset_index(drop=True)


def compute_skew_repeat_level(df):
    """
    One mean skew value for every Fraction x Repeat x Strand
    (one independent observation per repeat, used by Panel C).
    """
    records = []

    for (fraction_step, repeat_value), mask in df.groupby(["Fraction", "Repeat"]):
        for strand in STRAND_ORDER:
            skew = calculate_skew_from_arrays(
                mask[f"{strand}_A_raw"].to_numpy(dtype=float),
                mask[f"{strand}_T_raw"].to_numpy(dtype=float),
                mask[f"{strand}_C_raw"].to_numpy(dtype=float),
                mask[f"{strand}_G_raw"].to_numpy(dtype=float),
            )
            valid = skew[~np.isnan(skew)]
            if len(valid) == 0:
                continue

            records.append({
                "Fraction": fraction_step,
                "Repeat": repeat_value,
                "Strand": strand,
                "Skew": np.mean(valid),
                "N": len(valid),
            })

    return pd.DataFrame.from_records(records)


def summarize_panel_c_repeats(repeat_level_df):
    """Mean, SD, SEM across repeats for Group x Strand x Fraction."""
    if repeat_level_df.empty:
        return pd.DataFrame()

    summary = (
        repeat_level_df
        .groupby(["Group", "Strand", "Fraction"])["Skew"]
        .agg(["count", "mean", "std"])
        .reset_index()
        .rename(columns={"count": "N_repeats", "mean": "Mean", "std": "SD"})
    )

    summary["SEM"] = summary["SD"] / np.sqrt(summary["N_repeats"])

    if PANEL_C_ERROR_TYPE.lower() == "std":
        summary["Error"] = summary["SD"]
        summary["Error_Label"] = "SD"
    elif PANEL_C_ERROR_TYPE.lower() == "sem":
        summary["Error"] = summary["SEM"]
        summary["Error_Label"] = "SEM"
    else:
        raise ValueError("PANEL_C_ERROR_TYPE must be 'std' or 'sem'.")

    return summary


# ============================================================================
# PANEL A DATA
# ============================================================================

per_group_frames = []
wide_frames = []

for group_file in GROUP_FILES:
    path = os.path.join(DATA_DIR, f"{group_file}_sampling_results.csv")
    raw = pd.read_csv(path)

    df_group = compute_skew_long(raw)
    df_group["Group"] = group_file
    per_group_frames.append(df_group)

    wide_group = compute_skew_wide(raw)
    wide_group["Group"] = group_file
    wide_frames.append(wide_group)

pooled_df = pd.concat(per_group_frames, ignore_index=True)
wide_all = pd.concat(wide_frames, ignore_index=True)

fractions_by_group = {
    group: set(df["Fraction"].unique())
    for group, df in zip(GROUP_FILES, per_group_frames)
}

all_fractions = sorted(set().union(*fractions_by_group.values()))

if FRACTIONS_TO_SHOW_A is not None:
    all_fractions = [f for f in all_fractions if f in FRACTIONS_TO_SHOW_A]

for fraction in all_fractions:
    contributing = [g for g, fr in fractions_by_group.items() if fraction in fr]
    if len(contributing) < len(GROUP_FILES):
        missing = sorted(set(GROUP_FILES) - set(contributing))
        print(f"[warning] Panel A: Fraction={fraction} is missing data from: {missing}")

max_skew_a = pooled_df.loc[pooled_df["Fraction"].isin(all_fractions), "Skew"].max()
if pd.notna(max_skew_a) and max_skew_a > Y_LIM_A[1]:
    print(f"[warning] Panel A: max skew in data ({max_skew_a:.3f}) exceeds Y_LIM_A top ({Y_LIM_A[1]})")


# ============================================================================
# PANEL B DATA  (matched Mammals Gene+TE vs Gene-only)
# ============================================================================

raw_mammals = pd.read_csv(os.path.join(DATA_DIR, "Mammals_sampling_results.csv"))
raw_chordata = pd.read_csv(os.path.join(DATA_DIR, "Chordata_sampling_results.csv"))

if CHR_ID_COL not in raw_mammals.columns or CHR_ID_COL not in raw_chordata.columns:
    raise KeyError(f"Expected a '{CHR_ID_COL}' column in both Mammals and Chordata CSVs.")

raw_mammals[CHR_ID_COL] = raw_mammals[CHR_ID_COL].astype(str)
raw_chordata[CHR_ID_COL] = raw_chordata[CHR_ID_COL].astype(str)

mammal_chr_ids = set(raw_mammals[CHR_ID_COL].unique())
chordata_chr_ids = set(raw_chordata[CHR_ID_COL].unique())
common_chr_ids = mammal_chr_ids & chordata_chr_ids
missing_in_chordata = mammal_chr_ids - chordata_chr_ids

if missing_in_chordata:
    print(
        f"[warning] Panel B: {len(missing_in_chordata)} Chr_ID(s) from the Mammals "
        f"file were not found in Chordata and are excluded from Panel B."
    )

wide_te_all = compute_skew_wide(raw_mammals[raw_mammals[CHR_ID_COL].isin(common_chr_ids)])
wide_gene_all = compute_skew_wide(raw_chordata[raw_chordata[CHR_ID_COL].isin(common_chr_ids)])

for label, wide in (("Gene+TE", wide_te_all), ("Gene-only", wide_gene_all)):
    n_dup = int(wide.duplicated(subset=["Chr_ID", "Fraction"]).sum())
    if n_dup:
        print(
            f"[warning] Panel B: {label} has {n_dup} duplicated (Chr_ID, Fraction) rows "
            f"(same Chr_ID in more than one organism?). Check pairing."
        )

wide_te_matched, wide_gene_matched = restrict_to_common_rows(wide_te_all, wide_gene_all)

print(
    f"[info] Panel B: matched rows kept = {len(wide_te_matched):,} "
    f"(Gene+TE before matching: {len(wide_te_all):,}; "
    f"Gene-only before matching: {len(wide_gene_all):,})"
)

panel_b_df = pd.concat(
    [
        wide_to_long(wide_te_matched, "mammals_gene_te"),
        wide_to_long(wide_gene_matched, "mammals_gene"),
    ],
    ignore_index=True,
)

fractions_by_source = {
    source: set(df["Fraction"].unique())
    for source, df in panel_b_df.groupby("Source")
}

panel_b_fractions = sorted(set().union(*fractions_by_source.values()))

if FRACTIONS_TO_SHOW_B is not None:
    panel_b_fractions = [f for f in panel_b_fractions if f in FRACTIONS_TO_SHOW_B]

for fraction in panel_b_fractions:
    contributing = [s for s, fr in fractions_by_source.items() if fraction in fr]
    if len(contributing) < 2:
        missing = sorted(set(fractions_by_source) - set(contributing))
        print(f"[warning] Panel B: Fraction={fraction} is missing source(s): {missing}")

max_skew_b = panel_b_df.loc[panel_b_df["Fraction"].isin(panel_b_fractions), "Skew"].max()
if pd.notna(max_skew_b) and max_skew_b > Y_LIM_B[1]:
    print(f"[warning] Panel B: max skew in data ({max_skew_b:.4f}) exceeds Y_LIM_B top ({Y_LIM_B[1]})")


# ============================================================================
# PANEL C DATA
# ============================================================================

print_section("Preparing Panel C using ALL repeats")

panel_c_repeat_frames = []

for group_file in GROUP_FILES:
    path = os.path.join(DATA_DIR, f"{group_file}_sampling_results.csv")
    raw = pd.read_csv(path)

    if "Repeat" not in raw.columns:
        raise KeyError(f"Panel C requires a 'Repeat' column in {group_file}.")

    df_repeat = compute_skew_repeat_level(raw)
    df_repeat["Group"] = group_file
    panel_c_repeat_frames.append(df_repeat)

panel_c_repeat_level = pd.concat(panel_c_repeat_frames, ignore_index=True)

panel_c_repeat_level = panel_c_repeat_level[
    panel_c_repeat_level["Fraction"].isin(FRACTIONS_TO_SHOW_C)
].copy()

panel_c_summary = summarize_panel_c_repeats(panel_c_repeat_level)

print("\nPanel C repeat counts:")
repeat_count_table = (
    panel_c_summary
    .groupby(["Group", "Strand"])["N_repeats"]
    .agg(["min", "max", "mean"])
    .reset_index()
)
print(repeat_count_table.to_string(index=False))

insufficient = panel_c_summary[
    (panel_c_summary["N_repeats"] < PANEL_C_MIN_REPEATS)
    & (panel_c_summary["Fraction"] < 1.0)
    ]

if not insufficient.empty:
    print(
        f"\n[warning] Panel C: some Group x Strand x Fraction combinations below 100% "
        f"have fewer than {PANEL_C_MIN_REPEATS} repeats."
    )
else:
    print(
        "\n[info] Panel C: every fraction below 100% has enough repeats. The 100% "
        "fraction is one deterministic draw (the whole annotation), so it has no "
        "error bar by design."
    )

panel_c_fractions = sorted(panel_c_summary["Fraction"].unique())

for group in GROUP_ORDER:
    for strand in STRAND_ORDER:
        have = set(
            panel_c_summary.loc[
                (panel_c_summary["Group"] == group) & (panel_c_summary["Strand"] == strand),
                "Fraction",
            ]
        )
        missing = sorted(set(panel_c_fractions) - have)
        if missing:
            print(f"[warning] Panel C: {group}/{strand} is missing Fraction(s): {missing}")

max_skew_c = panel_c_summary["Mean"].max()
if pd.notna(max_skew_c) and max_skew_c > Y_LIM_C[1]:
    print(f"[warning] Panel C: max mean skew ({max_skew_c:.3f}) exceeds Y_LIM_C top ({Y_LIM_C[1]})")

print_section(f"PANEL C — mean ± {PANEL_C_ERROR_TYPE.upper()} across repeats")
print(
    panel_c_summary[
        ["Group", "Strand", "Fraction", "N_repeats", "Mean", "SD", "SEM", "Error"]
    ]
    .sort_values(["Group", "Strand", "Fraction"])
    .to_string(index=False)
)


# ============================================================================
# RAINCLOUD FUNCTIONS
# ============================================================================

def half_violin(ax, data, center, width, color, side="left"):
    data = np.asarray(data, dtype=float)

    if len(data) < 2:
        return
    if np.std(data) == 0:
        return

    kde = gaussian_kde(data)
    y_grid = np.linspace(data.min(), data.max(), 200)
    density = kde(y_grid)

    if density.max() == 0:
        return

    density = (density / density.max()) * width
    sign = -1 if side == "left" else 1

    ax.fill_betweenx(
        y_grid,
        center,
        center + sign * density,
        color=color,
        alpha=0.55,
        linewidth=0,
        zorder=2,
        )


def draw_group(
        ax,
        x0,
        skews,
        color,
        jitter_width=0.045,
        violin_width=0.16,
        box_width=0.05,
        box_offset=-0.03,
        jitter_offset=-0.10,
):
    skews = np.asarray(skews, dtype=float)

    if skews.size == 0:
        return None

    half_violin(ax, skews, x0, violin_width, color, side="right")

    ax.boxplot(
        skews,
        positions=[x0 + box_offset],
        widths=box_width,
        vert=True,
        patch_artist=True,
        showfliers=False,
        zorder=3,
        boxprops=dict(facecolor=color, alpha=0.85, linewidth=0.55, edgecolor="black"),
        medianprops=dict(color="black", linewidth=1.0),
        whiskerprops=dict(color="black", linewidth=0.65),
        capprops=dict(color="black", linewidth=0.65),
    )

    rng = np.random.default_rng(0)
    jitter = rng.uniform(-jitter_width, jitter_width, size=skews.size)

    ax.scatter(
        np.full(skews.size, x0 + jitter_offset) + jitter,
        skews,
        s=2.8,
        color=color,
        alpha=0.50,
        linewidths=0,
        zorder=1,
        )

    return (x0, float(skews.max()), float(skews.mean()), float(skews.std()))


def label_stats(ax, x0, top_y, mean, std, y_range, y_ceiling,
                fontsize=4.5, offset_frac=0.025):
    y_pos = min(top_y + offset_frac * y_range, y_ceiling)

    if y_range >= 0.1:
        fmt = f"{mean:.2f}±{std:.2f}"
    else:
        fmt = f"{mean:.3f}±{std:.3f}"

    ax.text(
        x0 + 0.05,
        y_pos,
        fmt,
        rotation=90,
        fontsize=fontsize,
        ha="center",
        va="bottom",
        zorder=6,
        )


# ============================================================================
# STATISTICS REPORT — ORIGINAL TABLES
# ============================================================================

def format_summary_table(df, group_cols, value_col="Skew", sort_cols=None):
    summary = (
        df.groupby(group_cols)[value_col]
        .agg(["count", "mean", "std", "median", "min", "max"])
        .reset_index()
    )
    if sort_cols is not None:
        summary = summary.sort_values(sort_cols)
    return summary


def print_figure_statistics():
    print_section("PANEL A — pooled all-groups skew statistics by Fraction and Strand")
    panel_a_stats = format_summary_table(
        pooled_df[pooled_df["Fraction"].isin(all_fractions)],
        ["Fraction", "Strand"],
        sort_cols=["Fraction", "Strand"],
    )
    print(panel_a_stats.to_string(index=False))
    save_table(panel_a_stats, "panelA_distribution_stats")

    print_section("PANEL B — matched Mammals Gene+TE vs Gene-only statistics")
    panel_b_stats = format_summary_table(
        panel_b_df[panel_b_df["Fraction"].isin(panel_b_fractions)],
        ["Fraction", "Source", "Strand"],
        sort_cols=["Fraction", "Source", "Strand"],
    )
    print(panel_b_stats.to_string(index=False))
    save_table(panel_b_stats, "panelB_distribution_stats")

    save_table(panel_c_summary, "panelC_group_means_sd")


# ============================================================================
# STATISTICS REPORT — EXTENDED TABLES
# ============================================================================

def equilibrium_table(wide, fractions=None):
    """
    Per fraction: N, means, medians, fold reduction, ratio to full annotation,
    % chromosomes where Mixed < both orientations, and % below thresholds.
    """
    wide = wide.dropna(subset=SKEW_COLS)

    full = wide[np.isclose(wide["Fraction"], 1.0)]
    full_plus = full["Plus_skew"].mean() if len(full) else np.nan
    full_minus = full["Minus_skew"].mean() if len(full) else np.nan

    rows = []

    for fraction, g in wide.groupby("Fraction"):
        if fractions is not None and not fraction_in(fraction, fractions):
            continue

        plus_mean = g["Plus_skew"].mean()
        minus_mean = g["Minus_skew"].mean()
        mixed_mean = g["Mixed_skew"].mean()

        plus_med = g["Plus_skew"].median()
        minus_med = g["Minus_skew"].median()
        mixed_med = g["Mixed_skew"].median()

        single_min = g[["Plus_skew", "Minus_skew"]].min(axis=1)

        row = {
            "Fraction": fraction,
            "N": len(g),
            "Plus_mean": plus_mean,
            "Minus_mean": minus_mean,
            "Mixed_mean": mixed_mean,
            "Plus_median": plus_med,
            "Minus_median": minus_med,
            "Mixed_median": mixed_med,
            "Fold_mean": ((plus_mean + minus_mean) / 2) / mixed_mean,
            "Fold_median": ((plus_med + minus_med) / 2) / mixed_med,
            "Plus_vs_full": plus_mean / full_plus,
            "Minus_vs_full": minus_mean / full_minus,
            "Pct_mixed_below_both": 100 * (g["Mixed_skew"] < single_min).mean(),
        }

        for thr in EQUILIBRIUM_THRESHOLDS:
            row[f"Pct_mixed_below_{thr}"] = 100 * (g["Mixed_skew"] < thr).mean()

        rows.append(row)

    return pd.DataFrame(rows)


def equilibrium_table_by_group(wide, fractions):
    frames = []
    for group in GROUP_ORDER:
        sub = wide[wide["Group"] == group]
        if sub.empty:
            continue
        table = equilibrium_table(sub, fractions)
        table.insert(0, "Group", group)
        frames.append(table)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def panel_c_group_table():
    """Group x fraction table from Panel C repeat means."""
    s = panel_c_summary

    mean = s.pivot_table(index=["Group", "Fraction"], columns="Strand", values="Mean")
    sd = s.pivot_table(index=["Group", "Fraction"], columns="Strand", values="SD")

    tbl = pd.DataFrame({
        "Plus": mean["Plus"],
        "Minus": mean["Minus"],
        "Mixed": mean["Mixed"],
    })

    tbl["Single_mean"] = (tbl["Plus"] + tbl["Minus"]) / 2
    tbl["Fold_single_over_mixed"] = tbl["Single_mean"] / tbl["Mixed"]
    tbl["Max_CV_pct"] = 100 * (sd / mean).max(axis=1)

    tbl = tbl.reset_index()

    full = tbl[np.isclose(tbl["Fraction"], 1.0)].set_index("Group")
    tbl["Plus_vs_full"] = tbl["Plus"] / tbl["Group"].map(full["Plus"])
    tbl["Minus_vs_full"] = tbl["Minus"] / tbl["Group"].map(full["Minus"])

    tbl["Group"] = pd.Categorical(tbl["Group"], GROUP_ORDER, ordered=True)
    return tbl.sort_values(["Group", "Fraction"]).reset_index(drop=True)


def clipped_points_report(df_long, ylim_top, fractions):
    d = df_long[df_long["Fraction"].isin(fractions)]
    table = (
        d.groupby(["Fraction", "Strand"])["Skew"]
        .agg(
            N="size",
            N_above_axis=lambda s: int((s > ylim_top).sum()),
            Max="max",
        )
        .reset_index()
    )
    return table


def panel_b_paired_table():
    cols = SKEW_COLS

    te = wide_te_matched.groupby(["Chr_ID", "Fraction"], as_index=False)[cols].mean()
    ge = wide_gene_matched.groupby(["Chr_ID", "Fraction"], as_index=False)[cols].mean()

    merged = te.merge(ge, on=["Chr_ID", "Fraction"], suffixes=("_TE", "_Gene"))

    rows = []
    for fraction, g in merged.groupby("Fraction"):
        if fraction not in panel_b_fractions:
            continue

        diff = g["Mixed_skew_TE"] - g["Mixed_skew_Gene"]

        try:
            p_value = wilcoxon(diff).pvalue
        except ValueError:
            p_value = np.nan

        te_mean = g["Mixed_skew_TE"].mean()
        gene_mean = g["Mixed_skew_Gene"].mean()

        rows.append({
            "Fraction": fraction,
            "N_pairs": len(g),
            "Mixed_mean_GeneTE": te_mean,
            "Mixed_mean_Gene": gene_mean,
            "Pct_change_mean": 100 * (te_mean - gene_mean) / gene_mean,
            "Mixed_median_GeneTE": g["Mixed_skew_TE"].median(),
            "Mixed_median_Gene": g["Mixed_skew_Gene"].median(),
            "Pct_chr_GeneTE_lower": 100 * (diff < 0).mean(),
            "Wilcoxon_p": p_value,
            "Plus_mean_GeneTE": g["Plus_skew_TE"].mean(),
            "Plus_mean_Gene": g["Plus_skew_Gene"].mean(),
            "Minus_mean_GeneTE": g["Minus_skew_TE"].mean(),
            "Minus_mean_Gene": g["Minus_skew_Gene"].mean(),
        })

    return pd.DataFrame(rows)


def pick_row(table, fraction):
    hit = table[np.isclose(table["Fraction"], fraction)]
    return hit.iloc[0] if len(hit) else None


def print_extended_statistics():
    os.makedirs(STATS_DIR, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Dataset size
    # ------------------------------------------------------------------
    print_section("EXTRA 1 — dataset size (chromosomes with valid values at 100%)")

    top = wide_all[np.isclose(wide_all["Fraction"], 1.0)].dropna(subset=SKEW_COLS)

    agg = {"N_chromosomes": ("Chr_ID", "size")}
    if "Organ_ID" in top.columns:
        agg["N_organisms"] = ("Organ_ID", "nunique")

    size_tbl = top.groupby("Group").agg(**agg).reindex(GROUP_ORDER).reset_index()

    repeats = (
        panel_c_summary.groupby("Group")["N_repeats"].max()
        .rename("Max_repeats_per_fraction")
        .reset_index()
    )
    size_tbl = size_tbl.merge(repeats, on="Group", how="left")

    total = {"Group": "SUM of groups", "N_chromosomes": int(size_tbl["N_chromosomes"].sum())}
    if "N_organisms" in size_tbl.columns:
        total["N_organisms"] = int(size_tbl["N_organisms"].sum())

    size_tbl = pd.concat([size_tbl, pd.DataFrame([total])], ignore_index=True)

    show(size_tbl)
    print("Note: Chordata may contain the mammal chromosomes, so the sum can count some twice.")
    save_table(size_tbl, "dataset_size")

    # ------------------------------------------------------------------
    # 2. Pooled equilibrium table
    # ------------------------------------------------------------------
    print_section("EXTRA 2 — pooled equilibrium table (all groups, Repeat "
                  f"{REPEAT_TO_USE})")

    eq_all = equilibrium_table(wide_all)
    show(eq_all)
    save_table(eq_all, "panelA_equilibrium_pooled")

    # ------------------------------------------------------------------
    # 3. Per-group equilibrium table
    # ------------------------------------------------------------------
    print_section("EXTRA 3 — per-group equilibrium (selected fractions)")

    eq_group = equilibrium_table_by_group(wide_all, REPORT_FRACTIONS)
    keep = ["Group", "Fraction", "N", "Fold_mean", "Pct_mixed_below_both"] + [
        f"Pct_mixed_below_{thr}" for thr in EQUILIBRIUM_THRESHOLDS
    ]
    show(eq_group[keep])
    save_table(eq_group, "per_group_equilibrium")

    # ------------------------------------------------------------------
    # 4. Panel C group table, stability and repeat variability
    # ------------------------------------------------------------------
    print_section("EXTRA 4 — Panel C group table (means across repeats)")

    c_tbl = panel_c_group_table()
    show(c_tbl[c_tbl["Fraction"].apply(lambda f: fraction_in(f, REPORT_FRACTIONS))])
    save_table(c_tbl, "panelC_group_table")

    below_100 = c_tbl[c_tbl["Fraction"] < 1.0]
    if len(below_100):
        worst = below_100.loc[below_100["Max_CV_pct"].idxmax()]
        print(
            f"\nLargest repeat-to-repeat variation (SD/mean) below 100%: "
            f"{worst['Max_CV_pct']:.2f}% ({worst['Group']}, fraction {worst['Fraction']})"
        )

    stab = c_tbl[c_tbl["Fraction"] >= STABILITY_FROM_FRACTION - 1e-9]
    stab_rows = []
    for group, g in stab.groupby("Group", observed=True):
        dev = pd.concat([(g["Plus_vs_full"] - 1).abs(), (g["Minus_vs_full"] - 1).abs()])
        stab_rows.append({
            "Group": group,
            "Max_abs_dev_vs_full_pct": 100 * dev.max(),
        })
    stab_tbl = pd.DataFrame(stab_rows)

    print(
        f"\nStability of Plus/Minus deviation from "
        f"{STABILITY_FROM_FRACTION * 100:g}% upward (max |value/full - 1|):"
    )
    show(stab_tbl)
    n_within = int((stab_tbl["Max_abs_dev_vs_full_pct"] <= 10).sum())
    print(f"Groups within 10% of the full-annotation value: {n_within} of {len(stab_tbl)}")
    save_table(stab_tbl, "panelC_stability_vs_full")

    full_tbl = c_tbl[np.isclose(c_tbl["Fraction"], 1.0)]
    if len(full_tbl):
        s_lo = full_tbl.loc[full_tbl["Single_mean"].idxmin()]
        s_hi = full_tbl.loc[full_tbl["Single_mean"].idxmax()]
        m_lo = full_tbl.loc[full_tbl["Mixed"].idxmin()]
        m_hi = full_tbl.loc[full_tbl["Mixed"].idxmax()]
        print(
            f"\nAt 100%: single-orientation deviation ranges "
            f"{s_lo['Single_mean']:.4f} ({s_lo['Group']}) to "
            f"{s_hi['Single_mean']:.4f} ({s_hi['Group']}) "
            f"= {s_hi['Single_mean'] / s_lo['Single_mean']:.1f}-fold; "
            f"mixed deviation ranges {m_lo['Mixed']:.4f} ({m_lo['Group']}) to "
            f"{m_hi['Mixed']:.4f} ({m_hi['Group']}) "
            f"= {m_hi['Mixed'] / m_lo['Mixed']:.1f}-fold."
        )

    # ------------------------------------------------------------------
    # 5. Points above the y-axis limits
    # ------------------------------------------------------------------
    print_section("EXTRA 5 — points above the y-axis limits")

    clip_a = clipped_points_report(pooled_df, Y_LIM_A[1], all_fractions)
    show(clip_a, f"Panel A (limit {Y_LIM_A[1]}):")
    print(f"Total above axis: {int(clip_a['N_above_axis'].sum())} of {int(clip_a['N'].sum())}")
    save_table(clip_a, "panelA_points_above_axis")

    clip_b = clipped_points_report(panel_b_df, Y_LIM_B[1], panel_b_fractions)
    show(clip_b, f"Panel B (limit {Y_LIM_B[1]}):")
    print(f"Total above axis: {int(clip_b['N_above_axis'].sum())} of {int(clip_b['N'].sum())}")
    save_table(clip_b, "panelB_points_above_axis")

    # ------------------------------------------------------------------
    # 6. Panel B matching check and paired statistics
    # ------------------------------------------------------------------
    print_section("EXTRA 6 — Panel B matching check and paired statistics")

    counts = (
        panel_b_df[panel_b_df["Fraction"].isin(panel_b_fractions)]
        .groupby(["Fraction", "Strand", "Source"]).size()
        .unstack("Source")
        .reset_index()
    )
    show(counts, "Chromosomes per source (must be identical after matching):")

    if {"mammals_gene", "mammals_gene_te"} <= set(counts.columns):
        same = bool((counts["mammals_gene"] == counts["mammals_gene_te"]).all())
        print(f"Counts identical across sources: {same}")

    paired = panel_b_paired_table()
    show(paired, "Paired comparison, Gene+TE versus Gene-only (same chromosomes):")
    save_table(paired, "panelB_paired_stats")

    # ------------------------------------------------------------------
    # 7. Sensitivity without the Chordata / Mammals overlap
    # ------------------------------------------------------------------
    print_section("EXTRA 7 — sensitivity: pooled Panel A without mammal chromosomes "
                  "inside Chordata")

    overlap_mask = (wide_all["Group"] == "Chordata") & wide_all["Chr_ID"].isin(mammal_chr_ids)
    n_overlap = int((overlap_mask & np.isclose(wide_all["Fraction"], 1.0)).sum())
    print(f"Chordata chromosomes that are also in the Mammals file (at 100%): {n_overlap}")

    sens = equilibrium_table(wide_all.loc[~overlap_mask], REPORT_FRACTIONS)
    show(sens)
    save_table(sens, "panelA_sensitivity_without_overlap")

    # ------------------------------------------------------------------
    # 8. Key numbers for the Results text
    # ------------------------------------------------------------------
    print_section("EXTRA 8 — key numbers for the Results text")

    thr_cols = [f"Pct_mixed_below_{thr}" for thr in EQUILIBRIUM_THRESHOLDS]

    for fraction in (0.01, 0.05, 0.10, 1.0):
        row = pick_row(eq_all, fraction)
        if row is None:
            continue
        thr_text = ", ".join(
            f"<{thr}: {row[col]:.1f}%" for thr, col in zip(EQUILIBRIUM_THRESHOLDS, thr_cols)
        )
        print(
            f"{fraction * 100:g}%: Plus {row['Plus_mean']:.4f}, Minus {row['Minus_mean']:.4f}, "
            f"Mixed {row['Mixed_mean']:.4f} | fold (mean) {row['Fold_mean']:.1f}, "
            f"fold (median) {row['Fold_median']:.1f} | Mixed lower than both: "
            f"{row['Pct_mixed_below_both']:.1f}% | {thr_text}"
        )

    for fraction in (0.01, 0.10, 1.0):
        row = pick_row(paired, fraction)
        if row is None:
            continue
        print(
            f"Panel B {fraction * 100:g}%: Mixed Gene+TE {row['Mixed_mean_GeneTE']:.4f} vs "
            f"Gene {row['Mixed_mean_Gene']:.4f} ({row['Pct_change_mean']:+.0f}%), "
            f"Gene+TE lower in {row['Pct_chr_GeneTE_lower']:.1f}% of chromosomes, "
            f"Wilcoxon p = {row['Wilcoxon_p']:.2e}, n = {int(row['N_pairs'])}"
        )

    print(f"\nTables saved in: {STATS_DIR}")


print_figure_statistics()
print_extended_statistics()


# ============================================================================
# FIGURE ASSEMBLY
# ============================================================================

fig = plt.figure(
    figsize=(FIG_WIDTH_MM * mm_to_in, FIG_HEIGHT_MM * mm_to_in),
    dpi=150,
)

outer_gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1], hspace=0.42)

ax_a = fig.add_subplot(outer_gs[0])
ax_b = fig.add_subplot(outer_gs[1])


# ============================================================================
# PANEL A
# ============================================================================

group_centers_a = [i * A_GROUP_SPACING for i in range(len(all_fractions))]

y_range_a = Y_LIM_A[1] - Y_LIM_A[0]
y_ceiling_a = Y_LIM_A[1] - 0.12 * y_range_a

for i, fraction in enumerate(all_fractions):
    base_x = group_centers_a[i]

    for j, strand in enumerate(STRAND_ORDER):
        x = base_x + (j - 1) * A_STRAND_SPACING

        subset = pooled_df[
            (pooled_df["Fraction"] == fraction) & (pooled_df["Strand"] == strand)
            ]["Skew"]

        result = draw_group(
            ax_a,
            x,
            subset.values,
            PALETTE[strand],
            jitter_width=A_JITTER_WIDTH,
            violin_width=A_VIOLIN_WIDTH,
            box_width=A_BOX_WIDTH,
            box_offset=A_BOX_OFFSET,
            jitter_offset=A_JITTER_OFFSET,
        )

        if result is not None:
            label_stats(
                ax_a,
                *result,
                y_range=y_range_a,
                y_ceiling=y_ceiling_a,
            )

if SHOW_CONNECTOR_LINES_A:
    for j, strand in enumerate(STRAND_ORDER):
        medians = []
        for fraction in all_fractions:
            subset = pooled_df[
                (pooled_df["Fraction"] == fraction) & (pooled_df["Strand"] == strand)
                ]["Skew"]
            medians.append(subset.median() if len(subset) else np.nan)

        xs = [
            group_centers_a[i] + (j - 1) * A_STRAND_SPACING
            for i in range(len(all_fractions))
        ]

        ax_a.plot(
            xs,
            medians,
            color=PALETTE[strand],
            linewidth=1.0,
            marker="o",
            markersize=3,
            markeredgecolor="black",
            markeredgewidth=0.4,
            zorder=5,
        )

ax_a.set_xticks(group_centers_a)
ax_a.tick_params(axis="both", labelsize=TICK_LABEL_SIZE_A)
ax_a.set_xticklabels([f"{fraction * 100:g}" for fraction in all_fractions])
ax_a.set_xlabel("Sampling fraction (%)", fontsize=6, labelpad=6)
ax_a.set_ylabel(
    r"Deviation from parity   $\sqrt{AT_{skew}^2+GC_{skew}^2}$",
    fontsize=5,
    labelpad=10,
)
ax_a.set_ylim(*Y_LIM_A)

ax_a.text(-0.075, 1.08, "A", transform=ax_a.transAxes,
          fontsize=9, fontweight="bold", va="top")

legend_a_handles = [
    Rectangle((0, 0), 1, 1, facecolor=PALETTE[strand], alpha=0.85, edgecolor="none")
    for strand in STRAND_ORDER
]

ax_a.legend(
    legend_a_handles,
    STRAND_ORDER,
    frameon=False,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.275),
    ncol=3,
    handlelength=1.0,
    handleheight=1.0,
    columnspacing=1.5,
    borderaxespad=0,
    fontsize=6,
)

sns.despine(ax=ax_a, top=True, right=True)
ax_a.spines["left"].set_position(("outward", 2))
ax_a.spines["bottom"].set_position(("outward", 2))


# ============================================================================
# PANEL B
# ============================================================================

group_centers_b = [i * B_GROUP_SPACING for i in range(len(panel_b_fractions))]

source_center_offset = {
    "mammals_gene_te": -B_SOURCE_SPACING / 2,
    "mammals_gene": B_SOURCE_SPACING / 2,
}

strand_sub_offset = {
    "Plus": -B_STRAND_SPACING,
    "Minus": 0.0,
    "Mixed": B_STRAND_SPACING,
}

source_shades = {
    "mammals_gene_te": MAMMALS_GENE_TE_SHADES,
    "mammals_gene": MAMMALS_GENE_SHADES,
}

y_range_b = Y_LIM_B[1] - Y_LIM_B[0]
y_ceiling_b = Y_LIM_B[1] - 0.12 * y_range_b

PANEL_B_SOURCES = ["mammals_gene_te", "mammals_gene"]

for i, fraction in enumerate(panel_b_fractions):
    base_x = group_centers_b[i]

    for source in PANEL_B_SOURCES:
        for strand in STRAND_ORDER:
            x = base_x + source_center_offset[source] + strand_sub_offset[strand]

            subset = panel_b_df[
                (panel_b_df["Fraction"] == fraction)
                & (panel_b_df["Source"] == source)
                & (panel_b_df["Strand"] == strand)
                ]["Skew"]

            result = draw_group(
                ax_b,
                x,
                subset.values,
                source_shades[source][strand],
                jitter_width=B_JITTER_WIDTH,
                violin_width=B_VIOLIN_WIDTH,
                box_width=B_BOX_WIDTH,
                box_offset=B_BOX_OFFSET,
                jitter_offset=B_JITTER_OFFSET,
            )

            if result is not None:
                label_stats(
                    ax_b,
                    *result,
                    y_range=y_range_b,
                    y_ceiling=y_ceiling_b,
                    fontsize=4.2,
                )

if SHOW_CONNECTOR_LINES_B:
    for source in PANEL_B_SOURCES:
        for strand in STRAND_ORDER:
            medians = []
            for fraction in panel_b_fractions:
                subset = panel_b_df[
                    (panel_b_df["Fraction"] == fraction)
                    & (panel_b_df["Source"] == source)
                    & (panel_b_df["Strand"] == strand)
                    ]["Skew"]
                medians.append(subset.median() if len(subset) else np.nan)

            xs = [
                group_centers_b[i] + source_center_offset[source] + strand_sub_offset[strand]
                for i in range(len(panel_b_fractions))
            ]

            ax_b.plot(
                xs,
                medians,
                color=source_shades[source][strand],
                linewidth=1.0,
                marker="o",
                markersize=2.5,
                markeredgecolor="black",
                markeredgewidth=0.35,
                zorder=5,
            )

ax_b.set_xticks(group_centers_b)
ax_b.tick_params(axis="both", labelsize=TICK_LABEL_SIZE_B)
ax_b.set_xticklabels([f"{fraction * 100:g}" for fraction in panel_b_fractions])
ax_b.set_xlabel("Sampling fraction (%)", fontsize=6, labelpad=6)
ax_b.set_ylabel(
    r"Deviation from parity   $\sqrt{AT_{skew}^2+GC_{skew}^2}$",
    fontsize=5,
    labelpad=10,
)
ax_b.set_ylim(*Y_LIM_B)

ax_b.text(-0.075, 1.08, "B", transform=ax_b.transAxes,
          fontsize=8, fontweight="bold", va="top")

legend_b_handles = []
legend_b_labels = []

for source in PANEL_B_SOURCES:
    for strand in STRAND_ORDER:
        legend_b_handles.append(
            Rectangle((0, 0), 1, 1,
                      facecolor=source_shades[source][strand],
                      alpha=0.9, edgecolor="none")
        )
        legend_b_labels.append(f"{SOURCE_LABELS_SHORT[source]}-{strand}")

fig.legend(
    legend_b_handles,
    legend_b_labels,
    frameon=False,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.335),
    ncol=6,
    handlelength=1.0,
    handleheight=1.0,
    columnspacing=1.2,
    fontsize=6,
)

sns.despine(ax=ax_b, top=True, right=True)
ax_b.spines["left"].set_position(("outward", 2))
ax_b.spines["bottom"].set_position(("outward", 2))


# ============================================================================
# PANEL C
# ============================================================================

Y_LABEL_C = r"Deviation from parity   $\sqrt{AT_{skew}^2+GC_{skew}^2}$"

if PANEL_C_MODE == "combined":

    ax_c = fig.add_subplot(outer_gs[2])
    ax_c.tick_params(axis="both", labelsize=TICK_LABEL_SIZE_C)

    for group in GROUP_ORDER:
        for strand in STRAND_ORDER:
            sub = panel_c_summary[
                (panel_c_summary["Group"] == group)
                & (panel_c_summary["Strand"] == strand)
                ].sort_values("Fraction")

            if sub.empty:
                continue

            x = sub["Fraction"].to_numpy()
            y = sub["Mean"].to_numpy()
            error = sub["Error"].fillna(0).to_numpy()

            ax_c.plot(
                x, y,
                color=GROUP_COLORS[group],
                linestyle=STRAND_LINESTYLE[strand],
                marker=STRAND_MARKER[strand],
                markersize=PANEL_C_MARKERSIZE,
                linewidth=PANEL_C_LINEWIDTH,
                markeredgecolor="none",
                zorder=4,
            )

            ax_c.errorbar(
                x, y,
                yerr=error,
                fmt="none",
                ecolor=GROUP_COLORS[group],
                elinewidth=PANEL_C_ERROR_LINEWIDTH,
                capsize=PANEL_C_ERROR_CAPSIZE,
                capthick=PANEL_C_ERROR_LINEWIDTH,
                alpha=0.75,
                zorder=3,
            )

    ax_c.set_xlabel("Sampling fraction of annotation(%) ", fontsize=6, labelpad=6)
    ax_c.set_ylabel(Y_LABEL_C, fontsize=5, labelpad=10)
    ax_c.set_ylim(*Y_LIM_C)

    ax_c.text(-0.075, 1.08, "C", transform=ax_c.transAxes,
              fontsize=8, fontweight="bold", va="top")

    sns.despine(ax=ax_c, top=True, right=True)
    ax_c.spines["left"].set_position(("outward", 2))
    ax_c.spines["bottom"].set_position(("outward", 2))

    group_legend_handles = [
        Line2D([0], [0], color=GROUP_COLORS[group], linewidth=1.6)
        for group in GROUP_ORDER
    ]

    group_legend = ax_c.legend(
        group_legend_handles,
        GROUP_ORDER,
        title="Group",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.5, 1.0),
        handlelength=1.5,
        fontsize=6.8,
        title_fontsize=7.2,
    )
    ax_c.add_artist(group_legend)

    strand_legend_handles = [
        Line2D([0], [0], color="black",
               linestyle=STRAND_LINESTYLE[strand],
               marker=STRAND_MARKER[strand],
               markersize=4, linewidth=1.2)
        for strand in STRAND_ORDER
    ]

    ax_c.legend(
        strand_legend_handles,
        STRAND_ORDER,
        title="Strand",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.5, 0.55),
        handlelength=2.0,
        fontsize=6.8,
        title_fontsize=7.2,
    )

elif PANEL_C_MODE == "split_1x3":

    inner_gs = outer_gs[2].subgridspec(1, 3, wspace=0.15)
    axes_c = [fig.add_subplot(inner_gs[i]) for i in range(3)]

    for k, strand in enumerate(STRAND_ORDER):
        ax_k = axes_c[k]

        for group in GROUP_ORDER:
            sub = panel_c_summary[
                (panel_c_summary["Group"] == group)
                & (panel_c_summary["Strand"] == strand)
                ].sort_values("Fraction")

            if sub.empty:
                continue

            x = sub["Fraction"].to_numpy()
            y = sub["Mean"].to_numpy()
            error = sub["Error"].fillna(0).to_numpy()

            ax_k.plot(
                x, y,
                color=GROUP_COLORS[group],
                linestyle="-",
                marker="o",
                markersize=2,
                linewidth=0.8,
                markeredgecolor="none",
                zorder=4,
            )

            ax_k.errorbar(
                x, y,
                yerr=error,
                fmt="none",
                ecolor=GROUP_COLORS[group],
                elinewidth=0.65,
                capsize=2.0,
                capthick=0.65,
                alpha=0.75,
                zorder=3,
            )

        ax_k.tick_params(axis="both", labelsize=TICK_LABEL_SIZE_C)
        ax_k.set_title(strand, fontsize=7)
        ax_k.set_xlabel("Fraction of annotation sampled (%)", fontsize=6)

        # Fraction values 0-1 shown as percentages (data unchanged).
        ax_k.set_xticks(panel_c_fractions)
        ax_k.set_xticklabels([f"{x * 100:g}" for x in panel_c_fractions])

        ax_k.set_ylim(*Y_LIM_C)

        if k == 0:
            ax_k.set_ylabel(Y_LABEL_C, fontsize=5, labelpad=10)
        else:
            ax_k.set_yticklabels([])

        sns.despine(ax=ax_k, top=True, right=True)
        ax_k.spines["left"].set_position(("outward", 2))
        ax_k.spines["bottom"].set_position(("outward", 2))

    axes_c[0].text(-0.20, 1.08, "C", transform=axes_c[0].transAxes,
                   fontsize=9, fontweight="bold", va="top")

    group_legend_handles = [
        Line2D([0], [0], color=GROUP_COLORS[group], marker="o",
               markersize=4, linewidth=1.4)
        for group in GROUP_ORDER
    ]

    fig.legend(
        group_legend_handles,
        GROUP_ORDER,
        title="Group",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.80, 0.30),
        handlelength=1.4,
        fontsize=6,
        title_fontsize=7,
    )

else:
    raise ValueError(
        f"Unknown PANEL_C_MODE: {PANEL_C_MODE!r}. Use 'combined' or 'split_1x3'."
    )


# ============================================================================
# FINAL LAYOUT AND EXPORT
# ============================================================================

fig.tight_layout(rect=(0, 0.02, 0.82, 0.98))

os.makedirs(OUT_DIR, exist_ok=True)

base = os.path.join(OUT_DIR, "Figure8_panelABC_std")

fig.savefig(f"{base}.pdf", bbox_inches="tight")
fig.savefig(f"{base}.svg", bbox_inches="tight")
fig.savefig(f"{base}-{DPI_RASTER}.png", dpi=DPI_RASTER, bbox_inches="tight")

print(f"\n[saved] figure files: {base}.pdf / .svg / -{DPI_RASTER}.png")

# Uncomment if you want the interactive plot to appear.
# plt.show()