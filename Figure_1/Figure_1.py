import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)
# ══════════════════════════════════════════════════════════════
#                      IMPORT STATEMENTS
# ══════════════════════════════════════════════════════════════

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import string
from matplotlib.lines import Line2D

# ══════════════════════════════════════════════════════════════
#                        USER CONFIG
# ══════════════════════════════════════════════════════════════

# --- Colors ---
COLOR_P = '#C62828'   # + annotation: red
COLOR_C = '#1565C0'   # − annotation: blue
COLOR_W = '#2E7D32'   # Combined ± annotation: green
COLOR_U = '#7D7D7D'   # Unannotated: gray

# --- Palettes (can be customized) ---
row1_palette = [COLOR_P, COLOR_C, COLOR_W]

# + strand nucleotide composition
row_plus_palette = ['#FFCDD2', '#EF9A9A', '#E57373', '#C62828']

# − strand nucleotide composition
row_minus_palette = ['#BBDEFB', '#90CAF9', '#64B5F6', '#1565C0']

# Combined ± nucleotide composition
row_combined_palette = ['#C8E6C9', '#A5D6A7', '#81C784', '#2E7D32']

# + and − ratio rows
row_plus_ratio_palette = ['#EF5350', '#C62828']
row_minus_ratio_palette = ['#42A5F5', '#1565C0']

# Combined ± ratio row
row_combined_ratio_palette = ['#66BB6A', '#2E7D32']

# Unannotated composition and ratio rows
row_unannotated_palette = ['#E0E0E0', '#BDBDBD', '#9E9E9E', '#616161']
row_unannotated_ratio_palette = ['#FFD700', '#B8860B']

# ──────────────────────────────────────────────────────────────
#   FIGURE MARGINS & SPACING
# ──────────────────────────────────────────────────────────────
FIG_SIZE = (42, 50)  # Width, Height (adjust as needed)

# Margins (adjust these like in Figure 2)
MARGIN_LEFT   = 0.04
MARGIN_RIGHT  = 0.99
MARGIN_TOP    = 0.98
MARGIN_BOTTOM = 0.04

# Spacing between subplots
WSPACE = 0.25  # Horizontal space between columns
HSPACE = 0.30  # Vertical space between rows

# Legend position
LEGEND_ANCHOR_X = 0.5
LEGEND_ANCHOR_Y = 0.0

# ──────────────────────────────────────────────────────────────
#   SPINE SETTINGS (Keep invisible for publication)
# ──────────────────────────────────────────────────────────────
SHOW_TOP_SPINE = False    # Set to True if you want top spine
SHOW_RIGHT_SPINE = False  # Set to True if you want right spine
SPINE_COLOR = 'black'
SPINE_LINEWIDTH = 1.0

# ──────────────────────────────────────────────────────────────
#   GROUP LABEL ORIENTATION
# ──────────────────────────────────────────────────────────────
GROUP_LABEL_ROTATION = 90    # 0 = horizontal, 90 = vertical, 45 = angled
GROUP_LABEL_PAD = 45        # Distance from plot
GROUP_LABEL_HA = 'center'   # Horizontal alignment: 'center', 'left', 'right'

# ──────────────────────────────────────────────────────────────
#   FONT SIZES
# ──────────────────────────────────────────────────────────────
FONT_TITLE        = 21  # Feature column titles (top)
FONT_GROUP_LABEL  = 25  # Group names on the left
FONT_PANEL_LABEL  = 18  # A1, B1, etc.
FONT_TICK_X       = 20  # X-axis tick labels
FONT_TICK_Y       = 20  # Y-axis tick labels
FONT_LEGEND       = 25  # Legend font size
FONT_STATS        = 18  # Mean±std text inside violin plots

# ──────────────────────────────────────────────────────────────
#   VIOLIN PLOT SETTINGS
# ──────────────────────────────────────────────────────────────
VIOLIN_BW_ADJUST = 0.7
VIOLIN_CUT = 0
VIOLIN_LINEWIDTH = 1

BOX_WIDTH = 0.25
BOX_WHISKER_LINEWIDTH = 1
BOX_MEDIAN_LINEWIDTH = 1

# Scatter dot for mean
MEAN_DOT_SIZE = 25
ERRORBAR_CAPSIZE = 3
ERRORBAR_LINEWIDTH = 1

# ──────────────────────────────────────────────────────────────
#   Y-AXIS LIMITS - SCIENTIFIC RECOMMENDATION: Use SAME limits
#   for fair comparison across all organisms
# ──────────────────────────────────────────────────────────────
YLIM_ANNOTATION_PORTION = (0, 100)      # Column 1 - consistent
YLIM_COMPOSITION       = (10, 50)        # Columns 2, 4, 6 - consistent
YLIM_UNANNOTATED_COMP  = (0, 60)        # Column 8 - consistent
YLIM_RATIOS            = (0.4, 1.6)     # Columns 3, 5, 7, 9 - consistent

# ──────────────────────────────────────────────────────────────
#   DASHED REFERENCE LINE COLOR
# ──────────────────────────────────────────────────────────────
REF_LINE_COLOR_P = COLOR_P      # + annotation reference line
REF_LINE_COLOR_C = COLOR_C      # − annotation reference line
REF_LINE_COLOR_W = COLOR_W      # ± annotation reference line
REF_LINE_COLOR_U = '#B8860B'    # Unannotated reference line

# ──────────────────────────────────────────────────────────────
#   FEATURE NAMES (Column headers)
# ──────────────────────────────────────────────────────────────
FEATURE_NAMES = [
    "Annotation fraction",
    "+annotation composition",
    "+annotation ratios",
    "−annotation composition",
    "−annotation ratios",
    "Combined composition",
    "Combined ratios",
    "Unannotated composition",
    "Unannotated ratios"
]

# ──────────────────────────────────────────────────────────────
#   X-AXIS LABELS
# ──────────────────────────────────────────────────────────────
XLABELS_ANNOTATION = ["+", "−", "±"]
XLABELS_BASES = ["A", "T", "C", "G"]
XLABELS_RATIOS = ["A/T", "G/C"]

# ──────────────────────────────────────────────────────────────
#   GROUP ORDERING
# ──────────────────────────────────────────────────────────────
ORDER_TOP_GROUP = "Chordata"
ORDER_BOTTOM_GROUP = "Plasmids"

# ──────────────────────────────────────────────────────────────
#   SAVE SETTINGS
# ──────────────────────────────────────────────────────────────
SAVE_NAME = str(config.output_dir(1) / 'Fig_1')
SAVE_DPI = 300

# ══════════════════════════════════════════════════════════════
#                     END OF CONFIG
# ══════════════════════════════════════════════════════════════

# ======================================
# Load master file
# ======================================
file_path = config.require(config.MASTER_CSV, 'Copy the master table into Dataframes/ (see README.md).')

df_raw = pd.read_csv(file_path, low_memory=False)
df_raw = df_raw[df_raw['Group'] != 'Mammals'].copy()

# Convert all non-group columns to numeric safely
for col in df_raw.columns:
    if col != "Group":
        df_raw[col] = pd.to_numeric(df_raw[col], errors="coerce")

# ======================================
# Data builders
# ======================================

def normalize_to_100(comp):
    """Rescale each row so A + T + C + G = 100."""
    total = comp[['A', 'T', 'C', 'G']].sum(axis=1).replace(0, np.nan)
    return comp[['A', 'T', 'C', 'G']].div(total, axis=0) * 100


def build_row1(df):
    """Annotation coverage / portion."""
    P_anno = df['Int_P%'] + df['UTR_P%'] + df['CDS_P%'] + df['Share_P%']
    C_anno = df['Int_C%'] + df['UTR_C%'] + df['CDS_C%'] + df['Share_C%']
    Wh_anno = df['Wh_ano']

    return pd.DataFrame({
        'P_anno': P_anno,
        'C_anno': C_anno,
        'Wh_anno': Wh_anno
    })


def build_plus_composition(df):
    """Nucleotide composition of + annotated regions (sums to 100)."""
    comp = df[['annot_P_A', 'annot_P_T', 'annot_P_C', 'annot_P_G']].rename(
        columns={
            'annot_P_A': 'A',
            'annot_P_T': 'T',
            'annot_P_C': 'C',
            'annot_P_G': 'G'
        }
    )
    return normalize_to_100(comp)


def build_minus_composition(df):
    """Nucleotide composition of − annotated regions (sums to 100)."""
    comp = df[['annot_C_A', 'annot_C_T', 'annot_C_C', 'annot_C_G']].rename(
        columns={
            'annot_C_A': 'A',
            'annot_C_T': 'T',
            'annot_C_C': 'C',
            'annot_C_G': 'G'
        }
    )
    return normalize_to_100(comp)


def build_combined_composition(df):
    """
    Size-weighted nucleotide composition after combining
    + and − annotated regions.
    Each strand is normalized to 100 first, then weighted by its
    annotated size, so the result sums to 100.
    """
    P_size = df['Int_P%'] + df['UTR_P%'] + df['CDS_P%'] + df['Share_P%']
    C_size = df['Int_C%'] + df['UTR_C%'] + df['CDS_C%'] + df['Share_C%']

    total_annotated = (P_size + C_size).replace(0, np.nan)

    # Prevent division by zero
    P_portion = P_size / total_annotated
    C_portion = C_size / total_annotated

    plus_comp = build_plus_composition(df)
    minus_comp = build_minus_composition(df)

    result = pd.DataFrame(index=df.index)

    for base in ['A', 'T', 'C', 'G']:
        result[base] = (
                plus_comp[base] * P_portion
                + minus_comp[base] * C_portion
        )

    return normalize_to_100(result)


def build_ratios(composition_df):
    """
    Calculate Chargaff ratios directly from a composition DataFrame.
    Output columns: A/T and G/C.
    """
    result = pd.DataFrame(index=composition_df.index)

    result['A/T'] = (
            composition_df['A']
            / composition_df['T'].replace(0, np.nan)
    )

    result['G/C'] = (
            composition_df['G']
            / composition_df['C'].replace(0, np.nan)
    )

    return result.replace([np.inf, -np.inf], np.nan)


def build_unannotated_composition(df):
    """Nucleotide composition in unannotated + regions."""
    return df[['Unanot_P_A_ano', 'Unanot_P_T_ano',
               'Unanot_P_C_ano', 'Unanot_P_G_ano']].rename(
        columns={
            'Unanot_P_A_ano': 'A',
            'Unanot_P_T_ano': 'T',
            'Unanot_P_C_ano': 'C',
            'Unanot_P_G_ano': 'G'
        }
    )


# ======================================
# Violin plotting function
# ======================================

def plot_violin_with_stats(ax, data, palette, ylim):
    """Plot violin with boxplot, mean dot, errorbar, and stats text."""

    # Remove invalid observations before plotting
    data = data.replace([np.inf, -np.inf], np.nan)
    data = data.clip(lower=0)

    sns.violinplot(
        data=data,
        ax=ax,
        palette=palette,
        inner=None,
        cut=VIOLIN_CUT,
        bw_adjust=VIOLIN_BW_ADJUST,
        density_norm='width',
        linewidth=VIOLIN_LINEWIDTH
    )

    sns.boxplot(
        data=data,
        ax=ax,
        width=BOX_WIDTH,
        boxprops={'facecolor': 'none', 'edgecolor': 'black'},
        showcaps=False,
        showfliers=False,
        whiskerprops={'linewidth': BOX_WHISKER_LINEWIDTH},
        medianprops={'color': 'black', 'linewidth': BOX_MEDIAN_LINEWIDTH}
    )

    means = data.mean()
    stds = data.std()

    for i, col in enumerate(data.columns):

        mean_val = means[col]
        std_val = stds[col]

        if pd.notna(mean_val):

            ax.scatter(
                i, mean_val,
                color='black',
                s=MEAN_DOT_SIZE,
                zorder=5
            )

            ax.errorbar(
                i,
                mean_val,
                yerr=std_val,
                fmt='none',
                ecolor='black',
                capsize=ERRORBAR_CAPSIZE,
                linewidth=ERRORBAR_LINEWIDTH,
                zorder=4
            )

            ax.text(
                i,
                ylim[1] * 0.97,
                f"{mean_val:.2f}±{std_val:.2f}",
                ha='center',
                va='top',
                fontsize=FONT_STATS,
                rotation=45,
                bbox=dict(
                    boxstyle="round,pad=0.20",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.70
                )
            )

    ax.set_ylim(ylim)
    ax.grid(False)

    # ── Set tick font sizes ──
    ax.tick_params(axis='x', labelsize=FONT_TICK_X)
    ax.tick_params(axis='y', labelsize=FONT_TICK_Y)

    # ── Control spines ──
    if not SHOW_TOP_SPINE:
        ax.spines['top'].set_visible(False)
    else:
        ax.spines['top'].set_visible(True)
        ax.spines['top'].set_color(SPINE_COLOR)
        ax.spines['top'].set_linewidth(SPINE_LINEWIDTH)

    if not SHOW_RIGHT_SPINE:
        ax.spines['right'].set_visible(False)
    else:
        ax.spines['right'].set_visible(True)
        ax.spines['right'].set_color(SPINE_COLOR)
        ax.spines['right'].set_linewidth(SPINE_LINEWIDTH)

    # Always show left and bottom spines
    ax.spines['left'].set_visible(True)
    ax.spines['left'].set_color(SPINE_COLOR)
    ax.spines['left'].set_linewidth(SPINE_LINEWIDTH)
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set_color(SPINE_COLOR)
    ax.spines['bottom'].set_linewidth(SPINE_LINEWIDTH)


# ======================================
# Plot setup: ORGANISMS AS ROWS, FEATURES AS COLUMNS
# ======================================
sns.set_theme(style="whitegrid", font_scale=0.9)

# Get ordered groups
groups_original = list(df_raw['Group'].dropna().unique())

ordered_groups = []
if ORDER_TOP_GROUP in groups_original:
    ordered_groups.append(ORDER_TOP_GROUP)

ordered_groups.extend(
    g for g in groups_original
    if g not in {ORDER_TOP_GROUP, ORDER_BOTTOM_GROUP, "Plasmid"}
)

if ORDER_BOTTOM_GROUP in groups_original:
    ordered_groups.append(ORDER_BOTTOM_GROUP)
elif "Plasmid" in groups_original:
    ordered_groups.append("Plasmid")

groups = ordered_groups
n_groups = len(groups)
n_features = len(FEATURE_NAMES)

# Create figure
fig, axes = plt.subplots(
    n_groups,
    n_features,
    figsize=FIG_SIZE,
    sharey=False
)

# Make axes consistently 2-D
if n_groups == 1:
    axes = axes.reshape(1, n_features)


# ======================================
# Main plotting loop
# ======================================
for row, group in enumerate(groups):

    df_group = df_raw[df_raw['Group'] == group].copy()

    # -------------------------------------------------
    # Column 1: Annotation portion
    # -------------------------------------------------
    ax = axes[row, 0]
    plot_violin_with_stats(
        ax,
        build_row1(df_group),
        row1_palette,
        YLIM_ANNOTATION_PORTION
    )
    ax.set_xticklabels(XLABELS_ANNOTATION, fontsize=FONT_TICK_X)

    # -------------------------------------------------
    # Column 2: + annotated composition
    # -------------------------------------------------
    plus_comp = build_plus_composition(df_group)

    ax = axes[row, 1]
    plot_violin_with_stats(
        ax,
        plus_comp,
        row_plus_palette,
        YLIM_COMPOSITION
    )
    ax.set_xticklabels(XLABELS_BASES, fontsize=FONT_TICK_X)

    # -------------------------------------------------
    # Column 3: + annotated ratios
    # -------------------------------------------------
    ax = axes[row, 2]
    plot_violin_with_stats(
        ax,
        build_ratios(plus_comp),
        row_plus_ratio_palette,
        YLIM_RATIOS
    )
    ax.axhline(1, color=REF_LINE_COLOR_P, linestyle='--', linewidth=1.2)
    ax.set_xticklabels(XLABELS_RATIOS, fontsize=FONT_TICK_X)

    # -------------------------------------------------
    # Column 4: − annotated composition
    # -------------------------------------------------
    minus_comp = build_minus_composition(df_group)

    ax = axes[row, 3]
    plot_violin_with_stats(
        ax,
        minus_comp,
        row_minus_palette,
        YLIM_COMPOSITION
    )
    ax.set_xticklabels(XLABELS_BASES, fontsize=FONT_TICK_X)

    # -------------------------------------------------
    # Column 5: − annotated ratios
    # -------------------------------------------------
    ax = axes[row, 4]
    plot_violin_with_stats(
        ax,
        build_ratios(minus_comp),
        row_minus_ratio_palette,
        YLIM_RATIOS
    )
    ax.axhline(1, color=REF_LINE_COLOR_C, linestyle='--', linewidth=1.2)
    ax.set_xticklabels(XLABELS_RATIOS, fontsize=FONT_TICK_X)

    # -------------------------------------------------
    # Column 6: ± combined composition
    # -------------------------------------------------
    combined_comp = build_combined_composition(df_group)

    ax = axes[row, 5]
    plot_violin_with_stats(
        ax,
        combined_comp,
        row_combined_palette,
        YLIM_COMPOSITION
    )
    ax.set_xticklabels(XLABELS_BASES, fontsize=FONT_TICK_X)

    # -------------------------------------------------
    # Column 7: ± combined ratios
    # -------------------------------------------------
    ax = axes[row, 6]
    plot_violin_with_stats(
        ax,
        build_ratios(combined_comp),
        row_combined_ratio_palette,
        YLIM_RATIOS
    )
    ax.axhline(1, color=REF_LINE_COLOR_W, linestyle='--', linewidth=1.2)
    ax.set_xticklabels(XLABELS_RATIOS, fontsize=FONT_TICK_X)

    # -------------------------------------------------
    # Column 8: Unannotated composition
    # -------------------------------------------------
    unannotated_comp = build_unannotated_composition(df_group)

    ax = axes[row, 7]
    plot_violin_with_stats(
        ax,
        unannotated_comp,
        row_unannotated_palette,
        YLIM_UNANNOTATED_COMP
    )
    ax.set_xticklabels(XLABELS_BASES, fontsize=FONT_TICK_X)

    # -------------------------------------------------
    # Column 9: Unannotated ratios
    # -------------------------------------------------
    ax = axes[row, 8]
    plot_violin_with_stats(
        ax,
        build_ratios(unannotated_comp),
        row_unannotated_ratio_palette,
        YLIM_RATIOS
    )
    ax.axhline(1, color=REF_LINE_COLOR_U, linestyle='--', linewidth=1.2)
    ax.set_xticklabels(XLABELS_RATIOS, fontsize=FONT_TICK_X)

    # ── Organism/group name on the left ──
    axes[row, 0].set_ylabel(
        group,
        fontsize=FONT_GROUP_LABEL,
        fontweight="bold",
        rotation=GROUP_LABEL_ROTATION,
        labelpad=GROUP_LABEL_PAD,
        va="center",
        ha=GROUP_LABEL_HA
    )


# ======================================
# Feature labels at the TOP
# ======================================
for col, feature in enumerate(FEATURE_NAMES):
    axes[0, col].set_title(
        feature,
        weight="bold",
        fontsize=FONT_TITLE,
        pad=12
    )


# ======================================
# Panel labels (A1, B1, etc.)
# ======================================
for i in range(n_groups):
    for j in range(n_features):
        panel_label = f"{i + 1}{string.ascii_uppercase[j]}"
        axes[i, j].text(
            0.02,
            0.02,
            panel_label,
            transform=axes[i, j].transAxes,
            fontsize=FONT_PANEL_LABEL,
            fontweight='bold',
            ha='left',
            va='bottom'
        )


# ======================================
# Layout with margins
# ======================================
plt.subplots_adjust(
    left=MARGIN_LEFT,
    right=MARGIN_RIGHT,
    top=MARGIN_TOP,
    bottom=MARGIN_BOTTOM,
    wspace=WSPACE,
    hspace=HSPACE
)

# ══════════════════════════════════════════════════════════════
#                      LEGEND
# ══════════════════════════════════════════════════════════════

legend_elements = [
    # Annotation portion (Column 1)
    Line2D([0], [0], marker='s', color='w', markerfacecolor=COLOR_P,
           markersize=12, label='+ Annotation'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor=COLOR_C,
           markersize=12, label='- Annotation'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor=COLOR_W,
           markersize=12, label='± Annotation (combined)'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor=COLOR_U,
           markersize=12, label='Unannotated'),
]

fig.legend(
    handles=legend_elements,
    loc='lower center',
    ncol=4,
    fontsize=FONT_LEGEND,
    bbox_to_anchor=(LEGEND_ANCHOR_X, LEGEND_ANCHOR_Y),
    frameon=True,
    edgecolor='gray'
)

# ======================================
# Export
# ======================================
plt.savefig(f'{SAVE_NAME}-{SAVE_DPI}.tiff',
            dpi=SAVE_DPI,
            format='tiff',
            bbox_inches=None,
            facecolor='white',
            edgecolor='none',
            pil_kwargs={"compression": "tiff_lzw"})

plt.savefig(f'{SAVE_NAME}-{SAVE_DPI}.svg',
            format='svg',
            bbox_inches=None,
            facecolor='white',
            edgecolor='none')

plt.savefig(f'{SAVE_NAME}-{SAVE_DPI}.png',
            dpi=SAVE_DPI,
            format='png',
            bbox_inches=None,
            facecolor='white',
            edgecolor='none')

# plt.show()