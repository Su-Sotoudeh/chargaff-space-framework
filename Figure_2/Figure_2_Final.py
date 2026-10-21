import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

# ══════════════════════════════════════════════════════════════
#                        USER CONFIG
# ══════════════════════════════════════════════════════════════

# --- Colors ---
COLOR_P = '#C62828'
COLOR_C = '#1565C0'
COLOR_A = '#2E7D32'
COLOR_W = '#000000'

# --- Scatter dots ---
DOT_SIZE = 0.9

# --- Mean lines ---
LINE_WIDTH = 1.8
LINE_ALPHA = 1.0

# --- Std shaded band ---
BAND_ALPHA = 0.22

# ──────────────────────────────────────────────────────────────
#   Y-AXIS LIMITS — edit each panel individually
#   Format: 'Name': {'A2T': (ymin, ymax), 'G2C': (ymin, ymax)}
# ──────────────────────────────────────────────────────────────
DEFAULT_YLIM_ORG = (0.80, 1.20)
DEFAULT_YLIM_VP  = (0.50, 1.80)

YLIM = {
    # ── Organisms (left two columns) ──────────────────────
    'Chordata'     : {'A2T': (0.80, 1.20), 'G2C': (0.80, 1.20)},
    'Viridiplantae': {'A2T': (0.80, 1.20), 'G2C': (0.80, 1.20)},
    'Arthropoda'   : {'A2T': (0.80, 1.20), 'G2C': (0.80, 1.20)},
    'Fungi'        : {'A2T': (0.80, 1.20), 'G2C': (0.80, 1.20)},
    'Archaea'      : {'A2T': (0.5, 1.5), 'G2C': (0.5, 1.5)},
    'Bacteria'     : {'A2T': (0.65, 1.5), 'G2C': (0.65, 1.5)},
    'Chloroplast'  : {'A2T': (0.7, 1.4), 'G2C': (0.7, 1.4)},
    'Mitochondria' : {'A2T': (0.2, 2), 'G2C': (0.2, 2)},

    # ── Viruses / Plasmids / Simulation (right two columns)
    '+ RNA Virus'  : {'A2T': (0.50, 1.80), 'G2C': (0.50, 1.80)},
    '- RNA Virus'  : {'A2T': (0.50, 1.2), 'G2C': (0.50, 1.2)},
    '+/- RNA Virus': {'A2T': (0.50, 1.4), 'G2C': (0.50, 1.4)},
    '+ DNA Virus'  : {'A2T': (0.50, 1.80), 'G2C': (0.50, 1.80)},
    '- DNA Virus'  : {'A2T': (0.2, 2), 'G2C': (0.2, 2)},
    '+/- DNA Virus': {'A2T': (0.50, 1.80), 'G2C': (0.50, 1.80)},
    'Plasmids'     : {'A2T': (0.4, 1.80), 'G2C': (0.4, 1.80)},
    'Simulation'   : {'A2T': (0.6, 1.7), 'G2C': (0.6, 1.7)},
}

# ──────────────────────────────────────────────────────────────
#   FIGURE MARGINS & SPACING
# ──────────────────────────────────────────────────────────────
FIG_SIZE      = (42, 50)

COL_WIDTHS    = [1, 1, 0.07, 1, 1]
ROW_HEIGHTS   = [1, 1, 1, 1, 1, 1, 1, 1]

HSPACE        = 0.5
WSPACE        = 0.4

MARGIN_LEFT   = 0.05
MARGIN_RIGHT  = 0.98
MARGIN_TOP    = 0.97
MARGIN_BOTTOM = 0.02

# Extra room under the last row for the x-axis labels + legend.
# (Raise it if the x labels touch the legend; lower it for less white space.)
GRIDSPEC_BOTTOM_PAD = 0.055
LEGEND_ANCHOR_Y     = 0.01   # legend sits just above bottom edge

# ──────────────────────────────────────────────────────────────
#   FONT SIZES
# ──────────────────────────────────────────────────────────────
FONT_PANEL_TITLE = 26   # group names (vertical labels on the left)
FONT_COL_HEADER  = 26
FONT_TICK        = 22
FONT_LEGEND      = 25
FONT_XLABEL      = 22   # x-axis label (bottom row only)

# ──────────────────────────────────────────────────────────────
#   AXIS LABEL SETTINGS
# ──────────────────────────────────────────────────────────────
YLABEL_PAD = 22    # distance of group names from the axis (was 8)
XLABEL_PAD = 10    # distance of x-axis label from the tick labels

# x-axis label text (shown on the bottom row only)
XLABEL_A2T = "Chromosomes (sorted by A/T(+) − A/T(−))"
XLABEL_G2C = "Chromosomes (sorted by G/C(+) − G/C(−))"

# ══════════════════════════════════════════════════════════════
#                     HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════

def get_ylim(name, ratio, is_org=True):
    default = DEFAULT_YLIM_ORG if is_org else DEFAULT_YLIM_VP
    return YLIM.get(name, {}).get(ratio, default)


def add_mean_band(ax, series, color):
    vals = series.dropna()
    if len(vals) == 0:
        return
    mean = vals.mean()
    std  = vals.std()
    ax.axhspan(mean - std, mean + std,
               color=color, alpha=BAND_ALPHA, linewidth=0, zorder=1)
    ax.axhline(mean, color=color, linewidth=LINE_WIDTH,
               linestyle='-', alpha=LINE_ALPHA, zorder=3)


def plot_panel(ax, df, cols_colors, sort_col,
               title, ylim, show_ylabel=True, xlabel=None):
    dfs = df.sort_values(sort_col).reset_index(drop=True) \
        if sort_col in df.columns else df.reset_index(drop=True)
    x = np.arange(len(dfs))

    # Reference line y=1
    ax.axhline(1, linestyle='--', linewidth=1.0,
               color='gray', alpha=0.6, zorder=1)

    # Mean bands first
    for col, color in cols_colors:
        if col in dfs.columns:
            add_mean_band(ax, dfs[col], color)

    # Scatter on top
    for col, color in cols_colors:
        if col in dfs.columns:
            ax.scatter(x, dfs[col], color=color,
                       s=DOT_SIZE, zorder=4, rasterized=True)

    ax.set_ylim(*ylim)
    ax.set_xlim(0, max(len(x) - 1, 1))

    # Group name (vertical label on the left)
    if show_ylabel:
        ax.set_ylabel(title, fontsize=FONT_PANEL_TITLE,
                      fontweight='bold', labelpad=YLABEL_PAD)
    else:
        ax.set_ylabel('')

    # x-axis label: only when a text is given (bottom row)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=FONT_XLABEL, labelpad=XLABEL_PAD)
    else:
        ax.set_xlabel('')

    ax.tick_params(axis='both', labelsize=FONT_TICK)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# ══════════════════════════════════════════════════════════════
#                       LOAD DATASETS
# ══════════════════════════════════════════════════════════════

df1 = pd.read_csv(
    config.require(config.MASTER_CSV, 'Copy the master table into Dataframes/ (see README.md).'),
    low_memory=False
)
df1 = df1[df1['Group'] != 'Mammals'].copy()
df1['A2T_P-A2T_C'] = df1['A2T_P'] - df1['A2T_C']
df1['G2C_P-G2C_C'] = df1['G2C_P'] - df1['G2C_C']

name_mapping = {
    'Chordata'     : 'Chordata',
    'Viride'       : 'Viridiplantae',
    'Plantae'      : 'Viridiplantae',
    'Viridiplantae': 'Viridiplantae',
    'Insect'       : 'Arthropoda',
    'Fungi'        : 'Fungi',
    'Archia'       : 'Archaea',
    'Bacteria'     : 'Bacteria',
    'Chloro_plast' : 'Chloroplast',
    'Mito'         : 'Mitochondria',
}
df1['Group'] = df1['Group'].replace(name_mapping)
df1['Prefix'] = df1['Chr_ID'].str[:2]

good_prefixes = [
    'AC','AE','AL','AM','AP','AY','BA','BX','CM','CP','CR','FN','FR','GK',
    'HE','HF','HG','L7','LN','LR','LT','NC','NZ','OU','OV','OW','OX','OY','OZ'
]
df1 = df1[df1['Prefix'].isin(good_prefixes)]
df1 = df1[(df1['A2T_W'] > 0.9) & (df1['A2T_W'] < 1.1)]
target_groups = df1['Group'].unique()[:8]

# Debug — verify group names from your CSV
print("Detected groups:", list(target_groups))

v_main = pd.read_csv(
    config.require(config.VIRUS_CSV),
    low_memory=False
)
v_ref = pd.read_csv(
    config.require(config.VIRUS_REF_CSV)
)[['Assembly', 'Molecule_type']]

df2 = v_main.merge(v_ref, left_on='Name', right_on='Assembly', how='inner')
df2['Genome'] = df2['Molecule_type']
df2['A2T_P-A2T_C'] = df2['A2T_P'] - df2['A2T_C']
df2['G2C_P-G2C_C'] = df2['G2C_P'] - df2['G2C_C']

conditions = {
    '+'  : (df2['annot_P_A'] != 0) & (df2['annot_C_A'] == 0),
    '-'  : (df2['annot_P_A'] == 0) & (df2['annot_C_A'] != 0),
    '+/-': (df2['annot_P_A'] != 0) & (df2['annot_C_A'] != 0),
}
virus_layout = [
    ('RNA', '+'), ('RNA', '-'), ('RNA', '+/-'),
    ('DNA', '+'), ('DNA', '-'), ('DNA', '+/-'),
]

df3 = pd.read_csv(
    config.require(config.PLASMID_CSV)
)
df3 = df3[(df3['annot_P_A'] != 0) & (df3['annot_C_A'] != 0)]
df3['A2T_P-A2T_C'] = df3['A2T_P'] - df3['A2T_C']
df3['G2C_P-G2C_C'] = df3['G2C_P'] - df3['G2C_C']

results_df = pd.read_csv(
    config.require(config.SIMULATION_CSV, 'Copy simulation_results.csv (artificial-chromosome simulation) into Dataframes/.')
)
results_df['A2T_P-A2T_C'] = results_df['A2T_P'] - results_df['A2T_C']
results_df['G2C_P-G2C_C'] = results_df['G2C_P'] - results_df['G2C_C']

# ══════════════════════════════════════════════════════════════
#                    BUILD PANEL LISTS
# ══════════════════════════════════════════════════════════════

right_datasets = []
for genome, cond in virus_layout:
    sub = df2[
        conditions[cond] & df2['Genome'].str.contains(genome, na=False)
        ].copy()
    right_datasets.append((f"{cond} {genome} Virus", sub))
right_datasets.append(('Plasmids',   df3))
right_datasets.append(('Simulation', results_df))

org_cols_A2T = [('A2T_P', COLOR_P), ('A2T_C', COLOR_C),
                ('A2T_W', COLOR_W), ('A2T_annot', COLOR_A)]
org_cols_G2C = [('G2C_P', COLOR_P), ('G2C_C', COLOR_C),
                ('G2C_W', COLOR_W), ('G2C_annot', COLOR_A)]
vp_cols_A2T  = [('A2T_P', COLOR_P), ('A2T_C', COLOR_C),
                ('A2T_W', COLOR_W), ('A2T_annot', COLOR_A)]
vp_cols_G2C  = [('G2C_P', COLOR_P), ('G2C_C', COLOR_C),
                ('G2C_W', COLOR_W), ('G2C_annot', COLOR_A)]
sim_cols_A2T = [('A2T_P', COLOR_P), ('A2T_C', COLOR_C),
                ('A2T_anno', COLOR_A)]
sim_cols_G2C = [('G2C_P', COLOR_P), ('G2C_C', COLOR_C),
                ('G2C_anno', COLOR_A)]

# ══════════════════════════════════════════════════════════════
#                      BUILD FIGURE
# ══════════════════════════════════════════════════════════════

GRID_BOTTOM = MARGIN_BOTTOM + GRIDSPEC_BOTTOM_PAD   # bottom edge of the last row
N_ROWS = 8
LAST_ROW = N_ROWS - 1

fig = plt.figure(figsize=FIG_SIZE)

fig.subplots_adjust(
    left   = MARGIN_LEFT,
    right  = MARGIN_RIGHT,
    bottom = MARGIN_BOTTOM,
    top    = MARGIN_TOP
)

gs = GridSpec(
    N_ROWS, 5,
    figure        = fig,
    width_ratios  = COL_WIDTHS,
    height_ratios = ROW_HEIGHTS,
    hspace        = HSPACE,
    wspace        = WSPACE,
    left          = MARGIN_LEFT,
    right         = MARGIN_RIGHT,
    bottom        = GRID_BOTTOM,   # room for x labels + legend
    top           = MARGIN_TOP
)

# ── 8 Rows ───────────────────────────────────────────────────
for row in range(N_ROWS):

    group = target_groups[row]
    dfg   = df1[df1['Group'] == group]

    # x-axis labels only on the bottom row
    xl_a2t = XLABEL_A2T if row == LAST_ROW else None
    xl_g2c = XLABEL_G2C if row == LAST_ROW else None

    # Org A2T — col 0
    ax = fig.add_subplot(gs[row, 0])
    plot_panel(ax, dfg, org_cols_A2T, 'A2T_P-A2T_C',
               group,
               get_ylim(group, 'A2T', is_org=True),
               show_ylabel=True, xlabel=xl_a2t)
    if row == 0:
        ax.set_title('A/T Ratio', fontsize=FONT_COL_HEADER,
                     fontweight='bold', pad=12)

    # Org G2C — col 1
    ax = fig.add_subplot(gs[row, 1])
    plot_panel(ax, dfg, org_cols_G2C, 'G2C_P-G2C_C',
               group,
               get_ylim(group, 'G2C', is_org=True),
               show_ylabel=False, xlabel=xl_g2c)
    if row == 0:
        ax.set_title('G/C Ratio', fontsize=FONT_COL_HEADER,
                     fontweight='bold', pad=12)

    # col 2 = spacer, intentionally empty

    # Right panels (viruses / plasmids / simulation)
    r_title, r_df = right_datasets[row]
    is_sim  = (r_title == 'Simulation')
    ra_cols = sim_cols_A2T if is_sim else vp_cols_A2T
    rg_cols = sim_cols_G2C if is_sim else vp_cols_G2C
    sa = 'A2T_P-A2T_C'
    sg = 'G2C_P-G2C_C'

    # Right A2T — col 3
    ax = fig.add_subplot(gs[row, 3])
    plot_panel(ax, r_df, ra_cols, sa,
               r_title,
               get_ylim(r_title, 'A2T', is_org=False),
               show_ylabel=True, xlabel=xl_a2t)
    if row == 0:
        ax.set_title('A2T Ratio', fontsize=FONT_COL_HEADER,
                     fontweight='bold', pad=12)

    # Right G2C — col 4
    ax = fig.add_subplot(gs[row, 4])
    plot_panel(ax, r_df, rg_cols, sg,
               r_title,
               get_ylim(r_title, 'G2C', is_org=False),
               show_ylabel=False, xlabel=xl_g2c)
    if row == 0:
        ax.set_title('G2C Ratio', fontsize=FONT_COL_HEADER,
                     fontweight='bold', pad=12)

# ── Vertical separator ────────────────────────────────────────
# Midpoint of spacer col 2 in figure coordinates
x_sep = (MARGIN_LEFT + (MARGIN_RIGHT - MARGIN_LEFT) *
         (COL_WIDTHS[0] + COL_WIDTHS[1] + COL_WIDTHS[2] * 0.5) /
         sum(COL_WIDTHS))

fig.add_artist(
    plt.Line2D(
        [x_sep, x_sep], [GRID_BOTTOM, MARGIN_TOP],
        transform=fig.transFigure,
        color='black', linewidth=2.0,
        linestyle='--', alpha=0.40
    )
)

# ══════════════════════════════════════════════════════════════
#                      GLOBAL LEGEND
# ══════════════════════════════════════════════════════════════

LEGEND_MARKER_SIZE = 25     # dot size in the legend
LEGEND_LINE_WIDTH  = 3.5     # mean-line width in the legend

legend_elements = [
    # ── COLUMN 1: Dots (Red and Blue) ──
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_P,
           markersize=LEGEND_MARKER_SIZE, label='+Annotation A/T or G/C'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_C,
           markersize=LEGEND_MARKER_SIZE, label='-Annotation A/T or G/C'),

    # ── COLUMN 2: Dots (Black and Green) ──
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_W,
           markersize=LEGEND_MARKER_SIZE, label='Whole Sequence A/T or G/C'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_A,
           markersize=LEGEND_MARKER_SIZE, label='±Annotation A/T or G/C'),

    # ── COLUMN 3: Mean Lines (Red and Blue) ──
    Line2D([0], [0], color=COLOR_P, linewidth=LEGEND_LINE_WIDTH,
           label='Mean +Annotation'),
    Line2D([0], [0], color=COLOR_C, linewidth=LEGEND_LINE_WIDTH,
           label='Mean -Annotation'),

    # ── COLUMN 4: Std Bands (Red and Blue) ──
    plt.Rectangle((0, 0), 1, 1, facecolor=COLOR_P, alpha=BAND_ALPHA,
                  edgecolor='none', label='Std +Annotation (A/T or G/C)'),
    plt.Rectangle((0, 0), 1, 1, facecolor=COLOR_C, alpha=BAND_ALPHA,
                  edgecolor='none', label='Std -Annotation (A/T or G/C)'),

    # ── COLUMN 5: Std Bands (Black and Green) ──
    plt.Rectangle((0, 0), 1, 1, facecolor=COLOR_W, alpha=BAND_ALPHA,
                  edgecolor='none', label='Std Whole Sequence (A/T or G/C)'),
    plt.Rectangle((0, 0), 1, 1, facecolor=COLOR_A, alpha=BAND_ALPHA,
                  edgecolor='none', label='Std ± Annotation (A/T or G/C)'),
]

fig.legend(
    handles        = legend_elements,
    loc            = 'lower center',
    ncol           = 5,              # 5 columns x 2 rows
    fontsize       = FONT_LEGEND,
    bbox_to_anchor = (0.5, LEGEND_ANCHOR_Y),
    frameon        = True,
    edgecolor      = 'gray',
    columnspacing  = 1.2,            # space between columns
    handlelength   = 2.5,            # length of line handles
    handletextpad  = 1.0,            # space between symbol and text
    labelspacing   = 0.8,            # space between rows
)



# ══════════════════════════════════════════════════════════════
#   SAVE FOR PUBLICATION
# ══════════════════════════════════════════════════════════════

SAVE_NAME = str(config.output_dir(2) / 'Fig_2')
DPI = 300

# ── 1. PDF — vector, best for journals (uncomment to use) ─────
'''plt.savefig(f'{SAVE_NAME}.pdf',
            dpi=600,                  # DPI only affects rasterized elements inside
            format='pdf',
            bbox_inches=None,         # respect your margins exactly
            facecolor='white',
            edgecolor='none')'''

# ── 2. TIFF — most journals accept this as primary format ─────
plt.savefig(f'{SAVE_NAME}-{DPI}.tiff',
            dpi=DPI,
            format='tiff',
            bbox_inches=None,
            facecolor='white',
            edgecolor='none',
            pil_kwargs={"compression": "tiff_lzw"})  # lossless compression

# ── 3. SVG — pure vector, no DPI needed ───────────────────────
plt.savefig(f'{SAVE_NAME}.svg',
            format='svg',
            bbox_inches=None,
            facecolor='white',
            edgecolor='none')

# ── 4. PNG — for presentations / supplementary / sharing ──────
plt.savefig(f'{SAVE_NAME}-{DPI}.png',
            dpi=DPI,
            format='png',
            bbox_inches=None,
            facecolor='white',
            edgecolor='none')

#plt.show()