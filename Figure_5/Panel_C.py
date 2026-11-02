# panel_c_heatmap.py
# Main heatmap with 8 rows: M, RC, C, R comparisons

import os
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# FILES
# ============================================================

MASTER_FILE = str(config.figure_input(5, '4k_kmers_normalised_MCRRC_PC0UW.csv'))
# Figure files are written into Figure_5/outputs/ (the inputs above are absolute paths)
os.chdir(config.output_dir(5))


# ============================================================
# CONDITIONS
# ============================================================

CONDITIONS = [
    "+M", "+C", "+R", "+RC",
    "-M", "-C", "-R", "-RC",
    "0M", "0C", "0R", "0RC",
    "UM", "UC", "UR", "URC",
    "WM", "WC", "WR", "WRC"
]

# ============================================================
# FUNCTION 1 - DETECT SEPARATOR
# ============================================================

def detect_separator(columns):
    first_col = str(columns[0])
    parts = first_col.split('_')

    if len(parts) >= 3:
        return '__'
    elif len(parts) == 2:
        return '_'
    else:
        return '_'

# ============================================================
# FUNCTION 2 - PREPARE FILE
# ============================================================

def prepare_file(file_path):
    print("\n" + "="*48)
    print("PROCESSING:", file_path)
    print("="*48)

    # Load data
    data = pd.read_csv(file_path, index_col=0)
    print("Original shape:", data.shape)

    # Remove Drosophila
    data = data[~data.index.astype(str).str.contains("Drosophila", case=False)]

    # Reset index
    data = data.reset_index(drop=True)
    print("After Drosophila removal:", data.shape)

    # Detect separator
    separator = detect_separator(data.columns)
    print("Separator:", repr(separator))

    # Get all unique k-mers
    columns = data.columns.tolist()
    unique_kmers = sorted(set(str(c).split(separator)[0] for c in columns))
    print("Unique k-mers:", len(unique_kmers))

    # Create chunks
    chunks = {condition: [] for condition in CONDITIONS}
    kmer_scores = []
    valid_kmers = []

    # Process each k-mer
    for kmer in unique_kmers:
        kmer_cols = [f"{kmer}{separator}{condition}" for condition in CONDITIONS]

        # Make sure all 20 columns exist
        if not all(col in data.columns for col in kmer_cols):
            continue

        # Extract this k-mer's 20 columns
        df_chunk = data[kmer_cols].copy()

        # Sort ALL ROWS by +M
        sort_col = f"{kmer}{separator}+M"
        df_chunk = df_chunk.sort_values(by=sort_col, ascending=True)

        # Calculate score for sorting k-mers
        score = df_chunk[sort_col].mean()
        kmer_scores.append(score)
        valid_kmers.append(kmer)

        # Reset index after sorting
        df_chunk = df_chunk.reset_index(drop=True)

        # Store each condition
        for condition in CONDITIONS:
            col = f"{kmer}{separator}{condition}"
            chunks[condition].append(df_chunk[[col]])

    print("Valid complete k-mers:", len(valid_kmers))

    # Sort k-mer chunks by +M score
    sorted_data = sorted(
        zip(kmer_scores, valid_kmers, *[chunks[condition] for condition in CONDITIONS]),
        key=lambda x: x[0]
    )

    # Extract sorted k-mer names
    sorted_kmers = [item[1] for item in sorted_data]

    # Extract each condition
    sorted_chunks = {}
    for condition_index, condition in enumerate(CONDITIONS):
        sorted_chunks[condition] = [item[condition_index + 2] for item in sorted_data]

    # Concatenate
    prepared = {}
    for condition in CONDITIONS:
        prepared[condition] = pd.concat(sorted_chunks[condition], axis=1)

    # Clean column names
    for condition in CONDITIONS:
        suffix = separator + condition
        prepared[condition].columns = [
            str(col)[:-len(suffix)] if str(col).endswith(suffix) else str(col)
            for col in prepared[condition].columns
        ]

    # Reorder columns according to sorted_kmers
    for condition in CONDITIONS:
        prepared[condition] = prepared[condition].reindex(columns=sorted_kmers)

    return {
        "data": prepared,
        "sorted_kmers": sorted_kmers,
        "n_rows": len(data)
    }

# ============================================================
# STEP 1 - PREPARE MASTER FILE
# ============================================================

master = prepare_file(MASTER_FILE)
master_data = master["data"]
master_sorted_kmers = master["sorted_kmers"]

print("\nMASTER SORTING")
print("Number of k-mers:", len(master_sorted_kmers))
print("First 10:", master_sorted_kmers[:10])

# ============================================================
# STEP 2 - CALCULATE MASTER PERCENTAGES
# ============================================================

df_plusM = master_data["+M"]
df_minusM = master_data["-M"]
df_plusRC = master_data["+RC"]
df_minusRC = master_data["-RC"]
df_plusC = master_data["+C"]
df_minusC = master_data["-C"]
df_plusR = master_data["+R"]
df_minusR = master_data["-R"]

# Use ALL rows
mean_plusM = df_plusM.mean()
mean_minusM = df_minusM.mean()
mean_plusRC = df_plusRC.mean()
mean_minusRC = df_minusRC.mean()
mean_plusC = df_plusC.mean()
mean_minusC = df_minusC.mean()
mean_plusR = df_plusR.mean()
mean_minusR = df_minusR.mean()

# Main 4-type normalization
total_sum = mean_plusM + mean_minusM + mean_plusRC + mean_minusRC
perc_plusM = (mean_plusM / total_sum) * 100
perc_minusM = (mean_minusM / total_sum) * 100
perc_plusRC = (mean_plusRC / total_sum) * 100
perc_minusRC = (mean_minusRC / total_sum) * 100

# C/R normalization
total_sum_C_R = mean_plusC + mean_minusC + mean_plusR + mean_minusR
perc_plusC = (mean_plusC / total_sum_C_R) * 100
perc_minusC = (mean_minusC / total_sum_C_R) * 100
perc_plusR = (mean_plusR / total_sum_C_R) * 100
perc_minusR = (mean_minusR / total_sum_C_R) * 100

# Master +M sorting - FINAL x-axis order
sorted_idx = perc_plusM.sort_values(ascending=False).index

# Reorder master values
perc_plusM = perc_plusM.reindex(sorted_idx)
perc_minusRC = perc_minusRC.reindex(sorted_idx)
perc_plusRC = perc_plusRC.reindex(sorted_idx)
perc_minusM = perc_minusM.reindex(sorted_idx)
perc_plusC = perc_plusC.reindex(sorted_idx)
perc_minusR = perc_minusR.reindex(sorted_idx)
perc_plusR = perc_plusR.reindex(sorted_idx)
perc_minusC = perc_minusC.reindex(sorted_idx)

# ============================================================
# STEP 3 - BUILD HEATMAP DATA
# ============================================================

heatmap_data = np.vstack([
    perc_plusM.values,
    perc_minusRC.values,
    perc_plusRC.values,
    perc_minusM.values,
    perc_plusC.values,
    perc_minusR.values,
    perc_plusR.values,
    perc_minusC.values
])

print("\n" + "="*40)
print("PANEL C - HEATMAP")
print("="*40)
print("Shape:", heatmap_data.shape)

# ============================================================
# STEP 4 - PLOT PANEL C
# ============================================================

n_rows, n_cols = heatmap_data.shape

# Gaps between rows
gap_heights = [0.04, 0.35, 0.04, 0.35, 0.04, 0.35, 0.04]
assert len(gap_heights) == n_rows - 1

# Build rows + gaps
row_heights = []
row_values = []

for i in range(n_rows):
    row_heights.append(1.0)
    row_values.append(heatmap_data[i, :])

    if i != n_rows - 1:
        row_heights.append(gap_heights[i])
        row_values.append(np.full(n_cols, np.nan))

# Matrix
data_stack = np.array(row_values)
masked = np.ma.masked_invalid(data_stack)

# Edges
x_edges = np.arange(n_cols + 1)
y_edges = np.concatenate([[0], np.cumsum(row_heights)])

# Figure
fig, ax = plt.subplots(figsize=(20, 4))
cmap = plt.get_cmap("viridis").copy()
cmap.set_bad(color='white')

mesh = ax.pcolormesh(x_edges, y_edges, masked, cmap=cmap, shading='flat', edgecolors='none')
ax.invert_yaxis()

# ============================================================
# COLORBAR WITH CUSTOM SIZE
# ============================================================

# Option 1: Use shrink parameter to control height (0.5 = 50% of default height)
cbar = fig.colorbar(mesh, ax=ax,
                    label='Percentage (%) in counting of 4 Type k-mer transition',
                    shrink=1,           # Reduce height to 60% of default
                    aspect=30,            # Control width (higher = thinner)
                    pad=0.02)             # Distance from the plot

# Option 2: If you want to control exact position with an axis
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="3%", pad=0.1)
# cbar = fig.colorbar(mesh, cax=cax, label='Percentage (%) in counting of 4 Type k-mer transition')

cbar.outline.set_visible(False)

# You can also adjust colorbar tick font size
cbar.ax.tick_params(labelsize=12)

# Adjust colorbar label font size
cbar.set_label('Percentage (%) in counting of 4 Type k-mer transition', fontsize=9)

# ============================================================

# Y tick labels
row_labels = ["+K", "-RC", "+RC", "-K", "+C", "-R", "+R", "-C"]
real_row_indices = list(range(0, len(row_heights), 2))
real_row_centers = [(y_edges[i] + y_edges[i+1]) / 2 for i in real_row_indices]

ax.set_yticks(real_row_centers)
ax.set_yticklabels(row_labels)

row_colors = ["red", "blue", "red", "blue", "red", "blue", "red", "blue"]
for label, color in zip(ax.get_yticklabels(), row_colors):
    label.set_color(color)
    label.set_fontweight('bold')

# Remove spines and ticks
ax.set_xticks([])
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(axis='y', length=0)
ax.tick_params(axis='x', length=0)

#plt.title("Panel C: K-mer Distribution - Main Comparison")
plt.xlabel("256 K-mers (sorted by +M)")



plt.tight_layout()
plt.savefig('Figure_5_Panel_C_300.png', dpi=300, bbox_inches='tight')
plt.savefig('Figure_5_Panel_C_600.png', dpi=600, bbox_inches='tight')
plt.savefig('Figure_5_Panel_C_600.svg', format = 'svg', bbox_inches='tight')
plt.savefig(f'Figure_5_Panel_C_300.tiff',dpi=300, format='tiff',facecolor='white', edgecolor='none',pil_kwargs={"compression": "tiff_lzw"})
plt.savefig('Figure_5_Panel_C.pdf', format='pdf', bbox_inches='tight')
plt.show()

print("\nPanel C saved as 'Panel_C_heatmap.png'")