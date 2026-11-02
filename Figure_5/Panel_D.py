# panel_d_heatmap.py
# Heatmap showing ONLY repeat element rows: CDS, UTR, INT, SINE, LINE, LTR, DNA

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

ELEMENT_FILES = {
    "+M-CDS":  '4k_kmers_normalised_MCRRC_PC0UW_CDS.csv',
    "+M-UTR":  '4k_kmers_normalised_MCRRC_PC0UW_UTR.csv',
    "+M-INT":  '4k_kmers_normalised_MCRRC_PC0UW_INT.csv',
    "+M-SINE": '4k_kmers_normalised_MCRRC_PC0UW_SINE.csv',
    "+M-LINE": '4k_kmers_normalised_MCRRC_PC0UW_LINE.csv',
    "+M-LTR":  '4k_kmers_normalised_MCRRC_PC0UW_LTR.csv',
    "+M-DNA":  '4k_kmers_normalised_MCRRC_PC0UW_DNA.csv'
}
# class tables are written by Figure_5_<class>.py into Figure_5/outputs/
ELEMENT_FILES = {k: str(config.figure_input(5, v)) for k, v in ELEMENT_FILES.items()}
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
# STEP 1 - PREPARE MASTER FILE (for sorting order)
# ============================================================

master = prepare_file(MASTER_FILE)
master_data = master["data"]
master_sorted_kmers = master["sorted_kmers"]

print("\nMASTER SORTING")
print("Number of k-mers:", len(master_sorted_kmers))
print("First 10:", master_sorted_kmers[:10])

# ============================================================
# STEP 2 - CALCULATE MASTER +M PERCENTAGES FOR SORTING
# ============================================================

df_plusM = master_data["+M"]
df_minusM = master_data["-M"]
df_plusRC = master_data["+RC"]
df_minusRC = master_data["-RC"]

mean_plusM = df_plusM.mean()
mean_minusM = df_minusM.mean()
mean_plusRC = df_plusRC.mean()
mean_minusRC = df_minusRC.mean()

total_sum = mean_plusM + mean_minusM + mean_plusRC + mean_minusRC
perc_plusM = (mean_plusM / total_sum) * 100

# Master +M sorting - FINAL x-axis order
sorted_idx = perc_plusM.sort_values(ascending=False).index

# ============================================================
# STEP 3 - PREPARE ELEMENT FILES
# ============================================================

element_results = {}
for label, file_path in ELEMENT_FILES.items():
    result = prepare_file(file_path)
    element_results[label] = result

# ============================================================
# STEP 4 - CALCULATE +M PERCENTAGE FOR EACH ELEMENT
# ============================================================

element_percentages = {}

for label, result in element_results.items():
    data = result["data"]

    plusM = data["+M"]
    minusM = data["-M"]
    plusRC = data["+RC"]
    minusRC = data["-RC"]

    mean_plusM_element = plusM.mean()
    mean_minusM_element = minusM.mean()
    mean_plusRC_element = plusRC.mean()
    mean_minusRC_element = minusRC.mean()

    total_element = mean_plusM_element + mean_minusM_element + mean_plusRC_element + mean_minusRC_element
    perc_element_plusM = (mean_plusM_element / total_element) * 100

    # Use MASTER +M ordering
    perc_element_plusM = perc_element_plusM.reindex(sorted_idx)
    element_percentages[label] = perc_element_plusM

    print(f"\n{label}")
    print("k-mers:", len(perc_element_plusM))
    print("non-NaN:", perc_element_plusM.notna().sum())
    print("sum:", perc_element_plusM.sum())

# ============================================================
# STEP 5 - BUILD HEATMAP DATA (ONLY ELEMENT ROWS)
# ============================================================

heatmap_data = np.vstack([
    element_percentages["+M-CDS"].values,
    element_percentages["+M-UTR"].values,
    element_percentages["+M-INT"].values,
    element_percentages["+M-SINE"].values,
    element_percentages["+M-LINE"].values,
    element_percentages["+M-LTR"].values,
    element_percentages["+M-DNA"].values
])

print("\n" + "="*40)
print("PANEL D - REPEAT ELEMENTS HEATMAP")
print("="*40)
print("Shape:", heatmap_data.shape)

# ============================================================
# STEP 6 - PLOT PANEL D (EACH ROW WITH ITS OWN COLOR SCALE)
# ============================================================

n_rows, n_cols = heatmap_data.shape

# No gaps between rows for panel D
gap_heights = [0.02, 0.02, 0.02, 0.02, 0.02, 0.02]
assert len(gap_heights) == n_rows - 1

# Build rows + gaps
row_heights = []
row_values = []
row_start_index = []

entry_idx = 0
for i in range(n_rows):
    row_start_index.append(entry_idx)
    row_heights.append(1.0)
    row_values.append(heatmap_data[i, :])
    entry_idx += 1

    if i != n_rows - 1:
        row_heights.append(gap_heights[i])
        row_values.append(np.full(n_cols, np.nan))
        entry_idx += 1

# Matrix
data_stack = np.array(row_values)
masked = np.ma.masked_invalid(data_stack)

# Edges
x_edges = np.arange(n_cols + 1)
y_edges = np.concatenate([[0], np.cumsum(row_heights)])
y_total = y_edges[-1]

# Figure - make room for the vertical label on the right
fig, ax = plt.subplots(figsize=(20, 5))
cmap = plt.get_cmap("viridis").copy()
cmap.set_bad(color='white')

# Plot each element row with its own color scale
element_meshes = {}
for i in range(n_rows):
    y_row_top = y_edges[row_start_index[i]]
    y_row_bottom = y_edges[row_start_index[i] + 1]
    row_masked = masked[row_start_index[i]:row_start_index[i] + 1, :]

    mesh_i = ax.pcolormesh(x_edges, [y_row_top, y_row_bottom], row_masked,
                           cmap=cmap, shading='flat', edgecolors='none')
    element_meshes[i] = mesh_i

ax.set_xlim(x_edges[0], x_edges[-1])
ax.set_ylim(y_edges[0], y_edges[-1])
ax.invert_yaxis()

# Y labels
row_labels = ["+K-CDS", "+K-UTR", "+K-INT", "+K-SINE", "+K-LINE", "+K-LTR", "+K-DNA"]

real_row_indices = list(range(0, len(row_heights), 2))
real_row_centers = [(y_edges[i] + y_edges[i+1]) / 2 for i in real_row_indices]

ax.set_yticks(real_row_centers)
ax.set_yticklabels(row_labels)

# All labels in black (or you can color-code them)
for label in ax.get_yticklabels():
    label.set_color('black')
    label.set_fontweight('bold')

# Remove spines and ticks
ax.set_xticks([])
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(axis='y', length=0)
ax.tick_params(axis='x', length=0)

plt.xlabel("256 K-mers (sorted by +M)")

plt.tight_layout()

# ============================================================
# STEP 7 - ADD INDEPENDENT COLORBARS FOR EACH ROW
# ============================================================

fig.canvas.draw()
pos = ax.get_position()

CBAR_WIDTH = 0.015
CBAR_PAD = 0.015
ROW_CBAR_SHRINK = 0.85

def data_y_to_fig_y(y_val):
    frac = 1 - (y_val / y_total)
    return pos.y0 + frac * pos.height

# One colorbar per row
for i in range(n_rows):
    fig_top_i = data_y_to_fig_y(y_edges[row_start_index[i]])
    fig_bottom_i = data_y_to_fig_y(y_edges[row_start_index[i] + 1])

    slot_height = fig_top_i - fig_bottom_i
    cbar_height = slot_height * ROW_CBAR_SHRINK
    margin = (slot_height - cbar_height) / 2

    cax_i = fig.add_axes([pos.x1 + CBAR_PAD, fig_bottom_i + margin, CBAR_WIDTH, cbar_height])
    cbar_i = fig.colorbar(element_meshes[i], cax=cax_i)
    cbar_i.outline.set_visible(False)
    cbar_i.ax.tick_params(labelsize=5, pad=1)

    # Show min and max values on each colorbar
    row_vmin = element_meshes[i].norm.vmin
    row_vmax = element_meshes[i].norm.vmax
    cbar_i.set_ticks([row_vmin, row_vmax])
    cbar_i.set_ticklabels([f"{row_vmin:.1f}", f"{row_vmax:.1f}"])

# ============================================================
# STEP 8 - ADD VERTICAL LABEL TO THE RIGHT OF COLORBARS
# ============================================================

# Calculate the position for the vertical label
# Place it to the right of all colorbars
label_x = pos.x1 + CBAR_PAD + CBAR_WIDTH + 0.02  # Slightly to the right of colorbars
label_y_center = (pos.y0 + pos.y1) / 2  # Center of the plot vertically

# Add vertical label
fig.text(label_x, label_y_center,
         'Percentage (%) of +K on 4-type transition',
         rotation=90,
         va='center',
         ha='center',
         fontsize=10)

# ============================================================
# SAVE
# ============================================================

plt.savefig('Figure_5_Panel_D_300.png', dpi=300, bbox_inches='tight')
plt.savefig('Figure_5_Panel_D_600.png', dpi=600, bbox_inches='tight')
plt.savefig('Figure_5_Panel_D.svg', format='svg', bbox_inches='tight')
plt.savefig('Figure_5_Panel_D.pdf', format='pdf', bbox_inches='tight')
plt.savefig(f'Figure_5_Panel_D_300.tiff',dpi=300, format='tiff',facecolor='white', edgecolor='none',pil_kwargs={"compression": "tiff_lzw"})

plt.show()

print("\nPanel D saved as 'Figure_5_Panel_D_300.png' and 'Figure_5_Panel_D_600.png'")











# ============================================================
# STEP 4b - PRINT + SAVE THE K-MER ORDER (paste right after STEP 4, before STEP 5)
# Needs: sorted_idx, perc_plusM (from STEP 2) and element_percentages (from STEP 4).
# ============================================================

N_SHOW = 10   # how many k-mers to print at each end of an order


def show_order(title, kmer_list):
    kmer_list = list(kmer_list)
    print(f"\n{title}")
    print(f"  first {N_SHOW} (highest): {kmer_list[:N_SHOW]}")
    print(f"  last  {N_SHOW} (lowest) : {kmer_list[-N_SHOW:]}")


master_order = list(sorted_idx)                       # the REAL x-axis order of Panel D (left -> right)
master_share = perc_plusM.reindex(sorted_idx)

print("\n" + "=" * 60)
print("K-MER ORDER USED IN PANEL D")
print("=" * 60)
print("Default (alphabetical) order?  ", master_order == sorted(master_order))
show_order("MASTER ORDER (used for ALL rows, x-axis left -> right)", master_order)

order_table = pd.DataFrame({
    "master_rank": range(1, len(master_order) + 1),
    "kmer": master_order,
    "master_+M_share": master_share.values,
})

for label, share in element_percentages.items():
    # share is already in master order; sort it by the class's OWN values (tie-break: alphabetical)
    own = share.dropna().sort_index().sort_values(ascending=False, kind="stable")
    own_order = list(own.index)

    # rank agreement between this class and the master (Spearman = Pearson on ranks, no scipy needed)
    both = pd.concat([master_share, share], axis=1, keys=["m", "c"]).dropna()
    rho = both["m"].rank().corr(both["c"].rank())
    top_overlap = len(set(own_order[:N_SHOW]) & set(master_order[:N_SHOW]))

    show_order(f"{label}: its OWN order (sorted by its own +M share)", own_order)
    print(f"  same as master order? {own_order == master_order}"
          f" | Spearman vs master = {rho:.3f}"
          f" | top-{N_SHOW} overlap with master = {top_overlap}/{N_SHOW}")

    order_table[f"{label}_share"] = share.values
    order_table[f"{label}_own_rank"] = share.rank(ascending=False, method="first").values

order_table.to_csv("Panel_D_kmer_order_per_class.csv", index=False)
print("\nSaved: Panel_D_kmer_order_per_class.csv  (master rank, master share, and each class's share + own rank per k-mer)")