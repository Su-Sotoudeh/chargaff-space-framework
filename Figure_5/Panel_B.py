import os
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# --- STEP 1: LOAD DATA ---
file_path = str(config.figure_input(5, '4k_kmers_normalised_MCRRC_PC0UW.csv'))
# Figure files are written into Figure_5/outputs/ (the inputs above are absolute paths)
os.chdir(config.output_dir(5))

df = pd.read_csv(file_path, index_col=0)
print('shape is',df.shape)


print(df.shape[1])

for i in range(0,df.shape[1],20):
    #print(i,i+1,i+2,i+3)
    A = df.iloc[:,i+0]
    B = df.iloc[:,i+1]
    C = df.iloc[:,i+2]
    D = df.iloc[:,i+3]
    E = df.iloc[:,i+4]
    F = df.iloc[:,i+5]
    G = df.iloc[:,i+6]
    H = df.iloc[:,i+7]
    J = df.iloc[:,i+8]
    K = df.iloc[:,i+9]
    M = df.iloc[:,i+10]
    N = df.iloc[:,i+11]
    O = df.iloc[:,i+12]
    P = df.iloc[:,i+13]
    Q = df.iloc[:,i+14]
    R = df.iloc[:,i+15]
    S = df.iloc[:,i+16]
    T = df.iloc[:,i+17]
    U = df.iloc[:,i+18]
    Y = df.iloc[:,i+19]




    ABCD =A+B+C+D
    EFGH =E+F+G+H
    JKMN =J+K+M+N
    OPQR =O+P+Q+R
    STUY =S+T+U+Y

    A = A/ABCD *100
    B = B/ABCD *100
    C = C/ABCD *100
    D = D/ABCD *100

    E = E/EFGH *100
    F = F/EFGH *100
    G = G/EFGH *100
    H = H/EFGH *100

    J = J/JKMN *100
    K = K/JKMN *100
    M = M/JKMN *100
    N = N/JKMN *100

    O = O/OPQR *100
    P = P/OPQR *100
    Q = Q/OPQR *100
    R = R/OPQR *100

    S = S/STUY *100
    T = T/STUY *100
    U = U/STUY *100
    Y = Y/STUY *100


    df.iloc[:,i+0] = A
    df.iloc[:,i+1] = B
    df.iloc[:,i+2] = C
    df.iloc[:,i+3] = D
    df.iloc[:,i+4] = E
    df.iloc[:,i+5] = F
    df.iloc[:,i+6] = G
    df.iloc[:,i+7] = H
    df.iloc[:,i+8] = J
    df.iloc[:,i+9] = K
    df.iloc[:,i+10] = M
    df.iloc[:,i+11] = N
    df.iloc[:,i+12] = O
    df.iloc[:,i+13] = P
    df.iloc[:,i+14] = Q
    df.iloc[:,i+15] = R
    df.iloc[:,i+16] = S
    df.iloc[:,i+17] = T
    df.iloc[:,i+18] = U
    df.iloc[:,i+19] = Y



# --- STEP 2: DROP DROSOPHILA ---
# We filter based on the index before we reset it
df = df[~df.index.astype(str).str.contains("Drosophila", case=False)]

# --- CRITICAL FIX: RESET INDEX ---
# This converts "Org1_Chr1", "Org1_Chr2" into 0, 1, 2, 3...
# This ensures every row is treated as an independent unit for sorting.
df = df.reset_index(drop=True)

# --- STEP 3: SORTING LOGIC ---

# 1. Get all column names
cols = df.columns.tolist()

# 2. Detect separator
first_col = str(cols[0])
parts = first_col.split('_')

if len(parts) >= 3:
    separator = '__'
elif len(parts) == 2:
    separator = '_'
else:
    separator = '_'

# 3. Extract unique k-mers
unique_kmers = sorted(list(set([c.split(separator)[0] for c in cols])))
#print(unique_kmers)
# 4. Define conditions
conditions = ["+M","+C","+R", "+RC", "-M","-C","-R", "-RC","0M","0C","0R","0RC","UM","UC","UR","URC","WM","WC","WR","WRC"]

# 5. Create lists to hold chunks for EACH condition separately
chunks_plusM = []
chunks_plusC = []
chunks_plusR = []
chunks_plusRC = []
chunks_minusM = []
chunks_minusC = []
chunks_minusR = []
chunks_minusRC = []
chunks_0M = []
chunks_0C = []
chunks_0R = []
chunks_0RC = []
chunks_UM = []
chunks_UC = []
chunks_UR = []
chunks_URC = []
chunks_WM = []
chunks_WC = []
chunks_WR = []
chunks_WRC = []

# We need a list to store the scores (calculated from +M) to sort ALL lists together
kmer_scores = []

for kmer in unique_kmers:
    kmer_cols = [f"{kmer}{separator}{cond}" for cond in conditions]

    if all(col in df.columns for col in kmer_cols):
        # Copy the data for this k-mer
        df_chunk = df[kmer_cols].copy()

        # A. Sort ROWS by +M (This determines the row order for ALL subplots)
        sort_col_name = f"{kmer}{separator}+M"
        df_chunk = df_chunk.sort_values(by=sort_col_name, ascending=True)
        #print(df_chunk)


        # B. Calculate Score for COLUMN SORTING (Mean of ALL values of +M)
        # This calculates the average of the +M column for this k-mer
        top_values = df_chunk[sort_col_name].iloc[:]
        score = top_values.mean()
        kmer_scores.append(score)

        # C. Reset index to make rows independent
        # This is technically redundant here since we reset the main DF,
        # but it ensures safety if the main DF reset is removed later.
        df_chunk = df_chunk.reset_index(drop=True)

        # D. Split the chunk into 4 separate DataFrames (one per condition)
        col_plusM  = f"{kmer}{separator}+M"
        col_plusC  = f"{kmer}{separator}+C"
        col_plusR  = f"{kmer}{separator}+R"
        col_plusRC = f"{kmer}{separator}+RC"
        col_minusM = f"{kmer}{separator}-M"
        col_minusC = f"{kmer}{separator}-C"
        col_minusR = f"{kmer}{separator}-R"
        col_minusRC= f"{kmer}{separator}-RC"
        col_0M = f"{kmer}{separator}0M"
        col_0C = f"{kmer}{separator}0C"
        col_0R = f"{kmer}{separator}0R"
        col_0RC= f"{kmer}{separator}0RC"
        col_UM = f"{kmer}{separator}UM"
        col_UC = f"{kmer}{separator}UC"
        col_UR = f"{kmer}{separator}UR"
        col_URC= f"{kmer}{separator}URC"
        col_WM = f"{kmer}{separator}WM"
        col_WC = f"{kmer}{separator}WC"
        col_WR = f"{kmer}{separator}WR"
        col_WRC= f"{kmer}{separator}WRC"



        chunks_plusM.append(df_chunk[[col_plusM]])
        chunks_plusC.append(df_chunk[[col_plusC]])
        chunks_plusR.append(df_chunk[[col_plusR]])
        chunks_plusRC.append(df_chunk[[col_plusRC]])
        chunks_minusM.append(df_chunk[[col_minusM]])
        chunks_minusC.append(df_chunk[[col_minusC]])
        chunks_minusR.append(df_chunk[[col_minusR]])
        chunks_minusRC.append(df_chunk[[col_minusRC]])
        chunks_0M.append(df_chunk[[col_0M]])
        chunks_0C.append(df_chunk[[col_0C]])
        chunks_0R.append(df_chunk[[col_0R]])
        chunks_0RC.append(df_chunk[[col_0RC]])
        chunks_UM.append(df_chunk[[col_UM]])
        chunks_UC.append(df_chunk[[col_UC]])
        chunks_UR.append(df_chunk[[col_UR]])
        chunks_URC.append(df_chunk[[col_URC]])
        chunks_WM.append(df_chunk[[col_WM]])
        chunks_WC.append(df_chunk[[col_WC]])
        chunks_WR.append(df_chunk[[col_WR]])
        chunks_WRC.append(df_chunk[[col_WRC]])

# 6. Sort the CHUNKS based on the SCORES
# We zip scores with chunks, sort by score, and then unzip them back
# key=lambda x: x[0] ensures we sort by the SCORE (the first element in the pair)
sorted_data = sorted(zip(kmer_scores, chunks_plusM,chunks_plusC,chunks_plusR, chunks_plusRC
                         , chunks_minusM,chunks_minusC,chunks_minusR, chunks_minusRC), key=lambda x: x[0])

print(sorted_data)

sorted_data = sorted(zip(kmer_scores, chunks_plusM,chunks_plusC,chunks_plusR, chunks_plusRC, chunks_minusM,chunks_minusC,chunks_minusR, chunks_minusRC,chunks_0M,chunks_0C,chunks_0R,chunks_0RC,chunks_UM,chunks_UC,chunks_UR,chunks_URC,chunks_WM,chunks_WC,chunks_WR,chunks_WRC), key=lambda x: x[0])
print('hi',sorted_data)
# Unzip the sorted data back into separate lists
# _ is used for the score since we don't need it after sorting
_, sorted_plusM,sorted_plusC,sorted_plusR, sorted_plusRC, sorted_minusM,sorted_minusC,sorted_minusR, sorted_minusRC,sorted_0M,sorted_0C,sorted_0R,sorted_0RC,sorted_UM,sorted_UC,sorted_UR,sorted_URC,sorted_WM,sorted_WC,sorted_WR,sorted_WRC = zip(*sorted_data)
#print(sorted_plusM,sorted_plusC,sorted_plusR, sorted_plusRC, sorted_minusM,sorted_minusC,sorted_minusR, sorted_minusRC, sorted_0M,sorted_0C,sorted_0R,sorted_0RC)

# 7. Concatenate the sorted chunks for each condition
df_plusM   = pd.concat(list(sorted_plusM), axis=1)
df_plusC   = pd.concat(list(sorted_plusC), axis=1)
df_plusR   = pd.concat(list(sorted_plusR), axis=1)
df_plusRC  = pd.concat(list(sorted_plusRC), axis=1)
df_minusM  = pd.concat(list(sorted_minusM), axis=1)
df_minusC  = pd.concat(list(sorted_minusC), axis=1)
df_minusR  = pd.concat(list(sorted_minusR), axis=1)
df_minusRC = pd.concat(list(sorted_minusRC), axis=1)
df_0M  = pd.concat(list(sorted_0M), axis=1)
df_0C  = pd.concat(list(sorted_0C), axis=1)
df_0R  = pd.concat(list(sorted_0R), axis=1)
df_0RC = pd.concat(list(sorted_0RC), axis=1)
df_UM  = pd.concat(list(sorted_UM), axis=1)
df_UC  = pd.concat(list(sorted_UC), axis=1)
df_UR  = pd.concat(list(sorted_UR), axis=1)
df_URC = pd.concat(list(sorted_URC), axis=1)
df_WM  = pd.concat(list(sorted_WM), axis=1)
df_WC  = pd.concat(list(sorted_WC), axis=1)
df_WR  = pd.concat(list(sorted_WR), axis=1)
df_WRC = pd.concat(list(sorted_WRC), axis=1)



# ---------------------------------------------------------------
#   FULL SPEARMAN SIMILARITY WITH VERTICAL AXIS REVERSED
# ---------------------------------------------------------------

from scipy.stats import spearmanr
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# 1. Put all 20 processed matrices into a dictionary
# ---------------------------------------------------------------

condition_matrices = {
    "+M": df_plusM,
    "-M": df_minusM,
    "0M": df_0M,
    "UM": df_UM,
    "WM": df_WM,

    "+C": df_plusC,
    "-C": df_minusC,
    "0C": df_0C,
    "UC": df_UC,
    "WC": df_WC,

    "+R": df_plusR,
    "-R": df_minusR,
    "0R": df_0R,
    "UR": df_UR,
    "WR": df_WR,

    "+RC": df_plusRC,
    "-RC": df_minusRC,
    "0RC": df_0RC,
    "URC": df_URC,
    "WRC": df_WRC
}

# ---------------------------------------------------------------
# 2. Spearman similarity function
# ---------------------------------------------------------------

def matrix_similarity(A, B):
    a = A.values.flatten()
    b = B.values.flatten()
    r, _ = spearmanr(a, b)
    return r

# ---------------------------------------------------------------
# 3. Compute similarity
# ---------------------------------------------------------------

names = list(condition_matrices.keys())
N = len(names)

sim = np.zeros((N, N))

for i, name1 in enumerate(names):
    for j, name2 in enumerate(names):
        sim[i, j] = matrix_similarity(condition_matrices[name1],
                                      condition_matrices[name2])

sim_df = pd.DataFrame(sim, index=names, columns=names)
print(sim_df.to_string(index=False))
# ---------------------------------------------------------------
# 4. Biological ordering (+,-,0,U,W) × (M,C,R,RC)
# ---------------------------------------------------------------

ordered_conditions = [
    "+M","-M","0M","UM","WM",
    "+C","-C","0C","UC","WC",
    "+R","-R","0R","UR","WR",
    "+RC","-RC","0RC","URC","WRC"
]

# Reverse only the vertical axis (index)
reversed_vertical = ordered_conditions[::-1]

# Apply ordering
sim_df = sim_df.loc[reversed_vertical, ordered_conditions]

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm



# ---------------------------------------------------------------
# 5. Plot final heatmap with vertical reversed and annotations
# ---------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, Normalize

# ---------------------------------------------------------------
# CUSTOM NONLINEAR COLOR MAP
#
# Negative side:
# -1.00  = dark blue
# -0.94  = medium blue
# -0.80  = light blue
#
# Central region:
# -0.80 → +0.80 = white
#
# Positive side:
# +0.80  = light red
# +0.94  = medium red
# +1.00  = dark red
# ---------------------------------------------------------------

# Important values are deliberately given more visual separation.
#
# Position in the colormap:
# -1.00 → 0.00
# -0.94 → 0.06
# -0.80 → 0.20
# +0.80 → 0.80
# +0.94 → 0.94
# +1.00 → 1.00

colors = [
    (0.00, "#0000A8"),   # -1.00  dark blue

    (0.06, "#4169E1"),   # -0.94  medium blue

    (0.20, "#D6E4FF"),   # -0.80  light blue

    (0.20, "#FFFFFF"),   # -0.80 → white
    (0.80, "#FFFFFF"),   # +0.80 → white

    (0.80, "#FFE0E0"),   # +0.80  light red

    (0.94, "#E94B4B"),   # +0.94  medium red

    (1.00, "#B00000")    # +1.00  dark red
]

custom_cmap = LinearSegmentedColormap.from_list(
    "Chargaff_Spearman_Custom",
    colors
)

# ---------------------------------------------------------------
# NORMALIZATION
# ---------------------------------------------------------------

norm = Normalize(
    vmin=-1,
    vmax=1
)

# ---------------------------------------------------------------
# CREATE FIGURE
# ---------------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(14, 12)
)

# ---------------------------------------------------------------
# CREATE MASK FOR ANNOTATIONS
# Only annotate the M box (bottom-left)
# ---------------------------------------------------------------

mask_annot = np.zeros_like(
    sim_df,
    dtype=bool
)

mask_annot[15:20, 0:5] = True

# ---------------------------------------------------------------
# ANNOTATION VALUES
# ---------------------------------------------------------------

annot_values = np.where(
    mask_annot,
    np.round(sim_df, 2),
    ""
).astype(str)

# ---------------------------------------------------------------
# DRAW HEATMAP
# ---------------------------------------------------------------

sns.heatmap(
    sim_df,

    cmap=custom_cmap,
    norm=norm,

    linecolor="black",
    ax=ax,

    annot=annot_values,
    fmt="",

    annot_kws={
        "size": 9
    },

    square=True,

    cbar_kws={
        "shrink": 0.8
    }
)

# ---------------------------------------------------------------
# THIN GRID LINES
# ---------------------------------------------------------------

for i in range(21):

    ax.axhline(
        i,
        color="black",
        linewidth=0.4,
        alpha=0.2
    )

    ax.axvline(
        i,
        color="black",
        linewidth=0.4,
        alpha=0.2
    )

# ---------------------------------------------------------------
# THICK BLOCK SEPARATORS
# ---------------------------------------------------------------

for i in range(0, 21, 5):

    ax.axhline(
        i,
        color="black",
        linewidth=2.5
    )

    ax.axvline(
        i,
        color="black",
        linewidth=2.5
    )

# ---------------------------------------------------------------
# DRAW BLACK LINES INSIDE THE BOTTOM-LEFT M BOX
# ---------------------------------------------------------------

# 1. Vertical line separating:
#    +/- (x = 0,1)
#    from 0/U/W (x = 2,3,4)

ax.plot(
    [2, 2],
    [15, 20],
    color="black",
    linewidth=1.5
)

# 2. Horizontal line separating:
#    +/- (y = 15,16)
#    from 0/U/W (y = 17,18,19)

ax.plot(
    [0, 5],
    [18, 18],
    color="black",
    linewidth=1.5
)

# ---------------------------------------------------------------
# TITLE
# ---------------------------------------------------------------

plt.title(
    "Spearman Similarity with Condition Grid",
    fontsize=18,
    pad=20
)

plt.tight_layout()

# ---------------------------------------------------------------
# SAVE THE IMAGE
# ---------------------------------------------------------------

plt.savefig(
    '20_type_4mer_heatmaplined_spearman_annotated-300.png',
    dpi=300,
    bbox_inches='tight')
plt.savefig(
    '20_type_4mer_heatmaplined_spearman_annotated-600.png',
    dpi=600,
    bbox_inches='tight')
plt.savefig(
    '20_type_4mer_heatmaplined_spearman_annotated.svg',
    bbox_inches='tight')


plt.savefig(f'20_type_4mer_heatmaplined_spearman_annotated.tiff',
            dpi=600,
            format='tiff',
            bbox_inches=None,
            facecolor='white',
            edgecolor='none',
            pil_kwargs={"compression": "tiff_lzw"})

plt.show()


# ---------------------------------------------------------------
# 6. OPTIONAL: Block borders
# ---------------------------------------------------------------

import numpy as np

# Create a mask for annotations
# This creates a boolean mask where True means "hide the annotation"
# We want to hide everything EXCEPT where value is 1 or -1

fig, ax = plt.subplots(figsize=(12, 10))
'''sns.heatmap(sim_df,
            cmap="coolwarm",
            vmin=-1,
            vmax=1,
            ax=ax,
            annot=np.where((sim_df >= 0.99) | (sim_df <= -0.99), np.round(sim_df, 2), ""),
            fmt="")'''


sns.heatmap(sim_df, cmap="coolwarm", vmin=-1, vmax=1, ax=ax,annot=True)


plt.title("Spearman Similarity with Condition Grid")
plt.tight_layout()
#plt.savefig('/media/genomics/OS/Genomics/Chargaff/Outputs/single_genome/20_type_4mer_heatmap_spearman.png',dpi= 900)
plt.show()






