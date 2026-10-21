import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)
# herei will do sumation CDSs as chromosme by chromsome for eukaerute
# for prokaryute was organims by organims

# and fraction poulation in percentage of plus annotation








import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize

name = 'Mammals'

# ------------------------------------------------------------
# 1. Load data & keep + strand
# ------------------------------------------------------------
INPUT_CSV = str(config.output_dir(8, "Single_annote_composition") / f"{name}_gene_by_chromosome.csv")   # written by Genes_info_extraction.py
OUTPUT_CSV = str(config.output_dir(8, "ML") / f"{name}_gene_fractions_compostion_chromosome_by_chromosome.csv")

df = pd.read_csv(INPUT_CSV)

# ------------------------------------------------------------
# Compatibility: if old file has seqid instead of Chr_ID
# ------------------------------------------------------------
if "Chr_ID" not in df.columns and "seqid" in df.columns:
    df = df.rename(columns={"seqid": "Chr_ID"})

if "Chr_ID" not in df.columns:
    raise ValueError("Missing chromosome column. Expected 'Chr_ID' or old column 'seqid'.")

# ------------------------------------------------------------
# If length column is missing, reconstruct it
# ------------------------------------------------------------
if "length" not in df.columns:
    df["length"] = df["A_raw"] + df["T_raw"] + df["G_raw"] + df["C_raw"]

print("gene count per Organ_ID:")
print(df.value_counts("Organ_ID"))

print("\ngene count per Organ_ID/chromosome:")
print(df.value_counts(["Organ_ID", "Chr_ID"]))

# Keep only + strand
df_plus = df[df["strand"] == "+"].copy()

# Remove invalid rows if any
df_plus = df_plus[df_plus["length"] > 0].copy()

# ------------------------------------------------------------
# 2. Subsampling fractions
# ------------------------------------------------------------
fractions = [1, 0.75, 0.50, 0.25, 0.10, 0.05, 0.01, 0.001]

results = []

for i, frac in enumerate(fractions):

    print(f"\nProcessing fraction: {frac}")

    if frac == 1.0 or frac == 1:
        sampled = df_plus.copy()

    else:
        # Sample CDS records separately within each organism/chromosome
        # This makes subsampling chromosome by chromosome
        sampled = (
            df_plus
            .groupby(["Organ_ID", "Chr_ID"], group_keys=False)
            .apply(
                lambda grp: grp.sample(
                    frac=frac,
                    random_state=42 + i
                ) if len(grp) > 0 else grp
            )
            .reset_index(drop=True)
        )

    if sampled.empty:
        print(f"Fraction {frac} produced empty sampled data. Skipping.")
        continue

    # ------------------------------------------------------------
    # Aggregate per organism/chromosome
    # ------------------------------------------------------------
    chr_sums = (
        sampled
        .groupby(["Organ_ID", "Chr_ID"], dropna=False)
        .sum(numeric_only=True)
    )

    # Remove empty chromosomes
    chr_sums = chr_sums[chr_sums["length"] > 0].copy()

    if chr_sums.empty:
        print(f"No valid chromosome sums for fraction {frac}. Skipping.")
        continue

    # ------------------------------------------------------------
    # 3. Compute chromosome-level metrics
    # ------------------------------------------------------------

    # Base percentages over CDS length
    chr_sums["A_perc"] = chr_sums["A_raw"] / chr_sums["length"] * 100
    chr_sums["T_perc"] = chr_sums["T_raw"] / chr_sums["length"] * 100
    chr_sums["G_perc"] = chr_sums["G_raw"] / chr_sums["length"] * 100
    chr_sums["C_perc"] = chr_sums["C_raw"] / chr_sums["length"] * 100

    # Safe ratios
    chr_sums["A_div_T"] = np.where(
        chr_sums["T_raw"] != 0,
        chr_sums["A_raw"] / chr_sums["T_raw"],
        np.nan
    )

    chr_sums["G_div_C"] = np.where(
        chr_sums["C_raw"] != 0,
        chr_sums["G_raw"] / chr_sums["C_raw"],
        np.nan
    )

    # Skews
    AT_sum = chr_sums["A_raw"] + chr_sums["T_raw"]
    GC_sum = chr_sums["G_raw"] + chr_sums["C_raw"]

    chr_sums["AT_skew"] = np.where(
        AT_sum != 0,
        (chr_sums["A_raw"] - chr_sums["T_raw"]) / AT_sum,
        np.nan
    )

    chr_sums["GC_skew"] = np.where(
        GC_sum != 0,
        (chr_sums["G_raw"] - chr_sums["C_raw"]) / GC_sum,
        np.nan
    )

    chr_sums["fraction"] = frac

    # Convert organism and Chr_ID from index back to columns
    chr_sums = chr_sums.reset_index()

    results.append(chr_sums)

# ------------------------------------------------------------
# 4. Combine all fractions
# ------------------------------------------------------------
if not results:
    raise ValueError("No results generated. Check input data and fractions.")

all_frac = pd.concat(results, ignore_index=True)

print("\nCombined chromosome-by-chromosome fraction table:")
print(all_frac.head().to_string(index=False))

# ------------------------------------------------------------
# 5. Baseline fraction = 1 and deviations
# ------------------------------------------------------------
baseline = (
    all_frac[all_frac["fraction"] == 1.0]
    .set_index(["Organ_ID", "Chr_ID"])
)

metrics = [
    "A_perc",
    "T_perc",
    "G_perc",
    "C_perc",
    "A_div_T",
    "G_div_C",
    "AT_skew",
    "GC_skew"
]

diff_data = []

# For every organism/chromosome and fraction,
# compute absolute difference from its own chromosome baseline
for (organism, chr_id), sub_data in all_frac.groupby(["Organ_ID", "Chr_ID"], dropna=False):

    key = (organism, chr_id)

    if key not in baseline.index:
        continue

    base_vals = baseline.loc[key, metrics]

    sub_data = sub_data.copy()

    for m in metrics:
        sub_data[f"diff_{m}"] = np.abs(sub_data[m] - base_vals[m])

    diff_data.append(sub_data)

if not diff_data:
    raise ValueError("No difference data generated. Baseline may be missing.")

df_diff = pd.concat(diff_data, ignore_index=True)

# ------------------------------------------------------------
# 6. Save output
# ------------------------------------------------------------
df_diff.to_csv(OUTPUT_CSV, index=False)

print("\nFinal chromosome-by-chromosome deviation table:")
print(df_diff.head().to_string(index=False))

print(f"\nSaved output to:")
print(OUTPUT_CSV)

print("\nSummary:")
print(f"Total rows in Organ_ID + strand data: {len(df_plus)}")
print(f"Total Organ_ID: {df_plus['Organ_ID'].nunique()}")
print(f"Total Organ_ID/chromosome pairs: {df_plus[['Organ_ID', 'Chr_ID']].drop_duplicates().shape[0]}")
print(f"Total output rows: {len(df_diff)}")
