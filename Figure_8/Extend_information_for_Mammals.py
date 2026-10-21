import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)
import pandas as pd
import numpy as np
import os

name = 'Mammals'
INPUT_CSV = str(config.output_dir(8, "Single_annote_composition") / f"{name}_gene_by_chromosome.csv")   # written by Genes_infor_extraction_Mammals.py
OUTPUT_CSV = str(config.output_dir(8, "ML") / f"{name}_gene_fractions_compostion_chromosome_by_chromosome.csv")

# ------------------------------------------------------------
# 1. Configuration & Metrics
# ------------------------------------------------------------
fractions = [1, 0.75, 0.50, 0.25, 0.10, 0.05, 0.01, 0.001]
metrics = ["A_perc", "T_perc", "G_perc", "C_perc", "A_div_T", "G_div_C", "AT_skew", "GC_skew"]

# Detect separator (Check if it's Tab or Comma)
with open(INPUT_CSV, 'r') as f:
    first_line = f.readline()
    sep = '\t' if '\t' in first_line else ','

def process_organism_group(df_org):
    """Perform the math for a single organism's data."""
    if df_org.empty: return None

    # Compatibility: Rename columns if necessary
    if "Chr_ID" not in df_org.columns and "seqid" in df_org.columns:
        df_org = df_org.rename(columns={"seqid": "Chr_ID"})

    if "length" not in df_org.columns:
        df_org["length"] = df_org["A_raw"] + df_org["T_raw"] + df_org["G_raw"] + df_org["C_raw"]

    # Only process plus strand
    df_plus = df_org[(df_org["strand"] == "+") & (df_org["length"] > 0)].copy()
    if df_plus.empty: return None

    org_results = []
    for i, frac in enumerate(fractions):
        if frac >= 1.0:
            sampled = df_plus.copy()
        else:
            # Group by Chr_ID to ensure subsampling is representative of each chromosome
            sampled = df_plus.groupby("Chr_ID", group_keys=False).apply(
                lambda grp: grp.sample(frac=frac, random_state=42 + i) if len(grp) > 0 else grp
            ).reset_index(drop=True)

        if sampled.empty: continue

        # Sum values per Chromosome
        chr_sums = sampled.groupby(["Organ_ID", "Chr_ID"]).sum(numeric_only=True)

        # Calculate Metrics
        chr_sums["A_perc"] = chr_sums["A_raw"] / chr_sums["length"] * 100
        chr_sums["T_perc"] = chr_sums["T_raw"] / chr_sums["length"] * 100
        chr_sums["G_perc"] = chr_sums["G_raw"] / chr_sums["length"] * 100
        chr_sums["C_perc"] = chr_sums["C_raw"] / chr_sums["length"] * 100
        chr_sums["A_div_T"] = np.where(chr_sums["T_raw"]!=0, chr_sums["A_raw"]/chr_sums["T_raw"], np.nan)
        chr_sums["G_div_C"] = np.where(chr_sums["C_raw"]!=0, chr_sums["G_raw"]/chr_sums["C_raw"], np.nan)

        at_sum = chr_sums["A_raw"] + chr_sums["T_raw"]
        gc_sum = chr_sums["G_raw"] + chr_sums["C_raw"]
        chr_sums["AT_skew"] = np.where(at_sum!=0, (chr_sums["A_raw"]-chr_sums["T_raw"])/at_sum, np.nan)
        chr_sums["GC_skew"] = np.where(gc_sum!=0, (chr_sums["G_raw"]-chr_sums["C_raw"])/gc_sum, np.nan)

        chr_sums["fraction"] = frac
        org_results.append(chr_sums.reset_index())

    if not org_results: return None

    all_org_frac = pd.concat(org_results, ignore_index=True)

    # Deviation from Baseline (Fraction 1.0)
    baseline = all_org_frac[all_org_frac["fraction"] == 1.0].set_index("Chr_ID")

    final_rows = []
    for chr_id, sub_data in all_org_frac.groupby("Chr_ID"):
        if chr_id not in baseline.index: continue
        base_vals = baseline.loc[chr_id, metrics]
        sub_data = sub_data.copy()
        for m in metrics:
            sub_data[f"diff_{m}"] = np.abs(sub_data[m] - base_vals[m])
        final_rows.append(sub_data)

    return pd.concat(final_rows, ignore_index=True)

# ------------------------------------------------------------
# 2. Memory-Efficient Loop (Grouped Streaming)
# ------------------------------------------------------------
_sep_name = 'tab' if sep == '\t' else 'comma'   # (backslash inside an f-string only works on Python >= 3.12)
print(f"Starting processing with separator: '{_sep_name}'")

if os.path.exists(OUTPUT_CSV):
    os.remove(OUTPUT_CSV)

first_write = True
# Using a larger chunksize to reduce overhead
reader = pd.read_csv(INPUT_CSV, sep=sep, chunksize=10**6)

# We will group by Organ_ID.
# Note: If Organ_IDs are sorted in your file, this is extremely efficient.
# If they are not sorted, this will still work but will process each appearance of an organism.
for chunk in reader:
    for organ_id, group_data in chunk.groupby("Organ_ID"):
        result_df = process_organism_group(group_data)

        if result_df is not None:
            result_df.to_csv(OUTPUT_CSV, mode='a', index=False, header=first_write)
            first_write = False

    print(f"Processed 1,000,000 rows...")

print(f"\nDone! Output saved to: {OUTPUT_CSV}")
