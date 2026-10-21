"""
01_load_validate.py
--------------------
Load the genome dataset, build Group2 taxonomic/viral subgroup labels,
and run defensive validation (missing columns, NaNs, dtype checks)
before any modeling happens.

NEW: Organ_ID (organism) and Chr_ID (chromosome) are now REQUIRED, because
train/test splitting is done by organism. Rows without an Organ_ID are
dropped (they cannot be assigned to a fold safely) and reported.

Output:
    - df_clean : validated DataFrame with Group2 column
    - validation report printed and saved to ../tables/00_validation_report.txt
"""

import pandas as pd
import numpy as np
import warnings
import importlib

warnings.filterwarnings("ignore", category=RuntimeWarning)
gs = importlib.import_module("grouped_split")
ORGAN_COL, CHR_COL = gs.ORGAN_COL, gs.CHR_COL

# ----------------------------------------------------------------------
# CONFIG -- input table. Default: Dataframes/Cleaned_datafram_whole_organims_2026_07.csv
# (see config.py; override with the CSPR_MASTER_CSV environment variable or --csv)
# ----------------------------------------------------------------------
_paths = importlib.import_module("paths")   # paths.py resolves everything through config.py
CSV_PATH = _paths.MASTER_CSV

ID_COLS = [ORGAN_COL, CHR_COL]
TEXT_COLS = ["Group", "Genome_Type", "Coding"] + ID_COLS
REQUIRED_COLS = TEXT_COLS + [
    "A_W", "T_W", "G_W", "C_W",
    "Z", "P_portion", "C_portion",
    "annot_P_A", "annot_P_T", "annot_P_G", "annot_P_C",
    "annot_C_A", "annot_C_T", "annot_C_G", "annot_C_C",
]

GROUP_ORDER = [
    "Chordata", "Arthropoda", "Plantae", "Fungi", "Archaea", "Bacteria",
    "Chloroplast", "Mitochondria", "Plasmids",
    "DNA Viruses(+/-)", "DNA Viruses(+)", "DNA Viruses(-)",
    "RNA Viruses(+/-)", "RNA Viruses(+)", "RNA Viruses(-)"
]

MIN_GROUP_SIZE = 3  # groups smaller than this are dropped (matches your original logic)


def virus_subgroup(row):
    if row["Group"] != "Viruses":
        return row["Group"]
    gt = str(row["Genome_Type"]).upper()
    cd = str(row["Coding"]).replace(" ", "")
    if "RNA" in gt:
        return "RNA Viruses(+)" if cd == "+" else "RNA Viruses(-)" if cd == "-" else "RNA Viruses(+/-)"
    if "DNA" in gt:
        return "DNA Viruses(+)" if cd == "+" else "DNA Viruses(-)" if cd == "-" else "DNA Viruses(+/-)"
    return "Unknown"


def load_and_validate(csv_path=CSV_PATH, verbose=True):
    report_lines = []

    def log(msg):
        if verbose:
            print(msg)
        report_lines.append(msg)

    log("=" * 70)
    log("DATA VALIDATION REPORT")
    log("=" * 70)

    # ---- 1. Load ----
    try:
        df = pd.read_csv(csv_path)
        df = df[df["Group"] != "Mammals"]
        log(f"[OK] Loaded {csv_path}")
        log(f"     Shape: {df.shape[0]} rows x {df.shape[1]} columns")
    except FileNotFoundError:
        raise FileNotFoundError(
            f"\n[FATAL] Could not find file at:\n  {csv_path}\n"
            f"Put the master table in the Dataframes/ folder (see README.md), set CSPR_MASTER_CSV, or pass --csv."
        )

    # ---- 2. Required columns present? ----
    missing_cols = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing_cols:
        log(f"[FATAL] Missing required columns: {missing_cols}")
        log(f"        Columns found in file: {list(df.columns)}")
        raise ValueError(f"Missing required columns: {missing_cols}")
    else:
        log(f"[OK] All {len(REQUIRED_COLS)} required columns present (incl. {ORGAN_COL}, {CHR_COL}).")

    # ---- 3. Coerce numerics defensively ----
    numeric_cols = [c for c in REQUIRED_COLS if c not in TEXT_COLS]
    for c in numeric_cols:
        before_na = df[c].isna().sum()
        df[c] = pd.to_numeric(df[c], errors="coerce")
        newly_bad = df[c].isna().sum() - before_na
        if newly_bad > 0:
            log(f"[WARN] Column '{c}': {newly_bad} non-numeric values coerced to NaN")

    # ---- 3b. IDs: organism id is mandatory for the organism-level split ----
    n_noorg = df[ORGAN_COL].isna().sum()
    if n_noorg:
        log(f"[WARN] {n_noorg} rows have no {ORGAN_COL}; dropped (cannot be assigned to a CV fold).")
        df = df[df[ORGAN_COL].notna()].copy()
    df[ORGAN_COL] = df[ORGAN_COL].astype(str)
    df[CHR_COL] = df[CHR_COL].astype(str)
    n_dup = df.duplicated(subset=[ORGAN_COL, CHR_COL]).sum()
    if n_dup:
        log(f"[WARN] {n_dup} duplicated ({ORGAN_COL}, {CHR_COL}) pairs -- check for repeated records.")

    # ---- 4. NaN report per column ----
    log("-" * 70)
    log("NaN counts per numeric column (after coercion):")
    for c in numeric_cols:
        n_na = df[c].isna().sum()
        pct = 100 * n_na / len(df)
        flag = "  <-- check this" if pct > 5 else ""
        log(f"  {c:15s}: {n_na:6d} NaN  ({pct:5.1f}%){flag}")

    # ---- 5. Build Group2 ----
    df["Group2"] = df.apply(virus_subgroup, axis=1)
    log("-" * 70)
    log("Group2 value counts  (rows / organisms / max chromosomes per organism):")
    for g in GROUP_ORDER:
        sub = df[df["Group2"] == g]
        n = len(sub)
        n_org = sub[ORGAN_COL].nunique()
        mx = sub.groupby(ORGAN_COL).size().max() if n else 0
        flag = "  <-- below MIN_GROUP_SIZE, will be dropped" if n < MIN_GROUP_SIZE else ""
        log(f"  {g:20s}: {n:6d} rows  {n_org:6d} organisms  max {mx} chr/organism{flag}")
    shared = df.groupby(ORGAN_COL)["Group2"].nunique()
    n_shared = int((shared > 1).sum())
    if n_shared:
        log(f"[INFO] {n_shared} {ORGAN_COL} values occur in more than one Group2 (e.g. nuclear + organelle "
            f"of the same organism). Splits and scores are done per Group2, so this is fine.")

    unknown_n = (df["Group2"] == "Unknown").sum()
    if unknown_n > 0:
        log(f"[WARN] {unknown_n} rows assigned Group2 == 'Unknown' "
            f"(Genome_Type neither DNA nor RNA, check these rows).")

    # ---- 6. Drop rows with NaN in the core columns needed for modeling ----
    core_cols = [
        "A_W", "T_W", "G_W", "C_W", "Z",
        "annot_P_A", "annot_P_T", "annot_P_G", "annot_P_C",
        "annot_C_A", "annot_C_T", "annot_C_G", "annot_C_C",
        "P_portion", "C_portion",
    ]
    n_before = len(df)
    df_clean = df.dropna(subset=core_cols).copy()
    n_after = len(df_clean)
    log("-" * 70)
    log(f"[INFO] Dropped {n_before - n_after} rows with NaN in core modeling columns.")
    log(f"[INFO] Final usable rows: {n_after}")

    # ---- 7. Sanity check: A+T+G+C ~ 100 for A_W..C_W ----
    comp_sum = df_clean["A_W"] + df_clean["T_W"] + df_clean["G_W"] + df_clean["C_W"]
    bad_sum = ((comp_sum < 95) | (comp_sum > 105)).sum()
    if bad_sum > 0:
        log(f"[WARN] {bad_sum} rows have A_W+T_W+G_W+C_W outside [95,105] "
            f"(expected ~100 if these are percentages). Check units.")
    else:
        log("[OK] A_W+T_W+G_W+C_W sums are consistent with percentage scale.")

    # ---- 8. Z range check ----
    z_bad = ((df_clean["Z"] < 0) | (df_clean["Z"] > 100)).sum()
    if z_bad > 0:
        log(f"[WARN] {z_bad} rows have Z outside [0,100]. Check Z is a percentage.")
    else:
        log("[OK] Z values within [0,100].")
    z_edge = ((df_clean["Z"] <= 0) | (df_clean["Z"] >= 100)).sum()
    if z_edge:
        log(f"[INFO] {z_edge} rows have Z == 0 or Z == 100 (single-stranded annotation). They are KEPT: "
            f"M0 and M4 and the ML models use the one annotated strand; M1/M2/M3/M5 may be undefined "
            f"on them and are scored only where defined (see the N column).")

    # ---- 9. Final group filter ----
    final_counts = df_clean["Group2"].value_counts()
    valid_groups = [g for g in GROUP_ORDER if final_counts.get(g, 0) >= MIN_GROUP_SIZE]
    dropped_groups = [g for g in GROUP_ORDER if final_counts.get(g, 0) < MIN_GROUP_SIZE]
    if dropped_groups:
        log(f"[INFO] Dropping groups below MIN_GROUP_SIZE={MIN_GROUP_SIZE}: {dropped_groups}")
    df_clean = df_clean[df_clean["Group2"].isin(valid_groups)].copy()

    log("=" * 70)
    log(f"VALIDATION COMPLETE. Usable dataset: {len(df_clean)} rows, "
        f"{df_clean['Group2'].nunique()} groups, {df_clean[ORGAN_COL].nunique()} organisms.")
    log("=" * 70)

    paths = importlib.import_module("paths")
    with open(paths.table_path("00_validation_report.txt"), "w") as f:
        f.write("\n".join(report_lines))

    return df_clean, valid_groups


if __name__ == "__main__":
    paths = importlib.import_module("paths")
    df_clean, valid_groups = load_and_validate()
    df_clean.to_csv(paths.table_path("df_clean.csv"), index=False)
    print(f"\nSaved cleaned dataframe to tables/df_clean.csv")
    print(f"Valid groups ({len(valid_groups)}): {valid_groups}")
