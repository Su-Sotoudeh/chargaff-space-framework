"""
transfer_CSV_files.py -- copy the per-chromosome result tables

    Fasta_GFF_inputs/Mammals/<organism>/<chromosome>.csv
 -> Dataframes/Mammals_CSV/<organism>/<chromosome>.csv        (input of Figure 3)

Existing files are overwritten with the new copy; nothing in the source is deleted.

Usage:
    python3 transfer_CSV_files.py
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import os
import shutil

input_dir = str(config.require(config.MAMMALS_WORK_DIR, "Run RepeatMasker and the annotation step first (RepeatMasker/README.md)."))
output_dir = str(config.MAMMALS_CHR_DIR)
os.makedirs(output_dir, exist_ok=True)

n = 0
for organism in sorted(os.listdir(input_dir)):
    org_path = os.path.join(input_dir, organism)
    if not os.path.isdir(org_path):
        continue
    csvs = [f for f in os.listdir(org_path) if f.endswith(".csv")]
    if not csvs:
        continue
    dest = os.path.join(output_dir, organism)
    os.makedirs(dest, exist_ok=True)
    for f in csvs:
        shutil.copy2(os.path.join(org_path, f), os.path.join(dest, f))
        n += 1
print(f"Copied {n} CSV file(s) to {output_dir}")
