"""
mix_repeatmasker_result.py -- progress check.  For every organism folder in
Fasta_GFF_inputs/Mammals/ list those where the number of <chr>.fna files differs from the
number of <chr>.fna.out files (i.e. RepeatMasker has not finished) and write the table to
RepeatMasker/outputs/<name>.csv.  Read-only for the genome folders.

Usage:
    python3 mix_repeatmasker_result.py                 # -> RepeatMasker/outputs/repeatmasker_progress.csv
    python3 mix_repeatmasker_result.py --dir /other/folder --name my_check
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import argparse
import os

import pandas as pd


def check_files(path, file_name):
    rows = []
    for organism in sorted(os.listdir(path)):
        folder = os.path.join(path, organism)
        if not os.path.isdir(folder):
            continue
        files = os.listdir(folder)
        n = {ext: sum(f.endswith(ext) for f in files) for ext in (".fna", ".out", ".gff.gz", ".csv", ".tbl")}
        if n[".fna"] != n[".out"]:
            rows.append({"Name": organism, ".FNA": n[".fna"], ".OUT": n[".out"], ".GFF": n[".gff.gz"], ".CSV": n[".csv"], ".TBL": n[".tbl"]})
    df = pd.DataFrame(rows, columns=["Name", ".FNA", ".OUT", ".GFF", ".CSV", ".TBL"])
    print(df.to_string(index=False) if len(df) else "All organisms have one .out per .fna.")
    out = config.output_dir("RepeatMasker") / f"{file_name}.csv"
    df.to_csv(out, index=False)
    print("Written:", out)
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=str(config.MAMMALS_WORK_DIR), help="folder with one sub-folder per organism")
    ap.add_argument("--name", default="repeatmasker_progress")
    a = ap.parse_args()
    check_files(str(config.require(a.dir, "Nothing to check yet: see RepeatMasker/README.md.")), a.name)
