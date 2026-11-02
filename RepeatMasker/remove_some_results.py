"""
remove_some_results.py -- delete leftover RepeatMasker files (*.cat.gz, *.masked) and *.csv
from every organism folder in Fasta_GFF_inputs/Mammals/.

!! The *.csv files are the per-chromosome results of eukaryotes_for_chargaff.  Run
!! transfer_CSV_files.py FIRST so they are copied to Dataframes/Mammals_CSV/.

Safe by default: it only LISTS what it would delete.  Add --delete to really delete.

Usage:
    python3 remove_some_results.py            # dry run
    python3 remove_some_results.py --delete   # delete
    python3 remove_some_results.py --delete --keep-csv   # delete only *.cat.gz and *.masked
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import argparse
import os

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--delete", action="store_true", help="really delete (default: dry run)")
ap.add_argument("--keep-csv", action="store_true", help="never delete *.csv files")
a = ap.parse_args()

root = str(config.require(config.MAMMALS_WORK_DIR, "Nothing to clean yet."))
exts = (".cat.gz", ".masked") if a.keep_csv else (".cat.gz", ".masked", ".csv")
count = 0
for organism in sorted(os.listdir(root)):
    folder = os.path.join(root, organism)
    if not os.path.isdir(folder):
        continue
    for f in sorted(os.listdir(folder)):
        if f.endswith(exts):
            count += 1
            path = os.path.join(folder, f)
            if a.delete:
                os.remove(path)
                print("deleted", path)
            else:
                print("would delete", path)
print(f"\n{count} file(s) {'deleted' if a.delete else 'would be deleted (dry run; use --delete)'}.")
