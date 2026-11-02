"""
linux_in_python_with_espicial_lib.py -- like Linux_in_python_RepeatMasker.py, but with a custom
repeat library (-lib) instead of -species.  Use it for clades that Dfam covers badly
(see filter_repeatmasker_library.py to inspect a library).

Input / output folder (same folder):  Fasta_GFF_inputs/Mammals/<organism>/

Usage:
    python3 linux_in_python_with_espicial_lib.py Neogale_vison --lib /path/to/RepeatMasker2.lib
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import argparse
import os
import subprocess

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("organism", help="folder name inside Fasta_GFF_inputs/Mammals/")
ap.add_argument("--lib", required=True, help="repeat library FASTA given to RepeatMasker -lib")
a = ap.parse_args()

input_dir = os.path.abspath(config.require(config.MAMMALS_WORK_DIR / a.organism, "Create it and put the per-chromosome .fna files in it."))
lib = os.path.abspath(config.require(a.lib, "Give the full path of your repeat library with --lib."))
fasta_files = sorted(f for f in os.listdir(input_dir) if f.endswith(".fna"))

for fasta_file in fasta_files:
    input_file = os.path.join(input_dir, fasta_file)
    if os.path.exists(os.path.join(input_dir, f"{fasta_file}.out")):
        print(f"Skipping {fasta_file} - .out file already exists.")
        continue
    print(f"Running RepeatMasker on {fasta_file}...")
    command = [
        config.REPEATMASKER_BIN,
        "-no_is", "-nolow", "-norna",
        "-pa", str(config.REPEATMASKER_THREADS), "-engine", "rmblast", "-lib", lib,
        input_file,
    ]
    subprocess.run(command, check=True)

for fasta_file in fasta_files:
    for ext in (".masked", ".cat.gz"):
        leftover = os.path.join(input_dir, f"{fasta_file}{ext}")
        if os.path.exists(leftover):
            os.remove(leftover)
            print(f"Removed: {leftover}")
