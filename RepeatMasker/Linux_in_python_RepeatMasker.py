"""
Linux_in_python_RepeatMasker.py -- run RepeatMasker on every <chromosome>.fna of one organism.

Input / output folder (same folder):  Fasta_GFF_inputs/Mammals/<organism>/
    <chr>.fna  ->  <chr>.fna.out  (+ .tbl, .masked, .cat.gz; the last two are removed afterwards)

Chromosomes that already have a .fna.out are skipped, so the script can be re-run to resume.
RepeatMasker program and threads come from config.py (env CSPR_REPEATMASKER, CSPR_RM_THREADS).
Install RepeatMasker first: see RepeatMasker/README.md.

Usage:
    python3 Linux_in_python_RepeatMasker.py Human_GRCh38.p14 9606
                                            ^ folder in Fasta_GFF_inputs/Mammals   ^ NCBI taxon ID (or species name)
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import argparse
import os
import subprocess


def RepetMasker(genome_folder_path, Taxon_ID):
    input_dir = os.path.abspath(genome_folder_path)
    fasta_files = sorted(f for f in os.listdir(input_dir) if f.endswith(".fna"))

    for fasta_file in fasta_files:
        input_file = os.path.join(input_dir, fasta_file)
        if os.path.exists(os.path.join(input_dir, f"{fasta_file}.out")):
            print(f"Skipping {fasta_file} - .out file already exists.")
            continue
        print(f"Running RepeatMasker on {fasta_file}...")
        command = [
            config.REPEATMASKER_BIN,
            "-species", str(Taxon_ID),
            "-no_is", "-nolow", "-norna",
            "-pa", str(config.REPEATMASKER_THREADS), "-engine", "rmblast",
            input_file,
        ]
        subprocess.run(command, check=True)

    # remove big files that are not needed any more
    for fasta_file in fasta_files:
        for ext in (".masked", ".cat.gz"):
            leftover = os.path.join(input_dir, f"{fasta_file}{ext}")
            if os.path.exists(leftover):
                os.remove(leftover)
                print(f"Removed: {leftover}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("organism", help="folder name inside Fasta_GFF_inputs/Mammals/ (e.g. Human_GRCh38.p14)")
    ap.add_argument("taxon", help="NCBI taxonomy ID or species name passed to RepeatMasker -species (e.g. 9606)")
    a = ap.parse_args()
    folder = config.require(config.MAMMALS_WORK_DIR / a.organism, "Create it and put the per-chromosome .fna files in it (Genome_file_dividing.py).")
    RepetMasker(folder, a.taxon)
