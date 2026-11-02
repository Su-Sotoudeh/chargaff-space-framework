"""
Genome_file_dividing.py -- split a whole-genome FASTA into one file per chromosome.

Input / output folder (same folder):  Fasta_GFF_inputs/Mammals/<organism>/
    <genome>.fna  or  <genome>.fna.gz   ->   <chromosome_id>.fna   (one per chromosome)

Scaffolds, unplaced contigs (NW_...) and patches are skipped.  The original whole-genome
file is kept unless --delete-original is given.

Usage:
    python3 Genome_file_dividing.py Human_GRCh38.p14                 # organism folder name
    python3 Genome_file_dividing.py Human_GRCh38.p14 --delete-original
    python3 Genome_file_dividing.py --file /path/genome.fna.gz --out /path/folder
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import argparse
import gzip
import os

from Bio import SeqIO


def _open_text(path):
    """Open plain or gzip-compressed FASTA as text (detected from the file's first bytes)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def Whole_genome_to_chr(file_path, output_path, delete_original=False):
    """Write every chromosome record of file_path to output_path/<record.id>.fna."""
    os.makedirs(output_path, exist_ok=True)
    try:
        handle = _open_text(file_path)
    except EOFError:
        print("EOFError: could not read", file_path)
        return None
    n_saved = 0
    with handle:
        try:
            for record in SeqIO.parse(handle, "fasta"):
                desc = record.description
                if "scaffold" in desc or "NW_" in desc or "patch" in desc:
                    continue
                if "chromosome" in desc.lower():
                    out_file = os.path.join(output_path, f"{record.id}.fna")
                    with open(out_file, "w") as out_f:
                        SeqIO.write(record, out_f, "fasta")
                    n_saved += 1
                    print(f"Saved: {out_file}")
                else:
                    print(f"Skipped (not a chromosome): {record.id}")
        except EOFError:
            print("EOFError: the compressed file is truncated; chromosomes saved so far are kept:", file_path)
            return None
    if delete_original and n_saved > 0:
        os.remove(file_path)
        print("Removed original:", file_path)
    return n_saved


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("organism", nargs="?", help="folder name inside Fasta_GFF_inputs/Mammals/")
    ap.add_argument("--file", help="whole-genome FASTA (.fna or .fna.gz); default: the only genome file in the organism folder")
    ap.add_argument("--out", help="output folder; default: the organism folder")
    ap.add_argument("--delete-original", action="store_true", help="delete the whole-genome file after splitting")
    a = ap.parse_args()
    if a.organism:
        folder = config.require(config.MAMMALS_WORK_DIR / a.organism, "Put the whole-genome .fna(.gz) into Fasta_GFF_inputs/Mammals/<organism>/")
    elif a.file:
        folder = _Path(a.file).resolve().parent
    else:
        ap.error("give an organism folder name or --file")
    if a.file:
        genome = _Path(a.file)
    else:
        found = sorted(f for f in folder.iterdir() if f.name.endswith((".fna", ".fna.gz")) and f.name.startswith("G"))
        if len(found) != 1:
            raise SystemExit(f"Expected exactly one whole-genome file starting with 'G' (GCF_/GCA_) in {folder}, found {len(found)}. Use --file.")
        genome = found[0]
    Whole_genome_to_chr(str(genome), a.out or str(folder), a.delete_original)
