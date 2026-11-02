"""
Automation.py -- the whole mammal chain for every organism in the reference list.

For each organism of the list:
  1. copy its genome (.fna.gz) and annotation (.gff.gz) from the downloaded NCBI package
         Fasta_GFF_inputs/Chordata/ncbi_dataset/data/<assembly accession>/
     to  Fasta_GFF_inputs/Mammals/<Organism_Name>/
  2. split the genome into one <chromosome>.fna per chromosome   (Genome_file_dividing.py)
  3. run RepeatMasker on every chromosome                        (Linux_in_python_RepeatMasker.py)
  4. annotate exons / introns / repeats per chromosome -> <chromosome>.csv
                                                                  (eukaryotes_for_chargaff.py)
Afterwards run, from this folder:
     python3 transfer_CSV_files.py       # CSV tables -> Dataframes/Mammals_CSV  (Figure 3)
     python3 ../Masking_GFF_OUT/produce_masked_fna.py   # -> Fasta_GFF_inputs/Masked_Fasta (Figures 4, 5)

Everything can be stopped and re-started: finished steps are skipped.

Usage:
    python3 Automation.py                                   # uses the list next to this script
    python3 Automation.py --list my_list.csv --only Mus_musculus
    python3 Automation.py --skip-annotation                 # stop after RepeatMasker

The reference list (.csv, or .tsv/.txt = tab separated) needs the columns
"Organism Name", "Organism Taxonomic ID" and the assembly accession (column
"Assembly Accession", or the first column when that name is absent).
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import argparse
import glob
import os
import shutil

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
_sys.path.insert(0, HERE)           # so the helper modules in this folder can be imported
from Genome_file_dividing import Whole_genome_to_chr            # noqa: E402
from Linux_in_python_RepeatMasker import RepetMasker            # noqa: E402

DEFAULT_LIST = os.path.join(HERE, "mamals_refrence_annotated_2025_7_12.csv")


def read_reference_list(path):
    sep = "\t" if str(path).lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    if "Assembly Accession" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Assembly Accession"})
    missing = {"Organism Name", "Organism Taxonomic ID"} - set(df.columns)
    if missing:
        raise SystemExit(f"Reference list {path} lacks column(s): {sorted(missing)}")
    return df


def copy_genome_and_gff(accession, source_root, dest):
    """Copy <accession>/G*.fna.gz and *.gff.gz into dest.  Returns the copied genome path or None."""
    pkg = os.path.join(source_root, accession)
    fna = glob.glob(os.path.join(pkg, "G*.fna.gz"))
    gff = glob.glob(os.path.join(pkg, "*.gff.gz"))
    if not fna or not gff:
        print(f"  ! genome or GFF missing in {pkg}")
        return None
    os.makedirs(dest, exist_ok=True)
    for f in (fna[0], gff[0]):
        target = os.path.join(dest, os.path.basename(f))
        if not os.path.exists(target):
            shutil.copy2(f, target)
    return os.path.join(dest, os.path.basename(fna[0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", default=DEFAULT_LIST, help="reference list of organisms")
    ap.add_argument("--only", help="process only this organism (name with underscores, e.g. Mus_musculus)")
    ap.add_argument("--skip-annotation", action="store_true", help="stop after RepeatMasker")
    a = ap.parse_args()

    source_root = str(config.require(
        config.ncbi_data_dir("Chordata"),
        "Download the mammal genomes first (Download_Data/README.md) into Fasta_GFF_inputs/Chordata/."))
    available = set(os.listdir(source_root))
    table = read_reference_list(a.list)
    bad_log = config.output_dir("RepeatMasker") / "bad_files.txt"

    for _, row in table.iterrows():
        organism = str(row["Organism Name"]).replace(" ", "_")
        if a.only and organism != a.only:
            continue
        accession, taxon = row["Assembly Accession"], row["Organism Taxonomic ID"]
        dest = str(config.MAMMALS_WORK_DIR / organism)
        print(f"\n=== {organism} ({accession})")

        if accession not in available:
            print("  we do not have this genome - please download it (Download_Data/README.md)")
            continue

        have_chr = os.path.isdir(dest) and any(f.endswith(".fna") and not f.startswith("G") for f in os.listdir(dest))
        if not have_chr:
            genome = copy_genome_and_gff(accession, source_root, dest)
            if genome is None:
                continue
            print("  splitting the genome into chromosomes")
            n = Whole_genome_to_chr(genome, dest, delete_original=True)   # the copy is removed, the download stays
            if not n:
                print("  ! splitting failed (corrupt or truncated .gz?)")
                with open(bad_log, "a") as z:
                    z.write(organism + "\n")
                continue

        print("  RepeatMasker")
        try:
            RepetMasker(dest, str(taxon))
        except Exception as e:                      # keep going with the next organism
            print("  ! RepeatMasker failed:", e)
            continue

        if a.skip_annotation:
            continue
        print("  annotation")
        try:
            from eukaryotes_for_chargaff import all_element_annotation_in_plus_and_minis
        except ImportError:
            raise SystemExit(
                "\nRepeatMasker/eukaryotes_for_chargaff.py is not in the repository.\n"
                "Add that file (it must define all_element_annotation_in_plus_and_minis(folder, organism_name)),\n"
                "or run with --skip-annotation.\n")
        all_element_annotation_in_plus_and_minis(dest, organism)


if __name__ == "__main__":
    main()
