import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)
# working as chromosome by chromosome
# resumable version: save each chromosome immediately
# =========================================================
# Resumable gene/CDS composition extraction by chromosome
#
# Behavior:
# - BATCH_SIZE = number of unfinished organisms to process this run
# - Reads existing CSV and skips completed (Organ_ID, Chr_ID)
# - Processes organisms in parallel
# - Appends new rows to one master CSV
# =========================================================

import os
import csv
import gzip
import re
from pathlib import Path
from multiprocessing import Pool, cpu_count

from Bio import SeqIO
from Bio.Seq import Seq


# =========================================================
# CONFIG
# =========================================================
NAME = "Archaea"

INPUT_ROOT = str(config.require(config.ncbi_data_dir(NAME), "Unzip the NCBI genome package of this group into Fasta_GFF_inputs/<Group>/ (see Fasta_GFF_inputs/README.md)."))   # Fasta_GFF_inputs/<NAME>/ncbi_dataset/data

OUTPUT_CSV = str(config.output_dir(8, "Single_annote_composition") / f"{NAME}_gene_by_chromosome.csv")

# IMPORTANT:
# This means "process next 10 unfinished organisms this run"
# It does NOT mean "always first 10 organisms"
BATCH_SIZE = 1000

# multiprocessing workers
MAX_WORKERS = 2

# chromosome filter(length size)
MIN_CHROMOSOME_BP = 10000000

# Choose feature type:
# Use "gene" for gene extraction
# Use "CDS" for CDS extraction
FEATURE_TYPE = "gene"

CSV_HEADER = [
    "Organ_ID",
    "Chr_ID",
    "A_raw",
    "T_raw",
    "C_raw",
    "G_raw",
    "length",
    "strand",
]


# =========================================================
# FILE HELPERS
# =========================================================
def open_maybe_gzip(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def normalize_text(text):
    return re.sub(r"\s+", " ", str(text).lower()).strip()


def find_genome_files(genome_dir):
    fasta_files = []
    gff_files = []

    for root, _, files in os.walk(genome_dir):
        for f in files:
            fl = f.lower()
            full = os.path.join(root, f)

            # Ignore extracted CDS FASTA files, because we want genome FASTA
            if fl.startswith("cds_"):
                continue

            if fl.endswith((".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")):
                fasta_files.append(full)

            elif fl.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
                gff_files.append(full)

    return fasta_files, gff_files


# =========================================================
# CHROMOSOME FILTER
# =========================================================
def is_valid_chromosome_name(seq_id, description, seq_len):
    text = normalize_text(f"{seq_id} {description}")

    # length filter
    if seq_len < MIN_CHROMOSOME_BP:
        return False

    # remove unwanted sequence types
    bad_keywords = [
        "scaffold",
        "contig",
        "unlocalized",
        "unplaced",
        "patch",
        "plasmid",
        "mitochond",
        "chloroplast",
    ]

    if any(word in text for word in bad_keywords):
        return False
    # Accept names like chr1, chr2, chrX, chrZ, etc.
    if re.search(r"\bchr[\s._-]*[0-9xyzwm]+\b", text) :
        return True
    if re.search(r">NC_", text) :
        return True
    if re.search(r">CM_", text) :
        return True
    if re.search(r">CP_", text) :
        return True

    if re.search(r"whole genome shotgun sequence", text) :
        return True

    if "chromosome un" in text:
        return False

    if re.search(r"\bchr\s*un\b", text):
        return False

    # Accept if description contains chromosome
    if "chromosome" in text:
        return True

    if 'complete sequence' in text:
        return True
    if 'complete genome' in text:
        return True
    if 'partial sequence' in text:
        return False
    print(text)
    return False


# =========================================================
# FASTA / GFF PARSING
# =========================================================
def load_chromosomes(fasta_files):
    """
    Load only valid chromosome sequences into memory for one organism.

    Returns:
        chrom_seqs = {
            chr_id: uppercase_sequence_string
        }
    """
    chrom_seqs = {}

    for fasta_path in fasta_files:
        with open_maybe_gzip(fasta_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_len = len(record.seq)

                if not is_valid_chromosome_name(record.id, record.description, seq_len):
                    continue

                chrom_seqs[str(record.id)] = str(record.seq).upper()

    return chrom_seqs


def parse_gff_features(gff_files, valid_chr_ids):
    """
    Collect selected GFF features grouped by chromosome.

    Returns:
        features_by_chr = {
            chr_id: [(start0, end0, strand), ...]
        }
    """
    features_by_chr = {chr_id: [] for chr_id in valid_chr_ids}

    target_feature = FEATURE_TYPE.lower()

    for gff_path in gff_files:
        with open_maybe_gzip(gff_path, "rt") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue

                parts = line.rstrip("\n").split("\t")

                if len(parts) != 9:
                    continue

                seqid, source, feature_type, start, end, score, strand, phase, attributes = parts

                if feature_type.lower() != target_feature:
                    continue

                if seqid not in features_by_chr:
                    continue

                try:
                    start0 = int(start) - 1
                    end0 = int(end)
                except ValueError:
                    continue

                if start0 < 0:
                    continue

                if end0 <= start0:
                    continue

                if strand not in {"+", "-"}:
                    strand = "."

                features_by_chr[seqid].append((start0, end0, strand))

    return features_by_chr


# =========================================================
# SEQUENCE EXTRACTION
# =========================================================
def count_bases(seq):
    seq = seq.upper()

    a = seq.count("A")
    t = seq.count("T")
    c = seq.count("C")
    g = seq.count("G")

    return a, t, c, g, len(seq)


def extract_feature_row(organ_id, chr_id, chrom_seq, start0, end0, strand):
    """
    Extract one feature sequence and count bases.
    """
    if end0 > len(chrom_seq):
        return None

    feature_seq = chrom_seq[start0:end0]

    if not feature_seq:
        return None

    a, t, c, g, length = count_bases(feature_seq)

    return [
        organ_id,
        chr_id,
        a,
        t,
        c,
        g,
        length,
        strand,
    ]


# =========================================================
# CSV OUTPUT / RESUME
# =========================================================
def ensure_output_csv(csv_path):
    """
    Create output CSV only if it does not already exist.
    Existing file is never overwritten.
    """
    if not os.path.exists(csv_path):
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(CSV_HEADER)

        print(f"[INIT] Created {csv_path}", flush=True)

    else:
        print(f"[INIT] Using existing file: {csv_path}", flush=True)


def append_rows(csv_path, rows):
    """
    Append rows safely to CSV.
    """
    if not rows:
        return

    with open(csv_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerows(rows)

        f.flush()
        os.fsync(f.fileno())


def load_completed_pairs(csv_path):
    """
    Read existing output CSV and return:

        completed_pairs = set of (Organ_ID, Chr_ID)
        completed_orgs = set of Organ_ID

    Important:
    If at least one row exists for a chromosome, this code treats
    that chromosome as completed.
    """
    completed_pairs = set()
    completed_orgs = set()

    if not os.path.exists(csv_path):
        return completed_pairs, completed_orgs

    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        if reader.fieldnames is None:
            return completed_pairs, completed_orgs

        if "Organ_ID" not in reader.fieldnames or "Chr_ID" not in reader.fieldnames:
            raise ValueError(
                f"Existing CSV must contain Organ_ID and Chr_ID columns. "
                f"Found columns: {reader.fieldnames}"
            )

        for row in reader:
            organ_id = str(row.get("Organ_ID", "")).strip()
            chr_id = str(row.get("Chr_ID", "")).strip()

            if organ_id and chr_id:
                completed_pairs.add((organ_id, chr_id))
                completed_orgs.add(organ_id)

    return completed_pairs, completed_orgs


# =========================================================
# STATUS SCAN
# =========================================================
def scan_organism_status(genome_dir, completed_pairs):
    """
    Scan one organism to know:
    - valid chromosomes
    - chromosomes already done
    - chromosomes remaining

    Returns dictionary:
        {
            organ_id,
            valid_chr,
            done_chr,
            remaining_chr,
            usable,
            genome_dir,
            error
        }
    """
    organ_id = os.path.basename(genome_dir.rstrip("/"))

    status = {
        "organ_id": organ_id,
        "genome_dir": genome_dir,
        "valid_chr": 0,
        "done_chr": 0,
        "remaining_chr": 0,
        "usable": False,
        "error": None,
    }

    fasta_files, gff_files = find_genome_files(genome_dir)

    if not fasta_files:
        status["error"] = "missing FASTA"
        return status

    if not gff_files:
        status["error"] = "missing GFF"
        return status

    try:
        chrom_seqs = load_chromosomes(fasta_files)

        if not chrom_seqs:
            status["error"] = "no valid chromosomes"
            return status

        valid_chr_ids = set(chrom_seqs.keys())

        valid_chr = len(valid_chr_ids)
        done_chr = sum((organ_id, chr_id) in completed_pairs for chr_id in valid_chr_ids)
        remaining_chr = valid_chr - done_chr

        status["valid_chr"] = valid_chr
        status["done_chr"] = done_chr
        status["remaining_chr"] = remaining_chr
        status["usable"] = True

        return status

    except Exception as e:
        status["error"] = str(e)
        return status


# =========================================================
# WORKER
# =========================================================
def process_one_organism(args):
    """
    Process one organism.

    It skips chromosomes that already exist in completed_pairs.

    Returns:
        organism_rows
    """
    genome_dir, completed_pairs = args

    organ_id = os.path.basename(genome_dir.rstrip("/"))

    fasta_files, gff_files = find_genome_files(genome_dir)

    if not fasta_files or not gff_files:
        print(f"[SKIP] {organ_id}: missing FASTA or GFF", flush=True)
        return []

    try:
        chrom_seqs = load_chromosomes(fasta_files)

        if not chrom_seqs:
            print(f"[SKIP] {organ_id}: no valid chromosomes found", flush=True)
            return []

        features_by_chr = parse_gff_features(gff_files, set(chrom_seqs.keys()))

        organism_rows = []

        analysed_chr = 0
        skipped_chr = 0
        no_feature_chr = 0

        for chr_id, chrom_seq in chrom_seqs.items():

            # Resume check by chromosome
            if (organ_id, chr_id) in completed_pairs:
                skipped_chr += 1
                print(
                    f"[SKIP-CHR] {organ_id} | {chr_id} already in CSV",
                    flush=True
                )
                continue

            feature_list = features_by_chr.get(chr_id, [])

            if not feature_list:
                no_feature_chr += 1
                print(
                    f"[NO-FEATURE] {organ_id} | {chr_id} has no {FEATURE_TYPE} features",
                    flush=True
                )
                continue

            chr_rows = []

            for start0, end0, strand in feature_list:
                row = extract_feature_row(
                    organ_id=organ_id,
                    chr_id=chr_id,
                    chrom_seq=chrom_seq,
                    start0=start0,
                    end0=end0,
                    strand=strand,
                )

                if row is not None:
                    chr_rows.append(row)

            if chr_rows:
                organism_rows.extend(chr_rows)
                analysed_chr += 1

                print(
                    f"[DONE-CHR] {organ_id} | {chr_id} | "
                    f"{FEATURE_TYPE} rows: {len(chr_rows)}",
                    flush=True
                )

            else:
                print(
                    f"[EMPTY-CHR] {organ_id} | {chr_id} produced no valid rows",
                    flush=True
                )

        print(
            f"[DONE-ORG] {organ_id} | "
            f"new chromosomes analysed: {analysed_chr} | "
            f"already done skipped: {skipped_chr} | "
            f"no-feature chromosomes: {no_feature_chr} | "
            f"new rows: {len(organism_rows)}",
            flush=True
        )

        return organism_rows

    except Exception as e:
        print(f"[ERROR] {organ_id}: {e}", flush=True)
        return []


# =========================================================
# MAIN
# =========================================================
def main():
    ensure_output_csv(OUTPUT_CSV)

    completed_pairs, completed_orgs = load_completed_pairs(OUTPUT_CSV)

    all_dirs = [str(p) for p in Path(INPUT_ROOT).iterdir() if p.is_dir()]
    all_dirs.sort()

    total_available_organisms = len(all_dirs)

    print("======================================", flush=True)
    print("[START]", flush=True)
    print(f"[INFO] Total available organisms in folder: {total_available_organisms}", flush=True)
    print(f"[INFO] Workers: {MAX_WORKERS}", flush=True)
    print(f"[INFO] Output: {OUTPUT_CSV}", flush=True)
    print(f"[INFO] Feature type: {FEATURE_TYPE}", flush=True)
    print(f"[INFO] Minimum chromosome size: {MIN_CHROMOSOME_BP:,}", flush=True)
    print(f"[INFO] Batch size requested this run: {BATCH_SIZE}", flush=True)
    print(f"[INFO] Already analysed chromosome pairs in CSV: {len(completed_pairs)}", flush=True)
    print(f"[INFO] Already analysed organisms in CSV: {len(completed_orgs)}", flush=True)
    print("======================================", flush=True)

    # -----------------------------------------------------
    # FAST SELECTION:
    # Do NOT scan all organisms.
    # Select next BATCH_SIZE organisms whose Organ_ID is not in CSV.
    # -----------------------------------------------------
    selected_dirs = []

    for genome_dir in all_dirs:
        organ_id = os.path.basename(genome_dir.rstrip("/"))

        # If this organism already has any saved chromosome in CSV,
        # consider it already processed and skip it.
        if organ_id in completed_orgs:
            continue

        selected_dirs.append(genome_dir)

        if len(selected_dirs) >= BATCH_SIZE:
            break

    print("--------------------------------------------------", flush=True)
    print(f"[INFO] Fast batch mode: not scanning all organisms", flush=True)
    print(f"[INFO] Organisms selected for this run: {len(selected_dirs)}", flush=True)

    if not selected_dirs:
        print("[INFO] No new organisms found to analyse.", flush=True)
        return

    print("[INFO] Selected organisms for this run:", flush=True)
    for d in selected_dirs:
        print(f"  - {os.path.basename(d.rstrip('/'))}", flush=True)

    print("--------------------------------------------------", flush=True)

    worker_args = [(genome_dir, completed_pairs) for genome_dir in selected_dirs]

    total_rows = 0
    finished_orgs = 0

    with Pool(processes=MAX_WORKERS) as pool:
        for i, rows in enumerate(pool.imap_unordered(process_one_organism, worker_args), start=1):

            append_rows(OUTPUT_CSV, rows)

            total_rows += len(rows)
            finished_orgs += 1

            print(
                f"[PROGRESS] {i}/{len(worker_args)} selected organisms finished | "
                f"new {FEATURE_TYPE} rows saved this run: {total_rows}",
                flush=True
            )

    print("======================================", flush=True)
    print("[FINISHED]", flush=True)
    print(f"[INFO] Batch size requested: {BATCH_SIZE}", flush=True)
    print(f"[INFO] Selected organisms this run: {len(selected_dirs)}", flush=True)
    print(f"[INFO] Finished organism jobs this run: {finished_orgs}", flush=True)
    print(f"[INFO] New {FEATURE_TYPE} rows saved this run: {total_rows}", flush=True)
    print(f"[INFO] Output file: {OUTPUT_CSV}", flush=True)
    print("======================================", flush=True)


if __name__ == "__main__":
    main()
