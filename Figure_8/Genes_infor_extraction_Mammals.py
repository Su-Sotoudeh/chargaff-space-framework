#!/usr/bin/env python3
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

# =========================================================
# Resumable Masked Element Composition Extraction by Chromosome
#
# Behavior:
# - Scans subfolders in INPUT_ROOT (each subfolder represents an Organism)
# - Processes next BATCH_SIZE unfinished organisms this run
# - Extracts contiguous target elements using DECODE_MAP (no GFF required)
# - Appends base count rows to a single master CSV:
#   Organ_ID, Chr_ID, A_raw, T_raw, C_raw, G_raw, length, strand
# =========================================================

import os
import csv
import gzip
import re
from pathlib import Path
from multiprocessing import Pool
import collections

# =========================================================
# CONFIG
# =========================================================
NAME = "Mammals"  # Group name matching directory structure

INPUT_ROOT = str(config.require(config.MASKED_FASTAS, "Put the masked FASTA files into Fasta_GFF_inputs/Masked_Fastas/<organism>/ (see Fasta_GFF_inputs/README.md)."))   # Fasta_GFF_inputs/Masked_Fastas/<organism>/

OUTPUT_CSV = str(config.output_dir(8, "Single_annote_composition") / f"{NAME}_elements_by_chromosome.csv")

# Number of unfinished organisms to process this run
BATCH_SIZE = 100

# Multiprocessing workers
MAX_WORKERS = 8

# Chromosome filter (minimum size in base pairs)
MIN_CHROMOSOME_BP = 5000000

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
# DECODE MAP & TARGETS
# =========================================================
DECODE_MAP = {
    '1': ('A', 'DNA', '+'), '2': ('T', 'DNA', '+'), '3': ('C', 'DNA', '+'), '4': ('G', 'DNA', '+'),
    '5': ('A', 'DNA', '-'), '6': ('T', 'DNA', '-'), '7': ('C', 'DNA', '-'), '8': ('G', 'DNA', '-'),
    'L': ('A', 'LTR', '+'), 'M': ('T', 'LTR', '+'), '@': ('C', 'LTR', '+'), '#': ('G', 'LTR', '+'),
    'l': ('A', 'LTR', '-'), 'm': ('T', 'LTR', '-'), '$': ('C', 'LTR', '-'), '%': ('G', 'LTR', '-'),
    'O': ('A', 'SINE', '+'), 'P': ('T', 'SINE', '+'), 'J': ('C', 'SINE', '+'), 'K': ('G', 'SINE', '+'),
    'o': ('A', 'SINE', '-'), 'p': ('T', 'SINE', '-'), 'j': ('C', 'SINE', '-'), 'k': ('G', 'SINE', '-'),
    '*': ('A', 'LINE', '+'), 'Y': ('T', 'LINE', '+'), 'U': ('C', 'LINE', '+'), 'I': ('G', 'LINE', '+'),
    '+': ('A', 'LINE', '-'), 'y': ('T', 'LINE', '-'), 'u': ('C', 'LINE', '-'), 'i': ('G', 'LINE', '-'),
    'Q': ('A', 'Int', '+'), 'W': ('T', 'Int', '+'), 'E': ('C', 'Int', '+'), 'R': ('G', 'Int', '+'),
    'q': ('A', 'Int', '-'), 'w': ('T', 'Int', '-'), 'e': ('C', 'Int', '-'), 'r': ('G', 'Int', '-'),
    'S': ('A', 'Exe', '+'), 'D': ('T', 'Exe', '+'), 'F': ('C', 'Exe', '+'), 'H': ('G', 'Exe', '+'),
    's': ('A', 'Exe', '-'), 'd': ('T', 'Exe', '-'), 'f': ('C', 'Exe', '-'), 'h': ('G', 'Exe', '-'),
    'Z': ('A', 'CDS', '+'), 'X': ('T', 'CDS', '+'), 'V': ('C', 'CDS', '+'), 'B': ('G', 'CDS', '+'),
    'z': ('A', 'CDS', '-'), 'x': ('T', 'CDS', '-'), 'v': ('C', 'CDS', '-'), 'b': ('G', 'CDS', '-'),
    'A': ('A', 'shre', '.'), 'T': ('T', 'shre', '.'), 'C': ('C', 'shre', '.'), 'G': ('G', 'shre', '.'),
    'a': ('A', 'un', '.'),   't': ('T', 'un', '.'),   'c': ('C', 'un', '.'),   'g': ('G', 'un', '.'),
}

TARGET_FEATURES = {'LINE', 'SINE', 'LTR', 'DNA', 'CDS', 'UTR', 'Int'}

def decode_char(char):
    return DECODE_MAP.get(char, (char, 'un', '.'))

# =========================================================
# FILE HELPERS
# =========================================================
def open_maybe_gzip(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")

def normalize_text(text):
    return re.sub(r"\s+", " ", str(text).lower()).strip()

def find_fasta_files(genome_dir):
    fasta_files = []
    for root, _, files in os.walk(genome_dir):
        for f in files:
            fl = f.lower()
            full = os.path.join(root, f)
            if fl.endswith((".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")):
                fasta_files.append(full)
    fasta_files.sort()
    return fasta_files

# =========================================================
# CHROMOSOME FILTER
# =========================================================
def is_valid_chromosome_name(seq_id, description, seq_len):
    text = normalize_text(f"{seq_id} {description}")

    # Length filter
    if seq_len < MIN_CHROMOSOME_BP:
        return False

    # Remove unwanted sequence types
    bad_keywords = [
        "scaffold", "contig", "unlocalized", "unplaced", "patch",
        "plasmid", "mitochond", "chloroplast",
    ]

    if any(word in text for word in bad_keywords):
        return False

    # Accept standard chromosomal tags
    if re.search(r"\bchr[\s._-]*[0-9xyzwm]+\b", text):
        return True
    if re.search(r">nc_", text) or seq_id.upper().startswith("NC_"):
        return True
    if re.search(r">cm_", text) or seq_id.upper().startswith("CM_"):
        return True
    if re.search(r">cp_", text) or seq_id.upper().startswith("CP_"):
        return True
    if re.search(r"whole genome shotgun sequence", text):
        return True
    if "chromosome un" in text or re.search(r"\bchr\s*un\b", text):
        return False
    if "chromosome" in text:
        return True
    if 'complete sequence' in text or 'complete genome' in text:
        return True
    if 'partial sequence' in text:
        return False

    return False

# =========================================================
# ELEMENT PARSING LOGIC
# =========================================================
def extract_elements_from_sequence(organ_id, chr_id, sequence):
    if not sequence:
        return []

    elements = []
    start = 0
    _, last_feat, last_strand = decode_char(sequence[0])

    def add_element_row(start_idx, end_idx, feat, strand):
        if feat in TARGET_FEATURES and strand in ('+', '-'):
            block = sequence[start_idx:end_idx]
            counts = collections.Counter()
            for ch in block:
                base, _, _ = decode_char(ch)
                counts[base.upper()] += 1

            elements.append([
                organ_id,
                chr_id,
                counts['A'],
                counts['T'],
                counts['C'],
                counts['G'],
                len(block),
                strand
            ])

    for i in range(1, len(sequence)):
        _, feat, strand = decode_char(sequence[i])

        if feat != last_feat or strand != last_strand:
            add_element_row(start, i, last_feat, last_strand)
            start = i
            last_feat = feat
            last_strand = strand

    add_element_row(start, len(sequence), last_feat, last_strand)
    return elements

# =========================================================
# CSV OUTPUT / RESUME
# =========================================================
def ensure_output_csv(csv_path):
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    if not os.path.exists(csv_path):
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(CSV_HEADER)
        print(f"[INIT] Created {csv_path}", flush=True)
    else:
        print(f"[INIT] Using existing file: {csv_path}", flush=True)

def append_rows(csv_path, rows):
    if not rows:
        return
    with open(csv_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(rows)
        f.flush()
        os.fsync(f.fileno())

def load_completed_pairs(csv_path):
    completed_pairs = set()
    completed_orgs = set()

    if not os.path.exists(csv_path):
        return completed_pairs, completed_orgs

    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        # Added delimiter='\t' to match output formatting so resume works correctly
        reader = csv.DictReader(f, delimiter='\t')
        if reader.fieldnames is None:
            return completed_pairs, completed_orgs

        for row in reader:
            organ_id = str(row.get("Organ_ID", "")).strip()
            chr_id = str(row.get("Chr_ID", "")).strip()
            if organ_id and chr_id:
                completed_pairs.add((organ_id, chr_id))
                completed_orgs.add(organ_id)

    return completed_pairs, completed_orgs

# =========================================================
# WORKER
# =========================================================
def process_one_organism(args):
    genome_dir, completed_pairs = args
    organ_id = os.path.basename(genome_dir.rstrip("/"))

    fasta_files = find_fasta_files(genome_dir)
    if not fasta_files:
        return organ_id, [], 0, 0, 0, [f"[SKIP] {organ_id}: missing FASTA files"]

    try:
        organism_rows = []
        analysed_chr = 0
        skipped_chr = 0
        empty_chr = 0
        chr_logs = []

        for fasta_path in fasta_files:
            current_seq = []
            current_header = ""

            with open_maybe_gzip(fasta_path, "rt") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith(">"):
                        if current_seq:
                            seq = "".join(current_seq)
                            seq_id = current_header.split()[0]

                            if is_valid_chromosome_name(seq_id, current_header, len(seq)):
                                if (organ_id, seq_id) in completed_pairs:
                                    skipped_chr += 1
                                else:
                                    chr_elements = extract_elements_from_sequence(organ_id, seq_id, seq)
                                    if chr_elements:
                                        organism_rows.extend(chr_elements)
                                        analysed_chr += 1
                                        chr_logs.append(f"[DONE-CHR] {organ_id} | {seq_id} | gene rows: {len(chr_elements)}")
                                    else:
                                        empty_chr += 1
                            current_seq = []
                        current_header = line[1:]
                    else:
                        current_seq.append(line)

                # Process final sequence in file
                if current_seq:
                    seq = "".join(current_seq)
                    seq_id = current_header.split()[0]
                    if is_valid_chromosome_name(seq_id, current_header, len(seq)):
                        if (organ_id, seq_id) in completed_pairs:
                            skipped_chr += 1
                        else:
                            chr_elements = extract_elements_from_sequence(organ_id, seq_id, seq)
                            if chr_elements:
                                organism_rows.extend(chr_elements)
                                analysed_chr += 1
                                chr_logs.append(f"[DONE-CHR] {organ_id} | {seq_id} | gene rows: {len(chr_elements)}")
                            else:
                                empty_chr += 1

        return organ_id, organism_rows, analysed_chr, skipped_chr, empty_chr, chr_logs

    except Exception as e:
        return organ_id, [], 0, 0, 0, [f"[ERROR] {organ_id}: {e}"]

# =========================================================
# MAIN
# =========================================================
def main():
    ensure_output_csv(OUTPUT_CSV)
    completed_pairs, completed_orgs = load_completed_pairs(OUTPUT_CSV)

    all_dirs = [str(p) for p in Path(INPUT_ROOT).iterdir() if p.is_dir()]
    all_dirs.sort()

    selected_dirs = []
    for genome_dir in all_dirs:
        organ_id = os.path.basename(genome_dir.rstrip("/"))
        if organ_id in completed_orgs:
            continue
        selected_dirs.append(genome_dir)
        if len(selected_dirs) >= BATCH_SIZE:
            break

    print("======================================", flush=True)
    print("[START]", flush=True)
    print(f"[INFO] Total available organisms in folder: {len(all_dirs)}", flush=True)
    print(f"[INFO] Workers: {MAX_WORKERS}", flush=True)
    print(f"[INFO] Output: {OUTPUT_CSV}", flush=True)
    print(f"[INFO] Feature type: target_features", flush=True)
    print(f"[INFO] Minimum chromosome size: {MIN_CHROMOSOME_BP:,}", flush=True)
    print(f"[INFO] Batch size requested this run: {BATCH_SIZE}", flush=True)
    print(f"[INFO] Already analysed chromosome pairs in CSV: {len(completed_pairs)}", flush=True)
    print(f"[INFO] Already analysed organisms in CSV: {len(completed_orgs)}", flush=True)
    print("======================================", flush=True)
    print("--------------------------------------------------", flush=True)
    print("[INFO] Fast batch mode: not scanning all organisms", flush=True)
    print(f"[INFO] Organisms selected for this run: {len(selected_dirs)}", flush=True)
    print("[INFO] Selected organisms for this run:", flush=True)
    for d in selected_dirs:
        print(f"  - {os.path.basename(d.rstrip('/'))}", flush=True)
    print("--------------------------------------------------", flush=True)

    if not selected_dirs:
        print("[INFO] No new organisms found to analyze.", flush=True)
        return

    worker_args = [(genome_dir, completed_pairs) for genome_dir in selected_dirs]
    total_rows = 0

    with Pool(processes=MAX_WORKERS) as pool:
        for i, result in enumerate(pool.imap_unordered(process_one_organism, worker_args), start=1):
            organ_id, rows, analysed_chr, skipped_chr, empty_chr, chr_logs = result

            # Print individual chromosome progress logs sequentially
            for log in chr_logs:
                print(log, flush=True)

            append_rows(OUTPUT_CSV, rows)
            total_rows += len(rows)

            print(
                f"[DONE-ORG] {organ_id} | "
                f"new chromosomes analysed: {analysed_chr} | "
                f"already done skipped: {skipped_chr} | "
                f"no-feature chromosomes: {empty_chr} | "
                f"new rows: {len(rows)}",
                flush=True
            )
            print(
                f"[PROGRESS] {i}/{len(worker_args)} selected organisms finished | "
                f"new gene rows saved this run: {total_rows}",
                flush=True
            )

    print("======================================", flush=True)
    print("[FINISHED]", flush=True)
    print(f"[INFO] Batch size requested: {BATCH_SIZE}", flush=True)
    print(f"[INFO] Selected organisms this run: {len(selected_dirs)}", flush=True)
    print(f"[INFO] Finished organism jobs this run: {len(selected_dirs)}", flush=True)
    print(f"[INFO] New gene rows saved this run: {total_rows}", flush=True)
    print(f"[INFO] Output file: {OUTPUT_CSV}", flush=True)
    print("======================================", flush=True)

if __name__ == "__main__":
    main()
