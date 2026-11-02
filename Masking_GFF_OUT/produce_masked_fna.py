import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import os
import pandas as pd
from Bio.Seq import Seq
import time
import sys
import gzip
from multiprocessing import cpu_count, Manager
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings
import gc
warnings.filterwarnings('ignore')

from all_functions import exon_and_intron_masking, opening_fasta_file
from all_functions import TE_masking_after_exon

# =====================================
# PATHS
# =====================================
# input : Fasta_GFF_inputs/Mammals/<organism>/  (<chr>.fna + <chr>.fna.out from RepeatMasker + *.gff.gz)
# output: Fasta_GFF_inputs/Masked_Fasta/<organism>/<chr>.fna.gz
input_root = str(config.require(config.MAMMALS_WORK_DIR, "Run the RepeatMasker steps first (see RepeatMasker/README.md): Fasta_GFF_inputs/Mammals/<organism>/ needs <chr>.fna, <chr>.fna.out and the *.gff.gz file."))
output_root = str(config.MASKED_FASTAS)

if not os.path.exists(output_root):
    os.makedirs(output_root)

FORCE_RESTART = False
NUM_WORKERS = max(1, cpu_count() - 2)
MAX_RETRIES = 3

# -------------------------------------
# PROGRESS TRACKER
# -------------------------------------
class ProgressTracker:
    def __init__(self):
        self.total = 0
        self.completed = 0
        self.failed = 0
        self.lock = Manager().Lock()

    def set_total(self, total):
        self.total = total

    def increment(self, success=True):
        with self.lock:
            if success:
                self.completed += 1
            else:
                self.failed += 1
            return self.completed, self.failed

progress = ProgressTracker()


# ==========================================================
# PROCESS ONE CHROMOSOME → MASK → SAVE FASTA.GZ
# ==========================================================
def process_chromosome(args):

    organism, f, o, gff_file, path, task_num, total_tasks, retry_count = args
    ch_id = f[:-4]
    row_index_name = f"{ch_id}"

    # Make organism subfolder in output
    org_out_dir = os.path.join(output_root, organism)
    if not os.path.exists(org_out_dir):
        os.makedirs(org_out_dir)

    out_fasta_file = os.path.join(org_out_dir, f"{row_index_name}.fna.gz")

    # Skip if exists and not restart
    if os.path.exists(out_fasta_file) and not FORCE_RESTART:
        print(f"[{task_num}/{total_tasks}] SKIP: {row_index_name} already done.")
        progress.increment(True)
        return row_index_name

    for attempt in range(MAX_RETRIES):
        try:
            print(f"[{task_num}/{total_tasks}] >>> START: {organism}/{row_index_name} (Attempt {attempt+1})")

            # Load sequence
            seq_obj = opening_fasta_file(path + f)

            # Mask exons/introns
            exon_mask_seq = exon_and_intron_masking(seq_obj, path + gff_file, ch_id)

            # Mask TEs
            final_masked_sequence = TE_masking_after_exon(exon_mask_seq, path + o, ch_id)


            # --- SAVE COMPRESSED FASTA (.gz) ---
            with gzip.open(out_fasta_file, "wt") as fw:
                fw.write(f">{row_index_name}\n")
                fw.write(final_masked_sequence + "\n")

            del seq_obj, exon_mask_seq, final_masked_sequence
            gc.collect()

            print(f"[{task_num}/{total_tasks}] <<< SAVED: {organism}/{row_index_name}.fna.gz")
            progress.increment(True)
            return row_index_name

        except Exception as e:
            print(f"[{task_num}/{total_tasks}] ❌ ERROR {row_index_name}: {e} (Attempt {attempt+1})")
            gc.collect()

            if attempt < MAX_RETRIES - 1:
                wait = (attempt + 1) * 5
                print(f"    Retrying in {wait} sec...")
                time.sleep(wait)
            else:
                progress.increment(False)
                return None

    return None


# ==========================================================
# MAIN
# ==========================================================
if __name__ == "__main__":

    start = time.time()

    print("="*70)
    print(" GENOME MASKING → SAVE STRUCTURED FASTA.GZ ")
    print("="*70)
    print(f"Workers: {NUM_WORKERS}")
    print("="*70)

    # Find previously processed chromosomes
    # Find previously processed chromosomes (skip empty files)
    already_done = set()

    for org in os.listdir(output_root):
        opath = os.path.join(output_root, org)
        if not os.path.isdir(opath):
            continue

        for f in os.listdir(opath):
            if not f.endswith(".fna.gz"):
                continue

            full_path = os.path.join(opath, f)

            # --- NEW: check if empty ---
            if os.path.getsize(full_path) == 0:
                print(f"EMPTY FILE FOUND → Deleting: {full_path}")
                os.remove(full_path)
                continue

            # If not empty, mark as done
            already_done.add(f.replace(".fna.gz", ""))


    print(f"Resume: {len(already_done)} existing chromosomes skipped.\n")

    organisms = os.listdir(input_root)
    tasks = []
    total_tasks = 0

    for organism in organisms:

        path = os.path.join(input_root, organism) + "/"
        if not os.path.exists(path):
            continue

        try:
            files = os.listdir(path)
        except:
            continue

        fna = [i for i in files if i.endswith(".fna")]
        out = [j for j in files if j.endswith(".out")]
        gff = [k for k in files if k.endswith(".gff.gz")]

        if not fna or not out or not gff:
            continue

        gff_file = gff[0]

        # Map .fna ↔ .out
        out_map = {o_file.replace(".fna.out", ""): o_file for o_file in out}

        for f in fna:
            ch = f.replace(".fna", "")
            if ch in already_done:
                continue
            if ch not in out_map:
                continue

            tasks.append((organism, f, out_map[ch], gff_file, path, 0, 0, 0))

    # Assign task numbers
    for i in range(len(tasks)):
        org, f, o, gff_file, path, _, _, _ = tasks[i]
        tasks[i] = (org, f, o, gff_file, path, i+1, len(tasks), 0)

    progress.set_total(len(tasks))

    print("="*60)
    print(f" WILL PROCESS: {len(tasks)} chromosomes")
    print("="*60)

    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as exe:
        future_map = {exe.submit(process_chromosome, t): t for t in tasks}

        for fut in as_completed(future_map):
            fut.result()

    print("\n" + "="*70)
    print(" DONE! Masked FASTA.GZ files written to:", output_root)
    print("="*70)
    print("Time:", round(time.time() - start, 2), "sec")
