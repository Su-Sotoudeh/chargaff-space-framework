import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import os
import pandas as pd
from Bio.Seq import Seq
import time
import sys
from multiprocessing import cpu_count, Manager
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings
import gc
warnings.filterwarnings('ignore')


#DNA  ATCG >> +1234 -5678
#LTR  ATCG >> +LM@# -lm$%
#SINE ATCG >> +OPJK -opjk
#LINE ATCG >> +*YUI -+yui
#Int  ATCG >> +QWER -qwer
#Exe  ATCG >> +SDFH -sdfh
#CDS  ATCG >> +ZXVB -zxvb
#shre ATCG >> ATCG (shared between genes)
#un   ATCG >> atcg




# --- CONFIGURATION ---
one_Trabite = str(config.require(config.MASKED_FASTAS, "Put the masked FASTA files into Fasta_GFF_inputs/Masked_Fastas/<organism>/ (see Fasta_GFF_inputs/README.md).")) + '/'   # Fasta_GFF_inputs/Masked_Fastas/
output_file = str(config.output_dir(5) / '4k_kmers_normalised_MCRRC_PC0UW_DNA.csv')
FORCE_RESTART = False
NUM_WORKERS = max(1, cpu_count()-1)  # Use fewer workers to avoid memory issues
MAX_RETRIES = 3  # Retry failed tasks up to 3 times

# --- K-MER GENERATION ---
A =  '1'
T =  '2'
C =  '3'
G =  '4'
a =  '5'
t =  '6'
c =  '7'
g =  '8'




k_mers = []
N = ['A', 'T', 'C', 'G']

for first in N:
    for secound in N:
        for third in N:
            for forth in N:
                k_mers.append(first+secound+third+forth)
main_k_mers = k_mers.copy() # our base wor is in capital  ATCG k-mer produced,next we can it lower or others
trans_tabel = str.maketrans({'A':A,'T':T,'C':C,'G':G,'a':a,'t':t,'c':c,'g':g})
Mps = [i.translate(trans_tabel) for i in k_mers]

C_trans_table = str.maketrans({T:A,A:T,C:G,G:C})
Cps = [i.translate(C_trans_table) for i in Mps]
Rps = []
for k in Mps:
    B = list(k)
    B.reverse()
    B = ''.join(B)
    Rps.append(B)

k_mers_RC = [str(Seq(k).reverse_complement()) for k in k_mers]
RCps = [i.translate(trans_tabel) for i in k_mers_RC]


Mcs = [str(i.lower().translate(trans_tabel)).replace(A,a) for i in k_mers]
C_trans_table = str.maketrans({t:a,a:t,c:g,g:c})
Ccs = [i.translate(C_trans_table) for i in Mcs]
Rcs = []
for k in Mcs:
    B = list(k)
    B.reverse()
    B = ''.join(B)
    Rcs.append(B)
RCcs = [str(i.lower().translate(trans_tabel)).replace(A,a) for i in k_mers_RC]



Mus =  [i.lower() for i in k_mers]
U_trans_table = str.maketrans({'a':'t','t':'a','c':'g','g':'c'})
Cus = [i.translate(U_trans_table) for i in Mus]
Rus = []
for k in Mus:
    B = list(k)
    B.reverse()
    B = ''.join(B)
    Rus.append(B)

RCus = []
for k in Cus:
    B = list(k)
    B.reverse()
    B = ''.join(B)
    RCus.append(B)





list_k_mers = [[Mps[i],Cps[i],Rps[i], RCps[i], Mcs[i],Ccs[i],Rcs[i], RCcs[i],Mus[i],Cus[i],Rus[i], RCus[i]] for i in range(len(k_mers))]

K_count_dic = {}
for i, sub in enumerate(list_k_mers):
    for j, sub2 in enumerate(sub):
        K_count_dic[sub2] = 0

print(list_k_mers)
column_names = []
for kmer in main_k_mers:
    for feat in ["+M","+C","+R", "+RC", "-M","-C","-R", "-RC","0M","0C","0R","0RC","UM","UC","UR","URC","WM","WC","WR","WRC"]:
        column_names.append(f"{kmer}_{feat}")


#b = both (+ and -)


# --- IMPORT FUNCTIONS ---
from all_functions import exon_and_intron_masking, opening_fasta_file
from all_functions import TE_masking_after_exon

# --- GLOBAL COUNTERS ---
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

# --- PROCESS FUNCTION WITH RETRY ---
def process_chromosome(args):
    """Worker function with retry logic"""
    organism, f, _, _, path, task_num, total_tasks, retry_count = args
    ch_id = f[:-7]
    row_index_name = f"{organism}_{ch_id}"

    for attempt in range(MAX_RETRIES):
        try:
            print(f"[{task_num}/{total_tasks}] >>> START: {organism} | {ch_id} (Attempt {attempt+1})", flush=True)

            # NEW: directly load pre‑masked fasta
            exon_TE_masked_sequecne = opening_fasta_file(path + f)



            #wtransaltiin_table = str.maketrans('1234LM@#OPJKQWERSDFHZXVB5678lm$%opjkqwersdfhzxvb','*YUI*YUI*YUI*YUI*YUI*YUI+yui+yui+yui+yui+yui+yui')
            #exon_TE_masked_sequecne = exon_TE_masked_sequecne.translate(wtransaltiin_table)
            L = len(exon_TE_masked_sequecne)


            # rolling k-mer window
            for i in range(L - 3):
                kmer = exon_TE_masked_sequecne[i:i+4]
                if kmer in K_count_dic:
                    K_count_dic[kmer] += 1

            row_data = {}
            for index, l in enumerate(main_k_mers):
                Mp  = K_count_dic[list_k_mers[index][0]]
                Cp  = K_count_dic[list_k_mers[index][1]]
                Rp  = K_count_dic[list_k_mers[index][2]]
                RCp = K_count_dic[list_k_mers[index][3]]
                Mc  = K_count_dic[list_k_mers[index][4]]
                Cc  = K_count_dic[list_k_mers[index][5]]
                Rc  = K_count_dic[list_k_mers[index][6]]
                RCc = K_count_dic[list_k_mers[index][7]]
                Mb  = Mp + Mc
                Cb = Cp + Cc
                Rb = Rp + Rc
                RCb = RCp + RCc
                Mu  = K_count_dic[list_k_mers[index][8]]
                Cu  = K_count_dic[list_k_mers[index][9]]
                Ru  = K_count_dic[list_k_mers[index][10]]
                RCu = K_count_dic[list_k_mers[index][11]]
                Mw  = Mp + Mc + Mu
                Cw  = Cp + Cc + Cu
                Rw  = Rp + Rc + Ru
                RCw = RCp + RCc + RCu



                total = Mp + Cp + Rp + RCp + Mc + Cc + Rc + RCc + Mb + Cb + Rb + RCb + Mu + Cu + Ru + RCu + Mw + Cw + Rw + RCw
                if total == 0: total = 1
                T_Mp = (Mp +Cp + Rp + RCp)
                pMp  = round(Mp / T_Mp * 100, 2)
                pCp  = round(Cp / T_Mp * 100, 2)
                pRp  = round(Rp / T_Mp * 100, 2)
                pRCp = round(RCp/ T_Mp * 100, 2)
                T_Mc = (Mc +Cc + Rc + RCc)
                pMc  = round(Mc / T_Mc * 100, 2)
                pCc  = round(Cc / T_Mc * 100, 2)
                pRc  = round(Rc / T_Mc * 100, 2)
                pRCc = round(RCc/ T_Mc * 100, 2)
                T_Mb = (Mb + Cb + Rb + RCb)
                pMb   = round(Mb / T_Mb * 100, 2)
                pCb   = round(Cb / T_Mb * 100, 2)
                pRb   = round(Rb / T_Mb * 100, 2)
                pRCb  = round(RCb/ T_Mb * 100, 2)
                T_Mu = (Mu + Cu + Ru + RCu)
                pMu   = round(Mu / T_Mu * 100, 2)
                pCu   = round(Cu / T_Mu * 100, 2)
                pRu   = round(Ru / T_Mu * 100, 2)
                pRCu  = round(RCu/ T_Mu * 100, 2)
                T_Mw = (Mw + Cw + Rw + RCw)
                pMw   = round(Mw / T_Mw * 100, 2)
                pCw   = round(Cw / T_Mw * 100, 2)
                pRw   = round(Rw / T_Mw * 100, 2)
                pRCw  = round(RCw/ T_Mw * 100, 2)


                row_data[f"{l}_+M"] = pMp
                row_data[f"{l}_+C"] = pCp
                row_data[f"{l}_+R"] = pRp
                row_data[f"{l}_+RC"] = pRCp
                row_data[f"{l}_-M"] = pMc
                row_data[f"{l}_-C"] = pCc
                row_data[f"{l}_-R"] = pRc
                row_data[f"{l}_-RC"] = pRCc
                row_data[f"{l}_0M"] = pMb
                row_data[f"{l}_0C"] = pCb
                row_data[f"{l}_0R"] = pRb
                row_data[f"{l}_0RC"] = pRCb
                row_data[f"{l}_UM"] = pMu
                row_data[f"{l}_UC"] = pCu
                row_data[f"{l}_UR"] = pRu
                row_data[f"{l}_URC"] = pRCu
                row_data[f"{l}_WM"] = pMw
                row_data[f"{l}_WC"] = pCw
                row_data[f"{l}_WR"] = pRw
                row_data[f"{l}_WRC"] = pRCw



            # Clean up memory
            del  exon_TE_masked_sequecne
            gc.collect()

            current, failed = progress.increment(success=True)
            print(f"[{task_num}/{total_tasks}] <<< DONE: {organism} | {ch_id} | Progress: {current}/{total_tasks}", flush=True)

            return (row_index_name, row_data, organism, ch_id)

        except Exception as e:
            error_msg = str(e)
            print(f"[{task_num}/{total_tasks}] ❌ ERROR {organism} | {ch_id}: {error_msg} (Attempt {attempt+1}/{MAX_RETRIES})", flush=True)

            # Clean up memory on failure
            gc.collect()

            if attempt < MAX_RETRIES - 1:
                wait_time = (attempt + 1) * 5  # Wait 5, 10, 15 seconds
                print(f"    ⏳ Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print(f"    ❌ All {MAX_RETRIES} attempts failed for {organism} | {ch_id}")
                progress.increment(success=False)

                # Print detailed error
                import traceback
                traceback.print_exc()

                return None

    return None

# --- SAVE FUNCTION ---
def save_to_csv(result, output_file):
    """Save single chromosome result to CSV"""
    if result is None:
        return False
    try:
        df_new = pd.DataFrame([result[1]], index=[result[0]])
        df_new.to_csv(output_file, mode='a', header=False)
        return True
    except Exception as e:
        print(f"Save Error: {e}", flush=True)
        return False

# --- MAIN ---
if __name__ == '__main__':
    start = time.time()

    print("="*70)
    print("PARALLEL K-MER ANALYSIS - WITH ERROR HANDLING & RETRY")
    print("="*70)
    print(f"Output: {os.path.abspath(output_file)}")
    print(f"Workers: {NUM_WORKERS} CPU cores (using CPU_COUNT - 2 for stability)")
    print(f"Max Retries: {MAX_RETRIES}")
    print("="*70)

    # =====================================================
    # RESUME LOGIC
    # =====================================================
    processed_chrs = set()
    existing_count = 0

    if FORCE_RESTART and os.path.exists(output_file):
        print("\n⚠ FORCE RESTART: Deleting old file...")
        os.remove(output_file)

    elif os.path.exists(output_file):
        try:
            df_existing = pd.read_csv(output_file, index_col=0)
            if len(df_existing) > 0:
                processed_chrs = set(df_existing.index)
                existing_count = len(processed_chrs)
                print(f"\n📂 RESUME MODE: Found existing file with {existing_count} processed chromosomes")
        except Exception as e:
            print(f"⚠ Warning: Could not read existing file: {e}")

    # List organisms
    organisms = os.listdir(one_Trabite)
    print(f"\n📁 Found {len(organisms)} organisms in folder.")

    # --- COLLECT ALL TASKS ---
    tasks = []
    organism_counts = {}

    for organism in organisms:
        path = one_Trabite + organism + '/'
        if not os.path.exists(path):
            continue
        try:
            files = os.listdir(path)
        except:
            continue

        # NEW: Load masked fasta files ONLY
        fna = [i for i in files if i.lower().endswith(('.fna', '.fa', '.fasta', '.fna.gz', '.fa.gz', '.fasta.gz'))]


        if len(fna) == 0:
            print(f"  ⚠ {organism}: No fasta files found - SKIPPED")
            continue




        org_total = 0
        org_skipped = 0
        org_to_process = 0

        for f in fna:
            # Remove extension safely (.fna / .fa / .fasta)
            if f.endswith('.gz') :
                base_id = f[:-7]
            elif f.endswith('.fna'):
                base_id = f[:-4]
            else:
                print('please chack type of file')
                print(f)

            full_name = f"{organism}_{base_id}"
            org_total += 1

            # Resume mode: skip if already processed
            if full_name in processed_chrs:
                org_skipped += 1
                continue

            # Every FASTA is a valid task (no .out or .gff needed)
            org_to_process += 1

            # Pass None for o and gff_file
            tasks.append((organism, f, None, None, path, 0, 0, 0))

        # Store summary
        if org_total > 0:
            organism_counts[organism] = {
                'total': org_total,
                'to_process': org_to_process,
                'skipped': org_skipped
            }


    # Add task numbers
    total_chr_to_process = len(tasks)
    for i, task in enumerate(tasks):
        tasks[i] = (task[0], task[1], task[2], task[3], task[4], i+1, total_chr_to_process, 0)

    progress.set_total(total_chr_to_process)

    # =====================================================
    # PRINT SUMMARY
    # =====================================================
    print("\n" + "="*70)
    print("PROCESSING SUMMARY")
    print("="*70)

    if existing_count > 0:
        print(f"📊 Already Processed: {existing_count} chromosomes")

    print(f"\n📊 Organism Breakdown:")
    for org, counts in organism_counts.items():
        total = counts['total']
        to_proc = counts['to_process']
        skipped = counts['skipped']
        print(f"  📦 {org}: {to_proc} to process, {skipped} already done, {total} total")

    print("="*70)
    print(f"  🚀 WILL PROCESS: {total_chr_to_process} chromosomes")
    print("="*70)

    if total_chr_to_process == 0:
        print("\n✅ All chromosomes already processed!")
        sys.exit()

    # Create output file with header
    if not os.path.exists(output_file):
        df_header = pd.DataFrame(columns=column_names)
        df_header.to_csv(output_file, index=True, index_label='Organism_Chr')

    processed_count = existing_count
    failed_tasks = []

    print(f"\n🚀 Starting parallel processing with {NUM_WORKERS} workers...")
    print("="*70)

    # --- PARALLEL EXECUTION ---
    completed_orgs = {}

    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        future_to_task = {executor.submit(process_chromosome, task): task for task in tasks}

        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                result = future.result()

                if save_to_csv(result, output_file):
                    processed_count += 1

                    if result:
                        org = result[2]
                        chr_name = result[3]
                        if org not in completed_orgs:
                            completed_orgs[org] = 0
                        completed_orgs[org] += 1
                        print(f"  💾 SAVED: {org} | {chr_name} | [{processed_count}/{existing_count + total_chr_to_process}]")

            except Exception as e:
                # This catches the "process terminated abruptly" error
                task_num = task[5]
                org = task[0]
                chr_id = task[1][:-4]
                print(f"❌ [{task_num}] Process failed for {org} | {chr_id}: {str(e)}")
                failed_tasks.append(task)
                progress.increment(success=False)

    # =====================================================
    # RETRY FAILED TASKS
    # =====================================================
    if failed_tasks:
        print("\n" + "="*70)
        print(f"🔄 RETRYING {len(failed_tasks)} FAILED TASKS...")
        print("="*70)

        # Update retry count in tasks
        retry_tasks = []
        for task in failed_tasks:
            # Increment retry count
            retry_tasks.append((task[0], task[1], task[2], task[3], task[4], task[5], task[6], task[7] + 1))

        with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
            future_to_task = {executor.submit(process_chromosome, task): task for task in retry_tasks}

            for future in as_completed(future_to_task):
                task = future_to_task[future]
                try:
                    result = future.result()

                    if save_to_csv(result, output_file):
                        processed_count += 1

                        if result:
                            org = result[2]
                            chr_name = result[3]
                            if org not in completed_orgs:
                                completed_orgs[org] = 0
                            completed_orgs[org] += 1
                            print(f"  💾 RETRY SUCCESS: {org} | {chr_name}")

                except Exception as e:
                    org = task[0]
                    chr_id = task[1][:-4]
                    print(f"❌ Retry failed for {org} | {chr_id}: {str(e)}")

    # --- FINAL SUMMARY ---
    print("\n" + "="*70)
    print("FINAL SUMMARY")
    print("="*70)

    print(f"\n📊 Organism Progress:")
    for org, counts in organism_counts.items():
        done = completed_orgs.get(org, 0)
        to_proc = counts['to_process']
        skipped = counts['skipped']
        total = counts['total']
        actual_done = done + skipped

        if actual_done == total:
            status = "✅ COMPLETE"
        else:
            status = f"⏳ {done}/{to_proc} processed"

        print(f"  {status} | {org}: {done} new + {skipped} skipped = {actual_done}/{total}")

    elapsed = time.time() - start
    print(f"\n⏱️ Total time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
    print(f"💾 Output: {os.path.abspath(output_file)}")
    print(f"📊 Total processed in file: {processed_count}")
    print("="*70)
    print("✅ ALL DONE!")
    print("="*70)
