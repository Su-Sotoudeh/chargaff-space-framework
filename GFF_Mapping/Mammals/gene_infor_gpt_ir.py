import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)
import os
import gzip
import pandas as pd
import numpy as np
import time
import re
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
import collections

# ----------------- Configuration & Decoding -----------------

# Minimum chromosome length: 5MB
MIN_CHROMOSOME_BP = 5000000

def is_valid_chromosome_name(description, seq_len):
    if seq_len < MIN_CHROMOSOME_BP:
        return False
    text = re.sub(r"\s+", " ", str(description).lower()).strip()
    bad_keywords = ["scaffold", "contig", "unlocalized", "unplaced", "patch", "plasmid", "mitochond", "chloroplast"]
    if any(word in text for word in bad_keywords):
        return False
    # Standard chromosome patterns
    if re.search(r"\bchr[\s._-]*[0-9xyzwm]+\b", text): return True
    if description.upper().startswith("NC_"): return True
    if "chromosome" in text: return True
    if 'complete sequence' in text or 'complete genome' in text: return True
    return False

def nucleotide_analysis_decoded(sequence):
    """
    Counts symbols from Masked FASTA and groups them.
    Rows: +Int, -Int, +UTR, -UTR, +CDS, -CDS, +TEs, -TEs, Shared, Unannotated
    """
    data = np.zeros((10, 4))
    sym_map = {
        'Q':(0,0), 'W':(0,1), 'E':(0,2), 'R':(0,3), # +Int (A,T,C,G)
        'q':(1,0), 'w':(1,1), 'e':(1,2), 'r':(1,3), # -Int
        'S':(2,0), 'D':(2,1), 'F':(2,2), 'H':(2,3), # +UTR
        's':(3,0), 'd':(3,1), 'f':(3,2), 'h':(3,3), # -UTR
        'Z':(4,0), 'X':(4,1), 'V':(4,2), 'B':(4,3), # +CDS
        'z':(5,0), 'x':(5,1), 'v':(5,2), 'b':(5,3), # -CDS
        # TEs group (DNA, LTR, SINE, LINE combined)
        '1':(6,0), '2':(6,1), '3':(6,2), '4':(6,3), 'L':(6,0), 'M':(6,1), '@':(6,2), '#':(6,3),
        'O':(6,0), 'P':(6,1), 'J':(6,2), 'K':(6,3), '*':(6,0), 'Y':(6,1), 'U':(6,2), 'I':(6,3),
        '5':(7,0), '6':(7,1), '7':(7,2), '8':(7,3), 'l':(7,0), 'm':(7,1), '$':(7,2), '%':(7,3),
        'o':(7,0), 'p':(7,1), 'j':(7,2), 'k':(7,3), '+':(7,0), 'y':(7,1), 'u':(7,2), 'i':(7,3),
        # Shared and Unannotated
        'A':(8,0), 'T':(8,1), 'C':(8,2), 'G':(8,3),
        'a':(9,0), 't':(9,1), 'c':(9,2), 'g':(9,3)
    }

    counts = collections.Counter(sequence)
    for sym, count in counts.items():
        if sym in sym_map:
            row, col = sym_map[sym]
            data[row, col] += count

    df = pd.DataFrame(data, columns=['A', 'T', 'C', 'G'])
    df.index = ['+Int', '-Int', '+UTR', '-UTR', '+CDS', '-CDS', '+TEs', '-TEs', 'Shared', 'Unannotated']
    return df

# ----------------- Parallel Worker -----------------

def process_one_organism(organism, input_root, output_csv, folder_path):
    """
    Processes one organism: Parses FASTA, calculates stats, saves individual CSV,
    and returns DataFrame for the master CSV.
    """
    try:
        # Time logging as per previous version
        start_time = time.time()
        print(f"[START] {organism}")

        organism_dir = os.path.join(input_root, organism)
        fasta_files = [f for f in os.listdir(organism_dir) if f.lower().endswith(('.fna.gz', '.fa.gz', '.fasta.gz'))]

        if not fasta_files:
            print(f"[SKIP] {organism}: No FASTA files found.")
            return None

        organism_rows = []

        for f_name in fasta_files:
            full_fasta_path = os.path.join(organism_dir, f_name)
            with gzip.open(full_fasta_path, 'rt') as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    seq = str(record.seq)
                    if not is_valid_chromosome_name(record.description, len(seq)):
                        continue

                    # Masking 'N's for length calculation
                    n_count = seq.count('n') + seq.count('N')
                    Length = len(seq) - n_count
                    if Length <= 0: continue

                    # Whole Genome Counts (Derived from counts of all specific symbols)
                    # This logic ensures we accurately count A,T,C,G across all masked categories
                    df_counts = nucleotide_analysis_decoded(seq)

                    # Totals for Whole Genome stats
                    total_counts = df_counts.sum() # Sums A, T, C, G columns across all rows
                    A_W, T_W, C_W, G_W = total_counts['A'], total_counts['T'], total_counts['C'], total_counts['G']

                    denom_W = (A_W + T_W + C_W + G_W)
                    gc_w = round((G_W + C_W) / denom_W * 100, 4) if denom_W != 0 else 0
                    a2t_w = round(A_W / T_W, 4) if T_W != 0 else 0
                    g2c_w = round(G_W / C_W, 4) if C_W != 0 else 0

                    elements = ['Int', 'UTR', 'CDS', 'TEs', 'Share', 'Unanot']
                    portions = []
                    main_dict = {
                        'Name': organism,
                        'Chr_ID': record.id,
                        'Length(M)': round(Length/1e6, 6),
                        'GC_W': gc_w,
                        'A2T_W': a2t_w,
                        'G2C_W': g2c_w
                    }

                    idx_ptr = 0
                    for el in elements:
                        if el in ("Share", "Unanot"):
                            vals = df_counts.iloc[idx_ptr].to_numpy() # 1 row
                            Ap, Tp, Cp, Gp = vals / 2
                            Ac, Tc, Cc, Gc = vals / 2
                            idx_ptr += 1
                        else:
                            Ap, Tp, Cp, Gp = df_counts.iloc[idx_ptr].to_numpy()   # Plus row
                            Ac, Tc, Cc, Gc = df_counts.iloc[idx_ptr+1].to_numpy() # Minus row
                            idx_ptr += 2

                        total_el = sum([Ap, Tp, Cp, Gp, Ac, Tc, Cc, Gc])
                        portions.append(total_el / Length * 100)

                        p_sum = sum([Ap, Tp, Cp, Gp])
                        c_sum = sum([Ac, Tc, Cc, Gc])

                        def get_stats(a, t, c, g, s):
                            if s == 0: return [0]*7
                            return [round((c+g)/s*100, 4), round(s/Length*100, 4), round(a/s*100, 4),
                                    round(t/s*100, 4), round(c/s*100, 4), round(g/s*100, 4), round(a/t, 4) if t!=0 else 0]

                        p_stats = get_stats(Ap, Tp, Cp, Gp, p_sum)
                        c_stats = get_stats(Ac, Tc, Cc, Gc, c_sum)

                        main_dict.update({
                            f'{el}_GC_P': p_stats[0], f'{el}_GC_C': c_stats[0],
                            f'{el}_P%': p_stats[1], f'{el}_C%': c_stats[1],
                            f'{el}_P2C': round(p_sum/c_sum, 4) if c_sum != 0 else 0,
                            f'{el}_A2T_ano': round((Ap+Ac)/(Tp+Tc), 4) if (Tp+Tc)!=0 else 0,
                            f'{el}_G2C_ano': round((Gp+Gc)/(Cp+Cc), 4) if (Cp+Cc)!=0 else 0,
                            f'{el}_P_A2T': p_stats[6], f'{el}_C_A2T': c_stats[6]
                        })

                    row_df = pd.DataFrame([main_dict])
                    # Insert derived summary columns
                    row_df.insert(3, 'Wh_ano', round(sum(portions[:4]), 4)) # Int + UTR + CDS + TEs
                    row_df.insert(4, 'Int_ano', round(portions[0], 4))
                    row_df.insert(5, 'UTR_ano', round(portions[1], 4))
                    row_df.insert(6, 'CDS_ano', round(portions[2], 4))
                    row_df.insert(7, 'TEs_ano', round(portions[3], 4))
                    row_df.insert(8, 'sum', round(sum(portions), 4))

                    organism_rows.append(row_df)

        if organism_rows:
            final_org_df = pd.concat(organism_rows, ignore_index=True)

            # Save individual organism CSV as per previous logic
            org_csv_path = os.path.join(folder_path, f"{organism}_masked_data.csv")
            final_org_df.to_csv(org_csv_path, index=False)

            elapsed = time.time() - start_time
            print(f"[DONE] {organism} (Time: {elapsed:.2f}s)")
            return final_org_df

    except Exception as e:
        print(f"[ERROR] {organism}: {e}")
    return None

# ----------------- Main Controller -----------------

def run_analysis(input_root, group_name):
    # Output Setup
    start_date = datetime.now().strftime("%Y-%m-%d")
    folder_name = os.path.join(str(config.output_dir("GFF_Mapping", group_name)), f"{group_name}_masked_analysis_{start_date}")
    os.makedirs(folder_name, exist_ok=True)

    master_csv = os.path.join(folder_name, f"{group_name}_masked_master.csv")

    # Discovery & Resume Logic
    all_orgs = sorted([d for d in os.listdir(input_root) if os.path.isdir(os.path.join(input_root, d))])

    processed_orgs = []
    if os.path.exists(master_csv):
        try:
            existing_df = pd.read_csv(master_csv, usecols=['Name'])
            processed_orgs = existing_df['Name'].unique().tolist()
        except:
            pass

    to_process = [o for o in all_orgs if o not in processed_orgs]

    print(f"Group: {group_name}")
    print(f"Total Organisms: {len(all_orgs)}")
    print(f"Already Processed: {len(processed_orgs)}")
    print(f"To Process: {len(to_process)}\n" + "-"*30)

    # Parallel Execution
    with ProcessPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(process_one_organism, org, input_root, master_csv, folder_name): org for org in to_process}

        for future in as_completed(futures):
            res_df = future.result()
            if res_df is not None:
                # Append to master CSV
                header = not os.path.exists(master_csv)
                res_df.to_csv(master_csv, mode='a', index=False, header=header)

if __name__ == "__main__":
    # Change these paths as needed
    INPUT_DIR = str(config.require(config.MASKED_FASTAS, "Put the masked FASTA files into Fasta_GFF_inputs/Masked_Fastas/<organism>/ (see Fasta_GFF_inputs/README.md)."))   # Fasta_GFF_inputs/Masked_Fastas/<organism>/
    GROUP = "Mammals"

    run_analysis(INPUT_DIR, GROUP)
