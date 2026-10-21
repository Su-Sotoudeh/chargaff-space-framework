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
import multiprocessing
from Bio import SeqIO
import sys

# ----------------- Helper Functions -----------------

def is_valid_chromosome_name(description):
    text = re.sub(r"\s+", " ", str(description).lower()).strip()
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
    if re.search(r"\bchr[\s._-]*[0-9xyzwm]+\b", text):
        return True
    if re.search(r">NC_", text):
        return True
    if re.search(r">CM_", text):
        return True
    if re.search(r">CP_", text):
        return True
    if re.search(r"whole genome shotgun sequence", text):
        return True
    if "chromosome un" in text:
        return False
    if re.search(r"\bchr\s*un\b", text):
        return False
    if "chromosome" in text:
        return True
    if 'complete sequence' in text:
        return True
    if 'complete genome' in text:
        return True
    if 'partial sequence' in text:
        return False
    return False

def nucleotide_analysis_for_exon(sequence_plus):
    A_entron_P = sequence_plus.count('Q')
    T_entron_P = sequence_plus.count('W')
    C_entron_P = sequence_plus.count('E')
    G_entron_P = sequence_plus.count('R')
    A_entron_C = sequence_plus.count('q')
    T_entron_C = sequence_plus.count('w')
    C_entron_C = sequence_plus.count('e')
    G_entron_C = sequence_plus.count('r')

    A_exon_P = sequence_plus.count('S')
    T_exon_P = sequence_plus.count('D')
    C_exon_P = sequence_plus.count('F')
    G_exon_P = sequence_plus.count('H')
    A_exon_C = sequence_plus.count('s')
    T_exon_C = sequence_plus.count('d')
    C_exon_C = sequence_plus.count('f')
    G_exon_C = sequence_plus.count('h')

    A_CDS_P = sequence_plus.count('Z')
    T_CDS_P = sequence_plus.count('X')
    C_CDS_P = sequence_plus.count('V')
    G_CDS_P = sequence_plus.count('B')
    A_CDS_C = sequence_plus.count('z')
    T_CDS_C = sequence_plus.count('x')
    C_CDS_C = sequence_plus.count('v')
    G_CDS_C = sequence_plus.count('b')

    A_shared = sequence_plus.count('A')
    T_shared = sequence_plus.count('T')
    C_shared = sequence_plus.count('C')
    G_shared = sequence_plus.count('G')
    A_unannotated = sequence_plus.count('a')
    T_unannotated = sequence_plus.count('t')
    C_unannotated = sequence_plus.count('c')
    G_unannotated = sequence_plus.count('g')

    A = [A_entron_P, A_entron_C, A_exon_P, A_exon_C, A_CDS_P, A_CDS_C, A_shared, A_unannotated]
    T = [T_entron_P, T_entron_C, T_exon_P, T_exon_C, T_CDS_P, T_CDS_C, T_shared, T_unannotated]
    C = [C_entron_P, C_entron_C, C_exon_P, C_exon_C, C_CDS_P, C_CDS_C, C_shared, C_unannotated]
    G = [G_entron_P, G_entron_C, G_exon_P, G_exon_C, G_CDS_P, G_CDS_C, G_shared, G_unannotated]
    dic = {'A': A, 'T': T, 'C': C, 'G': G}
    df = pd.DataFrame(dic)
    df.index = ['+entron', '-entron', '+UTR', '-UTR', '+CDS', '-CDS', 'Shared', 'Unannotated']
    return df

def exon_and_intron_masking(sequence, df_all, chr_id):
    ch_id_split = chr_id.split('.')[0]
    sequence = sequence.lower()
    sequence_plus = bytearray(sequence.encode('ascii'))
    Inreon_P_counter = 0
    Inreon_C_counter = 0
    UTR_P_counter = 0
    UTR_C_counter = 0
    CDS_P_counter = 0
    CDS_C_counter = 0

    df_all.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    filtered_df_all = df_all[df_all['seqid'].str.startswith(ch_id_split)]

    # First pass: genes
    filtered_df = filtered_df_all[filtered_df_all['type'].str.startswith('gene')]
    for annotation in filtered_df.to_numpy():
        start, end = int(annotation[3]), int(annotation[4])
        if end - start == 0:
            continue
        if annotation[6] == '+':
            Inreon_P_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('a'): sl[i] = ord('Q')
                elif ch == ord('t'): sl[i] = ord('W')
                elif ch == ord('c'): sl[i] = ord('E')
                elif ch == ord('g'): sl[i] = ord('R')
                elif ch == ord('q'): sl[i] = ord('A')
                elif ch == ord('w'): sl[i] = ord('T')
                elif ch == ord('e'): sl[i] = ord('C')
                elif ch == ord('r'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl
        elif annotation[6] == '-':
            Inreon_C_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('a'): sl[i] = ord('q')
                elif ch == ord('t'): sl[i] = ord('w')
                elif ch == ord('c'): sl[i] = ord('e')
                elif ch == ord('g'): sl[i] = ord('r')
                elif ch == ord('Q'): sl[i] = ord('A')
                elif ch == ord('W'): sl[i] = ord('T')
                elif ch == ord('E'): sl[i] = ord('C')
                elif ch == ord('R'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl

    # Second pass: exons
    filtered_df = filtered_df_all[filtered_df_all['type'].str.startswith('exon')]
    for annotation in filtered_df.to_numpy():
        start, end = int(annotation[3]), int(annotation[4])
        if end - start == 0:
            continue
        if annotation[6] == '+':
            UTR_P_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('Q'): sl[i] = ord('S')
                elif ch == ord('W'): sl[i] = ord('D')
                elif ch == ord('E'): sl[i] = ord('F')
                elif ch == ord('R'): sl[i] = ord('H')
                elif ch == ord('q'): sl[i] = ord('A')
                elif ch == ord('w'): sl[i] = ord('T')
                elif ch == ord('e'): sl[i] = ord('C')
                elif ch == ord('r'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl
        elif annotation[6] == '-':
            UTR_C_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('q'): sl[i] = ord('s')
                elif ch == ord('w'): sl[i] = ord('d')
                elif ch == ord('e'): sl[i] = ord('f')
                elif ch == ord('r'): sl[i] = ord('h')
                elif ch == ord('Q'): sl[i] = ord('A')
                elif ch == ord('W'): sl[i] = ord('T')
                elif ch == ord('E'): sl[i] = ord('C')
                elif ch == ord('R'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl

    # Third pass: CDS
    filtered_df = filtered_df_all[filtered_df_all['type'].str.startswith('CDS')]
    for annotation in filtered_df.to_numpy():
        start, end = int(annotation[3]), int(annotation[4])
        if end - start == 0:
            continue
        if annotation[6] == '+':
            CDS_P_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('S'): sl[i] = ord('Z')
                elif ch == ord('D'): sl[i] = ord('X')
                elif ch == ord('F'): sl[i] = ord('V')
                elif ch == ord('H'): sl[i] = ord('B')
                elif ch == ord('z'): sl[i] = ord('A')
                elif ch == ord('x'): sl[i] = ord('T')
                elif ch == ord('v'): sl[i] = ord('C')
                elif ch == ord('b'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl
        elif annotation[6] == '-':
            CDS_C_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('s'): sl[i] = ord('z')
                elif ch == ord('d'): sl[i] = ord('x')
                elif ch == ord('f'): sl[i] = ord('v')
                elif ch == ord('h'): sl[i] = ord('b')
                elif ch == ord('Z'): sl[i] = ord('A')
                elif ch == ord('X'): sl[i] = ord('T')
                elif ch == ord('V'): sl[i] = ord('C')
                elif ch == ord('B'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl

    return sequence_plus.decode('ascii'), Inreon_P_counter, Inreon_C_counter, UTR_P_counter, UTR_C_counter, CDS_P_counter, CDS_C_counter

# ----------------- Parallel Worker -----------------

def process_one_organism(final_number, organ_number, organism, all_genomes_folder_path, temp_csv_path):
    """
    Worker: Runs organism analysis and appends directly to the output CSV
    after completely analyzing all valid chromosomes in this organism.
    """
    try:
        if organism.endswith(".json") or organism.endswith(".jsonl"):
            return None

        organism_dir = os.path.join(all_genomes_folder_path, organism)
        files = os.listdir(organism_dir)
        fna = sorted([i for i in files if i.endswith('.fna.gz') and not i.startswith('cds')])
        gff = sorted([k for k in files if k.endswith('.gff.gz')])

        if not fna or not gff:
            return None

        try:
            fasta_file = gzip.open(os.path.join(organism_dir, fna[0]), 'rt')
        except (IOError, EOFError, IndexError):
            return None

        gff_path = os.path.join(organism_dir, gff[0])
        df_all = pd.read_csv(gff_path, sep='\t', comment='#')

        analysed_chr_counter = 0
        organism_rows = []

        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if 'mitochondrion' in seq_record.description:
                continue

            seq = str(seq_record.seq)
            l = len(seq) / 1000000
            if l < 1:
                continue
            if "chromosome" not in seq_record.description:
                continue

            ch_id = seq_record.id
            if ch_id.startswith('NW') or not is_valid_chromosome_name(seq_record.description):
                continue

            analysed_chr_counter += 1

            exon_mask_seq, Gene_P_counter, Gene_C_counter, Exon_P_counter, Exon_C_counter, CDS_P_counter, CDS_C_counter = \
                exon_and_intron_masking(seq, df_all, ch_id)
            if CDS_P_counter == 0 or CDS_C_counter == 0:
                continue

            A_W = seq.count('A') + seq.count("a")
            T_W = seq.count('T') + seq.count("t")
            C_W = seq.count('C') + seq.count("c")
            G_W = seq.count("G") + seq.count("g")
            N_W = seq.count('N') + seq.count("n")
            Length = len(seq) - N_W

            denom_W = (A_W + T_W + C_W + G_W)
            GC_countent_list_W = round((G_W + C_W) / denom_W * 100, 4) if denom_W != 0 else 0
            len_seq = round(Length / 1000000, 6)

            A2T = round(A_W / T_W, 4) if T_W != 0 else 0
            G2C = round(G_W / C_W, 4) if C_W != 0 else 0

            L = len(exon_mask_seq)
            df = nucleotide_analysis_for_exon(exon_mask_seq)
            df.loc["Chr_Len"] = {"A": L, 'T': L, "C": L, "G": L, "whole": L}
            n_count = exon_mask_seq.count('n') + exon_mask_seq.count('N')
            df.loc["Nn_number"] = {"A": n_count, 'T': n_count, "C": n_count, "G": n_count, "whole": n_count}
            df.loc["organism"] = {"A": organism, 'T': organism, "C": organism, "G": organism, "whole": organism}
            df.loc["chr_id"] = {"A": ch_id, 'T': ch_id, "C": ch_id, "G": ch_id, "whole": ch_id}

            gene_elements = ['Int', 'UTR', 'CDS', 'Share', 'Unanot']
            counter = 0
            portions = []
            main_dict = {
                'Name': organism,
                'Chr_ID': ch_id,
                'Length(M)': len_seq,
                'GC_W': GC_countent_list_W,
                'A2T_W': A2T,
                'G2C_W': G2C,
            }

            for element in gene_elements:
                if element in ("Share", "Unanot"):
                    Ap, Tp, Cp, Gp = df.iloc[counter]
                    Ap, Tp, Cp, Gp = Ap / 2, Tp / 2, Cp / 2, Gp / 2
                    Ac, Tc, Cc, Gc = df.iloc[counter]
                    Ac, Tc, Cc, Gc = Ac / 2, Tc / 2, Cc / 2, Gc / 2
                    counter += 1
                else:
                    Ap, Tp, Cp, Gp = df.iloc[counter]
                    Ac, Tc, Cc, Gc = df.iloc[counter + 1]
                    counter += 2
                try:
                    portions.append(sum([Ap, Tp, Cp, Gp, Ac, Tc, Cc, Gc]) / Length * 100)
                    P_portion = round(Ap + Tp + Cp + Gp, 4)
                    if P_portion == 0:
                        P_A_usage = P_T_usage = P_C_usage = P_G_usage = P_Portion = A2T_P = G2C_P = GC_countent_P = 0
                    else:
                        P_A_usage = round(Ap / P_portion * 100, 4)
                        P_T_usage = round(Tp / P_portion * 100, 4)
                        P_C_usage = round(Cp / P_portion * 100, 4)
                        P_G_usage = round(Gp / P_portion * 100, 4)
                        P_Portion = round(P_portion / Length * 100, 4)
                        A2T_P = round(Ap / Tp, 4) if Tp != 0 else 0
                        G2C_P = round(Gp / Cp, 4) if Cp != 0 else 0
                        GC_countent_P = round((Cp + Gp) / P_portion * 100, 4)

                    C_portion = round(Ac + Tc + Cc + Gc, 4)
                    if C_portion == 0:
                        C_A_usage = C_T_usage = C_C_usage = C_G_usage = C_portion = A2T_C = G2C_C = GC_countent_C = 0
                    else:
                        C_A_usage = round(Ac / C_portion * 100, 4)
                        C_T_usage = round(Tc / C_portion * 100, 4)
                        C_C_usage = round(Cc / C_portion * 100, 4)
                        C_G_usage = round(Gc / C_portion * 100, 4)
                        C_portion = round(C_portion / Length * 100, 4)
                        A2T_C = round(Ac / Tc, 4) if Tc != 0 else 0
                        G2C_C = round(Gc / Cc, 4) if Cc != 0 else 0
                        GC_countent_C = round((Cc + Gc) / C_portion * 100, 4)

                    A2T_anot = round((Ap + Ac) / (Tp + Tc) if (Tp + Tc) != 0 else 0, 4)
                    G2C_anot = round((Gp + Gc) / (Cp + Cc) if (Cp + Cc) != 0 else 0, 4)
                    P2C = round(P_portion / C_portion, 4) if C_portion != 0 else 0

                    PandC = sorted([P_portion, C_portion])
                    smaller_portion_2_Bigger = round(PandC[0] / PandC[1], 4) if PandC[1] != 0 else 0
                    dif2sum = round(abs(PandC[0] - PandC[1]) / ((PandC[0] + PandC[1]) * 2) * 100, 4) if (PandC[0] + PandC[1]) != 0 else 0

                except ZeroDivisionError:
                    portions.append(0)
                    P_A_usage = P_T_usage = P_C_usage = P_G_usage = P_Portion = A2T_P = G2C_P = GC_countent_P = 0
                    C_A_usage = C_T_usage = C_C_usage = C_G_usage = C_portion = A2T_C = G2C_C = GC_countent_C = 0
                    A2T_anot = G2C_anot = P2C = smaller_portion_2_Bigger = dif2sum = 0

                vaiavle_dic = {
                    f'{element}_GC_P': GC_countent_P,
                    f'{element}_GC_C': GC_countent_C,
                    f'{element}_P%': P_Portion,
                    f'{element}_C%': C_portion,
                    f'{element}_P2C': P2C,
                    f'{element}_S2B': smaller_portion_2_Bigger,
                    f'{element}_dif2sum': dif2sum,
                    f'{element}_A2T_ano': A2T_anot,
                    f'{element}_G2C_ano': G2C_anot,
                    f'{element}_P_A_ano': P_A_usage,
                    f'{element}_P_T_ano': P_T_usage,
                    f'{element}_P_C_ano': P_C_usage,
                    f'{element}_P_G_ano': P_G_usage,
                    f'{element}_P_A2T':   A2T_P,
                    f'{element}_P_G2C': G2C_P,
                    f'{element}_C_A_ano': C_A_usage,
                    f'{element}_C_T_ano': C_T_usage,
                    f'{element}_C_C_ano': C_C_usage,
                    f'{element}_C_G_ano': C_G_usage,
                    f'{element}_A2T_C': A2T_C,
                    f'{element}_G2C_C': G2C_C
                }
                main_dict.update(vaiavle_dic)

            # Insert metadata columns
            my_table = pd.DataFrame([main_dict])
            my_table.insert(4, 'Wh_ano', round(sum(portions[:4]), 4))
            my_table.insert(5, 'Int_ano', round(portions[0], 4))
            my_table.insert(6, 'UTR_ano', round(portions[1], 4))
            my_table.insert(7, 'CDS_ano', round(portions[2], 4))
            my_table.insert(8, 'share_ano', round(portions[3], 4))
            my_table.insert(9, 'Un_ano', round(portions[4], 4))
            my_table.insert(10, 'sum', round(sum(portions), 4))

            my_table['Gene_P_count'] = Gene_P_counter
            my_table['Gene_C_count'] = Gene_C_counter
            my_table['Exon_P_count'] = Exon_P_counter
            my_table['Exon_C_count'] = Exon_C_counter
            my_table['CDS_P_count'] = CDS_P_counter
            my_table['CDS_C_count'] = CDS_C_counter

            organism_rows.append(my_table)

        fasta_file.close()

        # If any valid chromosomes were analyzed, append this organism's results to CSV
        if organism_rows:
            combined_df = pd.concat(organism_rows, ignore_index=True)
            # Use file lock/safe append for multiprocessing
            write_header = not os.path.exists(temp_csv_path) or os.path.getsize(temp_csv_path) == 0
            combined_df.to_csv(temp_csv_path, mode='a', header=write_header, index=False)
            print(f"Finished & Saved organism: {organism} ({analysed_chr_counter} chromosomes processed)")
            return combined_df

        return None

    except Exception as e:
        print(f"Error in organism {organism}: {e}")
        return None

# ----------------- Parallel Resumable Controller -----------------

def gene_element_annotation_in_plus_and_minis(all_genomes_folder_path, organ_group):
    print('Analysis started...')
    start_time = time.time()
    start_run = datetime.now().strftime("%Y-%m-%d")
    print(f"Run started at: {start_run}")

    # Set up unique folders but simple output name
    results_folder = os.path.join(str(config.output_dir("GFF_Mapping", organ_group)), f"{organ_group}_refrence_Int_UTR_CDS_{start_run}")
    os.makedirs(results_folder, exist_ok=True)

    # The output filename remains a clean, static name
    output_csv = os.path.join(results_folder, f"{organ_group}_gene_info.csv")

    organism_list = os.listdir(all_genomes_folder_path)

    # Filter size list
    filtered_organism_list = []
    for i in organism_list:
        if i.endswith((".json", ".jsonl")):
            continue
        total_size = 0
        organism_dir = os.path.join(all_genomes_folder_path, i)
        if not os.path.isdir(organism_dir):
            continue
        for dirpath, _, filenames in os.walk(organism_dir):
            for filename in filenames:
                total_size += os.path.getsize(os.path.join(dirpath, filename))
        # (size filter disabled on purpose: process every organism)
        # if (total_size / 1024 ** 2) < 50000:
        filtered_organism_list.append(i)

    # --- RESUME PROGRESS CHECK ---
    processed_organisms = set()
    if os.path.exists(output_csv) and os.path.getsize(output_csv) > 0:
        try:
            existing_df = pd.read_csv(output_csv, usecols=['Name'])
            processed_organisms = set(existing_df['Name'].unique())
            print(f"Resuming run. Found {len(processed_organisms)} already processed organisms in existing results.")
        except Exception as e:
            print(f"Could not read existing results, starting fresh. Info: {e}")

    # Filter out already analyzed organisms
    run_list = [org for org in filtered_organism_list if org not in processed_organisms]
    print(f"Total organisms filtered by size: {len(filtered_organism_list)}")
    print(f"Remaining organisms to process: {len(run_list)}")

    if not run_list:
        print("All organisms have already been processed!")
        return

    # Run organisms in parallel
    max_workers = 4  # Adjusted as per original settings
    print(f"Starting parallel processing on {max_workers} workers...")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_org = {
            executor.submit(
                process_one_organism,
                len(run_list),
                index,
                organ,
                all_genomes_folder_path,
                output_csv
            ): organ for index, organ in enumerate(run_list)
        }

        for future in as_completed(future_to_org):
            org = future_to_org[future]
            try:
                future.result()
            except Exception as exc:
                print(f'{org} generated an exception: {exc}')

    end_time = time.time()
    print('Process time:', round((end_time - start_time) / 60, 2), 'minutes')
    print('Results saved to:', output_csv)
    print('Finished successfully!')

# ----------------- Execution -----------------
if __name__ == '__main__':
    one_Trabite = str(config.require(config.ncbi_data_dir("Chordata"), "Unzip the NCBI genome package of this group into Fasta_GFF_inputs/<Group>/ (see Fasta_GFF_inputs/README.md)."))   # Fasta_GFF_inputs/Chordata/ncbi_dataset/data
    gene_element_annotation_in_plus_and_minis(one_Trabite, 'Chordata')
