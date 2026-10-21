import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)
# parallelized_version_of_your_script.py
import os
import gzip
import pandas as pd
import numpy as np
import time
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from Bio import SeqIO
import sys  # For potential memory estimation if needed




# ----------------- your original functions (kept intact, minimal fix) -----------------
def nucleotide_analysis_for_exon(sequence_plus):
    A_entron_P = sequence_plus.count('Q')
    T_entron_P = sequence_plus.count('W')
    C_entron_P = sequence_plus.count('E')
    G_entron_P = sequence_plus.count('R')
    sum_entron_P = A_entron_P + T_entron_P + C_entron_P + G_entron_P
    A_entron_C = sequence_plus.count('q')
    T_entron_C = sequence_plus.count('w')
    C_entron_C = sequence_plus.count('e')
    G_entron_C = sequence_plus.count('r')
    sum_entron_C = A_entron_C + T_entron_C + C_entron_C + G_entron_C
    A_exon_P = sequence_plus.count('S')
    T_exon_P = sequence_plus.count('D')
    C_exon_P = sequence_plus.count('F')
    G_exon_P = sequence_plus.count('H')
    sum_exon_P = A_exon_P + T_exon_P + C_exon_P + G_exon_P
    A_exon_C = sequence_plus.count('s')
    T_exon_C = sequence_plus.count('d')
    C_exon_C = sequence_plus.count('f')
    G_exon_C = sequence_plus.count('h')
    sum_exon_C = A_exon_C + T_exon_C + C_exon_C + G_exon_C

    A_CDS_P = sequence_plus.count('Z')
    T_CDS_P = sequence_plus.count('X')
    C_CDS_P = sequence_plus.count('V')
    G_CDS_P = sequence_plus.count('B')
    sum_CDS_P = A_CDS_P + T_CDS_P + C_CDS_P + G_CDS_P

    A_CDS_C = sequence_plus.count('z')
    T_CDS_C = sequence_plus.count('x')
    C_CDS_C = sequence_plus.count('v')
    G_CDS_C = sequence_plus.count('b')
    sum_CDS_C = A_CDS_C + T_CDS_C + C_CDS_C + G_CDS_C

    A_shared = sequence_plus.count('A')
    T_shared = sequence_plus.count('T')
    C_shared = sequence_plus.count('C')
    G_shared = sequence_plus.count('G')
    A_unannotated = sequence_plus.count('a')
    T_unannotated = sequence_plus.count('t')
    C_unannotated = sequence_plus.count('c')
    G_unannotated = sequence_plus.count('g')
    sum_unannotated = A_unannotated + T_unannotated + C_unannotated + G_unannotated

    # pandas import kept local to match your style (but imported at file top too)
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
    sequence = sequence.lower()  # atcg all
    sequence_plus = bytearray(sequence.encode('ascii'))
    Inreon_P_counter = 0
    Inreon_C_counter = 0
    UTR_P_counter = 0
    UTR_C_counter = 0
    CDS_P_counter = 0
    CDS_C_counter = 0

    df_all.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    filtered_df_all = df_all[df_all['seqid'].str.startswith(ch_id_split)]

    # First pass: genes (introns)
    filtered_df = filtered_df_all[filtered_df_all['type'].str.startswith('gene')]
    for annotation in filtered_df.to_numpy():
        start = int(annotation[3])
        end = int(annotation[4])
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
        if annotation[6] == '-':
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

    # Second pass: exons (UTRs)
    filtered_df = filtered_df_all[filtered_df_all['type'].str.startswith('exon')]
    for annotation in filtered_df.to_numpy():
        start = int(annotation[3])
        end = int(annotation[4])
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
        if annotation[6] == '-':
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
        start = int(annotation[3])
        end = int(annotation[4])
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

                elif ch == ord('Q'): sl[i] = ord('Z')
                elif ch == ord('W'): sl[i] = ord('X')
                elif ch == ord('E'): sl[i] = ord('V')
                elif ch == ord('R'): sl[i] = ord('B')
                elif ch == ord('q'): sl[i] = ord('A')
                elif ch == ord('w'): sl[i] = ord('T')
                elif ch == ord('e'): sl[i] = ord('C')
                elif ch == ord('r'): sl[i] = ord('G')






            sequence_plus[(start - 1):end] = sl
        if annotation[6] == '-':
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

                elif ch == ord('q'): sl[i] = ord('z')
                elif ch == ord('w'): sl[i] = ord('x')
                elif ch == ord('e'): sl[i] = ord('v')
                elif ch == ord('r'): sl[i] = ord('b')
                elif ch == ord('Q'): sl[i] = ord('A')
                elif ch == ord('W'): sl[i] = ord('T')
                elif ch == ord('E'): sl[i] = ord('C')
                elif ch == ord('R'): sl[i] = ord('G')





            sequence_plus[(start - 1):end] = sl

    sequence_plus = sequence_plus.decode('ascii')
    return sequence_plus, Inreon_P_counter, Inreon_C_counter, UTR_P_counter, UTR_C_counter, CDS_P_counter, CDS_C_counter

# ----------------- parallel worker: exact per-organism body kept -----------------
def process_one_organism(final_number,organ_number,organism, all_genomes_folder_path):
    #print('start_procces')
    """
    Worker: runs the original per-organism loop body for one organism.
    Returns: the my_table DataFrame for that organism (or None).
    """
    try:
        if organism.endswith(".json") or organism.endswith(".jsonl"):
            return None

        organism_counter_local = None
        # replicate the local environment of original script
        files = os.listdir(all_genomes_folder_path + '/' + organism)
        fna = [i for i in files if i.endswith('.fna') and not i.startswith('cds')]
        gff = [k for k in files if k.endswith('.gff')]
        csv = [p for p in files if p.endswith('.csv')]
        fna.sort()
        gff.sort()

        try:
            fasta_file = open(all_genomes_folder_path + '/' + organism + '/' + fna[0], 'rt')
        except (IOError, EOFError, IndexError):
            return None

        # Load GFF once per organism to save memory (avoid reloading multiple times per chromosome)
        gff_path = all_genomes_folder_path + '/' + organism + '/' + gff[0]
        if gff_path.endswith('.gz'):
            with gzip.open(gff_path, 'rt') as gff_file:
                df_all = pd.read_csv(gff_file, sep='\t', comment='#')
        else:
            with open(gff_path, 'r') as gff_file:
                df_all = pd.read_csv(gff_file, sep='\t', comment='#')

        # reproduce all lists and counters exactly as in your original per-organism body
        A2T_list = []
        G2C_list = []
        P2C_number_list = [[], [], [], [], []]
        GC_countent_list_W = []
        GC_countent_list_P = [[], [], [], [], []]
        GC_countent_list_C = [[], [], [], [], []]
        len_seq_list = []
        plus_and_minis_strand_portion = [[], [], [], [], []]
        organims_name = []
        A2T_anot_list = [[], [], [], [], []]
        G2C_anot_list = [[], [], [], [], []]
        P_element_number_list = [[], [], [], [], []]
        C_element_number_list = [[], [], [], [], []]

        P_A_usage_list = [[], [], [], [], []]
        P_T_usage_list = [[], [], [], [], []]
        P_C_usage_list = [[], [], [], [], []]
        P_G_usage_list = [[], [], [], [], []]

        C_A_usage_list = [[], [], [], [], []]
        C_T_usage_list = [[], [], [], [], []]
        C_C_usage_list = [[], [], [], [], []]
        C_G_usage_list = [[], [], [], [], []]
        smaller_portion_2_Bigger = [[], [], [], [], []]
        dif2sum = [[], [], [], [], []]
        P_Portion_list = [[], [], [], [], []]
        C_portion_list = [[], [], [], [], []]
        P2C_list = [[], [], [], [], []]
        P2C_list_W = [[], [], [], [], []]

        A2T_P_list = [[], [], [], [], []]
        G2C_P_list = [[], [], [], [], []]
        A2T_C_list = [[], [], [], [], []]
        G2C_C_list = [[], [], [], [], []]

        chr_id_list = []
        All_annotation = []
        Int_annotaion = []
        UTR_annotation = []
        CDS_annotation = []
        Shared_annotation = []
        Un_annotation = []
        Sum_E = []

        Gene_P_counter_list = []
        Gene_C_counter_list = []
        Exon_P_counter_list = []
        Exon_C_counter_list = []
        CDS_P_counter_list = []
        CDS_C_counter_list = []

        analysed_chr_counter = 0
        lablename = 'refrence'
        organ_group = "fungis"

        # organism_counter is only for printing original-like info;
        # we don't have the global organism_counter here, so set a placeholder.
        # The original script used organism_counter incremented in outer loop; keep approximate info.
        # We'll set organism_counter_local to 0 so printed lines are still informative.
        organism_counter_local = 0

        my_table = None

        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if 'mitochondrion' in seq_record.description:
                continue
            analysed_chr_counter += 1
            seq = str(seq_record.seq)
            l = len(seq) / 1000000
            if l < 0.05:
                continue
            if "chromosome" not in seq_record.description:
                continue
            ch_id = seq_record.id
            if ch_id.startswith('NW'):
                continue

            chr_number = (seq_record.description.split())[0]

            # Optional memory check: Skip if sequence is too large (bytearray uses ~1 byte per base)
            estimated_mem_gb = len(seq) / (1024 ** 3)
            if estimated_mem_gb > 31:
                print(f"Skipping chromosome {ch_id} (organism: {organism}) due to high memory requirement: {estimated_mem_gb:.2f} GB estimated")
                continue

            exon_mask_seq, Gene_P_counter, Gene_C_counter, Exon_P_counter, Exon_C_counter, CDS_P_counter, CDS_C_counter = \
                exon_and_intron_masking(seq, df_all, ch_id)
            if CDS_P_counter == 0 or CDS_C_counter == 0:
                continue

            Gene_P_counter_list.append(Gene_P_counter)
            Gene_C_counter_list.append(Gene_C_counter)
            Exon_P_counter_list.append(Exon_P_counter)
            Exon_C_counter_list.append(Exon_C_counter)
            CDS_P_counter_list.append(CDS_P_counter)
            CDS_C_counter_list.append(CDS_C_counter)

            A_W = seq.count('A') + seq.count("a")
            T_W = seq.count('T') + seq.count("t")
            C_W = seq.count('C') + seq.count("c")
            G_W = seq.count("G") + seq.count("g")
            N_W = seq.count('N') + seq.count("n")
            Length = len(seq) - N_W
            organims_name.append(organism)
            # protect division by zero (original assumes non-zero denominator). Keep behavior but safe:
            denom_W = (A_W + T_W + C_W + G_W)
            if denom_W == 0:
                GC_countent_list_W.append(0)
            else:
                GC_countent_list_W.append(round((G_W + C_W) / denom_W * 100, 4))
            chr_id_list.append(ch_id)
            len_seq_list.append(round(Length / 1000000, 6))
            # avoid division by zero for T_W and C_W
            if T_W == 0:
                A2T_list.append(0)
            else:
                A2T_list.append(round((A_W) / (T_W), 4))
            if C_W == 0:
                G2C_list.append(0)
            else:
                G2C_list.append(round((G_W) / (C_W), 4))

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
            main_dict = {'Name': organims_name,
                         'Chr_ID': chr_id_list,
                         'Length(M)': len_seq_list,
                         'GC_W': GC_countent_list_W,
                         'A2T_W': A2T_list,
                         'G2C_W': G2C_list,}
            for e_number, element in enumerate(gene_elements):
                if element == "Share" or element == "Unanot":
                    Ap, Tp, Cp, Gp = df.iloc[counter + 0]
                    Ap, Tp, Cp, Gp = Ap / 2, Tp / 2, Cp / 2, Gp / 2
                    Ac, Tc, Cc, Gc = df.iloc[counter + 0]
                    Ac, Tc, Cc, Gc = Ac / 2, Tc / 2, Cc / 2, Gc / 2
                    counter += 1
                else:
                    Ap, Tp, Cp, Gp = df.iloc[counter + 0]
                    Ac, Tc, Cc, Gc = df.iloc[counter + 1]
                    counter += 2
                try:
                    u = 1 / Ap + 1 / Tp + 1 / Cp + 1 / Gp
                    u = 1 / Ac + 1 / Tc + 1 / Cc + 1 / Gc
                    portions.append(sum([Ap, Tp, Cp, Gp, Ac, Tc, Cc, Gc]) / Length * 100)

                    P_portion = round(Ap + Tp + Cp + Gp, 4)
                    # protect division by zero inside these calculations; original code would raise;
                    # we match behavior but return zeros when division by zero would occur
                    if P_portion == 0:
                        P_A_usage_list[e_number].append(0)
                        P_T_usage_list[e_number].append(0)
                        P_C_usage_list[e_number].append(0)
                        P_G_usage_list[e_number].append(0)
                        P_Portion_list[e_number].append(0)
                        A2T_P_list[e_number].append(0)
                        G2C_P_list[e_number].append(0)
                        GC_countent_list_P[e_number].append(0)
                    else:
                        P_A_usage_list[e_number].append(round(Ap / P_portion * 100, 4))
                        P_T_usage_list[e_number].append(round(Tp / P_portion * 100, 4))
                        P_C_usage_list[e_number].append(round(Cp / P_portion * 100, 4))
                        P_G_usage_list[e_number].append(round(Gp / P_portion * 100, 4))
                        P_Portion_list[e_number].append(round(P_portion / Length * 100, 4))
                        # protect Tp==0
                        if Tp == 0:
                            A2T_P_list[e_number].append(0)
                        else:
                            A2T_P_list[e_number].append(round(Ap / Tp, 4))
                        if Cp == 0:
                            G2C_P_list[e_number].append(0)
                        else:
                            G2C_P_list[e_number].append(round(Gp / Cp, 4))
                        P_element_number_list[e_number].append(CDS_P_counter)
                        GC_countent_list_P[e_number].append(round((Cp + Gp) / P_portion * 100, 4))

                    C_portion = round(Ac + Tc + Cc + Gc, 4)
                    if C_portion == 0:
                        C_A_usage_list[e_number].append(0)
                        C_T_usage_list[e_number].append(0)
                        C_C_usage_list[e_number].append(0)
                        C_G_usage_list[e_number].append(0)
                        C_portion_list[e_number].append(0)
                        A2T_C_list[e_number].append(0)
                        G2C_C_list[e_number].append(0)
                        GC_countent_list_C[e_number].append(0)
                    else:
                        C_A_usage_list[e_number].append(round(Ac / C_portion * 100, 4))
                        C_T_usage_list[e_number].append(round(Tc / C_portion * 100, 4))
                        C_C_usage_list[e_number].append(round(Cc / C_portion * 100, 4))
                        C_G_usage_list[e_number].append(round(Gc / C_portion * 100, 4))
                        C_portion_list[e_number].append(round(C_portion / Length * 100, 4))
                        if Tc == 0:
                            A2T_C_list[e_number].append(0)
                        else:
                            A2T_C_list[e_number].append(round(Ac / Tc, 4))
                        if Cc == 0:
                            G2C_C_list[e_number].append(0)
                        else:
                            G2C_C_list[e_number].append(round(Gc / Cc, 4))
                        C_element_number_list[e_number].append(CDS_C_counter)
                        GC_countent_list_C[e_number].append(round((Cc + Gc) / C_portion * 100, 4))

                    A2T_anot_list[e_number].append(round((Ap + Ac) / (Tp + Tc) if (Tp + Tc) != 0 else 0, 4))
                    G2C_anot_list[e_number].append(round((Gp + Gc) / (Cp + Cc) if (Cp + Cc) != 0 else 0, 4))

                    P2C_list[e_number].append(round(P_portion / C_portion, 4) if C_portion != 0 else 0)
                    P2C_number_list[e_number].append(round(CDS_P_counter / CDS_C_counter, 4) if CDS_C_counter != 0 else 0)
                    PandC = [P_portion, C_portion]
                    PandC.sort()
                    if PandC[1] == 0:
                        smaller_portion_2_Bigger[e_number].append(0)
                    else:
                        smaller_portion_2_Bigger[e_number].append(round(PandC[0] / PandC[1], 4))
                    dif2sum[e_number].append(round(abs(PandC[0] - PandC[1]) / ((PandC[0] + PandC[1]) * 2) * 100, 4) if (PandC[0] + PandC[1]) != 0 else 0)

                except ZeroDivisionError:
                    portions.append(0)
                    P_portion = round(Ap + Tp + Cp + Gp, 4)
                    P_A_usage_list[e_number].append(0)
                    P_T_usage_list[e_number].append(0)
                    P_C_usage_list[e_number].append(0)
                    P_G_usage_list[e_number].append(0)
                    P_Portion_list[e_number].append(0)
                    A2T_P_list[e_number].append(0)
                    G2C_P_list[e_number].append(0)
                    P_element_number_list[e_number].append(0)
                    GC_countent_list_P[e_number].append(0)

                    C_portion = round(Ac + Tc + Cc + Gc, 4)
                    C_A_usage_list[e_number].append(0)
                    C_T_usage_list[e_number].append(0)
                    C_C_usage_list[e_number].append(0)
                    C_G_usage_list[e_number].append(0)
                    C_portion_list[e_number].append(0)
                    A2T_C_list[e_number].append(0)
                    G2C_C_list[e_number].append(0)
                    C_element_number_list[e_number].append(0)
                    GC_countent_list_C[e_number].append(0)

                    A2T_anot_list[e_number].append(0)
                    G2C_anot_list[e_number].append(0)

                    P2C_list[e_number].append(0)
                    P2C_number_list[e_number].append(0)
                    PandC = [P_portion, C_portion]
                    PandC.sort()
                    smaller_portion_2_Bigger[e_number].append(0)
                    dif2sum[e_number].append(0)

                vaiavle_dic = {
                    f'{element}_GC_P': GC_countent_list_P[e_number],
                    f'{element}_GC_C': GC_countent_list_C[e_number],
                    f'{element}_P%': P_Portion_list[e_number],
                    f'{element}_C%': C_portion_list[e_number],
                    f'{element}_P2C': P2C_list[e_number],
                    f'{element}_S2B': smaller_portion_2_Bigger[e_number],
                    f'{element}_dif2sum': dif2sum[e_number],
                    f'{element}_A2T_ano': A2T_anot_list[e_number],
                    f'{element}_G2C_ano': G2C_anot_list[e_number],
                    f'{element}_P_E_num': P_element_number_list[e_number],
                    f'{element}_C_E_num': C_element_number_list[e_number],
                    f'{element}_P_A_ano': P_A_usage_list[e_number],
                    f'{element}_P_T_ano': P_T_usage_list[e_number],
                    f'{element}_P_C_ano': P_C_usage_list[e_number],
                    f'{element}_P_G_ano': P_G_usage_list[e_number],
                    f'{element}_P_A2T': A2T_P_list[e_number],
                    f'{element}_P_G2C': G2C_P_list[e_number],
                    f'{element}_C_A_ano': C_A_usage_list[e_number],
                    f'{element}_C_T_ano': C_T_usage_list[e_number],
                    f'{element}_C_C_ano': C_C_usage_list[e_number],
                    f'{element}_C_G_ano': C_G_usage_list[e_number],
                    f'{element}_A2T_C': A2T_C_list[e_number],
                    f'{element}_G2C_C': G2C_C_list[e_number]
                }

                main_dict = main_dict | vaiavle_dic
            All_annotation.append(round(sum(portions[:4]), 4))
            Int_annotaion.append(round(portions[0], 4))
            UTR_annotation.append(round(portions[1], 4))
            CDS_annotation.append(round(portions[2], 4))
            Shared_annotation.append(round(portions[3], 4))
            Un_annotation.append(round(portions[4], 4))
            Sum_E.append(round(sum(portions), 4))

            try:
                my_table = pd.DataFrame(main_dict)
            except ValueError:
                for i in main_dict:
                    print('hi', i, len(main_dict[i]))
            my_table.insert(4, 'Wh_ano', All_annotation)
            my_table.insert(5, 'Int_ano', Int_annotaion)
            my_table.insert(6, 'UTR_ano', UTR_annotation)
            my_table.insert(7, 'CDS_ano', CDS_annotation)
            my_table.insert(8, 'share_ano', Shared_annotation)
            my_table.insert(9, 'Un_ano', Un_annotation)
            my_table.insert(10, 'sum', Sum_E)
            my_table['Gene_P_count'] = Gene_P_counter_list
            my_table['Gene_C_count'] = Gene_C_counter_list
            my_table['Exon_P_count'] = Exon_P_counter_list
            my_table['Exon_C_count'] = Exon_C_counter_list
            my_table['CDS_P_count'] = CDS_P_counter_list
            my_table['CDS_C_count'] = CDS_C_counter_list

            column_width = 15
            header = " ".join(col.ljust(column_width) for col in my_table.columns)
            last_row_values = " ".join(str(val).ljust(column_width) for val in my_table.iloc[-1].values)
            #if analysed_chr_counter == 1:
                #print(header)
                #print("-" * len(header))
            print(final_number,organ_number, analysed_chr_counter)#, last_row_values)
            #here it should be saved----------------------------------------------------------------------
        # finished parsing fasta for this organism, return the DataFrame (or None if nothing)
        fasta_file.close()  # Explicitly close the fasta file
        return my_table

    except Exception as e:
        print(f"Error in organism {organism}: {e}")
        return None

# ----------------- parallel controller (keeps your outer behavior but parallel) -----------------
def gene_element_annotation_in_plus_and_minis(all_genomes_folder_path,organ_group):
    print('inalysis_started')
    start_time = time.time()
    current_time = datetime.now()
    start_run = current_time.strftime("%Y-%m-%d_%H-%M-%S")
    print(start_run)

    lablename = 'refrence'

    #refrence_list = list(refrence_list['Assembly Accession'])

    organism_list = os.listdir(all_genomes_folder_path)
    organism_number = len(organism_list)
    organism_counter = 0

    # run organisms in parallel
    results = []
    max_workers = 8    # multiprocessing.cpu_count()
    print(f"Starting parallel processing on {max_workers} workers...")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_org = {executor.submit(process_one_organism, organism_number,n,org, all_genomes_folder_path): org for n,org in enumerate(organism_list) if not (org.endswith(".json") or org.endswith(".jsonl"))}
        for future in as_completed(future_to_org):
            org = future_to_org[future]
            organism_counter += 1
            try:
                res = future.result()
                if isinstance(res, pd.DataFrame):
                    results.append(res)
                else:
                    # res may be None
                    pass
            except Exception as exc:
                print(f'{org} generated an exception: {exc}')

    end_time = time.time()
    print('procces time', round((end_time - start_time) / 60, 2), 'minit')

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    results_folder = os.path.join(str(config.output_dir("GFF_Mapping", organ_group)), f"{organ_group}_{lablename}_Int_UTR_CDS_{timestamp}")
    os.makedirs(results_folder, exist_ok=True)

    if results:
        try:
            final_df = pd.concat(results, ignore_index=True)
            output_filename = os.path.join(results_folder, f"{organ_group}_{timestamp}.csv")
            final_df.to_csv(output_filename, index=False)
            print(f"Results saved to: {output_filename}")
        except Exception as e:
            print("Error saving final Excel:", e)
    else:
        print("No results to save.")

    print('start running', start_run)
    print('end running  ', timestamp)
    print('procces finished')

# ----------------- run (keep your original variable name) -----------------
one_Trabite = str(config.require(config.ncbi_data_dir("Bacteria"), "Unzip the NCBI genome package of this group into Fasta_GFF_inputs/<Group>/ (see Fasta_GFF_inputs/README.md)."))   # Fasta_GFF_inputs/Bacteria/ncbi_dataset/data
if __name__ == '__main__':
    gene_element_annotation_in_plus_and_minis(one_Trabite,'Bacteria')